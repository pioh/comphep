/*
 * Lite LHAPDF mdl update - pure C, no LHAPDF library dependency.
 * Scans the filesystem to discover installed LHAPDF 6 PDF sets
 * and updates models/strfuns.mdl and models/beams.mdl.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#ifdef LHAPDF

#define MAX_PDFSETS 2048
#define MAX_LINES   512

typedef struct {
  char name[256];
  int nmembers;
} pdfset_info;

static char mdl_dir[1024];

static int cmp_pdfsets (const void *a, const void *b) {
  return strcmp (((const pdfset_info *) a)->name,
                ((const pdfset_info *) b)->name);
}

/* Read LHAPDF datadir path from .lhapdfpath */
static int read_datadir (char *buf, int bufsize) {
  const char *paths[] = {".lhapdfpath", "../.lhapdfpath",
                         "../../.lhapdfpath", NULL};
  int i;
  for (i = 0; paths[i]; i++) {
    FILE *f = fopen (paths[i], "r");
    if (f) {
      if (fscanf (f, "%1023s", buf) == 1) {
        fclose (f);
        return 1;
      }
      fclose (f);
    }
  }
  return 0;
}

/* Find the models directory (works from WDIR or WDIR/results) */
static int find_modelsdir (void) {
  FILE *f;
  f = fopen ("models/strfuns.mdl", "r");
  if (f) { fclose (f); strcpy (mdl_dir, "models"); return 1; }
  f = fopen ("../models/strfuns.mdl", "r");
  if (f) { fclose (f); strcpy (mdl_dir, "../models"); return 1; }
  return 0;
}

/* Scan LHAPDF datadir for installed PDF sets */
static int scan_pdfsets (const char *datadir, pdfset_info *sets, int maxsets) {
  DIR *dir;
  struct dirent *entry;
  int count = 0;

  dir = opendir (datadir);
  if (!dir) return 0;

  while ((entry = readdir (dir)) != NULL && count < maxsets) {
    char infopath[1024];
    FILE *f;

    if (entry->d_name[0] == '.') continue;

    snprintf (infopath, sizeof (infopath), "%s/%s/%s.info",
              datadir, entry->d_name, entry->d_name);
    f = fopen (infopath, "r");
    if (!f) continue;

    strncpy (sets[count].name, entry->d_name, 255);
    sets[count].name[255] = 0;
    sets[count].nmembers = 1;

    {
      char line[512];
      while (fgets (line, sizeof (line), f)) {
        int nm;
        if (sscanf (line, "NumMembers: %d", &nm) == 1) {
          sets[count].nmembers = nm;
          break;
        }
      }
    }

    fclose (f);
    count++;
  }

  closedir (dir);

  if (count > 1)
    qsort (sets, count, sizeof (pdfset_info), cmp_pdfsets);

  return count;
}

/* Check if setname is in available sets */
static int is_set_available (const char *setname,
                             const pdfset_info *sets, int nsets) {
  int i;
  for (i = 0; i < nsets; i++) {
    if (strcmp (sets[i].name, setname) == 0) return 1;
  }
  return 0;
}

/* Update strfuns.mdl with available PDF sets */
static int update_strfuns (const pdfset_info *sets, int nsets) {
  FILE *fin, *fout;
  char path[1024];
  char line[1024];
  char *kept_lines[MAX_LINES];
  int n_kept = 0;
  char header1[1024], header2[1024];
  int have_header1 = 0, have_header2 = 0;
  int i, num;

  header1[0] = 0;
  header2[0] = 0;

  snprintf (path, sizeof (path), "%s/strfuns.mdl", mdl_dir);
  fin = fopen (path, "r");
  if (!fin) return 0;

  while (fgets (line, sizeof (line), fin)) {
    if (line[0] == '=' && line[1] == '=') break;
    if (!have_header1 && strstr (line, "Strfuns")) {
      strcpy (header1, line);
      have_header1 = 1;
      continue;
    }
    if (!have_header2 && strstr (line, "Name") && strstr (line, "Number")) {
      strcpy (header2, line);
      have_header2 = 1;
      continue;
    }
    if (strstr (line, "LHA:")) continue;
    if (have_header2 && n_kept < MAX_LINES) {
      kept_lines[n_kept] = strdup (line);
      n_kept++;
    }
  }
  fclose (fin);

  fout = fopen (path, "w");
  if (!fout) {
    for (i = 0; i < n_kept; i++) free (kept_lines[i]);
    return 0;
  }

  if (header1[0]) fprintf (fout, "%s", header1);
  if (header2[0]) fprintf (fout, "%s", header2);

  for (i = 0; i < n_kept; i++) {
    fprintf (fout, "%s", kept_lines[i]);
    free (kept_lines[i]);
  }

  num = 100;
  for (i = 0; i < nsets; i++) {
    char name_buf[512];
    snprintf (name_buf, sizeof (name_buf), "LHA:%s(proton)", sets[i].name);
    fprintf (fout, "%-41s|%-14d|%-16d|%-18d\n", name_buf, num, 0, 0);
    num++;
    snprintf (name_buf, sizeof (name_buf), "LHA:%s(anti-proton)", sets[i].name);
    fprintf (fout, "%-41s|%-14d|%-16d|%-18d\n", name_buf, num, 0, 0);
    num++;
  }

  fprintf (fout, "==========================================================================================\n");
  fclose (fout);

  return 1;
}

/* Update beams.mdl - validate LHA: references against installed sets */
static int update_beams (const pdfset_info *sets, int nsets) {
  FILE *fin, *fout;
  char path[1024];
  char line[1024];
  char *lines[MAX_LINES];
  int nlines = 0;
  int i, modified = 0;

  if (nsets <= 0) return 0;

  snprintf (path, sizeof (path), "%s/beams.mdl", mdl_dir);
  fin = fopen (path, "r");
  if (!fin) return 0;

  while (fgets (line, sizeof (line), fin) && nlines < MAX_LINES) {
    lines[nlines] = strdup (line);
    nlines++;
  }
  fclose (fin);

  for (i = 0; i < nlines; i++) {
    char *p;
    int pipes = 0;

    /* Find 4th column (after 3rd pipe) */
    p = lines[i];
    while (*p && pipes < 3) {
      if (*p == '|') pipes++;
      p++;
    }
    if (pipes < 3) continue;

    /* Check if this column starts with LHA: */
    if (strncmp (p, "LHA:", 4) != 0) continue;

    /* Extract setname from LHA:<setname>(<particle>)... */
    {
      char setname[256];
      char *start = p + 4;
      char *paren = strchr (start, '(');
      int len;

      if (!paren) continue;
      len = paren - start;
      if (len <= 0 || len >= 256) continue;
      strncpy (setname, start, len);
      setname[len] = 0;

      if (!is_set_available (setname, sets, nsets)) {
        char newline[1024];
        int prefix_len = (p + 4) - lines[i];

        snprintf (newline, sizeof (newline), "%.*s%s%s",
                  prefix_len, lines[i], sets[0].name, paren);

        free (lines[i]);
        lines[i] = strdup (newline);
        modified = 1;
      }
    }
  }

  if (!modified) {
    for (i = 0; i < nlines; i++) free (lines[i]);
    return 1;
  }

  fout = fopen (path, "w");
  if (!fout) {
    for (i = 0; i < nlines; i++) free (lines[i]);
    return 0;
  }

  for (i = 0; i < nlines; i++) {
    fprintf (fout, "%s", lines[i]);
    free (lines[i]);
  }

  fclose (fout);
  return 1;
}

/* Main entry: scan LHAPDF datadir, update strfuns.mdl and beams.mdl */
int update_lhapdf_mdl (void) {
  char datadir[1024];
  pdfset_info *sets;
  int nsets;

  if (!read_datadir (datadir, sizeof (datadir))) return 0;
  if (!find_modelsdir ()) return 0;

  sets = (pdfset_info *) malloc (MAX_PDFSETS * sizeof (pdfset_info));
  if (!sets) return 0;

  nsets = scan_pdfsets (datadir, sets, MAX_PDFSETS);
  if (nsets <= 0) { free (sets); return 0; }

  update_strfuns (sets, nsets);
  update_beams (sets, nsets);

  free (sets);
  return 1;
}

#endif /* LHAPDF */
