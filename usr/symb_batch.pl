#!/usr/bin/perl -w
# Copyright (C) 2002-2009, CompHEP Collaboration
# Author: Alexander Sherstnev
# v. 0.1 created by Sherstnev A., date: 29/12/2002
# v. 0.2 
# v. 0.5 
# v. 1.0 
# v. 1.0  29.11.2007 1) renewed to work with the current CompHEP
#                    2) status option added
#                    3) multi-CPU calculations added
#                    4) version output
# v. 1.1  18.04.2008 1)  nice option added 
#                    2) ...

use Getopt::Long;
use IPC::Open3;
use File::Copy;
use POSIX qw(:sys_wait_h :signal_h :errno_h);

help() unless (@ARGV or -e "process.dat" );

#####################################################################################
# Take the input options
GetOptions(
           "n=s",       \$opt_n,                # 
           "d=s",       \$opt_d,                # set working directory name (insted of results)
           "recovery",  \$recovery,             # 
           "relink",    \$relink,               # 
           "show=s",    \$show,                 # show status/diagrams
           "mp=i",      \$mp,                   # set the number of parallel subprocesses
           "nice:s",    \$nnice,                # set priority for nice regime
           "h",         \$opt_h,                # short help
           "help",      \$longhelp,             # long help
           "version",   \$veropt,               # launch n_comphep in same mode
           "verbose",   \$verbose               # print more information 
           );

  my $filedat = "process.dat";
  my $rsltdir = "results";
  my $tmp_stor_dir=$rsltdir."_tmp_1";
  my $version = 1.1;
  my $versiondate = "18/04/2008";
  my $nice = "";

###################################################################################################
# print the script version
  version() if (defined $veropt);

###################################################################################################
# Print help message
  if($opt_h or $longhelp) {
    help ();
    exit (0);
  }

###################################################################################################
# check whether s_comphep.exe is available
  $chepname = ".compheppath";
  $comphep_dir = "tmp";
  open (DAT, "<$chepname") or die "\n*** there are no $chepname file";
  chomp ($comphep_dir = <DAT>);
  close DAT;
  $s_comphep = $comphep_dir."/bin/s_comphep.exe";
  if (!-e $s_comphep) {
    exit(0);
    print "\nsymb_batch.pl (error): there is no s_comphep.exe in $comphep_dir/bin. At first you must prepare\n";
    print "                         At first run ./configure and compile CompHEP!\n";
    exit (0);
  }

###################################################################################################
# set different names: $opt_n instead of process.dat & $opt_d instead of results
  $filedat = $opt_n if (defined $opt_n);
  if (defined $opt_d) {
    while(-e $tmp_stor_dir) {
      $i++;
      $tmp_stor_dir=$rsltdir."_tmp"."_$i";
    }
    rename ($rsltdir, $tmp_stor_dir);
    mkdir ($rsltdir, 0755);
    copy ("$tmp_stor_dir/Makefile", "$rsltdir/Makefile");
    copy ("$tmp_stor_dir/n_comphep","$rsltdir/n_comphep");
    copy ("$tmp_stor_dir/diag_view","$rsltdir/diag_view");
    chmod (0744,"$rsltdir/n_comphep");
    chmod (0744,"$rsltdir/diag_view");
  }

  if (-e "symb_batch.log") {
    rename("symb_batch.log", "symb_batch.log_old");
  }

###################################################################################################
# auxiliary operations
#1
  if (defined $recovery)
  {
    data_analysis();
    recovery();
    rename_results();
    exit(0);
  }

#2
  if (defined $relink) {
    relink();
    exit(0);
  }

#3
  if(defined $show) {
    show_check_names ();
  }

# 4
  if(defined $show and $show eq "stat") {
    show_status();
    exit(0);
  }

# 5
# Set nice priority for n_comphep
  if (defined $nnice) {
    if ($nnice eq "") {
      $nice = "nice -n 10";
    } else {
      if ($nnice gt "19") {
        print "symb_batch.pl (warning): $nnice is too high, I'll use 19";
        $nnice = 19;
      }
      $nice = "nice -n $nnice";
    }
  }

###################################################################################################
# check whether $rsltdir is empty. If not rename it.
  opendir(RST, $rsltdir) or die "\n*** there are no $rsltdir directory";
  $empty = 0;
  while ($f = readdir(RST)) { 
    $empty=1 unless ($f eq "Makefile" or $f eq "n_comphep" or $f eq "diag_view" or $f eq "." or $f eq "..");
  }
  if((-e $rsltdir) and ($empty == 1))
  {
    $i = 1;
    $tmp_rsltdir="$rsltdir"."_old_1";
    while(-e $tmp_rsltdir)
    {
      $i++;
      $tmp_rsltdir="$rsltdir"."_old"."_$i";
    }
    rename ($rsltdir, $tmp_rsltdir);
    mkdir ($rsltdir, 0755);
    copy ("$tmp_rsltdir/Makefile",  "$rsltdir/Makefile");
    copy ("$tmp_rsltdir/n_comphep", "$rsltdir/n_comphep");
    copy ("$tmp_rsltdir/diag_view", "$rsltdir/diag_view");
    chmod (0744,"$rsltdir/n_comphep");
    chmod (0744,"$rsltdir/diag_view");
  }

###################################################################################################
# Combine the results of previous calculations from the temporary PBS directories (this option
# is needed only if nocombine option was used before)
#  combine() if defined $opt_combine;

###################################################################################################
# main operation
  data_analysis ();
  form_safefile ();
  run ();
  rename_results ();
  exit(0);


###################################################################################################
###################################################################################################
#                                  subproutines
###################################################################################################
###################################################################################################

###################################################################################################
#
sub rename_results {
   if (defined $opt_d) {
     if (-e $opt_d) {
       $i=1;
       $old_dir=$opt_d."_old_1";
       while(-e $old_dir) {
         $i++;
         $old_dir=$opt_d."_old"."_$i";
       }
       rename ($opt_d, $old_dir);
     }
     rename ($rsltdir, $opt_d);
     rename ($tmp_stor_dir, $rsltdir);
#     rename($tmp_rsltdir, $rsltdir) if(defined $tmp_rsltdir);
   }
}

###################################################################################################
#
sub data_analysis {
  $excld = "";
  $sqexcld = "";
  open(DAT,"<$filedat") or die "\n*** there are no $filedat file";
  while($res=<DAT>)
  {
    next if( $res=~/^#/ or $res=~/^\s+/ );
    if ($res=~/[:]/)
    {
      $s=$';
      chomp $s;
      if ($`=~/(\s*)model(\s*)number(\s*)/)                             {$mdlnum =$s;}
      if ($`=~/(\s*)final(\s*)state(\s*)/)                              {$finalstate=$s;}
      if ($`=~/(\s*)beam(\s*)1(\s*)/)                                   {$beam1=$s;}
      if ($`=~/(\s*)beam(\s*)2(\s*)/)                                   {$beam2=$s;}
      if ($`=~/(\s*)strfun(\s*)1(\s*)/)                                 {$strfun1=$s;}
      if ($`=~/(\s*)strfun(\s*)2(\s*)/)                                 {$strfun2=$s;}
      if ($`=~/(\s*)beam(\s*)energy(\s*)1(\s*)/)                        {$energy1=$s;}
      if ($`=~/(\s*)beam(\s*)energy(\s*)2(\s*)/)                        {$energy2=$s;}
      if ($`=~/(\s*)exclude(\s*)diagrams(\s*)with(\s*)/)                {$exclde=$s;}
      if ($`=~/(\s*)keep(\s*)diagrams(\s*)with(\s*)/)                   {$keep=$s;}
      if ($`=~/(\s*)make(\s*)symbolic(\s*)calculations(\s*)\(yes\/no\)/){$smake =$s;}
      if ($`=~/(\s*)make(\s*)n_comphep(\s*)generator(\s*)\(yes\/no\)/)  {$nmake =$s;}
      if ($`=~/(\s*)excluded(\s*)diagrams(\s*)/)                        {$excld=$excld.$s}
      if ($`=~/(\s*)excluded(\s*)squared(\s*)diagrams(\s*)/)            {$sqexcld=$sqexcld.$s}
    }
  }
  close DAT;

#  print "excluded diagrams: $excld\n";
#  print "excluded squared diagrams: $sqexcld\n";

  if(!defined $mdlnum) {print "Error: model number is not defined\n;";exit(0);}
  if(!defined $finalstate){print "Error: process is not entered\n;";exit(0);}
  if(!defined $energy1 or !defined $energy2) {print "Error: energy is not defined\n;";exit(0);}
}

###################################################################################################
#
sub form_safefile {
  $sffl="safe_tmp";
  if(-e $sffl)
  {
    print "Warning! File $sffl exists, it is renamed to $sffl._old\n";
    rename($sffl,$sffl."_old");
  }
  open(SAF,"+>$sffl") or die "I can not open safe_open file ($!)";
  print SAF "#Model $mdlnum\n";
  print SAF "#nIn 0\n";
  print SAF "#nOut 0\n";
  print SAF "#beam_1$beam1\n";
  print SAF "#beam_2$beam2\n";
  print SAF "#beam_energy_1$energy1\n";
  print SAF "#beam_energy_2$energy2\n";
  if (defined $strfun1) {
    print SAF "#strfun_1 $strfun1\n";
  } else {
    print SAF "#strfun_1 OFF\n";
  }
  if (defined $strfun2) {
    print SAF "#strfun_2 $strfun1\n";
  } else {
    print SAF "#strfun_2 OFF\n";
  }
  print SAF "#Final_state$finalstate\n";
  print SAF "#Remove_Virtual $exclde\n";
  print SAF "#Keep_Virtual $keep\n";
  print SAF "#nSubproc(ampl) 0\n";
  print SAF "#nSubproc(squared) 0\n";
  print SAF "#ConservationLaw 0\n";
  print SAF "#Nc==inf  0\n";
  print SAF "#ExitCode 0\n";
  close SAF;
  rename($sffl,"tmp/safe");
  system("cp","tmp/safe","safe");

  $diagfl="digrams";
  open (DIAG,"+>$diagfl") or die "I can not open safe_open file ($!)";
  print DIAG "diagrams:$excld\n";
  print DIAG "squared diagrams:$sqexcld\n";
  close DIAG;
  rename ($diagfl,"tmp/.excluded_digrams");
#  system("cp","digrams","tmp/.excluded_digrams");

}

###################################################################################################
#
sub run {
  if (defined $mp) {
# parallel calculations
    $argv = "}}}}}}}}}}}}]}9";
    if (defined $show and $show eq "diag") {
      $argv = "}}}}}}}}}}}}9";
      $err = system ("./comphep -blind $argv > /tmp/.la_la_la_crazyfrog");
      $err = system ("./comphep -diagshow");
      print "err = $err\n" if $err;
      $argv = "]}9";
    }
    $err = system ("./comphep -blind $argv > /tmp/.la_la_la_crazyfrog");

# define ranges of diagrams
    @info = get_stat();
    $tot = $info[1] + $info[3];
    $nses = int($tot/$mp) + 1;
    my @diag;
    $diag[0] = 1;
    for ($i = 1; $i < $mp; $i++) {
      $diag[$i] = $nses * $i;
    }
    $diag[$mp] = $tot + 2;

    my %kids;
    my $wpid;
    my $pid;
    keys(%kids) = 1;
    $comb=" ";

# main loop
    for ($i = 0; $i < $mp; $i++) {
      $dir="dir_tmp_".$i;
      if (-e $dir) {
        $tmp_dir = "$dir"."_old";
        if (-e $tmp_dir) {
          system "rm -rf $tmp_dir";
        }
        rename ($dir, $tmp_dir);
      }
      mkdir $dir,0777;
      chdir $dir;
      $comb = $comb.$dir." ";
      system "ln -s  ../comphep ./comphep";
      system "ln -s  ../.compheppath ./.compheppath";
      system "ln -s  ../models ./models";
      system "ln -s  ../process.dat ./process.dat";
      system "cp -r ../tmp ./tmp";

      $argv = "-diagfirst ".$diag[$i]." -diaglast ".$diag[$i+1]." -blind "."]}9";
      $shell_scrpt=qq~#!/bin/sh
$nice ./comphep $argv
~;
      open(FH,">run.sh");
      print FH $shell_scrpt;
      close(FH);
      chmod 0755, 'run.sh';
      $pid = fork();
      die "$$:cannot fork: $!" unless defined $pid;
      if ($pid == 0) {
        exec("./run.sh") or die "WARNING: Couldn't calculate Session $i+1 of process: $!\n";
        exit 0;
      }
      else {
        $num =$i+1;
        print "process $num launched, pid = $pid\n" if defined $verbose;
        $kids{$pid} = $num;
      }
      chdir "../";
    }

# Wait until all of the subprocesses to be calculated
    $wpid = wait();
    while ($wpid ne "-1") {
      if (WIFEXITED($?)) {
        if ($?) {
          print "WARNING: check the session $kids{$wpid}, status: $? \n";
        }
        delete $kids{$wpid};
        if (keys %kids > 0) {
          print "we are waiting for @{[ %kids ]} \n";
        }
      } else {
        print "WARNING: alarm for the PID $wpid Session $kids{$wpid}, Status: $? \n";
      }
      $wpid = wait();
    }
    print "All of the subprocesses are finished. Combining... \n";
    $err = system "./archiv $comb";

    if ($nmake =~ /yes/) {
      $argv="-blind ]}}]}9";
      $err = system ("$nice ./comphep $argv");
    }
  } else {
# usual calculations
    $argv="}}}}}}}}}}}}]}";
    if (defined $show and $show eq "diag") {
      $argv="}}}}}}}}}}}}9";
      $err = system ("./comphep -blind $argv > /tmp/.la_la_la_crazyfrog");
      $err = system ("./comphep -diagshow");
      print "err = $err\n" if $err;
      $argv="]}";
    }

    if ($smake =~ /yes/) {
      $argv=$argv."]}";
      if ($nmake =~ /yes/) {
        $argv=$argv."]}}]}";
      }
    }
    $argv = $argv."9";
    $err = system ("$nice ./comphep -blind $argv");
    print "err = $err\n" if $err;
    print "\n*** n_comphep creation details have been written to symb_batch.log (make stdout/stderr)\n"
  }
}

###################################################################################################
#
sub show_diagrams {
  $err=system ("./comphep -blind }}}}}}}}}}}9");
  print "err = $err\n" if $err;
}

###################################################################################################
#
sub recovery {
  open(SAFE,"<tmp/safe") or die "I can not open tmp/safe file";
  while($str=<SAFE>)
  {
    next unless $str =~ /#ExitCode\s+(\d)\s+/;
    if($1 != 6) {
      print "recovery is imposible!\n";
      exit (0);
    }
  }
  close(SAFE);

  opendir(RST,$rsltdir) or die "I can not open $rsltdir dir ($!)";;
  $empty = 0;
  while($f = readdir(RST)) {
    $empty = 1 unless ($f eq "Makefile" or $f eq "n_comphep" or $f eq "diag_view" or $f eq "." or $f eq "..");
  }

  if($empty) {
    $argv = "]}}";
  } else {
    $argv = " ]}}}";
  }
  if ($nmake =~ /yes/) {
    $argv=$argv."]}";
  }
  $argv=$argv."9";
  $err = system ("./comphep -blind $argv");
  print "err = $err\n" if $err;
  print "\n*** n_comphep creation details have been written to symb_batch.log (make stdout/stderr)\n"
}

###################################################################################################
#
sub relink {
  open(FOUT,">symb_batch.log") or die "I can not open symb_batch.log";
  open(OUT,"make -C $rsltdir link 2>&1 |") or die " Can not run program: $!\n";
  while($str=<OUT>) {
    print FOUT $str;
  }
  close(FOUT);

  print "*** Error during link of n_comphep! err = $err.\n" if $err;
  print "\n*** operation details have been written to symb_batch.log (make stdout/stderr)\n"
}

###################################################################################################
#
sub show_check_names {
  foreach (split(/,/,$show)) {
    $show_p{$_} = 1;
  }

  unless (defined $show_p{"diag"} or $show_p{"stat"}) {
    print "\nUnknown command in -show! \n";
    print "Use ./symb_batch.pl -show diag|stat\n\n";
    exit (-1);
    return;
  }
}

###################################################################################################
#
sub show_status {
  @info = get_stat();
  $tot = $info[1] + $info[3];
  print "Diagram statistics: total = $tot, calculated = $info[2], deleted = $info[3]\n";
}

###################################################################################################
#
sub get_stat {
  $err = system ("./comphep -status > .diag_info.log");
  open(STAT,"<.diag_info.log") or die "I can not open .diag_info.log";
  while ($line = <STAT>) {
    if ($line =~ /Status:\scalculated\s(\d+)\sdiagrams\sfrom\s(\d+)\s\(deleted\s(\d+)\)/) {
      $proc_n[1] = $2;
      $proc_n[2] = $1;
      $proc_n[3] = $3;
    }
  }
  close(STAT);
  system "rm .diag_info.log";
  return (@proc_n);
}

###################################################################################################
#
sub help {
 print
" The symb_batch.pl usage:\n
$0 [-h] [--help]
$0 [-d dir] [-n file] [-show diag] [-nice NN]
$0 [-d dir] -recovery
$0 [-d dir] -relink
$0 [-d dir] -mp N [-nice NN]
$0 -show stat
$0 [--version]

-h              - print this message.
-help           - print a long help message.

-d dir          - use directory \"dir\" instead of directory \"results\"
                  in order to to prepare n_comphep.
-n filename     - use data from file \"filename\" (instead process.dat).
-recovery       - recovery script work in case of s_comphep
                  crash during symbolic calculation of squared diagrams 
                  (does not work with -mp).
-relink         - relink n_comphep.exe. It is useful if a user changes userFun.c
-mp N           - parallel calculation of N jobs on one machine
-nice nn        - run comphep with modified scheduling priority (nice), 
                  nn - is priority, if it is undefined nn = 10;
-show
      diag      - launch s_comphep in graphic mode for feynman
                  diagrams review and analysis.
      stat      - show current status of the script work 
                  (in not implemented yet).\n";

if (defined $longhelp)
{
  print 
"*******************************************************************
  The main goal of the script is to launch the symbolic part of the
  CompHEP generator in non-GUI mode. This allows the reuse of the script,
  avoiding typing mistakes, and also should facilitate using farms.

  First you have to install CompHEP on your system (e.g. in a directory
  {your_home_dir}/comphep_4.2.0).

  A default data file for the script is process.dat. Detailed instructions
  and explanations for the changes are contained in the file itself. 
  After that you can launch the script with no options. 
  The final result is a n_comphep.exe program in the directory \"results\".

  You can set another names for the data file or the results directory
  by the -n and -d options respectively. If the directory results
  is not empty before the script launching, it is renamed to results_old_0
  and the script creates a new results directory.
  
  Some auxiliary options:
  -recovery If s_comphep crashed during calculations of squared diagrams, 
            you can launch the script with the recovery option and s_comphep 
            resumes the computation from the last calculated diagram.

  -relink   If you have changed userFun.c file one should relink the n_comphep
            program. You can launch the script with the
            relink option and s_comphep relink the n_comphep program and
            save all details of the relinking to the symb_batch.log file.

  -mp N     If your machine has more then one CPU, it is possible to launch 
            several comphep jobs. This option prepares and runs N jobs in parallel. 
            Each jobs is bing run in directory \"dir_tmp_I\", where I=0,...,N-1.
            As soon as all jobs are over, archiv combines all data in ./tmp and 
            prepared n_comphep.

  -nice NN  This option says that comphep will be running with non-zero priority 
            in order to lower CPU demand for the task (nice is the standard tool 
            for that in Linux). Permitten range for nice is 1-19, the default 
            the status is equal to 10. 

  -show diag     - this option is applied if you'd like to exclude
                   by hand some diagrams from all those set by 
                   the program. The script launches s_comphep in GUI
                   mode in the feynman diagrams menu. After reviewing
                   the diagrams you have to finish the GUI session and
                   the script will go on.
                   
  -show stat     - show the current working status of the script. 
                 (i.e. how many diagrams have been calculated already)\n";}
 }

###################################################################################################
sub version {
  print "$0 (CompHEP symbolical tool) $version, date: $versiondate\n\n";
  print "Copyright (C) 2002-2008, CompHEP Collaboration.\n";
  print "This software is distributed under the CompHEP licence.\n";
  print "All details of the licence see at http://comphep.sinp.msu.ru\n";
  exit(0);
}
