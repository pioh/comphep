#!/usr/bin/perl -w
# Copyright (C) 2002-2009, CompHEP Collaboration
# Authors: Lev Dudko, Alexander Sherstnev
#
#####################################################################################
# CHANGELOG
# v. 0.1             1) created by Lev Dudko
# v. 0.2  19.02.2002
# v. 0.3  06.03.2002
# v. 0.5  23.06.2003 1) PBS/LSF parallel runs
# v. 1.0  08.11.2005 1) first production version
# v. 1.1  19.09.2006 1) report functionality added.
# v. 1.2  21.09.2006 1) a bug with the 'max' option fixed
# v. 1.3  06.10.2006 1) a bug with the 'max' option really fixed
# v. 1.4  01.12.2006 1) a bug of different cs errors in -lcs and -show cs fixed
#                    2) check wrong command in -run
#                    3) option -safe added (to prevent comphep launching if there is LOCK)
# v. 1.5  14.02.2007 1) a routine to print the script version added
#                    2) this version is included to comphep-4.4.4
#                    3) a bit more detailed warning if no dir results
#                    4) a bug in printing of CS error in -show cs fixed
#                    5) -b option descibed in help
# v  1.6  15.03.2007 1) add possbility for option -proc with option --add 
#                    2) add -lmix option
#                    3) remove automatic save of exist session.dat to batch.dat
# v. 1.7  06.06.2007 1) abolished...
# v  1.8  19.06.2007 1) add possbility generate events in the LHAef format (-lha)
#                    2) a small bug in -show (new option: model)
#                    3) -proc and -add work together properly...
# v  1.9  19.10.2007 1) add possbility of parallel run of N subprocesses (-mp N)
#                    2) fixed bug in -add -proc options together
#                    3) add option -verbose
# v  1.10 03.11.2007 1) now -show and -proc can work together
#                    2) long description revised
# v  1.11 20.11.2007 1) the script is compatible with comphep-4.4.96 or later
#                    2) rename -lha -> -lhaef
# v  1.12 18.04.2008 1) nice option added 
#                    2) improved warning/error messages
# v  1.13 17.11.2008 1) negative cross sections can processed on -nevnt
# v  1.14 22.03.2009 1) proper functioning of -kfactor with LHE files
# v  1.15 03.04.2009 1) copying several sections from session.dat to batch.dat
#                    2) extract session.dat from batch.dat for a concrete subproc
#                    3) merging prt files if old ones exist and -mp is used
# v  1.16 06.04.2009 1) internal cosmetic improvements
#                    2) new default format for events (LHE), and -cpyth1/-cpyth2 for olf formats
#                    3) TODO list
#                    4) standardized error/warning/info messages

#####################################################################################
# TODO
# 1) lmix for all formats
# 2) time reports and forecast
# 3) -auto option 
# 4) divede help sections, a-la ./num_batch.pl help run 
# 5) more reliable method to deal with random seeds, especially for -auto
# 6) for -auto: cleanstat if chi^2 > 2 only

use Getopt::Long;
use File::Path;
use File::Copy;
use POSIX qw(:sys_wait_h :signal_h :errno_h);

help() unless (@ARGV or !-e "results/batch.dat");

#####################################################################################
# Take the input options
GetOptions (
           "h",         \$opt_h,                # short help
           "help",      \$longhelp,             # long help
           "lcs",       \$lcs,                  # list cross sections
           "d=s",       \$opt_d,                # set working directory name (insted of results)
           "b=s",       \$opt_b,                # set batch file name (instead of batch.dat)
           "run:s",     \$run,                  # run production (main option)
           "auto=s",    \$autoreg,              # calculations are regulated by a given error
           "show:s",    \$show,                 # show different parts of batch.dat
           "proc=s",    \$proc,                 # set a partucular subprocess set
           "add:s" =>   \$ses2bat,              # add current session.dat to batch.dat
           "nses=i",    \$insession,            # set the initial number for prt
           "ginv",      \$ginv,                 # set gauge invarian width option
           "safe",      \$safeopt,              # launch n_comphep in same mode
           "version",   \$veropt,               # print out version of the script
           "cpyth1",    \$cpyth1,               # generate event in the cpyth1 format
           "cpyth2",    \$cpyth2,               # generate event in the cpyth2 format
           "geninfo",   \$geninfo,              # copy sev. sections from session.dat to batch.dat
           "extract=i", \$extract,              # extract session.dat from batch.dat for subproc i
           "nevnt=s",   \$nevnt,                # set the full number of events
           "call=s",    \$ncall,                # set the number of calls
           "kfactor=s", \$kfactor,              # scale cross section in files from "files"
           "files",     \$files,                # set all file names for the previous option
           "pbs:s",     \$pbs,                  # batch system options
           "lsf:s",     \$lsf,                  # batch system options
           "combine",   \$opt_combine,          # batch system options
           "nocombine", \$opt_nocombine,        # batch system options
           "lmix",      \$lmix,                 # list mixed subprocesses from mixed event file
           "mp=i",      \$mp,                   # set the number of parallel subprocesses
           "nice:s",    \$nnice,                # set priority for nice regime
           "verbose",   \$verbose,              # print more information 
           "silent",    \$silent                # print more information 
           );

###################################################################################################
# Set global variables - line numbers of blocks in session.dat
  my $number_subproc_num =  1;
  my $number_ses_num     =  2;
  my $number_model_num   =  3;
  my $number_inistate    =  4;
  my $number_params      =  5;
  my $number_width_sch   =  6;
  my $number_kin_sch     =  7;
  my $number_cuts        =  8;
  my $number_reguls      =  9;
  my $number_qcd_info    = 10;
  my $number_veg_call    = 11;
  my $number_veg_intl    = 12;
  my $number_disto       = 13;
  my $number_evnt        = 14;
  my $number_random      = 15;
  my $number_veg_grid    = 16;
  my $number_max         = 17;
  my $number_end         = 18;

###################################################################################################
# Set global variables
  my $fileses='session.dat';
  my $filebat='batch.dat';
  my $Nproc = 0;
  my $Nproc_current = 0;
  my $Nsession = 1;
  my $version = 1.16;
  my $datestamp = "06/04/2009";
  my $nice = "";
# variables for warning, error and info messages
  my $errr = "$0 (error):";
  my $warn = "$0 (warning):";
  my $info = "$0 (info):";

###################################################################################################
###################################################################################################
###################################################################################################
# print the script version
  version () if (defined $veropt);

###################################################################################################
# Print help message
  help () if ($opt_h or $longhelp);

###################################################################################################
# Change cross section values in event files according to kfactor.
# The file names are taken from the hash @ARGV
  change_cs (@ARGV) if (defined $kfactor && defined $files);

###################################################################################################
# List cross section values from protocol files (in the "prt" format).
# The file names are taken from the hash @ARGV
  list_prt (@ARGV) if defined $lcs;

###################################################################################################
# List cross section values from protocol files (in the "prt" format).
# The file names are taken from the hash @ARGV
  list_mix (@ARGV) if defined $lmix;

###################################################################################################
# Use "results" or another defined (-d dir) directory for the
# numerical calculations. The n_comphep file must be there.
  $opt_d = "results" unless defined $opt_d;
  chdir $opt_d or die "$errr There are no $opt_d dir: $!\nTry $0 --help\n";

###################################################################################################
# change default name of the batch file
  $filebat = $opt_b if defined $opt_b;

###################################################################################################
# Print the list of available processes in the
# n_comphep generator in the directory [dir] (or results).
  list_batch () if defined $show;

###################################################################################################
# Set nice priority for n_comphep
  set_nice () if defined $nnice;

###################################################################################################
# Set event numbers in each subprocesses according their contribution
# to the full cross section. The only input parameter is the full
# number of events requested
replace_nevnt (@ARGV) if defined $nevnt;

###################################################################################################
replace_calls (@ARGV) if defined $ncall;

###################################################################################################
replace_pars (@ARGV) if defined $geninfo;

###################################################################################################
extract_session (@ARGV) if defined $extract;

###################################################################################################
# remove LOCK file if it exists
  if (-e "LOCK") {
    print "\nThere is LOCK file in $opt_d\n";
    exit (0) if defined $safeopt;
    print "\n$warn LOCK file exists. It means n_comphep has not terminated properly or \n";
    print "                          it is running. the file is removed and n_comphep is launched\n";
    print "                          If you want to prevent the possibly unsafe behaviour use option -safe\n";
    unlink "LOCK";
  }

###################################################################################################
# check whether numerical comphep exists
  if (!-e "./n_comphep.exe") {
    print "$errr there is no n_comphep.exe. At first you must prepare\n";
    print "                        the n_comphep generator in directory $opt_d.\n";
    print "                        Use comphep symbolic part to do that.\n\n";
    exit (0);
  }

###################################################################################################
# Create batch.dat file in [dir] for the all available processes
# if it is absent.
  create_batch_dat () if !-s $filebat;

###################################################################################################
# With the option --add [file] include the parameters from session.dat (or [file])
# to the existed batch.dat file.
  if (defined $ses2bat) {
    $ses2bat=$fileses if $ses2bat eq "";
    old_ses ($ses2bat);
  }

###################################################################################################
# Combine the results of previous calculations from the temporary PBS directories (this option
# is needed only if nocombine option was used before)
  combine() if defined $opt_combine;

###################################################################################################
# two possible ways to do calculations
###################################################################################################
# Calculation are driven by a given error. -auto X% means calculations should stop as soon as 
# the total stat. error amounts to X%. If events are needed some more subtle condition is used
  if (defined $autoreg) {
    auto_calc () ;
    exit (0);
  }

###################################################################################################
# Calculation are driven by parameters set in batch.dat.
# Launch n_comphep in the standard butch regime
  run () if defined $run;

###################################################################################################
# end of the main program
  exit (0);

###################################################################################################
###################################################################################################
###################################################################################################
# working routines
sub old_ses {
  local $/="\n#";
  open(SES,"<$_[0]") or die "there are no $_[0] in $opt_d ($!)";
  @old_session = <SES>;

# Take the session number from the session.dat file
  if ($_[0] eq $fileses) {
    $old_session[2] =~/\s(\d+)\s/g;
    $Nsession=$1;
  }
  $old_session[1] =~/\s(\d+)\s/g;
  $Nproc_session = $1;
  close SES;

# Add parameters from session.dat (or --add [file]) file to batch.dat file
  open(BATCH, "+<$filebat") or die "Can't open $filebat for read/write: $!\n";
  local $/="++++++++++++++\n";
  my $tmp = "";
  while (<BATCH>) {
    @blocks = split(/\n#/);
    $blocks[1] =~/\s(\d+)\s/g;
    $Nproc_current = $1;

    if (defined $proc) {
       @proc_n = proc_nc();
    } else {
      @proc_n = ($Nproc_session);
    }

    if (grep { $Nproc_current eq $_ } @proc_n ) {
      $old_session[1] = $blocks[1]."\n#";
      $old_session[2] = $blocks[2]."\n#";
      $tmp .= join("",@old_session)."++++++++++++++\n";
      print "          The exist $ses2bat file is included to
      $filebat file for the process $Nproc_current.\n";
    } else {
      $tmp .= $_
    }
  }
  seek(BATCH,0,0);
  print BATCH $tmp;
  truncate (BATCH,tell(BATCH));
  close BATCH;
}

###################################################################################################
# set nice regime for n_comphep.exe
sub set_nice {
  $nice = "nice -n 10 ";
  if ($nnice ne "") {
    if ($nnice =~ m/d+/) {
      if ($nnice gt "19") {
        print "$warn $nnice is too high, 19 is used ";
        $nnice = 19;
      }
      $nice = "nice -n $nnice ";
    } else {
      print "$warn -nice argument should be a number, 10 is used";
    }
  }
}

###################################################################################################
# replace parameters, qcd info, init state, and width scheme.
# (can be extended to other session.dat sections)
sub replace_pars {
  open (SES, "<$fileses") or die "there are no $fileses in $opt_d ($!)";
  my $totses = "";
  while (<SES>) {
    $totses = $totses.$_;
  }
  close SES;
  @ses_blocks = split(/\n#/, $totses);

# replace sections #Initial_state, #Physical_Parameters, #Width_scheme, and #QCD
  open(BATCH, "+<$filebat") or die "$errr can't open $filebat in $opt_d: $!\n";
  local $/="++++++++++++++\n";
  my $tmp = "";
  while (<BATCH>) {
    @blocks = split(/\n#/);
    $blocks[$number_params]    = $ses_blocks[$number_params];
    $blocks[$number_inistate]  = $ses_blocks[$number_inistate];
    $blocks[$number_width_sch] = $ses_blocks[$number_width_sch];
    $blocks[$number_qcd_info]  = $ses_blocks[$number_qcd_info];
    $blocks[$number_veg_call]  = $ses_blocks[$number_veg_call];
    for ($i = 1; $i <= $#blocks; $i++) {
      $blocks[$i] = "\n#".$blocks[$i];
    }
    $tmp .= join ("",@blocks);
  }
  seek(BATCH,0,0);
  print BATCH $tmp;
  truncate (BATCH, tell(BATCH));
  close BATCH;
}

###################################################################################################
# extract session.dat from batch.dat for the i-th subprocess
sub extract_session {
  @sessions = get_sessions ();
  open (SES, ">$fileses") or die "$errr can't open $fileses in $opt_d: $!\n";
  $ncurses = $sessions[$extract-1];
  $ncurses =~s/\+\+\+\+\+\+\+\+\+\+\+\+\+\+\n//;
  print SES $ncurses;
  close SES;
}

###################################################################################################
sub create_batch_dat {
  open(BATCH, "+>$filebat") or die "$errr can't open $filebat in $opt_d: $!\n";
  system ("./n_comphep -blind }}9}");
  die "Can't run n_comphep:$@" if $@;
  my $string="}";
  print "You can calculate the processes as follow:\n";
  $/="\n#";
  for (;;) {
    open(SES,"<$fileses") or die "there are no session.dat in $opt_d ($!)";
    @block = <SES>;
    $block[1] =~/\s(\d+)\s/g;
    $Nproc_current = $1;
    if ($Nproc_current > $Nproc+1) {
      $string=$string.")";
      goto COMP;
    }
    last if $Nproc_current == $Nproc;

    $Nproc=$Nproc_current;
    if ($old_session[1] and $old_session[1] eq $block[1]) {
      print BATCH @old_session;
      print "
   The exist session.dat file parameters are included to new
   batch.dat file for the process $Nproc_current.\n"
    } else {
      print BATCH @block
    }
    print BATCH "++++++++++++++\n";
    $block[1]=~ s/(\n\#)/\n/;
    print $block[1];
    close(SES);
    $string=$string."]";
COMP: system ("./n_comphep", "-blind", "$string.'}9}'");
    die "Can't run n_comphep:$@" if $@;
  }
  print "File ${opt_d}/batch.dat is created.\n";
  close BATCH;
}

###################################################################################################
# the main script function. Corresponds to -run and defines different targets of the option 
sub run {
# Fill the steps to run from input options and check the commands
  my %run_p;
  keys(%run_p) = 0;
  if ($run eq "") {
    $run_p{"vegas"} = 1;
    $run_p{"max"}   = 1;
    $run_p{"evnt"}  = 1;
  }
  foreach (split(/,/,$run)) {
    $run_p{$_} = 1;
    if (($_ ne "cleangrid") and ($_ ne "cleanstat") and ($_ ne "clean") and 
        ($_ ne "vegas") and ($_ ne "max") and ($_ ne "evnt")) {
      print "\n$errr unknown command($_) in -run! \n";
      print "        possible commands: vegas,max,evnt,cleangrid,cleanstat,clean or their combinations divided by commas\n";
      return;
    }
  }

  my $nginv="";
  $nginv="}}{" if defined $ginv;
  my $ndifffrmt="]]";
  $ndifffrmt="]" if defined $cpyth1;
  $ndifffrmt="" if defined $cpyth2;

  local $/="++++++++++++++\n";
  open(BATCH, "+<$filebat") or die "$errr can't open $filebat for read: $!\n";
  my @batch = <BATCH>;
  close BATCH;

  my $nmp = defined $mp ? $mp : $#batch+1;
  my %kids;
  my $wpid;
  my $pid;
  keys(%kids) = 1;

# Set number of subprocesses to run from input options
  if (defined $proc) {
    @proc_n = proc_nc ();
    print "$warn There are only ", $#batch+1," processes, not $proc_n[$#proc_n].\n" unless $proc_n[$#proc_n] <= $#batch+1;
  } else {
    @proc_n = 1 .. $#batch + 1;
  }

  while (-e "prt_".$Nsession and ! defined $insession) {
    $Nsession++;
  }
  $Nsession = $insession if defined $insession;

# Run numerical CompHEP generator, loop by subprocess 
  foreach $prc (@proc_n) {
    $Nsession = $prc if defined $autoreg;
    $Nprocess = $prc - 1;
    prepare_session ($Nsession, $Nprocess);

    my $max_string = (length( $blocks[17]) > 50 ) ? "]]]]]$nginv]]]]}]]]]]]]]}]]]]}]]]}9}" : 
                                                    "]]]]]$nginv]]]]}]]]]]]]]}]]]}9}";
    my $vegas_string = "]]]]]$nginv]]]]}]]]}9}"                  ;
    my $event_string = "]]]]]$nginv]]]]}]]]]]]]]}]]$ndifffrmt}9}";
    my $cleangrid_string= "]]]]]]]]]}]]]]]]]}}9}"                ;
    my $cleanstat_string= "]]]]]]]]]}]]]]]]}}9}"                 ;
    my $clean_string= "]]]]]]]]]}]]]]]]}}]}}]}{9}"               ;
    if (defined $silent) {
      $max_string       .= " 1>stdout 2> stderr";
      $vegas_string     .= " 1>stdout 2> stderr";
      $event_string     .= " 1>stdout 2> stderr";
      $cleangrid_string .= " 1>stdout 2> stderr";
      $cleanstat_string .= " 1>stdout 2> stderr";
      $clean_string     .= " 1>stdout 2> stderr";
    }

# Parallel or lsf/pbs calculation mode
# Create subdirectories for each subprocess, setup there the session.dat
# and run run.sh in PBS or LSF batch system
    if (defined $pbs or defined $lsf or defined $mp) {
      my ($vegas,$max,$event) = ('','','');
      my $comphep='';
      my $batch_options='';
      my $b_prefix='';
      $vegas = '#' unless defined $run_p{"vegas"};
      $max   = '#' unless defined $run_p{"max"};
      $event = '#' unless defined $run_p{"evnt"};
      $comphep="COMPHEP=".$ENV{"COMPHEP"}." \n export COMPHEP\n" if defined $ENV{"COMPHEP"};
      if (defined $pbs) {
        $batch_options = qq~
#PBS -I
##PBS -l nodes=1:farm
~;
        $pbs eq '' ? $b_prefix = "qsub" : $b_prefix = $pbs;
      } elsif (defined $lsf) {
        $batch_options = qq~
#LSF -I
~;
        $lsf eq '' ? $b_prefix = "bsub" : $b_prefix = $lsf;
      }
      $dirname = 'res_'.$Nsession;
      mkdir ($dirname,0777);
      chdir $dirname;
      system ("ln -s  ../LHAIndex-comphep.txt ./LHAIndex-comphep.txt");
      system ("ln -s  ../n_comphep ./n_comphep");
      system ("ln -s  ../n_comphep.exe ./n_comphep.exe") if -s "../n_comphep.exe";
      system ("mv  ../session.dat .");
      $shell_scrpt=qq~#!/bin/sh
$batch_options
$comphep
${vegas}$nice./n_comphep -blind ${vegas_string}
${max}$nice./n_comphep -blind ${max_string}
${event}$nice./n_comphep -blind ${event_string}
~;
      open (FH,">run.sh");
      print FH $shell_scrpt;
      close (FH);
      chmod 0755, 'run.sh';
      $pid = fork();
      die "$$:cannot fork: $!" unless defined $pid;
      if ($pid == 0) {
        exec("$b_prefix ./run.sh") or die "WARNING: Couldn't calculate Session $Nsession of process $Nprocess: $!\n";
        print "CHILD finished, Session: $Nsession with process: $Nprocess \n" if defined $verbose;
        exit 0;
      } else {
        $kids{$pid}=$Nsession;
      }
      chdir "../";
      $Nsession++;
      foreach $ccpid (keys %kids) { 
        if (kill(0, $ccpid)) {
          print "alive pid, session: $ccpid $kids{$ccpid} \n" if defined $verbose;
        } else {
          print " MISSED PID removed (pid, session): $ccpid $kids{$ccpid} \n ";
          delete $kids{$ccpid};
        }
      }
      print scalar(keys %kids)," Active PID, Session: @{[ %kids ]} \n" if defined $verbose;

      while (keys %kids>$nmp-1) {
        $wpid=wait();
        if (WIFEXITED($?)) {
          print "PID $wpid Session $kids{$wpid} finished, Status: $? \n" if defined $verbose;
          if($?){print "WARNING: check the session $kids{$wpid}, status: $? \n";}
          delete $kids{$wpid};
        }
      }
    } else {
# usual batch calculation mode
      if (defined $run_p{"vegas"}) {
        $status=system ("$nice./n_comphep -blind $vegas_string");
        unless (defined $silent) {
          open (PRT,"<prt_$Nsession") or die "$errr Can't open prt_$Nsession in $opt_d ($!)";
          $/="\n";
          while (<PRT>) {
            print if /\s\#IT|\s\<\s\>/
          };
          print "\n";
          close PRT;
        }
      }
      $status = system ("$nice./n_comphep -blind $max_string")       if defined $run_p{"max"};
      $status = system ("$nice./n_comphep -blind $event_string")     if defined $run_p{"evnt"};
      $status = system ("$nice./n_comphep -blind $cleangrid_string") if defined $run_p{"cleangrid"};
      $status = system ("$nice./n_comphep -blind $cleanstat_string") if defined $run_p{"cleanstat"};
      $status = system ("$nice./n_comphep -blind $clean_string")     if defined $run_p{"clean"};
      if (0 == $status) {
        open(SES,"<$fileses") or die "Can't read session.dat in $opt_d ($!)";
        local undef $/;
        $batch[$Nprocess] = <SES>;
        $batch[$Nprocess] .= "++++++++++++++\n";
        close SES;
        open(BATCH, "+>$filebat") or die "Can't open $filebat for read/write: $!\n";
        print BATCH @batch;
        close BATCH;
        $Nsession++;
      } else {
        print "$errr Something wrong in calculations: $status \n" unless defined $status;
        exit (-1);
      }
    }
  }

# Wait until all of the subprocesses to be calculated
  $wpid = wait ();
  while ( $wpid ne "-1" ) {
    if (WIFEXITED($?)) {
      print "PID $wpid Session $kids{$wpid} finished, Status: $? \n" if defined $verbose;
      if($?) {print "WARNING: check the session $kids{$wpid}, status: $? \n";}
      delete $kids{$wpid};
      if (keys %kids > 0 ) {print "we are waiting for @{[ %kids ]} \n";}
    } else {
      print "WARNING: alarm for the PID $wpid Session $kids{$wpid}, Status: $? \n";
    }
    $wpid=wait();
  }
  print "All subprocesses calculated\n" unless defined $silent;

# clean cubic maxima
  if (defined $run_p{"clean"} or defined $run_p{"cleanstat"} or defined $run_p{"cleangrid"}) {
    replace_max ();
  }

# include the new session.dat parameters to batch.dat, move prt_ and events_
# files to main results dir and remove temporary res_ directories
  combine () if !defined $opt_nocombine;
}

###################################################################################################
sub prepare_session {
  my $Nproc = $_[1];
  local $/="++++++++++++++\n";
  open(BATCH, "+<$filebat") or die "Can't open $filebat for read: $!\n";
  my @batch = <BATCH>;
  close BATCH;

  print "The current session is number $_[0]\n" unless defined $silent;
  @blocks = split(/\n#/,$batch[$Nproc]);
  $blocks[2] =~s/\d+/$_[0]/g;
  print "$blocks[1] \n" unless defined $silent;
  foreach (@blocks) {
    $_="\n#".$_ unless $_ eq "\n" or $_ eq "";
  }
  open (SES,"+>$fileses") or die "$errr Can't create session.dat in $opt_d: $!";
  print SES @blocks;
  close SES;
}

###################################################################################################
sub auto_calc {
# check and move old prt/batch/session to old_results_X
  my @oldprts = <prt_*>;
  if (0 < $#oldprts) {
    $num=1;
    $rsltdir="old_results"."_$num";
    while (-e $rsltdir) {
      $num++;
      $rsltdir="old_results"."_$num";
    }
    mkdir ($rsltdir, 0755) or die "$errr can't create  $rsltdir directory: $!\n";
    print "$warn there are prt files in $opt_d. These files, session.dat, and batch.dat are saved in $rsltdir\n";
    move ("session.dat",  "$rsltdir/session.dat");
    copy ("batch.dat",    "$rsltdir/batch.dat");
    system ("mv prt_* $rsltdir/.");
  }

# main loop
  my $cleanstat=1;
  if ($autoreg < 0.000001) {
    print "$warn Requested error is too small. It may require too much time to reach the accuracy.\n";
  }
  $silent=1;
  print "Start culculations!\n";
  if (defined $nevnt) {
    one_cycle ();
  } else {
    one_cycle ($cleanstat);
    $nrun = 1;
    $calcerr=mid_cs_report ();
    while (0. < $calcerr - $autoreg and $nrun < 10) {
       print "Run $nrun: accomulated err = $calcerr, requested error = $autoreg. Runs continue\n\n";
      replace_ncall_niter ();
      one_cycle ();
      $calcerr=mid_cs_report ();
      $nrun++;
    }
    print "Run $nrun: accomulated err = $calcerr, requested error = $autoreg. Runs over\n\n";
  }
}

###################################################################################################
sub one_cycle {
  $run = "vegas";
  run ();
  print "Vegas estimation is done\n";

  if (defined $_[0]) {
    $clr_subpr = mid_cs_report ("ksi");
    if (defined $proc) {
      $proc_bkup = $proc;
    } else {
      $remove_proc=1;
    }
    $proc = $clr_subpr;
    $run = "cleanstat";
    run ();
    print "Clean up statistics for subprocesses: $proc\n";
    if (defined $remove_proc) {
      undef $proc;
    } else {
      $proc=$proc_bkup;
    }
    $run = "vegas";
    run ();
    print "Vegas calculations is done\n";
  }
}

###################################################################################################
sub vegas_calc_duration {
  my $Nprocess = $_[0];
  print "$warn The function still not implemented...\n";
}

###################################################################################################
sub proc_nc {
  $proc =~ s/(\d+)-(\d+)/join(',',$1 .. $2)/eg;
  my @list = sort {$a <=> $b} split(/,/,$proc);
  my %temp= ();
  @proc_n = grep { !$temp{$_}++ } @list;
  return (@proc_n);
}

###################################################################################################
# combine results obtained in several res_N: copy corresponding parts to batch.dat 
# and move events_N.txt and prt_N to results
sub combine {
  my @list = <res_*>;
  my @oldprts = <prt_*>;
  %prthere = ();
  foreach $prt (@oldprts) {$prthere{$prt} = 1;}

  foreach $resdir (@list) {
    $resdir=~/res\_(\d+)/g;
    my $resn=$1;
    my $prt="prt_$1";
    my @keepfiles;
    $ses2bat="$resdir/$fileses";
    old_ses ("$resdir/$fileses") if -e "$resdir/$fileses" ;
    if (-e "$resdir/$prt") {
      $keepfiles[0] = "$resdir/$prt";
      if (-e "$prt") {
        mergeprt ("$prt", "$resdir/$prt");
      } else {
        copy ("$resdir/$prt","$prt");
      }
    }

    if (-e "$resdir/events_$resn.txt") {
      $evnt="events_$resn.txt";
      $keepfiles[1]="$resdir/$evnt";
      move ("$resdir/$evnt","$evnt");
    }
    eval {rmtree("$resdir")};
    if ($@) {
      print "$0 (error): can't delete directory $resdir, reason:$@";
    }
    else {
      print "$0 (info): the $resdir directory is removed, @keepfiles files are kept\n";
    }
  }
}

####################################################################################################
# add new section appeared in (new) prt file $_[1] to (old) prt file $_[0]
# if subprocess names in the file are different, mv $_[0] _[0].old and mv $_[1] _[0]
sub mergeprt {
  local $/="\n#";
  open (OLDPRT,"<$_[0]") or die "$0 (error):there are no $_[0] in $opt_d ($!)";
   $totses = "";
  while (<OLDPRT>) {$totses = $totses.$_;}
  @old_prt = split (/\n#/, $totses);
  open (NEWPRT,"<$_[1]") or die "$0 (error):there are no $_[1] in $opt_d ($!)";
  $totses = "";
  while (<NEWPRT>) {$totses = $totses.$_;}
  @new_prt = split (/\n#/, $totses);
  close (OLDPRT);
  close (NEWPRT);

  if ($old_prt[1] ne $new_prt[1]) {
    print "$warn $_[1] and $_[0] contains information about different subprocesses.\n";
    print "              So, the old prt file is kept\n";
    move ("$_[0]","$_[0].old");
    copy ("$_[1]","$_[0]");
    return;
  }

  %blkthere = ();
  foreach $blk (@old_prt) {$blkthere{$blk} = 1;}

  open (OLDPRT,">$_[0]") or die "$0 (error): there are no $_[0] in $opt_d ($!)";
  print OLDPRT $old_prt[0];
  for ($i = 1; $i <= $#old_prt; $i++) {
    print OLDPRT "\n#".$old_prt[$i];
  }
  my $strt = 0;
  for ($i = 1; $i <= $#new_prt; $i++) {
    if ($new_prt[$i] =~ m/END\s+/) {$strt = $i;}
  }
  for ($i = $strt; $i <= $#new_prt; $i++) {
    unless (defined $blkthere{$new_prt[$i]}) {print OLDPRT "\n#".$new_prt[$i];}
  }
  close (OLDPRT);
  close (NEWPRT);
}

###################################################################################################
sub list_batch {
  foreach (split(/,/,$show)) {
    $show_p{$_} = 1;
  }
  unless (defined $show_p{"proc"} or $show_p{"model"} or $show_p{"inistat"} or $show_p{"params"}  or
                  $show_p{"cuts"}  or $show_p{"regs"}    or $show_p{"qcd"}     or
                  $show_p{"vegas"} or $show_p{"events"}  or $show_p{"cs"}      ) {
    print "\nUnknown command in -show! \n";
    print "Use ./num_batch.pl -show proc|inistat|params|cuts|regs|qcd|vegas|events|cs\n\n";
    return;
  } else {
    if (defined $proc) {
       @proc_n = proc_nc ();
    }
    my $cs_tot     = 0.0;
    my $cs_err_tot = 0.0;
    print "List of available subprocesses:\n";

    @sessions = get_sessions ();
    @proc_n = 1..$#sessions+1;
    @proc_n = proc_nc () if defined $proc;
    foreach (@sessions) {
      @blocks = split(/\n#/,$_) ;
      $blocks[1] =~/\s(\d+)\s/g;
      $Nproc_current = $1;
      if (grep { $Nproc_current eq $_ } @proc_n) {
        print "-----------------------------------------------------------\n" if defined $show_p{"inistat"} or $show_p{"params"} or
                                                                                         $show_p{"cuts"}    or $show_p{"regs"};
        print "$blocks[1] \n" if defined $show_p{"proc"};
        print "$blocks[1]:\n$blocks[3]  \n" if defined $show_p{"model"};
        print "$blocks[1]:\n$blocks[4]  \n" if defined $show_p{"inistat"};
        print "$blocks[1]:\n$blocks[5]  \n" if defined $show_p{"params"};
        print "$blocks[1]:\n$blocks[8]  \n" if defined $show_p{"cuts"};
        print "$blocks[1]:\n$blocks[9]  \n" if defined $show_p{"regs"};
        print "$blocks[1]:  $blocks[10] \n" if defined $show_p{"qcd"};
        print "$blocks[1]:  $blocks[11] \n" if defined $show_p{"vegas"};
        print "$blocks[1]:  $blocks[14] \n" if defined $show_p{"events"};
        if ($show_p{"cs"}) {
          if ($blocks[12] =~ /Vegas_integral\s+([-]*\d\.\d+E[+-]\d+)\s+([-]*\d\.\d+E[+-]\d+)\s+([-]*\d\.\d+E[+-]\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
            $cs = 0.0;
            $cs_err = 0.0;
            $ksi_test = 0;
            if (0.0 != $1) {
              $cs = $2/$1;
              if (1 < $4) {$ksi_test = ($3 - $2 * $2 / $1) / ($4 - 1);}
              if ($1 > 0.0) {
                $cs_err = sqrt($1)/abs($2);
              } else {
                $cs_err = sqrt(-$1)/abs($2);
              }
            }
            $cs_tot     = $cs_tot     + $cs;
            $cs_err_tot = $cs_err_tot + $cs_err * $cs_err * $cs * $cs;
            printf "%s: cross section [pb] = %1.4e   +\/- %1.2e ( %1.2e %% ), chi2 = %1.2e\n", 
            $blocks[1], $cs, $cs_err*abs($cs), 100.*$cs_err,$ksi_test;
          } else {
            printf "$errr can't decode the Vegas_integral string for %s\n", $blocks[1];
          }
        }
      }
    }
    if ($show_p{"cs"}) {
      my $cs_err_tot_p = sqrt($cs_err_tot) * 100/$cs_tot;
      printf "\n Total CS [pb] = %1.4e   +\/- %1.2e ( %1.2e %% )\n",$cs_tot,sqrt($cs_err_tot),$cs_err_tot_p;
    }
  }
  exit (0);
}

###################################################################################################
sub get_sessions {
  local $/="++++++++++++++\n";
  open (BATCH, "+<$filebat") or die "$errr can't open $filebat: $!\n";
  @ses = <BATCH>;
  close (BATCH);
  return (@ses);
}

###################################################################################################
sub mid_cs_report {
  my $ksi_gt_1 = "";
  @sessions = get_sessions ();
  @proc_n = 1..$#sessions+1;
  @proc_n = proc_nc () if defined $proc;
  my $cs_tot = 0.;
  my $cs_err_tot = 0.;
  foreach (@sessions) {
    @blocks=split(/\n#/,$_);
    $blocks[$number_subproc_num] =~/\s(\d+)\s/g;
    $Nproc_current = $1;
    if (grep { $Nproc_current eq $_ } @proc_n) {
      if ($blocks[$number_veg_intl] =~ /Vegas_integral\s+([-]*\d\.\d+E[+-]\d+)\s+([-]*\d\.\d+E[+-]\d+)\s+([-]*\d\.\d+E[+-]\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
        if (0.0 != $1) {
          $cs = $2/$1;
          if ($1 > 0.0) {
            $cs_err = sqrt($1)/abs($2);
          } else {
            $cs_err = sqrt(-$1)/abs($2);
          }
        } else {
          $cs = 0.0;
          $cs_err = 0.0;
        }
        $cs_tot     = $cs_tot     + $cs;
        $cs_err_tot = $cs_err_tot + $cs_err * $cs_err * $cs * $cs;
        $ksi_test   = ($3 - $2 * $2 / $1) / ($4 - 1);
        if ($ksi_test > 2.) {
	  $ksi_gt_1 .= "$Nproc_current,";
	}
	printf "%s: cross section [pb] = %1.4e   +\/- %1.2e ( %1.2e %% )\n", $blocks[$number_subproc_num], $cs, $cs_err * abs($cs), 100. * $cs_err if defined $verbose;
      }
    }
  }

  my $cs_err_tot_p = sqrt ($cs_err_tot)/$cs_tot;
  printf "Total CS [pb] = %1.4e   +\/- %1.2e ( %1.2e %% )\n",$cs_tot,sqrt($cs_err_tot),100*$cs_err_tot_p if defined $verbose;
  if ($_[0] eq "ksi") {
    return $ksi_gt_1;
  }
  return $cs_err_tot_p;
}

###################################################################################################
sub list_prt {
  print "Listed prt files: @ARGV \n\n" ;
  my $tot_cs   = 0;
  my $tot_err  = 0;
  my $tot_evnt = 0;
  my ($cs,$cs_err);

  print "\n   Cross section [pb]  Error %   nCall     chi**2     N of events\n\n";
  $i = 0;
  foreach (@ARGV) {
    open(PRT,"<$_") or die "Can't open $_ ($!)";
    $/="\n";
    $cs = 0.0;
    while ($line = <PRT>) {
      if ($line =~ /\s\<\s\>\s+([-]*\d+\.\d+E[+-]\d+)\s+([-]*\d+\.\d+E[+-]\d+)\s+\d+\s+.+/) {
        $cs = $1;
      }
    }
    close(PRT);
    $tot_cs  = $tot_cs + $cs;
    $i++;
  }

  $i = 0;
  foreach (@ARGV) {
    open(PRT,"<$_") or die "Can't open $_ ($!)";
    $lline = 0.0;
    while ($line = <PRT>) {
      if ($line =~ /\#Subprocess/) {
        $nline = $line;
      }
      if ($line =~ /\s\<\s\>\s+([-]*\d+\.\d+E[+-]\d+)\s+([-]*\d+\.\d+E[+-]\d+)\s+\d+\s+.+/) {
        $cs = $1;
        $cs_err =abs($1 * $2) / 100.;
        $revnt = 10000;
        $negwarn = "";
        if (defined $nevnt) {
          if ($cs > 0.0) {
            $revnt = sprintf("%.0f", $nevnt * $cs / $tot_cs);
          } else {
            $negwarn = "The subprocess cross section is negative! no events...";
            $revnt = sprintf("%.0f", 0);
          }
        }
        $line =~ /(\s\<\s\>\s+[-]*\d+\.\d+E[+-]\d+\s+[-]*\d+\.\d+E[+-]\d+\s+\d+\s+.+)/;
        $lline = $1;
      }
    }
    close(PRT);
    print " $nline $lline $revnt $negwarn\n";
    $tot_evnt = $tot_evnt + $revnt;
    $tot_err = $tot_err + $cs_err * $cs_err;
    $i++;
  }
  my $tot_err_p = sqrt($tot_err) * 100/$tot_cs;
  printf "\n Total CS [pb] = %1.4e   +\/- %1.2e ( %1.2e %% )\n",$tot_cs,sqrt($tot_err),$tot_err_p;
  print " Total Number of events = $tot_evnt\n\n";
  exit(0);
}

###################################################################################################
# list subprocess names and cross sections in an event file with multiple subprocesses. 
# It's working for cpyth1 only!
sub list_mix {
  print "\n";
  foreach (@ARGV) {
    open(PRT,"<$_") or die "Can't open $_ ($!)";
    print "List of subprocesses in mixed event file: $_ \n\n" ;
    $/="\n"; 
    while (<PRT>) {
    if (/(#PROCESS|\#Cross_section|\#Number_of_subprocesses|\#Total_cross|\#Events_mixed)/g) {
      $_=~ s/\#// ; $_=~ s/\s{3,4}/ /g;
      $_=~ s/Cross_section\(Width\)/    Cross_section/g;
      print;
    }
    last if /\#Nproc ================== Events/;}
    print "\n\n";
    close PRT;
 }
 exit(0);
}

###################################################################################################
# changes the total number of events for each subprocess which should be prepared 
# so that the the number sum is equal to a requested number of events
sub replace_nevnt {
  my $oldfilebat;
  my @cs;
  my @cs_err;

  $i=0;
  $tot_cs = 0;
  $tot_nevt = 0;
  $tot_cs_err = 0;
  $oldfilebat = $filebat.".BAK";
  move ($filebat, $oldfilebat);
  open(BATCH, "+<$oldfilebat") or die "Can't open $oldfilebat for read/write: $!\n";
  while ($line = <BATCH>) {
    if ($line =~ /\#Vegas_integral\s+([-]*\d\.\d+E[+-]\d+)\s+([-]*\d\.\d+E[+-]\d+)\s+([-]*\d\.\d+E[+-]\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
      $cs[$i] = 0.0;
      if ($2/$1 > 0.0) {
        $cs[$i] = $2/$1;
      }
      $cs_err[$i] = sqrt(abs($1))/abs($2);
      $tot_cs_err = $tot_cs_err + $cs_err[$i] * $cs_err[$i];
      $tot_cs = $tot_cs + $cs[$i];
      $i++;
    }
  }
  seek(BATCH,0,0);

  $i=0;
  my $tmp="";
  while ($line = <BATCH>) {
    if ($line =~ /\#Subprocess/) {
      print $line;
    }
    if ($line =~ /(\#Events\s+\d+\s+\d+\s+\d\.\d+\s+\d\.\d+\s+)(\d+)/) {
      $new_nevnt = sprintf("%.0f", $nevnt * $cs[$i] / $tot_cs);
      $line = $1.$new_nevnt."\n";
      printf " Cross section [pb] = %1.4e +\/- %1.2e (Nevnt = %i)\n", $cs[$i], $cs_err[$i], $new_nevnt;
      $tot_nevt = $tot_nevt + $new_nevnt;
      $i++;
    }
    $tmp .= $line;
  }
  close(BATCH);

  open(BATCH, ">$filebat") or die "Can't open $filebat for read/write: $!\n";
  print BATCH $tmp;
  close(BATCH);
  print "\n";
  printf " Total cross section [pb]  = %1.4e +\/- %1.2e (Ntotal = %i)\n", $tot_cs, sqrt($tot_cs_err), $tot_nevt;
  exit(0);
}

###################################################################################################
sub replace_ncall_niter {
  $niter_encfactor = 1.0;
  $ncall_encfactor = 1.5;
  $ncall_encfactor = $_[0] if (defined $_[0]);
  $niter_encfactor = $_[1] if (defined $_[1]);
  @sessions = get_sessions ();
  @proc_n = 1..$#sessions+1;
  @proc_n = proc_nc () if defined $proc;
  $Nproc_current=1;
  $tmp="";
  foreach (@sessions) {
    if (grep {$Nproc_current eq $_} @proc_n) {
      @blocks=split(/\n#/,$_);
      foreach $blk (@blocks) {
        if ($blk =~/(Vegas_calls )(\d+)x(\d+)/) {
          print "replace_ncall_niter: 2 -> $2, 3 -> $3\n";
          $new_ncalls = int($ncall_encfactor*$2);
          $new_niters = int($niter_encfactor*$3);
          print "replace_ncall_niter: new_ncalls -> $new_ncalls, new_niters -> $new_niters\n";
          $blk = $1.$new_ncalls."x".$new_niters;
        }
        if ($blk =~/\+\+\+\+\+\+\+\+\+\+\+\+\+\+\n/) {
          $tmp .= $blk;
        } else {
          $tmp .= $blk."\n#";
        }
      }
    } else {
      $tmp .= $_;
    }
    $Nproc_current++;
  }
  move ($filebat, $filebat.".old");
  open (BATCH, ">$filebat") or die "$errr can't open $filebat in $opt_d: $!\n";
  print BATCH $tmp;
  close (BATCH);
}

###################################################################################################
sub check_chi_squared {
#  @sessions = get_sessions ();
#  @proc_n = 1..$#sessions+1;
#  if (defined $proc) {
#    @proc_n_orig = proc_nc ();
#    @proc_n = proc_nc ();
#  } else {
#  foreach (@sessions) {
#    if (grep {$Nproc_current eq $_} @proc_n) {
#      @blocks=split(/\n#/,$_);
#      foreach $blk (@blocks) {
#        if ($blk =~/(Vegas_calls )(\d+)x(\d+)/) {
#          print "replace_ncall_niter: 2 -> $2, 3 -> $3\n";
#	  $new_ncalls = int($ncall_encfactor*$2);
#          $new_niters = int($niter_encfactor*$3);
#          print "replace_ncall_niter: new_ncalls -> $new_ncalls, new_niters -> $new_niters\n";
#          $blk = $1.$new_ncalls."x".$new_niters;
#        }
#        if ($blk =~/\+\+\+\+\+\+\+\+\+\+\+\+\+\+\n/) {
#          $tmp .= $blk;
#        } else {
#          $tmp .= $blk."\n#";
#        }
#      }
#    } else {
#      $tmp .= $_;
#    }
#    $Nproc_current++;
#  }
#  move ($filebat, $filebat.".old");
#  open (BATCH, ">$filebat") or die "$errr can't open $filebat in $opt_d: $!\n";
#  print BATCH $tmp;
#  close (BATCH);
}

###################################################################################################
# changes the total cross sections in event files (kept either in cpyth1 or LHE formats)
sub change_cs {
  local $/="\n";

  foreach(@ARGV) {
    my $oldcs;
    my $newcs;
    my $file=$_;
    my $tmp = "";
    open(EVTFILE,"+<$file") or die "Can't open $file for read/write: $!\n";
    while (<EVTFILE>) {
      if(/(.+<crossSection.+errorMinus=.+)(\d+\.\d+E[+-]\d+)(.+errorPlus.+)(\d+\.\d+E[+-]\d+)(.+)(\d+\.\d+E[+-]\d+)\s*<\/crossSection>/) {
        $tmp .= "$1".sprintf("%E",$2*$kfactor)."$3".sprintf("%E",$4*$kfactor)."$5".sprintf("%E",$6*$kfactor)."</crossSection>\n";
        $oldcs = $6;
        $newcs = $6*$kfactor;
      }elsif(/(.+<crossSection.+)(\d+\.\d+E[+-]\d+)\s*<\/crossSection>/) {
         $tmp .= "$1 ".sprintf("%E",$2*$kfactor)."</crossSection>\n";
         $oldcs = $2;
         $newcs = $2*$kfactor;
      }elsif( /(.+section.+)\s(\d+\.\d+E[+-]\d+)\s+/ ){
        $tmp .= "$1 ".sprintf("%E",$2*$kfactor)."\n";
        $oldcs = $2;
        $newcs = $2*$kfactor;
      } else {
        $tmp .=$_;
      }
    }

    seek(EVTFILE,0,0);
    print EVTFILE $tmp;
    truncate (BATCH,tell(EVTFILE));
    close EVTFILE;
    print "file: $file, old cs: $oldcs, new cs: ".sprintf("%E",$newcs)."\n";
  }
  exit(0);
}

###################################################################################################
sub replace_calls {
  my $oldfilebat = $filebat.".BAK";
  my $tmp = "";
  move ($filebat, $oldfilebat);

  open(BATCH, "+<$oldfilebat") or die "Can't open $oldfilebat for read/write: $!\n";
  while ($line = <BATCH>) {
    if ($line =~ /(\#Vegas_calls\s+)(\d+)(x\d+)/) {
      $line = $1.$ncall.$3."\n";
    }
    $tmp .= $line;
  }
  close(BATCH);

  open(BATCH, ">$filebat") or die "Can't open $filebat for read/write: $!\n";
  print BATCH $tmp;
  close(BATCH);

  print "Numbers of calls have been replaced\n";
  exit(0);
}

###################################################################################################
sub replace_max {
  my $oldfilebat = $filebat.".BAK";
  my $tmp = "";
  move ($filebat, $oldfilebat);

  open(BATCH, "+<$oldfilebat") or die "Can't open $oldfilebat for read/write: $!\n";
  while ($line = <BATCH>) {
    $line =~ s/Max:/!Max/;
    $tmp .= $line;
  }
  close(BATCH);

  open (BATCH, ">$filebat") or die "Can't open $filebat for read/write: $!\n";
  print BATCH $tmp;
  close (BATCH);
}

###################################################################################################
sub help
{
  print
  "Usage: \n
$0 [-d dir] [-b file] [-run [vegas,max,evnt,cleanstat,cleangrid,clean]] [-proc N1,N2-N3,N4] [-nses N] [-cpyth1] [-cpyth2] [-ginv] [-mp N] [-nice NN] [-verbose]
$0 [-d dir] [-b file] [-proc N1,N2-N3,N4] [--add [file]]
$0 [-d dir] [-b file] [-auto X] [-nevnt Nevents]
$0 [-d dir] [-b file] [-nevnt Nevents]
$0 [-d dir] [-b file] [-call Ncalls]
$0 [-d dir] [-b file] [-geninfo]
$0 [-d dir] [-b file] [-extract N1]
$0 [-d dir] [-b file] [-show [proc,inistat,params,cuts,regs,qcd,vegas,events,cs]]
$0 [-lcs [prt file names]]
$0 [-lmix [names of files with mixed events]]
$0 [-kfactor [factor] -files [event file names]]
$0 [-h] [--help]
$0 [--version]

/********************************************************************************************/
/***                    Launch numerical CompHEP in batch (single machine)                ***/
-run vegas         - run vegas calculations;
-run max           - run maxima search in cubes;
-run evnt          - generate events;
-run               - run all three steps: vegas/maxima search/ event generation;

-run cleanstat     - clean statistics only (batch.dat remembers the vegas grid);
-run cleangrid     - clean grid only (batch.dat remembers statistics);
-run clean         - clean both statistics and grid;
-proc n1,n2-n3,n4  - run calculation for these processes only;
-nses n            - n is a number of the first session, it is 
                     calculated automatically by default; 
-ginv              - gauge invariance ON (See CompHEP manual); 
-safe              - prevent to launch n_comphep if it is locked by the LOCK file; 
-nice nn           - run CompHEP with modified scheduling priority (nice), nn - is priority, 
                     if it is undefined nn = 10;
-cpyth1            - generate events in the old CompHEP format cpyth1, 
                     compatible with cpyth-1.*;
-cpyth2            - generate events in the old CompHEP format cpyth2, 
                     compatible with cpyth-2.*;
-d dir_name        - use the directory new_dir instead of results for calculations;
-b file_name       - use the file file_name instead of batch.dat for calculations;
-auto X            - set needed calculation accuracy. Be default, it is an accuracy of the 
                     total cross sectoin. If it is combined with -nevnt it is stat. accuracy
                     for each subprocc.

/********************************************************************************************/
/***                        Batch system or Parallel regime (computer farms)              ***/
-mp N              - parallel calculation of N subprocesses on one machine
-pbs [pbs prefix]  - run calculations in parallel with the PBS batch system (e.g. qsub);
-lsf [lsf prefix]  - run calculations in parallel with the LSF batch system (e.g. bsub -I);
-nocombine         - do not combine automatically the results of
                     PBS (LSF) calculations from the temporary subdirectories
                     it can be done with the -combine option later;
-combine           - combine the results of previous calculations from the
                     temporary PBS directories (you need it only if you used
                     -nocombine option before);

/********************************************************************************************/
/***                          Statistics and file manipulations                           ***/
-show ...          - print different parts of batch file
   [proc]          - subprocess list
   [inistat]       - initial state: beams, sqrt(S), PDF (of each subprocess)
   [params]        - print parameter list (for each subprocess)
   [cuts]          - print cuts applied list (for each subprocess)
   [regs]          - print regularisation list (for each subprocess)
   [qcd]           - print qcd information (for each subprocess)
   [vegas]         - print vegas information (for each subprocess)
   [events]        - print events information (for each subprocess)
   [cs]            - print cross section of each subprocess and the full one;

-lcs  LIST         - report all cross sections in the results from prt_* files,
                     LIST is the list of prt files (e.g. prt_* or prt_1 prt_2 ...);
-lmix LIST         - list of subprocesses and cross sections from the mixed event 
                     files from the LIST;
--add [file]       - add parameters from exist session.dat
                     file (or [file])  to batch.dat file in [dir] for the
                     subprocesses -proc n1,n2-n3,n4  or one subprocess from
                     session.dat (or [file]) without option -proc ...
-kfactor {factor}  - change cross sections in event files, the event file
                     names are set by the command -files
-files {filenames} - set event file names to change cross sections in the
                     files (by an option -files)
-nevnt Nevents     - set the full number of events. The script sets the number
                     of event for each subprocess according to its contribution
                     to the full cross section;
-call Ncalls       - set a number of calls in vegas sessions;
-geninfo           - copy several sections from session.dat to batch.dat for all subprocesses
-extract N1        - extract session.dat from batch.dat for the subprocess N1

/********************************************************************************************/
/***                             Help and other options                                   ***/
-h                 - print this message
--help             - print the long help message
--version|         - print version of the script\n";

  if (defined $longhelp) {
    print
"

/********************************************************************************************/
1. Introduction.

  Numerical batch mode is a useful tool to calculate cross sections and generate 
  events in tasks with the large number of subprocesses or where large computer 
  resources are required. It helps to perform the large scale calculations, even 
  the calculations on computer farms (with the PBS and LSF batch systems installed). 

  At first, a user should prepare a numerical CompHEP program - event generator -
  in GUI mode (see details in the CompHEP manual). This script uses special 'blind' 
  keys to launch the program in the batch regime. The numerical program can also 
  be prepared by the symb_batch.pl script. 

/********************************************************************************************/
2. Configuration of the calculation job.

  At first, the Monte-Carlo generator should be configured. All parameters for 
  numerical sessions are stored in the file batch.dat. It is built with session.dat. 
  So, the user should prepare session.dat with the generator GUI interface. 

  cd results
  ./n_comphep&

  The user sets all required parameters (PDF, QCD scale, kinematic cuts,
  regularizations, etc.) for the first process, check the vegas menu and 
  terminate the program. The session.dat file will be created. 

  After that, the user runs $0 without options in order to create the configuration
  file results/batch.dat with the same parameters for all subprocesses
  as (s)he sets for the first subprocess. The script must run from the CompHEP
  working directory (CWDIR):

  $0

  By default, the directory 'results' is used by this script to find and
  run the generator n_comphep.exe. The user can change the default name by
  the option -d {new_dir}. Use the option for all manipulation where you need
  a reference to the n_comphep.exe, e.g.:

  $0 -d new_dir

  If you want to change parameters of a particular subprocess (or a set of subprocesses) 
  you can edit the results/batch.dat file manually or use GUI: choose the required 
  subprocess in GUI, change the parameters, check the vegas menu and return to the main 
  GUI menu (this operation saves the parameters to session.dat). After that, the
  script can insert the new parameters to the results/batch.dat file (add session.dat
  to batch.dat) if you launch it with the command:

  $0 --add

  If you use the option with -proc n1,n2-n3,n4 it adds exactly the same parameters as in 
  session.dat (or FILE from --add [FILE]) to results/batch.dat (or to NEW_DIR/BATCH_FILE 
  with options -d NEW_DIR -b BATCH_FILE) for the the n1,n2-n3,n4 subprocesses.

/********************************************************************************************/
3. How to start the calculations.

  If the configuration process of the program has been finished you are ready
  to start calculations. It can be done with the command -run:

  $0 -run

  with this command n_comphep.exe will calculate cross sections, find maxima
  and generate events for all of the subprocesses.

  The full process can be divided to separate steps:
  $0 -run vegas   : calculation of cross sections
  $0 -run max     : search for maxima
  $0 -run evnt    : event generation, can be used after vegas
                    and max calculations only!

  There are several additional options, useful in practical calculations:
  $0 -run cleanstat   : clean the vegas statistic
  $0 -run cleangrid   : clean the vegas grid
  $0 -run clean       : clean the vegas grid and statistic

  The option -proc allows a specific set of subprocesses be calculated.
  For instance,

  $0 -run -proc 1,3-5,17,2

  means calculations for 1st, 3rd, 4th, 5th, 17th, 8th subprocesses.
  Be default, CompHEP uses the LHE format for event files (see details in 
  hep-ph/0609017), and we ecourage the format since it is a de facto 
  standard in the comunity of authors of Monte-Carlo codes and it is useful 
  as a simple and independent of software format for software environments 
  of experiments. But the old CompHEP formats are still supported via 
  additional options -cpyth1   (compatible with the cpyth-1.* ComHEP-PYTHIA 
  interface) and -cpyth2 (compatible with the cpyth-2.* ComHEP-PYTHIA 
  interface). 

  CompHEP stores all session details in protocol files, prt_XYZ, where XYZ is
  equal to the session number. By default, the number of a new session will 
  be XYZ+1, where XYZ is the last session number. The user can change the 
  behaviour with the option -nses. E.g. 

  $0 -d new_dir -run vegas -proc 2-10,12 -nses 20 -nice 19

  means calculations for the subprocesses 2-10,12 and the first session (for the
  2nd subprocess) will have the number 20. n_comphep.exe located in the directory
  new_dir will do thew calculations. Option -nice says that n_comphep.exe will be 
  running with non-zero priority in order to lower CPU demand for the task (nice 
  is the standard tool for that in Linux). Permitten range for nice is 1-19, the 
  default the status is equal to 10.

  WARNING! if prt_20 exists already, n_comphep.exe will append the details of the
  2th subprocess to the file, even the file corresponds to another subprocess.

/********************************************************************************************/
4. Parallel calculations. 

  The commands described above provide a simple method of calculation of several 
  subprocess one after another. This way is appropriate on a machine with one CPU. 
  If there is a machine with several CPU available, several jobs can be run at 
  the same time. For the parallel calculations of subprocesses a special option 
  should be used: 

  $0 -run vegas -mp N [-nocombine] [-proc n1,n2-n3,n4] [-nice NN] [-verbose]

  where N corresponds to the number of CPUs/cores on your machine (N subprocesses 
  will be run at the same time). As soon as one subprocess has been completed 
  the next subprocess starts. This command will create temporary subdirectories 
  for each subprocess and run each subprocess in the corresponding temporary 
  directory. Since all of the subprocesses have been calculated $0 combines results 
  automatically and removes the temporary directories. The option -nocombine helps 
  to prevent the behaviour, in this case all the subdirectories are kept. One can 
  be possible to combine the results later with the option -combine. -nice can be 
  used in the regime too. 

  $0 -combine

/********************************************************************************************/
5. Calculations with LSF or PBS systems.

  If a user has an access to a BATCH system (LSF or PBS), there is a way 
  to launch CompHEP in parallel regime. There are two options for launching 
  CompHEP jobs by means of the systems: -lsf and -pbs with almost the same 
  functionality. 

  $0 -run vegas -pbs

  This command creates temporary subdirectories for each subprocess, prepares 
  PBS batch jobs, and runs the jobs for all subprocess with the PBS command: 
  qsub -I run.sh

  or for LSF system
  $0 -run vegas -lsf
  will submit
  bsub -I run.sh

  As soon as the jobs have been completed $0 saves the results in 
  the usual CompHEP way and removes the temporary subdirectories automatically. 

  If your PBS or LSF system does not provide you the interactive option ( -I ) 
  you must add the option -nocombine: 

  $0 -run vegas -lsf  -nocombine

  In this case $0 submits the LSF jobs and stops immediately. As soon as 
  all of the submitted jobs have been finished the user can combine results 
  collected in temporary subdirectories by an option: 

  $0 -combine

  If the user wants to change the default method to submit the jobs it is possible
  to run something like:

  $0 -run vegas -lsf 'bsub -I -q large -N -B -o log.out'

  Interactive parallel calculations can be launched with the empty prefix of the
  command:

  $0 -run vegas -lsf ' '

  In this case $0 starts calculations of the subprocesses in parallel and 
  combines the results at the end. 

  WARNING 1!
  Since each batch system has its own machinery and configuration subtleties, 
  we can not guarantee that CompHEP can start on your batch system even if 
  you use LSF or PBS. 

  WARNING 2! 
  This way of calculating presumes a direct access to n_comphep.exe and the 
  CompHEP installation area (to PDF tables). So, if your computer farm nodes use 
  separate file systems, you must edit the prepared scripts run.sh (in res_i, where 
  i=1,2,...) and add a code which copies all necessary files (the subdirectory
  and the PDF tables from the CompHEP area) to the farm nodes. Certainly, the
  code should copy all results back, to the corresponding subdirectories.

/********************************************************************************************/
6. Presentation of results.

  There are two special options to summarise information in batch.dat or
  in protocol files.

  The option -lcs prints the cross sections and statistical errors from the
  protocol files. Examples of usage of the command:

  $0 -lcs prt_1 prt_2 prt_5
  $0 -lcs results/prt_*

  The option -lmix prints the mixed subprocesses and cross sections from the
  mixed event files

  $0 -lmix mix_events*.pev

  The command -show prints separate parts of batch.dat. For example,
  -show proc    : prints the full subprocess list,
  -show inistat : prints initial states: beams, sqrt(S), PDF
  -show params  : prints the list of physical parameters
  -show cuts    : prints the cuts applied list
  -show regs    : prints the regularization list
  -show qcd     : prints qcd information
  -show vegas   : prints vegas information
  -show events  : prints events information
  -show cs      : prints the the total one and cross sections of each subprocess

  Examples:

  $0 -show inistat
  $0 -d new_dir -show cs

/********************************************************************************************/
7. Other manipulations in batch mode.

  By default, CompHEP sets the numbers of events for the event generation equal each
  other. Since subprocesses have different cross sections, not all of the events will
  be used after the generation. The command -nevnt sets the total number of events
  for all subprocesses according to their contributions to the total cross section.
  For example,

  $0 -nevnt 100000

  means the full number of events generated by n_comphep.exe in all event files will
  be 100000. Since the further mixing of the event files is a statistical process,
  the real number of events in the final file will be between 100000 - sqrt(100000)
  and 100000.

  The command -call changes the number of calls in vegas sessions. Usage:

  $0 -call 100000

  The commands -kfactor and -files are useful if a user wants to change cross
  sections of subprocesses in the event files prepared, so batch.dat is not changed.
  Sometimes the option is needed. Usage:

  $0 -kfactor 2.0 -files results/events_*

  It means all cross section in the files (results/events_*) will be doubled.

  \n";

  }
  exit (0);
}

###################################################################################################
sub version {
  print "$0 (CompHEP numerial tool) $version, date: $datestamp\n\n";
  print "Copyright (C) 2002-2009, CompHEP Collaboration.\n";
  print "This software is distributed under the CompHEP licence.\n";
  print "All details of the licence see at http://comphep.sinp.msu.ru\n";
  exit (0);
}
