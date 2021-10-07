#!/bin/bash
pref=running_code_to_generate_expected_IBD_and_ROH_distribution
memory=1000
# bsub -J "$pref[1-72]" -n 1 -o cluster_output/$pref.%J.%I.out -e cluster_output/$pref.%J.%I.err -q normal -G team281 -R "select[mem>$memory] rusage[mem=$memory]" -M $memory /lustre/scratch115/projects/autozyg/BiB/code_for_popgen_paper/run_code_to_generate_expected_IBD_and_ROH_distribution.sh

mylist=parameters_to_run.unique.txt
whichjob=$LSB_JOBINDEX

max_T_exact_model=`head -n $whichjob  $mylist|tail -n 1|cut -f 1`
oldkinship=`head -n $whichjob  $mylist|tail -n 1|cut -f 2`
whichRecentConsang=`head -n $whichjob  $mylist|tail -n 1|cut -f 3`
whichNeTrajectory=`head -n $whichjob  $mylist|tail -n 1|cut -f 4`
whichpop=`head -n $whichjob  $mylist|tail -n 1|cut -f 5`

/software/team163/mf13/R-3.4.0/bin/Rscript run_code_to_generate_expected_IBD_and_ROH_distribution.flexible_options.R $max_T_exact_model $oldkinship $whichRecentConsang $whichNeTrajectory $whichpop
