#! /usr/bin/env bash

RUN=$1
OF=scripts/run_${RUN}.sh

echo > $OF

echo "#PBS -l walltime=24:00:00" >> $OF
echo "cd /user1/pwaller/master-a4/bumps/bat" >> $OF
echo "build/diphoton2012-limitsetting $RUN 1 1000 $RUN 1 1 1 1 0 0 1 output-supersyst IsoTemplate_c0.01_SuperSyst &> log_${RUN}" >> $OF

chmod u+x scripts/run_${RUN}.sh

echo "echo Completed run | mail -s 'Completed run ${RUN}' pwaller" >> $OF
