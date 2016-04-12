# test a simple job array script
bsub -n 1 -J test[1-10] -We 30 -R "rusage[mem=1]" -eo err.%I -oo log.%I 'echo
$LSB_JOBINDEX >> test.txt'

# use job array to call another script
#bsub -n 1 -R "rusage[mem=8]" -We 60 -J myruns[1-10] -oo log.%I.out -eo
#log.%I.err 'Rscript /home/rathik/scripts/bsub_testscripts/test_jobarray.R
#$LSB_JOBINDEX'
