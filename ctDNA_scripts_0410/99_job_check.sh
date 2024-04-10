#!/bin/bash 

path=$1 

cd $path

# first, create a list of all the slurm output files
ls slurm-* > slurm_outputs.txt

# then extract the job info from each slurm output file
echo "all_jobs" > job_status.txt
while read -r LINE; do
    slurm_id=$(echo $LINE | cut -d'.' -f1 | cut -d'-' -f2)
    sacct --format JobId,JobName,NNodes,Partition,NCPUs,State,ReqMem,MaxRSS,Elapsed,CPUTime,TimeLimit,ExitCode,Start,End -j $slurm_id >> job_status.txt
done < slurm_outputs.txt
