#!/bin/bash
#PBS -l walltime=71:00:00,select=1:ncpus=1:mem=16gb
#PBS -N cofracQC
#PBS -o /scratch/st-ljfoster-1/logs/cofracQC/cofracQC.out
#PBS -e /scratch/st-ljfoster-1/logs/trinary/cofracQC.err
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS -A st-ljfoster-1

if [ -e ~/environments/ppicluster_venv/bin/activate]
then
  source ~/environments/ppicluster_venv/bin/activate
else  
  echo "File not found: ~/environments/ppicluster_venv/bin/activate"
fi
module load gcc/9.1.0
module load openmpi/3.1.5
module load netcdf/4.7.3
module load r/3.6.2-py3.7.3

cd /scratch/st-ljfoster-1/staceyri/trinary/

Rscript /scratch/st-ljfoster-1/staceyri/cofracQC/R/measure_all_datasets.R 
