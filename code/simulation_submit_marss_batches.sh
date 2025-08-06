#!/bin/bash

#SBATCH --job-name=submit_marss
#SBATCH --account=your_account_name ######!!!! Please fill in your slurm account name
#SBATCH --partition=standard
#SBATCH --output=/home/%u/log/submit_marss.out   
#SBATCH --error=/home/%u/log/submit_marss.err 
#SBATCH --mail-type=END 
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
module load R/4.3.2-mkl ##### modify to accommodate your R version
for i in {1..16}  
do
  echo "Starting batch $i at $(date '+%Y-%m-%d %H:%M:%S')"
  Rscript code/simulation_submit_marss_batches.R $i
  sleep 2h
done