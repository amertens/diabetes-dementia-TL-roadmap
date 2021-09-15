#!/bin/bash
# Job name:
#SBATCH --job-name=test
#
# Account:
#SBATCH --account=co_biostat
#
# Partition:
#SBATCH --partition=savio3
#SBATCH --qos=biostat_savio3_normal
#
# Specify one task:
#SBATCH --ntasks=1
#
# Number of processors for threading:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#
# Wall clock limit:
#SBATCH --time=72:00:00
#SBATCH --error=test_job_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wangzeyi@berkeley.edu
#
## Command(s) to run (example):
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load r/3.6.3
R CMD BATCH --no-save /global/home/users/wangzeyi/exps/2021/0915/parallel_correctcsr-nospeedglm-noenv-meanltmle_202109_5k_SL.glm/1.R /global/home/users/wangzeyi/exps/2021/0915/parallel_correctcsr-nospeedglm-noenv-meanltmle_202109_5k_SL.glm/1.Rout
