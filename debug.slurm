#!/bin/bash
#SBATCH --job-name=radialprofs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:02:00
#SBATCH --partition=debug
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fmcaporaso@mi.unc.edu.ar


. /etc/profile

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun julia jl_flagship/PERFILES.jl