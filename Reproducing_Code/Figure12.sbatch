#!/bin/bash
#SBATCH --job-name fir_deri1        # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition short         # Slurm partition to use
#SBATCH --array=1-100               # How many tasks in the array
#SBATCH --time 0-11:40             # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=1000           # Memory limit for each tasks (in MB)
#SBATCH -o out/fir_deri1_%j.out        # File to which STDOUT will be written
#SBATCH -e out/fir_deri1_%j.err        # File to which STDERR will be written
#SBATCH --mail-type=ALL            # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yikun@uw.edu # Email to which notifications will be sent

module load R
i=${SLURM_ARRAY_TASK_ID}
Rscript ./Figure12.R ${i}
