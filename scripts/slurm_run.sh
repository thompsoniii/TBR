#!/bin/bash
#SBATCH --job-name=senior_design_testing           # Job name
#SBATCH --mail-user=kthompson309@gatech.edu # E-mail address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail preferences
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=00:30:00
#SBATCH --output=output_logs/testing.out

## MODIFY FOLLOWING SECTION ##

# define sbatch nodes, ntasks-per-node, mem-per-cpu, time, output
module load anaconda3
conda init
conda activate openmc
# pip install -e .
date
srun -np 24 python arc-2.py 14
# python3 arc-2.py 13
# python3 arc-2.py 12
date
# load appropiate modules
# run python code
# print date
