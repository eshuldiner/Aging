#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 2:00:00
#SBATCH --array=1-4
#SBATCH --mem=25g
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=FAIL


#######################

arg="${3}tubaseq_inp_files/${1}_${2}.inp"
parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${arg})
python3 convert_to_cells.py $parameters
