#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 2:00:00
#SBATCH --mem=8000
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=FAIL

#######################
# Example call below:
# sh run_tumor_processing.sh Example 2 /labs/mwinslow/Emily/tubaseq_pipeline_upload/ 4
#######################

#ml python/3.6.4
#module load miniconda/3
python3 process_tumors.py "--project=${1}" "--parameter=${2}" "--root=${3}"
