#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 6:00:00
#SBATCH --mem=200
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=NONE

#######################
# Example call:
# sh tubaseq.sh Example 2 /labs/mwinslow/Emily/tubaseq_pipeline_upload/ 4
#######################

PROJECT_ID=$1
PARAMETER_ID=$2
ROOT=$3
NSAMPLES=$4
array_lab="1-$4"

python3 project_set_up.py "--project=${PROJECT_ID}" "--parameter=${PARAMETER_ID}" "--root=${ROOT}"

#Processes gzipped fastq files in parallel, returns files with raw read counts for tumors.
jid1=$(sbatch --array=$array_lab --parsable run_count_reads_array.sh ${PROJECT_ID} ${PARAMETER_ID} ${ROOT})

#Filtering and clustering or barcodes
jid2=$(sbatch --dependency=afterany:${jid1} --array=$array_lab --parsable run_filtering.sh ${PROJECT_ID} ${PARAMETER_ID} ${ROOT}) 

#Convert reads to cells
jid3=$(sbatch --dependency=afterany:${jid2} --array=$array_lab --parsable run_convert_to_cells.sh ${PROJECT_ID} ${PARAMETER_ID} ${ROOT})

sbatch --dependency=afterany:${jid3} --parsable run_status_check.sh ${PROJECT_ID} ${PARAMETER_ID} ${ROOT}

#Process tumors
jid4=$(sbatch --dependency=afterany:${jid3} --parsable run_tumor_processing.sh ${PROJECT_ID} ${PARAMETER_ID} ${ROOT})
