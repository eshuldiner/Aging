#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 5:00:00
#SBATCH --array=1-24
#SBATCH --mem=10g
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=FAIL

#######################

# Sample call:
# sh run_count_reads_array.sh HCTS102 3 /labs/mwinslow/Emily/ 63

#ml python/3.6.4
#module load miniconda/3

arg="${3}tubaseq_inp_files/${1}_${2}.inp"
echo $arg

#Parameters
parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${arg})

python3 count_reads.py $parameters
#python3 count_reads.py --sample=SL664 --sample_path=/oak/stanford/groups/dpetrov/emilys/SingleGene3TubaSeq/01.RawData/LA66_ES_17_L1_1.fq.gz --barcode_length=20 --spikein_barcode_lengths=20,20,20 --sgRNAs=NT1,NT2,Safe19,Safe34,Safe36,Spi --sgids=AGTTGCTC,ACGTCGAA,AAGAGGTC,AAGGCTAG,CATAGCTC,ACGCTAGA --bc_dist=0 --R1_regex_looseness=2 --R2_regex_looseness=3 --root=/labs/mwinslow/Emily/ --project_name=AgingSingleGene3 --parameter_name=2c

