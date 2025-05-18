#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 4:00:00
#SBATCH --array=1-1
#SBATCH --mem=25g
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=FAIL

# sh run_filtering.sh JDH9_L1_2 2c /labs/mwinslow/Emily/ 1
# sbatch run_filtering.sh JDH9_L3_4 2c /labs/mwinslow/Emily/ 4
#sh run_filtering.sh
#######################


#ml python/3.6.4
#module load miniconda/3

#arg=/labs/mwinslow/Emily/tubaseq_inp_files/testerarray.inp
arg="${3}tubaseq_inp_files/${1}_${2}.inp"
#arg=/labs/mwinslow/Emily/tubaseq_inp_files/JDH9_L1_2_2c_LS03_69_L2_redo.inp
#arg=/labs/mwinslow/Emily/tubaseq_inp_files/JDH9_L3_4_redos.inp
#arg=/labs/mwinslow/Emily/tubaseq_inp_files/EA_Variants_2.inp

echo $arg
#Parameters
parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${arg})
echo $parameters

python3 filter_tumors.py $parameters

