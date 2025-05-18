#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 2:00:00
#SBATCH --array=1-4
#SBATCH --mem=25g
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=FAIL


#######################
#sh run_convert_to_cells.sh Jackie1 0 /labs/mwinslow/Emily/ 13
#sh run_convert_to_cells.sh AgingScreen2KPTC_AllLanes 2 /labs/mwinslow/Emily/ 92
#sh run_convert_to_cells.sh AgingScreen1FiltV2 2c /labs/mwinslow/Emily/ 1
# sh run_convert_to_cells.sh JDH9_L1_2 2c /labs/mwinslow/Emily/ 1
# sbatch run_convert_to_cells.sh JDH9_L3_4 2c /labs/mwinslow/Emily/ 4

# sh run_convert_to_cells.sh AgingInertTransplantV2 2 /labs/mwinslow/Emily/ 49

#ml python/3.6.4
#module load miniconda/3

arg="${3}tubaseq_inp_files/${1}_${2}.inp"
#arg=/labs/mwinslow/Emily/tubaseq_inp_files/JDH9_L1_2_2c_LS03_69_L2_redo.inp
#arg=/labs/mwinslow/Emily/tubaseq_inp_files/JDH9_L3_4_redos.inp


echo $arg
#Parameters
parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${arg})

#echo $parameters

#python3 convert_to_cells.py --sample=MT501_merged --sample_path=/oak/stanford/groups/dpetrov/emilys/AgingScreen1/MT501_merged_R1_001.fastq.gz --barcode_length=20 --spikein_barcode_lengths=20,20,20 --sgRNAs=Apc_1,Apc_2,Arid1a_1,Arid1a_2,Arid1b_1,Arid1b_2,Atm_1,Atm_2,Brca2_1,Brca2_2,Cdkn2a_1,Cdkn2a_2,Cdkn2c_1,Cdkn2c_3,Cmtr2_1,Cmtr2_2,Dnmt3a_1,Dnmt3a_2,Fbxw7_1,Fbxw7_2,Kdm6a_1,Kdm6a_2,Keap1_1,Keap1_2,Kmt2c_1,Kmt2c_2,Lkb1_1,Lkb1_4,NT3,NT4,Nf1_1,Nf1_2,Nf2_1,Nf2_2,NT1,NT2,p53_1,p53_4,Pcna_1,Pcna_2,Pole_1,Pole_2,Pten_1,Pten_2,Rb1_1,Rb1_4,Rbm10_1,Rbm10_4,Rnf43_1,Rnf43_2,Setd2_1,Setd2_V4,Smad4_1,Smad4_2,Smarca4_1,Smarca4_2,Stag2_1,Stag2_2,Spi --sgids=AGGAGTCC,GATTCTCG,ACAATGGC,TCGGATCT,TTAGCCTG,AATCCGGT,CCTCAACA,TGACGTTC,CGCTACTT,CGACTAGT,GACCATAG,CACTGCTA,CTGAACCT,TGAGTGTC,ACTGCATG,CTAGTGAC,ACCGATTC,TCTGCCTT,GTCCATTC,GTTATGGC,GCAGCTAT,GCTATACG,GAGAGGTA,TTCCGAAG,CAAGCTCA,GTTGCAAG,GCTCTAGT,ACCTCTAC,CTAGGCTA,AATGCTGG,ACGTGTGT,AATTGCCG,ACTCTACC,CGATAGCT,AGTTGCTC,ACGTCGAA,CTACTCTC,ATACAGGC,ACAGTCGT,TCGCTATG,GAGCTTGA,ACGTACAG,ATCGTTGG,ATCACGTG,GAGAACTC,TCATGAGC,ACCTTAGG,GCTCACTT,AACACAGG,ACGGAACA,TGGTCCTT,TTGCATCG,GCAATGCT,GAGGTGTT,ACAATCCG,TGGTTCCA,GCAACAAG,CCTGAATC,ACGCTAGA --bc_dist=0 --R1_regex_looseness=2 --R2_regex_looseness=3 --root=/labs/mwinslow/Emily/ --project_name=AgingScreen1FiltV2 --parameter_name=2c

python3 convert_to_cells.py $parameters
