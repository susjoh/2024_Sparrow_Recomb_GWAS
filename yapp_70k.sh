#!/bin/bash

#$ -V
#$ -cwd
#$ -N yapp_70k
#$ -o o_files/
#$ -e e_files/
#$ -l h="c6"


# Array job run as  qsub -t 1-29 yapp_70k.sh

# activate genedrop conda environment

source /ceph/software/conda/etc/profile.d/conda.sh &&
conda activate sjohnston &&

# make directory in /scratch for this job

SCRATCH=/scratch/$USER/$JOB_ID/yapp &&
mkdir -p $SCRATCH  &&

# run command
#yapp mendel 70K_200K_maf_geno_mind_v5
yapp phase --reg $SGE_TASK_ID 70K_200K_maf_geno_mind_v5
yapp recomb 70K_200K_maf_geno_mind_v5_$SGE_TASK_ID

#yapp recomb 70K_200K_maf_geno_mind_v5_1
#yapp recomb 70K_200K_maf_geno_mind_v5_2
#yapp recomb 70K_200K_maf_geno_mind_v5_3
#yapp recomb 70K_200K_maf_geno_mind_v5_4
#yapp recomb 70K_200K_maf_geno_mind_v5_5
#yapp recomb 70K_200K_maf_geno_mind_v5_6
#yapp recomb 70K_200K_maf_geno_mind_v5_7
#yapp recomb 70K_200K_maf_geno_mind_v5_8
#yapp recomb 70K_200K_maf_geno_mind_v5_9
#yapp recomb 70K_200K_maf_geno_mind_v5_10
#yapp recomb 70K_200K_maf_geno_mind_v5_11
#yapp recomb 70K_200K_maf_geno_mind_v5_12
#yapp recomb 70K_200K_maf_geno_mind_v5_13
#yapp recomb 70K_200K_maf_geno_mind_v5_14
#yapp recomb 70K_200K_maf_geno_mind_v5_15
#yapp recomb 70K_200K_maf_geno_mind_v5_17
#yapp recomb 70K_200K_maf_geno_mind_v5_18
#yapp recomb 70K_200K_maf_geno_mind_v5_19
#yapp recomb 70K_200K_maf_geno_mind_v5_20
#yapp recomb 70K_200K_maf_geno_mind_v5_21
#yapp recomb 70K_200K_maf_geno_mind_v5_22
#yapp recomb 70K_200K_maf_geno_mind_v5_23
#yapp recomb 70K_200K_maf_geno_mind_v5_24
#yapp recomb 70K_200K_maf_geno_mind_v5_25
#yapp recomb 70K_200K_maf_geno_mind_v5_26
#yapp recomb 70K_200K_maf_geno_mind_v5_27
#yapp recomb 70K_200K_maf_geno_mind_v5_28
#yapp recomb 70K_200K_maf_geno_mind_v5_29

# move scratch back to directory

rsync -av $SCRATCH ./

rm -rf $SCRATCH

#deactivate conda environment
conda deactivate 