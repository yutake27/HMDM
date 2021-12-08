#!/bin/bash
#$ -cwd
#$ -l s_core=1
#$ -l h_rt=0:30:00

source /home/4/16B09097/.bashrc
conda activate mypython
SF_ID=$1
BLASTDB=../../../../../blastdb/pdbaa_20210407/pdbaa
scop_sf_fasta=../../../scop/20210330/scop_sf_represeq_lib20210330.fa
python psiblast4sf.py $SF_ID --blastdb ${BLASTDB} --scop_sf_fasta ${scop_sf_fasta}
