#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=1:20:00
#$ -j y
#$ -m a
#$ -M takei@cb.cs.titech.ac.jp


fasta=${1}
data_dir=${2}

source /home/4/16B09097/.bashrc
conda activate mypython


python ../../../../P3CMQA/src/preprocess.py -f ../../../fasta/${data_dir}/${fasta}.fasta -o ../../../profile/${data_dir} -d $uniref90 -n 14
