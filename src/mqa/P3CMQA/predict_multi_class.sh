#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=1:30:00
#$ -j y

target=${1}
dataset=${2}
dir=MQA_dataset

source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load intel/17.0.4.196
module load cuda/9.2.148 cudnn/7.6
conda activate mypython


cd ../../../../P3CMQA/src

python predict_multi_class.py -d ../../${dir}/scwrl_pdb/${dataset}/${target} -f ../../${dir}/fasta/${dataset}/${target}.fasta -o ../../${dir}/score/${dataset}/3dcnn_multi -g 0 -p ../../${dir}/profile/${dataset}/${target}
