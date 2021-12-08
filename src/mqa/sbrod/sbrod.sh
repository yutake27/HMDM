#!/bin/bash
#$ -cwd
#$ -l s_core=1
#$ -l h_rt=0:10:00
#$ -j y

target_pdb_dir=${1}
dir=MQA_dataset

source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
conda activate mypython

python sbrod.py ${target_pdb_dir}
