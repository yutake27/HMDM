#!/bin/bash

#$ -cwd
#$ -l s_core=1
#$ -l h_rt=3:00:00

source /home/4/16B09097/.bashrc
conda activate mypython

target=$1
dataset=$2
python pisces_modeling.py $target -d $dataset
