#!/bin/bash

#$ -cwd
#$ -l s_core=1
#$ -l h_rt=1:30:00

source /home/4/16B09097/.bashrc
conda activate mypython

sf_id=$1
dataset=$2
python scop_modeling.py $sf_id -d $dataset
