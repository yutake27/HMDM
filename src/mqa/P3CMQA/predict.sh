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

mkdir -p ../../../scwrl_pdb/${dataset}/${target}

cd ../../../pdb/${dataset}/${target}/sampling
scwrl_dir=../../../../scwrl_pdb/${dataset}/${target}

for i in *;do
    if [ ! -e ${scwrl_dir}/${i} ];then
        Scwrl4 -i ${i} -o ${scwrl_dir}/${i}
    fi
done

cd ../../../../../P3CMQA/src

python predict.py -d ../../${dir}/scwrl_pdb/${dataset}/${target} -f ../../${dir}/fasta/${dataset}/${target}.fasta -o ../../${dir}/score/${dataset}/P3CMQA -g 0 -p ../../${dir}/profile/${dataset}/${target}

echo finish
