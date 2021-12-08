#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=1:00:00
#$ -j y

dataset=${1}
target=${2}
dir=MQA_dataset

source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load cuda/11.0.194 cudnn/7.6
exec_dir=`pwd`

cd ../../../../bin/DeepAccNet
DeepAcc_dir=`pwd`
source venv/bin/activate

input_dir=../../${dir}/pdb/${dataset}/${target}/'sampling'
output_dir=../../${dir}/score/${dataset}/DeepAccNet/${target}
mkdir -p ${output_dir}

# generate bert features
python extractBert.py ${input_dir} ${output_dir} --modelpath ProtBert-BFD
# DeepAccNet (standard)
python DeepAccNet.py -v -lt ${input_dir} ${output_dir}
cd ${exec_dir}
python mv_local_score.py ../../../score/${dataset}/DeepAccNet/${target}/DeepAccNet

# DeepAccNet-MSA
cd ${DeepAcc_dir}
python DeepAccNet.py -v --bert ${input_dir} ${output_dir}
cd ${exec_dir}
python mv_local_score.py ../../../score/${dataset}/DeepAccNet/${target}/DeepAccNet-Bert

python get_global_score.py ../../../score/${dataset}/DeepAccNet/${target}

echo finish
