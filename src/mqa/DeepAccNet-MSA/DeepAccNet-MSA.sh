#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=1:30:00
#$ -j y

dataset=${1}
target=${2}

source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load cuda/11.0.194 cudnn/7.6
conda activate mypython
exec_dir=`pwd`

# path
MQA_dataset_dir=`readlink -f ../../..`
tool_dir=`readlink -f ../../../../bin`
uniclust_dir=`readlink -f ../../../../../hhsuitedb/uniclust`
db_name=UniRef30_2020_06
uniclust_db_dir=${uniclust_dir}/${db_name}
uniclust_db_tar=${db_name}_hhsuite.tar.gz

DeepAccNet_dir=${tool_dir}/DeepAccNet
DeepAccNet_MSA_dir=${tool_dir}/DeepAccNet-MSA
trRosetta_dir=${tool_dir}/trRosetta

fasta=${MQA_dataset_dir}/fasta/${dataset}/${target}.fasta
input_dir=${MQA_dataset_dir}/pdb/${dataset}/${target}/'sampling'
output_dir=${MQA_dataset_dir}/score/${dataset}/DeepAccNet-MSA/${target}
mkdir -p ${output_dir}

profile_dir=${MQA_dataset_dir}/profile/hhblits/${db_name}/${dataset}/${target}
mkdir -p ${profile_dir}
output_a3m=${profile_dir}/${target}.a3m
output_hhr=${profile_dir}/${target}.hhr
output_hhm=${profile_dir}/${target}.hhm

feature_dir=${output_dir}/feature
mkdir ${feature_dir}
output_dist_npz=${feature_dir}/${target}_distPred.npz

# runtime envinronment
source ${DeepAccNet_dir}/venv/bin/activate

# generate MSA by hhblits
core=14
if [ ! -e ${output_a3m} ];then
    cp_cmd="cp ${uniclust_db_dir}/${uniclust_db_tar} ${TMPDIR}"
    echo ${cp_cmd}
    eval ${cp_cmd}
    cd ${TMPDIR}
    time tar -xvf ${uniclust_db_tar} --use-compress-prog=pigz
    uniclust_db=${TMPDIR}/${db_name}
    hhblits -i ${fasta} -d ${uniclust_db} -o ${output_hhr} -ohhm ${output_hhm} -oa3m ${output_a3m} -cpu ${core}
fi

# predict ContactMap by trRosetta
cd ${trRosetta_dir}
if [ ! -e ${output_dist_npz} ];then
    python network/predict.py -m model2019_07 ${output_a3m} ${output_dist_npz}
fi

# run DeepAccNet-MSA
cd ${DeepAccNet_MSA_dir}
python ErrorPredictorMSA.py --verbose ${output_dist_npz} ${input_dir} ${output_dir}
cd ${exec_dir}
python get_global_score.py ${output_dir} 
echo finish
