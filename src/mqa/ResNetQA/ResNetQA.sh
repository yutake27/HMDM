#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=1:30:00
#$ -j y

dataset=${1}
target=${2}

source /home/4/16B09097/.bashrc
. /etc/profile.d/modules.sh
module load cuda cudnn/7.6
exec_dir=`pwd`

input_dir=`readlink -f ../../../pdb/${dataset}/${target}/sampling`
output_dir_rel=../../../score/${dataset}/ResNetQA/${target}
mkdir -p ${output_dir_rel}
output_dir=`readlink -f ${output_dir_rel}`
fasta_path=`readlink -f ../../../fasta/${dataset}/${target}.fasta`

bin_path=`readlink -f ../../../../bin`

features_dir=${output_dir}/${target}_OUT
seq_features=${features_dir}/${target}_contact/feat_${target}_uce3/${target}.inputFeatures.pkl
dist_features=${features_dir}/DistancePred/${target}.pairPotential.DFIRE16.pkl
output_pkl_path=${output_dir}/${target}.QA.pkl
quality_type=GDTTS


# generate features by RaptorX-3DModeling
if [ ! -e ${seq_features} ] || [ ! -e ${dist_features} ]; then
    cd ${bin_path}/RaptorX-3DModeling
    conda activate RaptorX
    export CUDA_ROOT=$CUDA_HOME
    ./Server/RaptorXFolder.sh -n 0 -o ${output_dir} ${fasta_path}
fi

# prediction by ResNetQA
cd ${bin_path}/ResNetQA/main
conda activate mypython

features_dir=${output_dir}/${target}_OUT
seq_features=${features_dir}/${target}_contact/feat_${target}_uce3/${target}.inputFeatures.pkl
dist_features=${features_dir}/DistancePred/${target}.pairPotential.DFIRE16.pkl

# using GDT_TS model
gdtts_model_output_pkl_path=${output_dir}/${target}.QA.GDTTS.pkl
quality_type=GDTTS
python ResNetQA.py ${seq_features} ${dist_features} ${input_dir}/ ${gdtts_model_output_pkl_path} ${quality_type} -device_id 0

# using GDT_TS rank model
gdtts_ranking_model_output_pkl_path=${output_dir}/${target}.QA.GDTTS_Ranking.pkl
quality_type=GDTTS_Ranking
python ResNetQA.py ${seq_features} ${dist_features} ${input_dir}/ ${gdtts_ranking_model_output_pkl_path} ${quality_type} -device_id 0


# compress features
cd ${output_dir}
tar -czf ${target}_feature.tar.gz ${target}_OUT
rm -rf ${features_dir}

# convert pkl to csv
cd ${exec_dir}
python get_global_score.py ${gdtts_model_output_pkl_path} ${gdtts_ranking_model_output_pkl_path}

echo finish
