#!/bin/bash
#$ -cwd
#$ -l q_node=1
#$ -l h_rt=24:00:00

. /etc/profile.d/modules.sh
module load cuda/10.1.105 cudnn/7.6 openmpi r

source /home/4/16B09097/.bashrc
conda activate keras
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"

target=$1
datadir=$2
fasta=`readlink -f ../../../fasta/$datadir/$target.fasta`
pdbdir=`readlink -f ../../../pdb/$datadir/$target/sampling`
outpath=`readlink -f ../../../score/$datadir/proq3/${target}`

if [ -e ${outpath}.csv ]; then
    echo proq3 for $target have already been completed!
    exit
fi

mkdir $outpath
txtpath=$outpath/$target.txt

if [ -e $txtpath ]; then
    rm $txtpath
fi

for i in $pdbdir/*;do
    echo $i >> $txtpath
done

run_proq3.sh -fasta $fasta -l $txtpath -outpath $outpath -ncores 7

python make_proq3_df_for_a_target.py $outpath

cd ${outpath}/..
tar -cvzf ${target}.tar.gz ${target}
rm -rf ${target}
