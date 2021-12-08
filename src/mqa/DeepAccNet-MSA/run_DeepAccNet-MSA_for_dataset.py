from pathlib import Path
import subprocess
import argparse


def predict(dataset_name: str, target_name: str, qsub: bool):
    cmd = ['./DeepAccNet-MSA.sh', str(dataset_name), str(target_name)]
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab','-N', 'DANM_' + target_name] + cmd
        subprocess.run(cmd)
    else:
        print(' '.join(cmd))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_dir', type=str, help='pdb dir of dataset. for example ../../../pdb/target_50')
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    dataset_name = pdb_dir.stem
    out_score_dir = Path('../../../score') / dataset_name / 'DeepAccNet-MSA'
    for target in pdb_dir.glob('*'):
        out_score_path = (out_score_dir / target.stem).with_suffix('.csv')
        if not out_score_path.exists():  
            predict(dataset_name, target.stem, args.qsub)
