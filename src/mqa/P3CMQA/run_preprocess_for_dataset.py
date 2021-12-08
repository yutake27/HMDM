from pathlib import Path
import subprocess
import argparse


def qsub_predict(target_name, dataset_name, qsub):
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'preprocess_' + target_name, './preprocess.sh', target_name, dataset_name]
        subprocess.run(cmd)
    else:
        cmd = ['./preprocess.sh', target_name, dataset_name]
        print(' '.join(cmd))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_dir', type=str, help='pdb dir of dataset. for example ../../../pdb/target_50')
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    profile_dir = Path('../../../profile') / pdb_dir.stem

    for target in pdb_dir.glob('*'):
        if not ((profile_dir / target.stem / target.stem).with_suffix('.acc20')).exists():
            qsub_predict(target.stem, pdb_dir.stem, args.qsub)
