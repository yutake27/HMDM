from pathlib import Path
import subprocess
import argparse


def qsub_predict(target_name, dataset_name, qsub=False):
    if not qsub:
        cmd = ['./predict.sh', target_name, dataset_name]
        print(' '.join(cmd))
    else:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'predict_' + target_name, './predict.sh', target_name, dataset_name]
        subprocess.run(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_dir', type=str, help='pdb dir of dataset. for example ../../../pdb/target_50')
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    out_score_dir = Path('../../../score') / pdb_dir.stem / 'P3CMQA'

    for target in pdb_dir.glob('*'):
        out_score_path = (out_score_dir / target.stem).with_suffix('.csv')
        if not out_score_path.exists():
            qsub_predict(target.stem, pdb_dir.stem, args.qsub)
