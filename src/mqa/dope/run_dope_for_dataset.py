from pathlib import Path
import subprocess
import argparse


def qsub_predict(target_dir_path, qsub):
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'dope_' + target_dir_path.stem, './dope.sh', str(target_dir_path)]
        subprocess.run(cmd)
    else:
        cmd = ['./dope.sh', str(target_dir_path)]
        print(' '.join(cmd))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_dir', type=str, help='pdb dir of dataset. for example ../../../pdb/target_50')
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    out_score_dir = Path('../../../score') / pdb_dir.stem / 'dope'
    for target in pdb_dir.glob('*'):
        out_score_path = (out_score_dir / target.stem).with_suffix('.csv')
        if not out_score_path.exists():
            qsub_predict(target, args.qsub)
