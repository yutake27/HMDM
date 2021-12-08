from pathlib import Path
import subprocess
import argparse


def qsub_predict(target_name, dataset_name):
    cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'predict_' + target_name,'./predict_multi_class.sh', target_name, dataset_name]
    subprocess.run(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_dir', type=str, help='pdb dir of dataset. for example ../../../pdb/target_50')
    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)
    out_score_dir = Path('../../../score') / pdb_dir.stem / 'P3CMQA_multi'

    for target in pdb_dir.glob('*'):
        out_score_path = (out_score_dir / target.stem).with_suffix('.csv')
        if not out_score_path.exists():
            qsub_predict(target.stem, pdb_dir.stem)
