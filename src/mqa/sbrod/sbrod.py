import argparse
import subprocess
import pandas as pd
from pathlib import Path


def run_sbrod(model_dir):
    output = subprocess.check_output(['sbrod', model_dir])
    return output.decode()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dir', type=str, help='target pdb dir of a dataset. ex) ../../../pdb/target_10/1A3Q')
    args = parser.parse_args()

    target_dir = Path(args.target_dir)
    dataset_name = target_dir.parent.stem
    output_csv_path = (Path('../../../score') / dataset_name / 'sbrod' /
                       f'{target_dir.stem}').with_suffix('.csv')

    model_name_array = []
    sbrod_score_array = []
    output = run_sbrod(target_dir / 'sampling' / '*.pdb')
    output_list = output.split()
    df = pd.DataFrame({'sbrod': output_list[1::2]}, index=output_list[::2])
    df.index = list(map(lambda x: Path(x).stem, df.index))
    df.to_csv(output_csv_path)
