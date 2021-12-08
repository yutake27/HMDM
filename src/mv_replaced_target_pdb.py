import argparse
import os
import pandas as pd
from pathlib import Path
import shutil

"""
mv targets that does not meet the criteria
"""



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('final_tmscore_csv_path', type=str, help='Path to the csv of final sampled tmscore')
    args = parser.parse_args()

    tmscore_csv_path = Path(args.final_tmscore_csv_path)
    dataset_name = tmscore_csv_path.parent.stem
    from_pdb_dir = Path('../pdb') / dataset_name
    to_pdb_dir = Path('../pdb') / (dataset_name + '_replaced')
    to_pdb_dir.mkdir(exist_ok=True)

    df = pd.read_csv(tmscore_csv_path)
    target_list = df['target'].unique()

    for target in from_pdb_dir.glob('*'):
        if target.stem not in target_list:
            print(target, to_pdb_dir)
            to_pdb_target_dir = to_pdb_dir / target.stem
            if to_pdb_target_dir.exists():
                for d in target.glob('*'):
                    shutil.move(str(d), str(to_pdb_target_dir))
                os.rmdir(str(target))
            else:
                shutil.move(str(target), str(to_pdb_dir))


if __name__ == '__main__':
    main()

