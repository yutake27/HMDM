import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def get_global_score(npz_path: str):
    x = np.load(npz_path)
    global_score = np.mean(x['lddt'])
    return global_score


def get_df(dir_path: Path) -> pd.DataFrame:
    model_name_list = []
    global_score_list = []


    for npz_path in dir_path.glob('*.npz'):
        global_score = get_global_score(str(npz_path))
        model_name_list.append(npz_path.stem)
        global_score_list.append(global_score)
    method_name = dir_path.stem
    df = pd.DataFrame({'Score': global_score_list}, index=model_name_list)
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('score_dir', type=str, 
                        help='directory of local score. ex) ../../../score/dataset/method/target')
    
    args = parser.parse_args()

    score_dir = Path(args.score_dir)
    msa_score_dir = score_dir
    df = get_df(msa_score_dir).rename({'Score': 'DeepAccNet-MSA'}, axis=1)
    output_path = score_dir.with_suffix('.csv')
    df.to_csv(output_path)
