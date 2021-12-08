import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd


def pkl2df(pkl_path: str) -> pd.DataFrame:
    model_list = []
    global_score_list = []
    with open(pkl_path, 'rb') as f:
        score_dict = pickle.load(f)
    for model, values in score_dict.items():
        global_score = values['global']
        model_name = Path(model).stem
        model_list.append(model_name)
        global_score_list.append(global_score)
    df = pd.DataFrame({'Score': global_score_list}, index=model_list)
    return df



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gdtts_model_score_pkl_path', type=str, 
                        help='Path to the pickle file that contains MQA scores (GDTTS model).')
    parser.add_argument('gdtts_ranking_model_score_pkl_path', type=str,
                        help='Path to the pickle file that contains MQA scores (GDTTS Ranking model)')
    args = parser.parse_args()

    gdtts_pkl_path = Path(args.gdtts_model_score_pkl_path)
    gdtts_model_df= pkl2df(gdtts_pkl_path).rename({'Score': 'ResNetQA'}, axis=1)
    gdtts_ranking_model_df = pkl2df(args.gdtts_ranking_model_score_pkl_path).rename({'Score': 'ResNetQA-R'}, axis=1)
    df = pd.merge(gdtts_model_df, gdtts_ranking_model_df, left_index=True, right_index=True, how='inner')
    output_path = gdtts_pkl_path.parent.with_suffix('.csv')
    df.to_csv(output_path)
