import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sampling_csv', type=str,
                        help='csv of models after sampling. e.g. ../../../score/dataset/dataset_sampling.csv')
    parser.add_argument('target_list_csv', type=str,
                        help='csv of final target list. e.g. ../../../scop/date/target_list_final.csv')
    args = parser.parse_args()

    # load tmscore
    tmscore_df = pd.read_csv(args.sampling_csv)
    # load target list
    target_df = pd.read_csv(args.target_list_csv, index_col=0)
    class_dict = {'1000000': 'All alpha', '1000001': 'All beta', '1000002': 'alpha / beta', '1000003': 'alpha + beta'}
    target_df['Class'] = target_df['CL'].apply(lambda x: class_dict[str(x)])

    # merge dataframe
    df = pd.merge(tmscore_df, target_df, on='SF-DOMID', how='inner')
    sampling_csv_path = Path(args.sampling_csv)
    dataset_name = sampling_csv_path.parent.stem
    output_path = sampling_csv_path.with_name(f'{dataset_name}_final.csv')
    df.to_csv(output_path)


if __name__ == '__main__':
    main()
