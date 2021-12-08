import argparse
from pathlib import Path

import pandas as pd


def get_final_target_list(df: pd.DataFrame,
                          max_gdtts_min_threshold: float = 0.7,
                          gdtts_min_threshold: float = 0.4,
                          num_model_min_threshold: int = 50,
                          num_target: int = 100):
    fil_df = df.query('GDT_TS >= @gdtts_min_threshold')
    sel_df = fil_df.groupby('target').filter(
        lambda x: x['GDT_TS'].max() >= max_gdtts_min_threshold
        and len(x) > num_model_min_threshold)
    target_sorted_by_hitnum = sel_df.groupby('target').head(1).sort_values(['3', '2', '1'], ascending=False)
    final_target_list = list(target_sorted_by_hitnum[: num_target]['target'])
    return final_target_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sampling_csv', type=str,
                        help='csv of sampling models. e.g. ../../../score/dataset/dataset_sampling.csv')
    parser.add_argument('--output_target_list_path', type=str,
                        help='output path of final targets. ',
                        default='../../../pisces/20210225/' +
                        'multidomain_target_100_coverage60_gdtts40_num_model_50_cullpdb_pc20_res2.0_R0.25.csv')
    parser.add_argument('--pisces_multidomain_csv', type=str,
                        default='../../../pisces/20210225/multidomain_cullpdb_pc20_res2.0_R0.25.csv',
                        help='multi-domain target list of pisces.')
    parser.add_argument('--pisces_multidomain_hitnum_csv', type=str,
                        default='../../../pisces/20210225/multidomain_hit_num_cullpdb_pc20_res2.0_R0.25.csv',
                        help='hit num csv of multi-domain target of pisces')
    args = parser.parse_args()

    # load tmscore
    tmscore_df = pd.read_csv(args.sampling_csv)
    # load target info
    pisces_multidomain_df = pd.read_csv(args.pisces_multidomain_csv, index_col=0)
    pisces_multidomain_df['target'] \
        = [row['PDB_ID'] + '_' + row['Chain'] for _, row in pisces_multidomain_df.iterrows()]
    pisces_hitnum_df = pd.read_csv(args.pisces_multidomain_hitnum_csv, index_col=0)
    target_df = pd.merge(pisces_multidomain_df, pisces_hitnum_df,
                         left_on='target', right_index=True, how='left')

    # merge dataframe
    df = pd.merge(tmscore_df, target_df, on='target', how='left')

    # get final dataset
    final_target_list = get_final_target_list(df)
    final_df = df.query('target in @final_target_list')
    print(final_df.groupby('target').agg({'GDT_TS': max, 'model': len}).sort_values('GDT_TS')[: 20])
    print(final_df.groupby('target').agg({'GDT_TS': max, 'model': len}).sort_values('model')[: 20])

    # save final target list
    final_target_df = target_df.query('target in @final_target_list').sort_values('target').reset_index(drop=True)
    print(final_target_df)
    final_target_df.to_csv(args.output_target_list_path)

    sampling_csv_path = Path(args.sampling_csv)
    dataset_name = sampling_csv_path.parent.stem
    output_path = sampling_csv_path.with_name(f'{dataset_name}_final.csv')
    final_df.to_csv(output_path)


if __name__ == '__main__':
    main()
