import pandas as pd
from pathlib import Path
import argparse
import numpy as np


"""
Replace the targets whose modeling results do not meet the criteria.
When replacing a target, select a different target from the same superfamily.
"""


def select_target_from_sf(csv_path, replace_target_sfdomid):
    sf_id = csv_path.stem
    df = pd.read_csv(csv_path, index_col=0)
    sorted_df = df.sort_values(['3', '2', '1'], ascending=False)
    for index, row in sorted_df.iterrows():
        if index not in replace_target_sfdomid:
            target = pd.DataFrame(sorted_df.loc[index]).T
            break

    len_sf = len(sorted_df)
    target['SF_ID'] = sf_id
    target['len_SF'] = len_sf
    target = target.rename_axis('SF-DOMID').reset_index()
    return target


def select_replacing_targets(tmscore_df, less_threshold: int = 50, max_gdtts_threshold=0.7):
    group = tmscore_df.groupby('target')
    less_target = group.filter(lambda x: len(x) < less_threshold)
    less_target['Reason'] = 'less'
    low_target = group.filter(lambda x: x['GDT_TS'].max() < max_gdtts_threshold)
    low_target['Reason'] = 'low'
    replace_target_df = pd.concat([less_target, low_target]).groupby(['target', 'SF-DOMID', 'Reason'], as_index=False).agg({'GDT_TS': 'max', 'model': 'count'})
    replace_target_df = replace_target_df.drop_duplicates(subset='target')
    return replace_target_df


def main(sampling_score_csv_path, scop_cla_csv_path, sf_csv_dir, output_path, save):
    tmscore_df = pd.read_csv(sampling_score_csv_path)
    replace_target_df = select_replacing_targets(tmscore_df)
    replace_target_sfdomid = replace_target_df['SF-DOMID'].unique()
    scop_cla_df = pd.read_csv(scop_cla_csv_path)
    tmscore_cla_df = pd.merge(tmscore_df, scop_cla_df, on='SF-DOMID', how='left')
    replace_target_info_df = pd.merge(replace_target_df, scop_cla_df, on='SF-DOMID', how='left')
    print('Targets to be replaced (including targets that have already been replaced)')
    print(replace_target_info_df.sort_values('SF')[
          ['target', 'SF', 'SF-DOMID', 'Reason', 'GDT_TS', 'model']].reset_index(), '\n')
    for name, group in replace_target_info_df.groupby('SF'):
        tmscore_target_list = tmscore_cla_df.query('SF == @name')['target'].unique()
        replace_target_list = group['target'].unique()
        new_target_list = set(tmscore_target_list) - set(replace_target_list)
        if len(new_target_list) == 0:
            print('SF to replace the target:', name)
            print(group[['target', 'SF', 'SF-DOMID', 'Reason', 'GDT_TS', 'model']].reset_index(drop=True), '\n')

    blast_csv_dir = Path(sf_csv_dir)
    target_df_list = []
    for csv in blast_csv_dir.glob('*.csv'):
        target = select_target_from_sf(csv, replace_target_sfdomid)
        target_df_list.append(target)
    target_df = pd.concat(target_df_list)
    concat_df = pd.merge(target_df, scop_cla_df, on='SF-DOMID', how='left')
    if save:
        concat_df.to_csv(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampling_score_csv_path', type=str,
                        help='Path to the csv with the tmscore of sampled model.',
                        default='../../../score/scop_cl_equal_globular100_identity95/'
                                + 'scop_cl_equal_globular100_identity95_sampling.csv')
    parser.add_argument('--scop_cla_csv_path', type=str,
                        help='Path to the csv with the number of SF',
                        default='../../../scop/20210330/scop-cla-latest.csv')
    parser.add_argument('--blast_xml_dir', type=str,
                        help='Path to the directory of the blast_xml',
                        default='../../../blast_xml/pdbaa_20210407/scop_sf_check')
    parser.add_argument('--out_csv_path', type=str,
                        help='output path',
                        default='../../../scop/20210330/globular_cl_equal_target_100_replace.csv')
    parser.add_argument('--save', '-s', action='store_true', help='save target list')
    args = parser.parse_args()

    main(args.sampling_score_csv_path, args.scop_cla_csv_path, args.blast_xml_dir, args.out_csv_path, args.save)
