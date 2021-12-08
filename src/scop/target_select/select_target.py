import pandas as pd
from pathlib import Path
import argparse


def select_target_from_sf(csv_path):
    sf_id = csv_path.stem
    df = pd.read_csv(csv_path, index_col=0)
    sorted_df = df.sort_values(['3', '2', '1'], ascending=False)
    target = pd.DataFrame(sorted_df.iloc[0]).T
    len_sf = len(sorted_df)
    target['SF_ID'] = sf_id
    target['len_SF'] = len_sf
    target = target.rename_axis('SF-DOMID').reset_index()
    return target


def main(scop_cla_csv_path, csv_dir, output_path):
    scop_cla_df = pd.read_csv(scop_cla_csv_path)
    blast_csv_dir = Path(args.blast_xml_dir)
    target_df_list = []
    for csv in blast_csv_dir.glob('*.csv'):
        target = select_target_from_sf(csv)
        target_df_list.append(target)
    target_df = pd.concat(target_df_list)
    concat_df = pd.merge(target_df, scop_cla_df, on='SF-DOMID', how='left')
    concat_df.to_csv(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--scop_cla_csv_path', type=str,
                        help='Path to the csv with the number of SF',
                        default='../../../scop/20210330/scop-cla-latest.csv')
    parser.add_argument('--blast_xml_dir', type=str,
                        help='Path to the directory of the blast_xml',
                        default='../../../blast_xml/pdbaa_20210407/scop_sf_check')
    parser.add_argument('--out_csv_path', type=str,
                        help='output path',
                        default='../../../scop/20210330/globular_cl_equal_target_100.csv')
    args = parser.parse_args()

    main(args.scop_cla_csv_path, args.blast_xml_dir, args.out_csv_path)

