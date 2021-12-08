import argparse
import subprocess
from pathlib import Path

import pandas as pd


def scop_modeling(sf_domid: int, dataset_name: str, qsub: bool):
    cmd = ['./scop_modeling.sh', str(sf_domid), dataset_name]
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'modeling_' + str(sf_domid)] + cmd
        subprocess.run(cmd)
    else:
        print(' '.join(cmd))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    args = parser.parse_args()

    target_csv_path = '../scop/20210330/globular_cl_equal_target_100_minnum_50_final.csv'
    dataset_name = 'scop_cl_equal_globular100_identity95_coverage60'
    df = pd.read_csv(target_csv_path).sort_values('len_SF', ascending=False)

    out_gdt_dir = Path('../tmscore') / dataset_name
    for index, row in df.iterrows():
        sf_domid = row['SF-DOMID']
        pdb_id = row['SF-PDBID']
        chain = row['SF-PDBREG'].split(':')[0]
        out_gdt_path = (out_gdt_dir / (pdb_id + '_' + chain)).with_suffix('.csv')
        if not out_gdt_path.exists():
            scop_modeling(sf_domid, dataset_name, args.qsub)


if __name__ == '__main__':
    main()
