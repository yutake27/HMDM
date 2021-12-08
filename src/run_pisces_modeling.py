import argparse
import subprocess
from pathlib import Path

import pandas as pd


def run_modeling(target: int, dataset_name: str, qsub: bool):
    cmd = ['./pisces_modeling.sh', str(target), dataset_name]
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'modeling_' + str(target)] + cmd
        subprocess.run(cmd)
    else:
        print(' '.join(cmd))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    args = parser.parse_args()

    target_csv_path = '../pisces/20210225/multidomain_hit_num_cullpdb_pc20_res2.0_R0.25.csv'
    dataset_name = 'pisces_multidomain100_identity95_coverage60'
    target_num = 150
    df = pd.read_csv(target_csv_path, index_col=0)
    df = df.sort_values(['3', '2', '1'], ascending=False)[: target_num]

    out_gdt_dir = Path('../tmscore') / dataset_name
    for index, row in df.iterrows():
        target = index
        out_gdt_path = (out_gdt_dir / target).with_suffix('.csv')
        if not out_gdt_path.exists():
            run_modeling(target, dataset_name, args.qsub)


if __name__ == '__main__':
    main()

