import argparse
import pandas as pd
import subprocess
from pathlib import Path


def run_psiblast4sf(sf_id: int, qsub: bool = False):
    cmd = ['./psiblast4sf.sh', str(sf_id)]
    if qsub:
        cmd = ['qsub', '-g', 'tga-ishidalab', '-N', 'SF_' + str(sf_id)] + cmd
        subprocess.run(cmd)
    else:
        print(' '.join(cmd))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--qsub', '-q', action='store_true', help='qsub mode')
    parser.add_argument('--sf_csv_path', type=str,
                        help='Path to the csv with the number of SF',
                        default='../../../scop/20210330/SF_globular_list.csv')
    parser.add_argument('--blast_xml_dir', type=str,
                        help='Path to the directory of the blast_xml',
                        default='../../../blast_xml/pdbaa_20210407/scop_sf_check')
    args = parser.parse_args()

    df = pd.read_csv(args.sf_csv_path)
    df = df.query('CL != "1000004"')
    sf_list = df.reset_index().groupby('CL').apply(lambda x: x.sort_values('len_SF', ascending=False)[: 25]).reset_index(drop=True)
    print(sf_list.value_counts('CL'))
    blast_xml_dir = Path(args.blast_xml_dir)

    for sf in sf_list['SF']:
        blast_csv_path = (blast_xml_dir / str(sf)).with_suffix('.csv')
        if not blast_csv_path.exists():
            run_psiblast4sf(sf, args.qsub)
