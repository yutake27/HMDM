import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


def get_HitNum_from_xml(xml_path: str):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    hit_num_list = [None] * 3
    for i, iteration in enumerate(root[8]):
        hit_num = len(iteration[4])
        hit_num_list[i] = hit_num
    return hit_num_list


def main(blast_xml_dir, out_path):
    blast_xml_dir = Path(args.blast_xml_dir)
    entry_list = []
    hit_num_list = []
    for xml in blast_xml_dir.glob('*.xml'):
        print(xml.stem)
        hit_num = get_HitNum_from_xml(xml)
        entry_list.append(xml.stem)
        hit_num_list.append(hit_num)
    df = pd.DataFrame(hit_num_list, columns=[1, 2, 3], index=entry_list)
    df = df.sort_values([3, 2, 1], ascending=False)
    df.to_csv(out_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast_xml_dir', type=str,
                        help='Path to the directory of the blast_xml',
                        default='../../../pisces/20210225/multi-domain_blast')
    parser.add_argument('--out_csv_path', type=str,
                        help='output path',
                        default='../../../pisces/20210225/multidomain_hit_num_cullpdb_pc20_res2.0_R0.25.csv')
    args = parser.parse_args()

    main(args.blast_xml_dir, args.out_csv_path)

