import os
import subprocess
from pathlib import Path

import pandas as pd


def run_cdhit(input_fasta_path: str,
              output_fasta_path: str,
              threshold: float = 0.95) -> None:
    cmd = ['cd-hit', '-i', input_fasta_path,
           '-o', output_fasta_path, '-c', str(threshold)]
    subprocess.run(cmd)


def parse_clstr_line(line: str) -> (str, int, int):
    hit_info = line.split()[2][1: -3]
    hit_id, iteration, hit_index = hit_info.split('-')
    return hit_id, int(iteration), int(hit_index)


def parse_clstr_file(clstr_file_path: str) -> list:
    entry_list = []
    with open(clstr_file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == '>':
                clstr_id = int(line.split()[1])
            else:
                hit_id, iteration, hit_index = parse_clstr_line(line)
                entry_list.append([hit_id, iteration, hit_index, clstr_id])
    return entry_list


def make_df_from_clstr_file(clstr_file_path: str) -> pd.DataFrame:
    entry_list = parse_clstr_file(clstr_file_path)
    df = pd.DataFrame(entry_list,
                      columns=['hit_id', 'iteration', 'hit_index', 'clstr_id'])
    sort_df = df.sort_values(['iteration', 'hit_index'])
    return sort_df


def clstr_blast_hit(fasta_path: str) -> pd.DataFrame:
    output_fasta_path = fasta_path[: -6] + '_cdhit.fasta'
    output_clstr_path = output_fasta_path + '.clstr'
    if not Path(output_clstr_path).exists():
        run_cdhit(fasta_path, output_fasta_path)
        os.remove(output_fasta_path)
    df = make_df_from_clstr_file(output_clstr_path)
    return df
