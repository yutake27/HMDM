import subprocess
import argparse
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_TMscore(result):
    lines = result.split('\n')
    for line in lines:
        line_split = line.split()
        if len(line_split) == 0:
            continue
        elif line_split[0] == 'TM-score':
            tmscore = float(line_split[2])
        elif line_split[0] == 'GDT-TS-score=':
            gdtts = line_split[1]
        elif line_split[0] == 'GDT-HA-score=':
            gdtha = line_split[1]
    return tmscore, gdtts, gdtha


def run_TMscore(native_pdb, model_pdb):
    cmd = ['TMscore', model_pdb, native_pdb, '-outfmt', '-1']
    result = subprocess.check_output(cmd)
    return result.decode('utf-8')


def get_gdt(native_pdb, model_pdb):
    result = run_TMscore(native_pdb, model_pdb)
    tmscore, gdtts, gdtha = parse_TMscore(result)
    return tmscore, gdtts, gdtha


def get_gdt_for_target(native_pdb_path, model_pdb_dir, blast_xml_csv_path, out_gdt_path):
    model_array = []
    tmscore_array = []
    gdtts_array = []
    gdtha_array = []

    for model in model_pdb_dir.iterdir():
        model_array.append(model.stem)
        tmscore, gdtts, gdtha = get_gdt(native_pdb_path, model)
        tmscore_array.append(tmscore)
        gdtts_array.append(gdtts)
        gdtha_array.append(gdtha)

    df = pd.DataFrame({'TMscore': tmscore_array, 'GDT_TS': gdtts_array, 'GDT_HA': gdtha_array}, index=model_array)
    df = df.astype('float')
    df = df.sort_index()
    df['target'] = [index.rsplit('_', 4)[0] for index in df.index]
    df['template'] = [index.split('_', 2)[2].rsplit('_', 1)[0] for index in df.index]
    df_template = pd.read_csv(blast_xml_csv_path, index_col=0)
    df = pd.merge(df, df_template, left_on='template', right_index=True, how='left')
    df.to_csv(out_gdt_path)


def get_gdt_for_target_df(native_pdb_path, model_pdb_dir, blast_xml_csv_path) -> pd.DataFrame:
    model_array = []
    tmscore_array = []
    gdtts_array = []
    gdtha_array = []

    for model in model_pdb_dir.iterdir():
        model_array.append(model.stem)
        tmscore, gdtts, gdtha = get_gdt(native_pdb_path, model)
        tmscore_array.append(tmscore)
        gdtts_array.append(gdtts)
        gdtha_array.append(gdtha)

    df = pd.DataFrame({'TMscore': tmscore_array, 'GDT_TS': gdtts_array, 'GDT_HA': gdtha_array}, index=model_array)
    df = df.astype('float')
    df = df.sort_index()
    df['target'] = [index.rsplit('_', 4)[0] for index in df.index]
    df['template'] = [index.split('_', 2)[2].rsplit('_', 1)[0] for index in df.index]
    df_template = pd.read_csv(blast_xml_csv_path, index_col=0)
    df = pd.merge(df, df_template, left_on='template', right_index=True, how='left')
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('target', type=str, help='target name')
    parser.add_argument('--blastdb', '-b', type=str, help='blastdb name', default='pdbaa_20200712')
    parser.add_argument('--dataset_name', '-d', type=str, help='name of the dataset', default='target_10')
    args = parser.parse_args()
    native_pdb = (Path('../native_pdb') / args.dataset_name / args.target).with_suffix('.pdb')
    model_pdb_dir = Path('../pdb')/args.dataset_name/args.target
    df_template_path = (Path('../blast_xml') / args.blastdb / args.dataset_name / args.target).with_suffix('.csv')
    out_dir = Path('../tmscore') / args.dataset_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = (out_dir/args.target).with_suffix('.csv')
    get_gdt_for_target(native_pdb, model_pdb_dir, df_template_path, out_path)
