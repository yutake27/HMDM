import argparse
import tarfile
from pathlib import Path

import tqdm
import pandas as pd

route_path = Path('..')

"""
# Compress necessary files of the dataset
Source
* pdb   / dataset_name / target / sampling
* fasta / dataset_name
* native_pdb / dataset_name
* score


Output
* README.md
    write about dataset details
* pdb
* native_pdb
* fasta
* score
"""


def compress_pdb(dataset_name: str, output_dir: Path) -> Path:
    print('Compress pdb')
    output_tar_pdb_path = output_dir / 'pdb.tar.gz'
    if output_tar_pdb_path.exists():
        return output_tar_pdb_path
    pdb_dir = route_path / 'pdb' / dataset_name
    with tarfile.open(output_tar_pdb_path, 'w:gz') as f:
        for target in tqdm.tqdm(pdb_dir.glob('*'), total=100):
            for pdb in (target / 'sampling').glob('*.pdb'):
                f.add(str(pdb), 'pdb/' + target.stem + '/' + pdb.name)
    return output_tar_pdb_path


def compress_native_pdb(tf: tarfile.TarFile, dataset_name: str, arcdir: str) -> None:
    print('Compress native pdb')
    native_pdb_dir = route_path / 'native_pdb' / dataset_name
    for native_pdb in native_pdb_dir.glob('*.pdb'):
        pdb_dir = route_path / 'pdb' / dataset_name / native_pdb.stem
        if pdb_dir.exists():
            tf.add(str(native_pdb), arcname=arcdir + '/native_pdb/' + native_pdb.name)


def compress_fasta(tf: tarfile.TarFile, dataset_name: str, arcdir: str) -> None:
    print('Compress fasta')
    fasta_dir = route_path / 'fasta' / dataset_name
    for fasta in fasta_dir.glob('*.fasta'):
        pdb_dir = route_path / 'pdb' / dataset_name / fasta.stem
        if pdb_dir.exists():
            tf.add(str(fasta), arcname=arcdir + '/fasta/' + fasta.name)


rename_columns_dict = {'model': 'Model', 'target': 'Target', 'template': 'Template', 'seq_len': 'SeqLength'}

label_columns = [
    'Model', 'Target', 'Template', 'GDT_TS', 'GDT_HA',
    'SeqLength', 'identity', 'positive', 'coverage',
    'identity(-misres)', 'positive(-misres)', 'coverage(-misres)', 'num_misres'
]

score_columns = [
    'Model', 'Target', 'identity(%)', 'positive(%)', 'coverage(%)',
    'identity(-misres)(%)', 'positive(-misres)(%)', 'coverage(-misres)(%)',
    'DOPE', 'SOAP', 'SBROD', 'ProQ2D', 'ProQRosCenD', 'ProQRosFAD', 'ProQ3D',
    'P3CMQA', 'DeepAccNet', 'DeepAccNet-Bert'
]


def make_scop_score(dataset_name: str, output_dir: Path) -> (str, str, str):
    """load scop final score and split it into target, label, and mqa score.

    Args:
        dataset_name (str): Created dataset name
        output_dir (Path): Output directory path
    Return:
        (str): path to target.csv
        (str): path to label.csv
        (str): path to score.csv
    """
    csv_path = route_path / 'score' / dataset_name / (dataset_name + '_final_all_score.csv')
    df = pd.read_csv(csv_path, index_col=0)
    output_score_dir = output_dir / 'score'
    output_score_dir.mkdir(exist_ok=True)
    # Rename columns
    df = df.rename(rename_columns_dict, axis=1)
    # Drop columns
    label_df = df[label_columns]
    label_output_path = output_score_dir / 'label.csv'
    label_df.to_csv(label_output_path)
    if dataset_name[: 4] == 'scop':
        target_df = df[[
            'Target', 'SeqLength', 'SF-DOMID', 'SF', 'len_SF',
            'FA-DOMID', 'FA-PDBID', 'FA-PDBREG', 'FA-UNIID', 'FA-UNIREG', 'SF-PDBID',
            'SF-PDBREG', 'SF-UNIID', 'SF-UNIREG', 'TP', 'CL', 'CF', 'FA', 'Class'
        ]]
    elif dataset_name[: 6] == 'pisces':
        target_df = df[[
            'Target', 'SeqLength', 'IDs', 'Exptl.', 'resolution', 'R-factor',
            'FreeRvalue', 'PDB_ID', 'Chain', 'Domain_num'
        ]]
        target_df = target_df.rename({'Domain_num': 'DomainNum'}, axis=1)
    else:
        raise ValueError()
    target_df = target_df.groupby('Target').head(1).reset_index(drop=True)
    target_output_path = output_score_dir / 'target.csv'
    target_df.to_csv(target_output_path)
    score_df = df[score_columns]
    score_output_path = output_score_dir / 'mqa_score.csv'
    score_df.to_csv(score_output_path)
    return target_output_path, label_output_path, score_output_path


def compress_score(tf: tarfile.TarFile, dataset_name: str, arcdir: str, output_dir: Path) -> None:
    print('Compress score file')
    target_path, label_path, score_path = make_scop_score(dataset_name, output_dir)
    arcdir_prefix = arcdir + '/data/'
    tf.add(target_path, arcdir_prefix + 'target.csv')
    tf.add(label_path, arcdir_prefix + 'label.csv')
    tf.add(score_path, arcdir_prefix + 'mqa_score.csv')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_pdb_path', type=str)
    parser.add_argument('output_dataset_name', type=str, help='[Single-domain, Multi-domain]')
    args = parser.parse_args()
    dataset_name = Path(args.dataset_pdb_path).stem
    assert dataset_name[:4] == 'scop' or dataset_name[:6] == 'pisces'
    output_dataset_name = args.output_dataset_name

    output_dir = route_path / 'dataset' / output_dataset_name
    output_dir.mkdir(parents=True, exist_ok=True)

    # Compress pdb
    output_tar_pdb_path = compress_pdb(dataset_name, output_dir)

    # Compress other files
    output_tar_path = (output_dir / output_dataset_name).with_suffix('.tar.gz')
    with tarfile.open(output_tar_path, 'w:gz') as tf:
        arcdir: str = output_dataset_name
        # pdb
        tf.add(output_tar_pdb_path, arcname=arcdir + '/pdb.tar.gz')
        # native pdb
        compress_native_pdb(tf, dataset_name, arcdir)
        # fasta
        compress_fasta(tf, dataset_name, arcdir)
        # score
        compress_score(tf, dataset_name, arcdir, output_dir)
        # README
        dataset_dir = route_path / 'dataset'
        if dataset_name[:4] == 'scop':
            readme_file = 'README-single.md'
        else:  # dataset_name[: 6] == 'pisces'
            readme_file = 'README-multi.md'
        readme_path = dataset_dir / readme_file
        tf.add(readme_path, arcname=arcdir + '/README.md')


if __name__ == '__main__':
    main()

