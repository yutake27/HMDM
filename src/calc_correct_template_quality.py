import argparse
import numpy as np
import pandas as pd
from Bio import SearchIO, SeqIO
from pathlib import Path

"""
    calculate identity, positive, and coverage considering missing residues
"""


def calc_coverage_from_pir(pir_file):
    query_record, template_record = SeqIO.parse(pir_file, 'pir')
    query_seq, template_seq = query_record.seq, template_record.seq
    template_not_gap_indices = np.where(np.array(template_seq) != '-')[0]
    aln_start, aln_end = template_not_gap_indices[0], template_not_gap_indices[-1]
    query_aln_seq_array = np.array(query_seq)[aln_start: aln_end + 1]
    aln_len = len(np.where(query_aln_seq_array != '-')[0])
    return aln_len


def calc_correct_template_quality(hit, pir_file: Path, query_seq_len: int):
    assert pir_file.exists()
    aln_annotation = hit[0].aln_annotation['similarity']
    hit_seq = str(hit[0].hit.seq)
    ident_num = hit[0].ident_num
    positive_num = hit[0].pos_num
    start = hit[0].query_start
    end = hit[0].query_end
    coverage = end - start
    hit_seq_end = - (query_seq_len - end) if query_seq_len != end else None
    pir_template_seq = str(list(SeqIO.parse(pir_file, 'pir'))[1].seq)[start: hit_seq_end]
    missing_num = 0
    not_ident_num = 0
    not_pos_num = 0
    assert len(hit_seq) == len(pir_template_seq)
    for i, (xr, pr) in enumerate(zip(hit_seq, pir_template_seq)):
        if xr != pr:
            missing_num += 1
            if aln_annotation[i] == ' ':
                continue
            elif aln_annotation[i] == '+':
                not_pos_num += 1
            else:
                not_pos_num += 1
                not_ident_num += 1
    return missing_num, ident_num - not_ident_num, positive_num - not_pos_num, ident_num, positive_num, coverage


def calc_target(dataset_name: str, target_name: str, xml_path: Path, num_iteration: int = 3):
    qresult_list = list(SearchIO.parse(xml_path, 'blast-xml'))
    pir_dir = Path('../pir') / dataset_name / target_name
    quality_list = []
    for iteration in range(1, num_iteration + 1):
        qresult = qresult_list[iteration - 1]
        query_seq_len = qresult.seq_len
        for hit in qresult:
            pir_name = '_'.join(hit.id.split('|')[1:]) + '_' + str(iteration) + '.pir'
            pir_fil_file = pir_dir / 'filtering' / pir_name
            pir_nofil_file = pir_dir / 'nofiltering' / pir_name
            if pir_fil_file.exists():
                pir_file = pir_fil_file
            elif pir_nofil_file.exists():
                pir_file = pir_nofil_file
            else:
                continue
            missing_num, identity, positive, identity_search, positive_search, coverage_search\
                = calc_correct_template_quality(hit, pir_file, query_seq_len)
            coverage = calc_coverage_from_pir(pir_file)
            quality_list.append([pir_file.stem, missing_num, identity, positive, coverage,
                                 identity_search, positive_search, coverage_search, query_seq_len])
    df = pd.DataFrame(quality_list,
                      columns=['template', 'num_missing_residues', 'identity(-misres)',
                               'positive(-misres)', 'coverage(-misres)', 'identity', 'positive',
                               'coverage', 'target_seq_len']).sort_values('template')
    df['target'] = target_name
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-d', type=str, help='Dataset name. e.g. scop_dataset',
                        default='scop_cl_equal_globular100_identity95')
    args = parser.parse_args()
    dataset = args.dataset
    pdbaa = 'pdbaa_20210407'
    template_csv_dir = Path('../blast_xml') / pdbaa / dataset

    result_df_list = []
    for csv_path in template_csv_dir.glob('*.csv'):
        target_name = csv_path.stem
        xml_path = csv_path.with_suffix('.xml')
        target_df = pd.read_csv(csv_path, index_col=0)
        template_quality_target_df = calc_target(dataset, target_name, xml_path)
        print(template_quality_target_df)
        result_df_list.append(template_quality_target_df)
    result_df = pd.concat(result_df_list).sort_values(['target', 'template']).reset_index(drop=True)
    output_path = Path('../score') / dataset / 'template_quality_correct.csv'
    result_df.to_csv(output_path)


if __name__ == '__main__':
    main()

