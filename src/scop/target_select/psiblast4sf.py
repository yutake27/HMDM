import argparse
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def psiblast_xml(fasta_path, out_xml_path, blastdb, iteration='3', max_sequence=10000):
    cmd = ['psiblast', '-query', fasta_path, '-db', blastdb, '-out', out_xml_path,
           '-num_iterations', iteration, '-outfmt', '5', '-evalue', '1.0e-3',
           '-max_target_seqs', str(max_sequence)]
    subprocess.run(cmd)


def getSeq(sf_domid, scop_sf_fasta):
    with open(scop_sf_fasta, 'r') as f:
        lines = f.readlines()
    index = [i for i, line in enumerate(lines) if line[1: 8] == sf_domid][0]
    seq = lines[index + 1][: -1]
    return seq


def getFASTA(scop_text, scop_sf_fasta, out_fasta_path):
    if not out_fasta_path.exists():
        pdb_id = scop_text.split()[1]
        chain = scop_text.split()[2][0]
        seq = getSeq(scop_text.split()[5], scop_sf_fasta)
        seq = Seq(seq)
        seq_r = SeqRecord(seq, id=pdb_id + '_' + chain, description=scop_text)
        SeqIO.write(seq_r, out_fasta_path, 'fasta')


def getSCOP_text(scop_cla_text, sf_domid):
    with open(scop_cla_text, 'r') as f:
        lines = f.readlines()
    sf_text_array = list(filter(lambda x: 'SF=' + sf_domid in x, lines))
    return sf_text_array


def get_HitNum_from_xml(xml_path: str):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    hit_num_list = [None] * 3
    for i, iteration in enumerate(root[8]):
        hit_num = len(iteration[4])
        hit_num_list[i] = hit_num
    return hit_num_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('sf_id', type=str, help='super family id of scop')
    parser.add_argument('--blastdb', '-b', type=str,
                        default='../../../../../blastdb/pdbaa_20200712/pdbaa')
    parser.add_argument('--scop_cla_text', '-s', type=str,
                        default='../../../scop/20210330/scop-cla-latest.txt')
    parser.add_argument('--scop_sf_fasta', '-f', type=str,
                        default='../../../scop/20210330/scop_sf_represeq_lib20210330.fa')
    args = parser.parse_args()

    dataset_name = 'scop_sf_check'

    scop_text_array = getSCOP_text(args.scop_cla_text, args.sf_id)
    out_fasta_dir = Path('../../../fasta') / dataset_name / args.sf_id
    out_fasta_dir.mkdir(parents=True, exist_ok=True)
    out_xml_dir = Path('../../../blast_xml') / Path(args.blastdb).parent.name / dataset_name / args.sf_id
    out_xml_dir.mkdir(parents=True, exist_ok=True)

    sf_domid_list = []
    hit_num_list = []
    for scop_text in scop_text_array:
        sf_domid = scop_text.split()[5]
        print('SF DOMID: {}'.format(sf_domid))
        sf_domid_list.append(sf_domid)
        fasta_path = (out_fasta_dir / sf_domid).with_suffix('.fasta')
        if not fasta_path.exists():
            getFASTA(scop_text, args.scop_sf_fasta, fasta_path)
        out_xml_path = (out_xml_dir / sf_domid).with_suffix('.xml')
        if not out_xml_path.exists():
            psiblast_xml(fasta_path, out_xml_path, args.blastdb)
        try:
            hit_num = get_HitNum_from_xml(out_xml_path)
        except Exception:
            psiblast_xml(fasta_path, out_xml_path, args.blastdb)
            hit_num = get_HitNum_from_xml(out_xml_path)
        hit_num_list.append(hit_num)
    out_csv_path = out_xml_dir.with_suffix('.csv')
    if not out_csv_path.exists():
        df = pd.DataFrame(hit_num_list, index=sf_domid_list, columns=['1', '2', '3'])
        df.to_csv(out_csv_path)
