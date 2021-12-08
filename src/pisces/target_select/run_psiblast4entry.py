import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def psiblast(fasta_path, out_xml_path, blast_db, iteration='3', max_sequence=10000):
    cmd = ['psiblast', '-query', fasta_path, '-db', blast_db, '-out', out_xml_path,
           '-num_iterations', iteration, '-outfmt', '5', '-evalue', '1.0e-3',
           '-max_target_seqs', str(max_sequence)]
    subprocess.run(cmd)


def main():
    multidomain_csv_path = '../../../pisces/20210225/multidomain_cullpdb_pc20_res2.0_R0.25.csv'
    multidomain_df = pd.read_csv(multidomain_csv_path)
    pdb_fasta_path = '../../../PDBseq/pdb_seqres.txt'
    out_fasta_dir = Path('../../../pisces/20210225/multi-domain_fasta')
    out_fasta_dir.mkdir(parents=True, exist_ok=True)
    out_xml_dir = Path('../../../pisces/20210225/multi-domain_blast')
    out_xml_dir.mkdir(parents=True, exist_ok=True)
    blast_db_path = '../../../../../blastdb/pdbaa_20210407/pdbaa'
    records_dict = SeqIO.to_dict(SeqIO.parse(pdb_fasta_path, 'fasta'))

    for index, rows in multidomain_df.iterrows():
        seq_id = rows['PDB_ID'].lower() + '_' + rows['Chain']
        out_fasta_path = (out_fasta_dir / seq_id.upper()).with_suffix('.fasta')
        out_xml_path = (out_xml_dir / seq_id.upper()).with_suffix('.xml')
        if out_fasta_path.exists() and out_xml_path.exists():
            continue
        try:
            print(seq_id)
            seq = records_dict[seq_id]
        except KeyError:
            print(seq_id)
        else:
            SeqIO.write(seq, out_fasta_path, 'fasta')
            psiblast(out_fasta_path, blast_db_path, out_xml_path)


if __name__ == '__main__':
    main()
