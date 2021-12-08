import argparse

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from make_dataset import MakeDataset


class scop_sf_dom(MakeDataset):
    def __init__(self, sf_domid: str, dataset_name: str, blast_db: str, scop_sf_fasta: str, scop_cla_text: str):
        """modeling for scop SF_DOMID

        Args:
            sf_domid (str): SF_DOMID
            dataset_name (str): name of the dataset
            blast_db (str): path to the blast database of the pdbaa
            scop_sf_fasta (str): path to the scop_sf_represeq
            scop_cla_text (str): path to the scop-cla-latest
        """

        self.sf_domid = sf_domid
        self.scop_sf_fasta = scop_sf_fasta
        self.scop_cla_text = scop_cla_text
        self.scop_text = self._search_SCOP(self.sf_domid)
        pdb_id, chain, resnum_start, resnum_end = self._parse_SCOP(self.scop_text)
        target_name = pdb_id + '_' + chain
        super().__init__(target_name, dataset_name, blast_db, resnum_start=resnum_start, resnum_end=resnum_end)
        print(target_name)

    def _search_SCOP(self, sf_domid: str) -> str:
        """search information for sf_domid from scop cla

        Args:
            sf_domid (str): SF_DOMID

        Returns:
            str: sf_domid line of the scop-cla-latest
        """
        with open(self.scop_cla_text, 'r') as f:
            lines = f.readlines()
            lines = list(filter(lambda line: sf_domid in line, lines))
            scop_text = lines[0]
        return scop_text

    @staticmethod
    def _parse_SCOP(scop_text: str) -> (str, str, str, str):
        """parse line of the scop-cla-latest

        Args:
            scop_text (str): line of the scop-cla-latest
        Returns:
            str: pdb id.
            str: chain name.
            str: residue number of the start.
            str: residue number of the end.
        """
        scop_split = scop_text.split()
        pdb_id = scop_split[6]
        chain = scop_split[7].split(':')[0]
        resnum_start, resnum_end = scop_split[7].split(':')[1].split('-')
        return pdb_id, chain, resnum_start, resnum_end

    def _get_seq_from_SCOP(self) -> str:
        """get sequence from scop_sf_represeq

        Returns:
            str: sequence (not including header)
        """

        with open(self.scop_sf_fasta, 'r') as f:
            lines = f.readlines()
        index = [i for i, line in enumerate(lines) if line[1: 8] == self.sf_domid][0]
        seq = lines[index + 1][: -1]
        return seq

    def _get_fasta(self) -> None:
        """generate fasta of the SF_DOMID
        """
        if not self.native_fasta_path.exists():
            seq = self._get_seq_from_SCOP()
            seq = Seq(seq)
            seq_r = SeqRecord(seq, id=self.pdb_id + '_' + self.chain, description=self.scop_text)
            SeqIO.write(seq_r, self.native_fasta_path, 'fasta')

    def _xml2pir(self) -> None:
        """Save SF-DOMID into the csv.
        """
        super()._xml2pir()
        df = pd.read_csv(self.xml_csv_path, index_col=0)
        df['SF-DOMID'] = self.sf_domid
        df.to_csv(self.xml_csv_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Download pdb and align to scop text.')
    parser.add_argument('sf_domid', type=str, help='super family id of scop-cla-latest.txt')
    parser.add_argument('--dataset_name', '-d', type=str, help='name of dataset', default='target_50')
    parser.add_argument('--blast_db', '-b', type=str, help='blastdb path',
                        default='../../../blastdb/pdbaa_20210407/pdbaa')
    parser.add_argument('--scop_sf_fasta', type=str, help='Path to scop_sf_represeq_lib fasta',
                        default='../scop/20210330/scop_sf_represeq_lib20210330.fa')
    parser.add_argument('--scop_cla_text', type=str, help='Path to scop-cla-latext.txt',
                        default='../scop/20210330/scop-cla-latest.txt')
    args = parser.parse_args()

    sd = scop_sf_dom(args.sf_domid, args.dataset_name, args.blast_db, args.scop_sf_fasta, args.scop_cla_text)
    sd.make_dataset()

