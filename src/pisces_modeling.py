import argparse

from Bio import SeqIO

from make_dataset import MakeDataset


class pisces_entry(MakeDataset):
    def __init__(self, target_name: str, dataset_name: str, blast_db: str, pdb_fasta: str):
        """modeling for pisces entry

        Args:
            target_name (str): target name. PDBID_CHAIN. For example, 1ADB_A.
            dataset_name (str): name of the dataset
            blast_db (str): path to the blast database of the pdbaa
            pdb_fasta (str): path to the fasta file which contains all entries in the PDB archive.
        """

        super().__init__(target_name, dataset_name, blast_db)
        self.pdb_fasta = pdb_fasta
        print(self.target_name)

    def _get_fasta(self) -> None:
        """generate fasta of the SF_DOMID
        """
        if not self.native_fasta_path.exists():
            records_dict = SeqIO.to_dict(SeqIO.parse(self.pdb_fasta, 'fasta'))
            seq_id = self.pdb_id.lower() + '_' + self.chain
            seq = records_dict[seq_id]
            SeqIO.write(seq, self.native_fasta_path, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Download pdb and align to scop text.')
    parser.add_argument('target_name', type=str, help='PISCES target name. For example, 1ADB_A.')
    parser.add_argument('--dataset_name', '-d', type=str, required=True, help='name of dataset')
    parser.add_argument('--blast_db', '-b', type=str, help='blastdb path',
                        default='../../../blastdb/pdbaa_20210407/pdbaa')
    parser.add_argument('--pdb_fasta', '-p', type=str, help='All fasta of pdb.',
                        default='../PDBseq/pdb_seqres.txt')
    args = parser.parse_args()

    pe = pisces_entry(args.target_name, args.dataset_name, args.blast_db, args.pdb_fasta)
    pe.make_dataset()

