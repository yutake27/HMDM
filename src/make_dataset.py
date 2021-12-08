import os
import subprocess
from abc import ABCMeta, abstractmethod
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd
import prody.atomic.atomgroup
from prody import parsePDB, writePDB

from download_protein import DownloadProtein
from get_gdt import get_gdt_for_target_df
from modeller_modeling import modeller_modeling
from seq import AlignSeq, ReadSeq
from xml2pir import blast_xml


class MakeDataset(metaclass=ABCMeta):
    def __init__(self, target_name: str, dataset_name: str, blast_db: str,
                 resnum_start: int = None, resnum_end: int = None):
        """Abstract modeling class

        Args:
            target_name (str): target name. PDBID_CHAIN. For example, 1ADB_A.
            dataset_name (str): name of the dataset
            blast_db (str): path to the blast database of the pdbaa
        """
        self.target_name = target_name
        self.pdb_id, self.chain = self.target_name.split('_')
        self.dataset_name = dataset_name
        self.blast_db = blast_db
        self.resnum_start = resnum_start
        self.resnum_end = resnum_end

        out_fasta_dir = Path('../fasta') / self.dataset_name
        out_fasta_dir.mkdir(parents=True, exist_ok=True)
        self.native_fasta_path = (out_fasta_dir / (self.target_name)).with_suffix('.fasta')
        native_pdb_dir = Path('../native_pdb') / self.dataset_name
        native_pdb_dir.mkdir(parents=True, exist_ok=True)
        self.native_pdb_path = (native_pdb_dir / self.target_name).with_suffix('.pdb')
        self.template_pdb_dir = Path('../template_pdb/')

        self.xml_dir = Path('../blast_xml') / Path(self.blast_db).parent.name / self.dataset_name
        self.xml_dir.mkdir(parents=True, exist_ok=True)
        self.xml_path = (self.xml_dir / self.target_name).with_suffix('.xml')
        self.xml_csv_path = (self.xml_dir / self.target_name).with_suffix('.csv')

        self.out_pir_dir = Path('../pir') / self.dataset_name / self.target_name
        self.out_pir_dir.mkdir(parents=True, exist_ok=True)
        self.out_pir_dir_nofiltering = self.out_pir_dir / 'nofiltering'
        self.out_pir_dir_nofiltering.mkdir(exist_ok=True)
        self.out_pir_dir_filtering = self.out_pir_dir / 'filtering'
        self.out_pir_dir_filtering.mkdir(exist_ok=True)

        self.out_pdb_dir = Path('../pdb') / self.dataset_name / self.target_name
        self.out_pdb_dir.mkdir(parents=True, exist_ok=True)
        self.out_pdb_dir_nofiltering = self.out_pdb_dir / 'nofiltering'
        self.out_pdb_dir_nofiltering.mkdir(exist_ok=True)
        self.out_pdb_dir_filtering = self.out_pdb_dir / 'filtering'
        self.out_pdb_dir_filtering.mkdir(exist_ok=True)

    @abstractmethod
    def _get_fasta(self) -> str:
        """Get fasta sequence.

        Returns:
            str: sequence (not including header)
        """
        pass

    @staticmethod
    def _test_match(fasta_seq: str, mol: prody.atomic.atomgroup.AtomGroup) -> None:
        """test that the fasta sequence matches to the pdb sequence.

        Args:
            fasta_seq (str): Sequence of the fasta.
            mol (prody.atomic.atomgroup.AtomGroup): PDB object read by ProDy.
        """
        pdb_seq, pdb_resnum = ReadSeq.mol2seq(mol, insert_gap=False)
        fasta_seq_array = np.array(list(fasta_seq))
        pdb_seq_array = np.copy(fasta_seq_array)
        pdb_seq_array[pdb_resnum - 1] = list(pdb_seq)
        num_diff = np.count_nonzero(fasta_seq_array != pdb_seq_array)
        num_missing = len(fasta_seq) - len(pdb_seq)
        assert num_diff < len(fasta_seq) * 0.05
        print('length:', len(fasta_seq))
        print('num different residues between pdb and fasta:', num_diff)
        print('num missing residues:', num_missing)

    def _get_pdb(self) -> None:
        """download the pdb and fix residue numbers.
        """
        if not self.native_pdb_path.exists():
            tmp_pdb_path = self.pdb_id + '.pdb'
            DownloadProtein.download_native_pdb(self.pdb_id, self.chain, tmp_pdb_path)
            fasta_seq = ReadSeq.fasta2seq(self.native_fasta_path)
            mol = parsePDB(tmp_pdb_path, chain=self.chain)
            # if the first residue number is negative, correct the residue numbers
            if mol.getResnums()[0] < 0:
                resnums = mol.getResnums()
                new_resnums = resnums - resnums[0] + 1
                mol.setResnums(new_resnums)
            new_resnums = DownloadProtein.correct_resnums(mol)
            mol.setResnums(new_resnums)
            pdb_seq, pdb_resnums = ReadSeq.mol2seq(mol, insert_gap=True)
            align_pseq, align_fseq, align_findices, align_pindices = AlignSeq.align_fasta_and_pdb(fasta_seq, pdb_seq)
            sel_mol_resnums = pdb_resnums[align_pindices]
            sel_mol = mol.select('resnum {}'.format(reduce(lambda a, b: str(a) + ' ' + str(b), sel_mol_resnums)))
            assert sel_mol is not None
            if self.resnum_start is not None and self.resnum_end is not None:
                pdb_resnum_start, pdb_resnum_end = str(sel_mol_resnums[0]), str(sel_mol_resnums[-1])
                if pdb_resnum_start != self.resnum_start or pdb_resnum_end != self.resnum_end:
                    print('The residue number of pdb and the specified number are different')
                    print('pdb start: {}, specified start: {}, pdb end: {}, specified end:{}'
                          .format(pdb_resnum_start, self.resnum_start, pdb_resnum_end, self.resnum_end))
            fasta_resnum = align_findices + 1
            convert_resnum_dict = dict(zip(sel_mol_resnums, fasta_resnum))
            new_resnum = [convert_resnum_dict[resnum] for resnum in sel_mol.getResnums()]
            sel_mol.setResnums(new_resnum)
            self._test_match(fasta_seq, sel_mol)
            writePDB(str(self.native_pdb_path), sel_mol)
            os.remove(tmp_pdb_path)

    def _psiblast_xml(self, iteration: str = '3', evalue_threshold: str = '1.0e-3') -> None:
        """psiblast against pdb.

        Args:
            iteration (str, optional): the number of the iteration. Defaults to '3'.
        """
        cmd = ['psiblast', '-query', self.native_fasta_path, '-db', self.blast_db, '-out', self.xml_path,
               '-num_iterations', iteration, '-outfmt', '5', '-evalue', evalue_threshold]
        subprocess.run(cmd)

    def _xml2pir(self) -> None:
        """convert blast xml to the pir format.
        """
        bx = blast_xml(self.xml_path, self.native_fasta_path, self.template_pdb_dir)
        nofiltering_template_df = bx.convert_xml_to_pir_nofiltering(self.out_pir_dir_nofiltering)
        filtering_template_df = bx.convert_xml_to_pir_filtering(self.out_pir_dir_filtering)
        nofiltering_template_df['Method'] = 'nofiltering'
        filtering_template_df['Method'] = 'filtering'
        concat_df = pd.concat([nofiltering_template_df, filtering_template_df])
        concat_df.to_csv(self.xml_csv_path)

    def preprocess(self) -> None:
        """Generate fasta, download pdb, psiblast, and generate pir files.
        """
        self._get_fasta()
        self._get_pdb()
        self._psiblast_xml()
        self._xml2pir()

    def check_standard(self):
        df = pd.read_csv(self.xml_csv_path)
        df['positive/seq_len'] = df['positive'] / df['seq_len']
        if len(df[df['positive/seq_len'] > 0.6]) < 6:
            print('does not meet standard')
            exit()

    def modeling(self) -> None:
        """modeling using modeller.
        """
        mm = modeller_modeling(self.target_name)
        # modeling for no-filtering templates
        mm.modeling_dir(self.out_pir_dir_nofiltering, self.template_pdb_dir, self.out_pdb_dir_nofiltering)
        # modeling for filtering templates
        mm.modeling_dir(self.out_pir_dir_filtering, self.template_pdb_dir, self.out_pdb_dir_filtering)

    def get_gdt(self) -> None:
        """get GDT score using TM-Score.
        """
        out_gdt_dir = Path('../tmscore') / self.dataset_name
        out_gdt_dir.mkdir(parents=True, exist_ok=True)
        out_gdt_path = (out_gdt_dir / self.target_name).with_suffix('.csv')
        # get GDT for no-filtering dataset
        nofiltering_gdt_df = get_gdt_for_target_df(self.native_pdb_path,
                                                   self.out_pdb_dir_nofiltering,
                                                   self.xml_csv_path)
        nofiltering_gdt_df = nofiltering_gdt_df.query('Method == "nofiltering"').sort_index()
        # get GDT for filtering dataset
        filtering_gdt_df = get_gdt_for_target_df(self.native_pdb_path,
                                                 self.out_pdb_dir_filtering,
                                                 self.xml_csv_path)
        filtering_gdt_df = filtering_gdt_df.query('Method == "filtering"').sort_index()
        gdt_df = pd.concat([nofiltering_gdt_df, filtering_gdt_df])
        gdt_df.to_csv(out_gdt_path)

    def make_dataset(self) -> None:
        """Do all process.
        """
        self.preprocess()
        self.modeling()
        self.get_gdt()

