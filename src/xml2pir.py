import copy
import shutil
import time
from pathlib import Path
from typing import List, Dict

import numpy as np
import pandas as pd
from Bio import SearchIO, SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.PDB.PDBExceptions import PDBConstructionException
from prody import parsePDB, writePDB

import clstr_blast_hit
from seq import AlignSeq, ReadSeq
from download_protein import DownloadProtein


class blast_hit:
    def __init__(self, hit, target_fasta, template_pdb_dir):
        self.target_fasta = target_fasta
        self.query = copy.deepcopy(hit[0][0].query)
        self.qseq_start = hit[0][0].query_start
        self.qseq_end = hit[0][0].query_end
        self.hit = copy.deepcopy(hit[0][0].hit)
        print(self.hit.id)
        self.hit_chain = self.hit.id.split('_')[1]
        self.hit_pdbid = self.hit.id.split('_')[0]
        self.hit.id = self.hit.id.split('_')[0]
        self.template_pdb_dir = template_pdb_dir
        self.__convert_hit_to_pir()

    def __correct_pdb_resnums(self, mol):
        """correct pdb residue number such as 1A,1B, and 1C

        Args:
            mol (prody.atomic.atomgroup.AtomGroup): pdb mol

        Returns:
            np.ndarray : correct residue number
        """
        new_resnum_list = []
        plus = 0
        resnum_before = float('inf')
        resid_before = float('inf')
        for resnum, resid in zip(mol.getResnums(), mol.getResindices()):
            if resnum_before == resnum and resid_before != resid:
                plus += 1
            new_resnum_list.append(resnum + plus)
            resnum_before = resnum
            resid_before = resid
        return np.array(new_resnum_list)

    def __remove_pdb_resnum_alphabet(self, pdb_file: str):
        """remove the alphabet from the residue number in pdb file

        Args:
            pdb_file (str) : pdb file path
        """
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
            lines = list(filter(lambda line: len(line) > 0, lines))
            for i, line in enumerate(lines):
                if i > 0:
                    line = line[:26] + ' ' + line[27:]
                    lines[i] = line
        with open(pdb_file, 'w') as f:
            f.writelines(lines)

    def __download_template_pdb(self) -> str:
        """download template pdb files by prody.
        In some cases, the residue numbers are the same but the residues are different,
        so they are corrected.　Also, HETATM is excluded.
        Some residue numbers contain alphabetic characters (e.g., 1A), which cannot be processed,
        so the alphabetic characters are deleted.

        Args:
            pdb_id (str): PDB ID.

        Returns:
            str: output path of the pdb.
        """
        template_pdb_dir = Path(self.template_pdb_dir)
        template_pdb_file = (template_pdb_dir / self.hit_pdbid).with_suffix('.pdb')
        template_pdb_chain_file = (template_pdb_dir / (self.hit_pdbid + '_' + self.hit_chain)) .with_suffix('.pdb')
        if template_pdb_file.exists():
            return template_pdb_file
        elif template_pdb_chain_file.exists():
            return template_pdb_chain_file
        else:
            pdb_file = DownloadProtein.download_template_pdb(self.hit_pdbid, self.hit_chain, '.')
            mol = parsePDB(pdb_file)
            correct_resnums = self.__correct_pdb_resnums(mol)
            mol.setResnums(correct_resnums)
            mol = mol.select('not hetatm')
            writePDB(pdb_file, mol)
            self.__remove_pdb_resnum_alphabet(pdb_file)
            output_template_path = template_pdb_dir / Path(pdb_file).name
            shutil.move(pdb_file, output_template_path)
            time.sleep(3)
            return output_template_path

    def __convert_qseq(self) -> None:
        """Convert query sequence.
        """
        self.query.id = 'TARGET'
        self.query.name = None
        self.query.description = 'sequence:TARGET:1::' + str(len(self.target_fasta)) + ':::::'
        self.query.seq = self.target_fasta.seq[: self.qseq_start] + \
            self.query.seq + self.target_fasta.seq[self.qseq_end:]

    def __convert_hit_residue_x(self) -> None:
        """Convert 'X' residue in fasta to 'M'.
        """
        seq_array = np.array(list(self.hit.seq))
        np.place(seq_array, seq_array == 'X', 'M')
        self.hit.seq = ''.join(seq_array)

    def __convert_hseq(self, align_tseq: np.ndarray, align_hindices: np.ndarray, hseq_start: int, hseq_end: int):
        """Convert hit sequence according to the result of alignment.
        Insert gaps where there is no structure in the template pdb.

        Args:
            align_tseq (np.ndarray): [description]
            align_hindices (np.ndarray): [description]
            hseq_start (int): [description]
            hseq_end (int): [description]
        """
        self.hit.name = None
        hit_not_missing_indices = np.where(np.array(list(self.hit.seq)) != '-')[0]
        missing_list = hit_not_missing_indices[np.setdiff1d(np.arange(len(align_tseq)), align_hindices)]
        hseq_list = np.array(list(self.hit.seq))
        hseq_list[missing_list] = '-'
        hseq = ''.join(hseq_list)
        self.hit.seq = '-' * self.qseq_start + hseq + '-' * (len(self.target_fasta) - self.qseq_end)
        self.hit.description = 'structure:' + self.hit.id + ':' + str(hseq_start) + ':' + self.hit_chain[0]\
            + ':' + str(hseq_end) + ':' + self.hit_chain[0] + '::::'

    def __convert_hit_to_pir(self) -> None:
        """convert blast hit to a pir file.
        """
        self.__convert_qseq()
        template_pdb_path = self.__download_template_pdb()
        self.hit.id = Path(template_pdb_path).stem
        template_seq, template_resnums = ReadSeq.pdb2seq(template_pdb_path, self.hit_chain[0], insert_gap=True)
        hit_seq_fil = ''.join(list(filter(lambda x: x != '-', self.hit.seq)))
        align_hseq, align_tseq, align_hindices, align_tindices = AlignSeq.align_fasta_and_pdb(hit_seq_fil, template_seq)
        hseq_start, hseq_end = template_resnums[align_tindices[0]], template_resnums[align_tindices[-1]]
        self.__convert_hseq(align_tseq, align_hindices, hseq_start, hseq_end)

    def save_pir(self, out_pir_path: str) -> None:
        """save pir file.

        Args:
            out_pir_path (str): output path of the pir file.
        """
        self.query.seq = Seq(str(self.query.seq))
        self.hit.seq = Seq(str(self.hit.seq))
        self.query.annotations = {'molecule_type': 'protein'}
        self.hit.annotations = {'molecule_type': 'protein'}
        SeqIO.write((self.query, self.hit), out_pir_path, 'pir')


class blast_xml:
    def __init__(self, xml_path: str, fasta_path: str, out_template_pdb_dir: str):
        """Generate pir file from blast result in xml format.

        Args:
            xml_path (str): Path to the blast result in xml format.
            fasta_path (str): Path to the target fasta file.
            out_template_pdb_dir (str): Path to the directory of template pdbs.
        """
        self.xml_path = Path(xml_path)
        self.qresult_list = self.read_xml(str(xml_path))
        self.target_name = self.qresult_list[0].id
        self.fasta_path = Path(fasta_path)
        self.target_fasta = SeqIO.read(self.fasta_path, 'fasta')
        self.out_template_pdb_dir = out_template_pdb_dir

    def read_xml(self, xml_path: str) -> List[SearchIO._model.query.QueryResult]:
        """Parse the result of blast in xml format.

        Args:
            xml_path (str): path of the blast result in xml format.

        Returns:
            List[SearchIO._model.query.QueryResult]: List of the QueryResult of each iteration.
        """
        qresult_list = list(SearchIO.parse(xml_path, 'blast-xml'))
        # convert hit.id to PDBID_CHAIN
        for qresult in qresult_list:
            for hit in qresult:
                hit_pdbid = hit.id.split('|')[1]
                hit_chain = hit.id.split('|')[2]
                hit.id = hit_pdbid + '_' + hit_chain
        return qresult_list

    def filter_gap_from_seq(self, seq: Seq) -> Seq:
        """delete gaps ('-') from sequence.

        Args:
            seq (Seq): a sequence that includes gaps.

        Returns:
            Seq: Sequence without gaps.
        """
        filtered_seq = ''.join(list(filter(lambda x: x != '-', seq)))
        return Seq(filtered_seq)

    def generate_seqrecord_from_hit(self,
                                    hit: SeqRecord.SeqRecord,
                                    iteration: int,
                                    hit_index: int) -> SeqRecord.SeqRecord:
        """From blast Hit, generate SeqRecord whose id includes iteration and hit number information,
        and whose sequence has no gap ('-').

        Args:
            hit (SeqRecord.SeqRecord): Hit of QueryResult
            iteration (int): nubmer of blast iteration (Iteration_iter_num)
            hit_index (int): index of hit (hit_num - 1)

        Returns:
            SeqRecord.SeqRecord: new SeqRecord
        """
        hit_seqrecord = copy.deepcopy(hit[0][0].hit)
        hit_id = hit_seqrecord.id
        new_id = hit_id + '-' + str(iteration) + '-' + str(hit_index)
        hit_seqrecord.id = new_id
        hit_seqrecord.seq = self.filter_gap_from_seq(hit_seqrecord.seq)
        return hit_seqrecord

    def get_seqrecord_list_from_qresult(self,
                                        qresult: SearchIO._model.query.QueryResult,
                                        iteration: int) -> List[SeqRecord.SeqRecord]:
        """get seqrecord list from qresult of a iteration.

        Args:
            qresult (SearchIO._model.query.QueryResult): QueryResult of a iteration
            iteration (int): the number of the iteration.

        Returns:
            List[SeqRecord.SeqRecord]: list of the SeqRecord.
        """
        seqrecord_list = []
        for hit_index, hit in enumerate(qresult):
            hit_seqrecord = self.generate_seqrecord_from_hit(hit, iteration, hit_index)
            seqrecord_list.append(hit_seqrecord)
        return seqrecord_list

    def make_fasta_file_from_seqrecord_list(self, seqrecord_list: List, out_fasta_path: str) -> None:
        """make fasta file from the list of the SeqRecord.

        Args:
            seqrecord_list (list): the list of the SeqRecord.
            out_fasta_path (str): output path of the fasta.
        """
        SeqIO.write(seqrecord_list, out_fasta_path, 'fasta')

    def xml2fasta(self, out_fasta_path: str) -> None:
        """Convert the hit sequences in a single iteration into a single fasta.

        Args:
            out_fasta_path (str): output path of the fasta.
        """
        seqrecord_list = []
        for iteration, qresult in enumerate(self.qresult_list, start=1):
            seqrecord_list_iteration = self.get_seqrecord_list_from_qresult(qresult, iteration)
            seqrecord_list.extend(seqrecord_list_iteration)
        self.make_fasta_file_from_seqrecord_list(seqrecord_list, out_fasta_path)

    def __calc_coverage_considering_missing_residue_from_pir(self, pir_file: str) -> int:
        """calc coverage from pir file considering missing residues.

        Args:
            pir_file (str): Path to pir file.

        Returns:
            int: Coverage considering missing residue.
        """
        query_record, template_record = SeqIO.parse(pir_file, 'pir')
        query_seq, template_seq = query_record.seq, template_record.seq
        template_not_gap_indices = np.where(np.array(template_seq) != '-')[0]
        aln_start, aln_end = template_not_gap_indices[0], template_not_gap_indices[-1]
        query_aln_seq_array = np.array(query_seq)[aln_start: aln_end + 1]
        coverage_misres = len(np.where(query_aln_seq_array != '-')[0])
        return coverage_misres

    def __get_template_quality(self, seq_len: int, hit: SearchIO._model.hit.Hit, pir_path: str) -> Dict:
        """Get template quality from hit.
        Get the following labels as template quality.
        * identity
        * positive
        * coverage
        * identity - missing residues
        * positive - missing residues
        * coverage - missing residues
        * missing residues
        * sequence length

        Args:
            seq_len (int): Sequence length.
            hit (SearchIO._model.hit.Hit): hit of blast.
            pir_path (str): Path to pir file.

        Returns:
            Dict: Dict of template quality.
        """
        aln_annotation = hit[0].aln_annotation['similarity']
        hit_seq = str(hit[0].hit.seq)
        identity = hit[0].ident_num
        positive = hit[0].pos_num
        start = hit[0].query_start
        end = hit[0].query_end
        coverage = end - start
        hit_seq_end = - (seq_len - end) if seq_len != end else None
        pir_template_seq = str(list(SeqIO.parse(pir_path, 'pir'))[1].seq)[start: hit_seq_end]
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
        coverage_misres = self.__calc_coverage_considering_missing_residue_from_pir(pir_path)

        return {'seq_len': seq_len, 'identity': identity, 'positive': positive, 'coverage': coverage,
                'identity(-misres)': identity - not_ident_num, 'positive(-misres)': positive - not_pos_num,
                'coverage(-misres)': coverage_misres, 'num_misres': missing_num}

    def convert_hit_to_pir(self, seq_len: int, hit: SearchIO._model.hit.Hit, iteration: int, out_pir_dir: Path) -> Dict:
        """convert the hit to the pir file.

        Args:
            seq_len (int): Length of sequence.
            hit (SearchIO._model.hit.Hit): hit read by biopython.
            iteration (int): the number of the iteration.
            out_pir_dir (Path): path of the output pir directory.

        Return:
            Dict: Dict of template quality.
        """

        template_name = hit.id + '_' + str(iteration)
        out_pir_path = (out_pir_dir / template_name).with_suffix('.pir')
        if not out_pir_path.exists():
            bh = blast_hit(hit, self.target_fasta, self.out_template_pdb_dir)
            bh.save_pir(out_pir_path)
        template_quality_dict = self.__get_template_quality(seq_len, hit, out_pir_path)
        template_quality_dict['template'] = template_name
        return template_quality_dict

    def __make_hit_info_df(self, hit_info_array: List, hit_id_array: List) -> pd.DataFrame:
        """make DataFrame about hit information (identity, positive, evalue, length of the sequence).

        Args:
            hit_info_array (List): the list of the hit information.
            hit_id_array (List): the list of the hit ids.

        Returns:
            pd.DataFrame: DataFrame about hit information.
        """
        df = pd.DataFrame(hit_info_array,
                          columns=['seq_len', 'evalue', 'identity', 'positive', 'align_len'],
                          index=hit_id_array)
        df = df.sort_index()
        return df

    def convert_xml_to_pir_nofiltering(self,
                                       out_pir_dir: Path,
                                       iteration_array: List = [1, 2, 3],
                                       num_hit: int = 10,
                                       identity_threshold: float = 0.95,
                                       coverage_minimum_threshold: float = 0.6) -> pd.DataFrame:
        """Convert blast xml to pir file.
        Select a specific number of templates for each iteration.
        Similarity between hit sequences is not taken into account when selecting templates.
        The pdb of the target itself is not used as a template.
        Do not select the exact same template (same template pdb id and alignment)
        even if iterations are different.

        Args:
            out_pir_dir (Path): Path to the output directory of the pir files.
            iteration_array (list, optional): List of iterations to generate pir. Defaults to [1, 2, 3].
            num_hit (int, optional): The number of pir files to generate in each iteration. Defaults to 10.
            identity_threshold (float, optional): Identity threshold for template selection.
            If the identity of the hit is greater than this threshold, the hit is not selected as a template.
            coverage_minimum_threshold (float, optional): Minimum threshold of coverage.
            Models with coverage lower than this value are excluded.

        Returns:
            pd.DataFrame: DataFrame about hit information.
        """
        target_pdb_id = self.target_name[: 4].upper()
        hseq_dic = {}
        template_quality_dict_array = []

        for iteration in iteration_array:
            hit_num = 0
            try:
                qresult = self.qresult_list[iteration - 1]
            except Exception as e:
                print(e)
                break
            seq_len = qresult.seq_len
            pdb_id_array = []
            hseq_array = []
            for hit in qresult:
                hseq = str(hit[0][0].hit.seq)
                hit_pdb_id = hit.id[: 4]
                coverage = hit[0].query_end - hit[0].query_start
                if target_pdb_id != hit_pdb_id \
                        and hit[0].ident_num / seq_len < identity_threshold \
                        and hit_num < num_hit \
                        and not (hit.id in hseq_dic and hseq_dic[hit.id] == hseq) \
                        and coverage / seq_len >= coverage_minimum_threshold:
                    try:
                        template_quality_dict = self.convert_hit_to_pir(seq_len, hit, iteration, out_pir_dir)
                    except (FileNotFoundError, AttributeError,
                            ValueError, KeyError, PDBConstructionException, AssertionError) as e:
                        print(e)
                        continue
                    hit_num += 1
                    pdb_id_array.append(hit_pdb_id)
                    hseq_array.append(hseq)
                    hseq_dic[hit.id] = hseq
                    template_quality_dict_array.append(template_quality_dict)
        template_df = pd.DataFrame(template_quality_dict_array).set_index('template')
        return template_df

    def convert_xml_to_pir_from_clstr_df(self,
                                         out_pir_dir: Path,
                                         df: pd.DataFrame,
                                         num_hit: int,
                                         identity_threshold: float = 0.95,
                                         coverage_minimum_threshold: float = 0.6) -> pd.DataFrame:
        """select templates from clustering DataFrame and convert hit to pir.

        Args:
            out_pir_dir (Path): Path to the output directory of the pir files.
            df (pd.DataFrame): A data frame containing the results of clustering.
            num_hit (int): Number of templates to select for each iteration.
            identity_threshold (float, optional): Identity threshold for template selection.
            If the identity of the hit is greater than this threshold, the hit is not selected as a template.
            coverage_minimum_threshold (float, optional): Minimum threshold of coverage.
            Models with coverage lower than this value are excluded.

        Returns:
            pd.DataFrame: DataFrame about hit information.
        """
        target_pdb_id = self.target_name[: 4].upper()
        selected_clstr_list = []
        template_quality_dict_array = []

        for iteration, group in df.groupby('iteration'):
            hit_num = 0
            qresult = self.qresult_list[iteration - 1]
            seq_len = qresult.seq_len
            for index, row in group.iterrows():
                clstr_id = row['clstr_id']
                hit = qresult[row['hit_index']]
                hit_pdb_id = hit.id[: 4].upper()
                coverage = hit[0].query_end - hit[0].query_start
                if hit_pdb_id != target_pdb_id \
                        and clstr_id not in selected_clstr_list \
                        and hit[0].ident_num / seq_len < identity_threshold \
                        and hit_num < num_hit \
                        and coverage / seq_len >= coverage_minimum_threshold:
                    try:
                        template_quality_dict = self.convert_hit_to_pir(seq_len, hit, iteration, out_pir_dir)
                    except (FileNotFoundError, AttributeError,
                            ValueError, KeyError, PDBConstructionException, AssertionError) as e:
                        print(e)
                        continue
                    hit_num += 1
                    template_quality_dict_array.append(template_quality_dict)
                    selected_clstr_list.append(clstr_id)
        template_df = pd.DataFrame(template_quality_dict_array).set_index('template')
        return template_df

    def convert_xml_to_pir_filtering(self,
                                     out_pir_dir: Path,
                                     iteration_array=[1, 2, 3],
                                     num_hit=10) -> pd.DataFrame:
        """Convert blast xml to pir file. Filtering similar alignment using cd-hit.
        Select a specific number of templates for each iteration.
        Cluster the hit sequences using cd-hit, and select only one template from each cluster.
        The pdb of the target itself is not used as a template.

        Args:
            out_pir_dir (Path): Path to the output directory of the pir files.
            iteration_array (list, optional): List of iterations to generate pir. Defaults to [1, 2, 3].
            num_hit (int, optional): The number of pir files to generate in each iteration. Defaults to 10.

        Returns:
            pd.DataFrame: DataFrame about hit information.
        """
        out_fasta_path = str(self.xml_path.with_suffix('.fasta'))
        if not Path(out_fasta_path).exists():
            self.xml2fasta(out_fasta_path)
        df = clstr_blast_hit.clstr_blast_hit(out_fasta_path)
        hit_df = self.convert_xml_to_pir_from_clstr_df(out_pir_dir, df, num_hit=num_hit)
        return hit_df
