import argparse
import re
import shutil
import subprocess
import time
from functools import reduce
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SearchIO, SeqIO, pairwise2
from Bio.Alphabet import generic_protein
from Bio.PDB import PDBList, parse_pdb_header
from Bio.Seq import Seq
from prody import parsePDB, writePDB


def convert_pdb_to_seq(pdb_path):
    mol = parsePDB(pdb_path)
    resindices = mol.getResindices()
    resindices_diff = np.diff(resindices, prepend=float('inf'))
    res_start_indices = np.where(resindices_diff != 0)[0]
    seq = ''.join(np.array(list(mol.getSequence()))[res_start_indices])
    resnums = mol.getResnums()[res_start_indices]
    return mol, seq, resnums


def read_fasta(fasta_path):
    fasta = SeqIO.read(fasta_path, 'fasta')
    return fasta.seq


def correct_misalignment_by_X(aseq_array: np.ndarray, bseq_array:np.ndarray):
    adel_list = []
    bdel_list = []
    for i, bres in enumerate(bseq_array):
        if i != len(aseq_array) - 1 and bres == '-' and aseq_array[i] == 'X' and aseq_array[i + 1] == '-':
            adel_list.append(i + 1)
            bdel_list.append(i)
    aseq_del = np.delete(aseq_array, adel_list, 0)
    bseq_del = np.delete(bseq_array, bdel_list, 0)
    return aseq_del, bseq_del


def align_seq(seqa, seqb):
    alignment_list = pairwise2.align.globalms(seqa, seqb, 5, -10, -2, -1)
    alignment = None
    var_indices_min = float('inf')
    for a in alignment_list:
        var_indices = np.var(np.where(np.array(list(a.seqA)) != '-')[0])
        if var_indices < var_indices_min:
            alignment = a
            var_indices_min = var_indices
    align_seqa, align_seqb = np.array(list(alignment.seqA)), np.array(list(alignment.seqB))
    align_seqa, align_seqb = correct_misalignment_by_X(align_seqa, align_seqb)
    align_seqa, align_seqb = align_seqa[np.where(align_seqb != '-')], align_seqb[np.where(align_seqa != '-')]
    align_aindices = np.where(align_seqb != '-')[0]
    align_bindices = np.where(align_seqa != '-')[0]
    return align_seqa, align_seqb, align_aindices, align_bindices



def correct_pdb_resnum(pdb_path, fasta_path):
    """align pdb with fasta, and correct residue numbers in pdb.

    Args:
        pdb_path (str): pdb file path. correct the residue number, and overwrite it.
        fasta_path (str): fasta file path
    """
    mol, pdb_seq, resnums = convert_pdb_to_seq(pdb_path)
    fasta_seq = read_fasta(fasta_path)
    align_pdb_seq, align_fasta_seq, align_pdb_indices, align_fasta_indices = align_seq(pdb_seq, fasta_seq)
    align_mol = mol.select('resindex {}'.format(reduce(lambda a, b: str(a) + ' ' + str(b), align_pdb_indices)))
    convert_dic = {pdb_index: fasta_index for pdb_index, fasta_index in zip(align_pdb_indices, align_fasta_indices + 1)}
    correct_resnums = np.array([convert_dic[resindex] for resindex in align_mol.getResindices()])
    align_mol.setResnums(correct_resnums)
    writePDB(pdb_path, align_mol)


