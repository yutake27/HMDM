import gzip
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import prody.atomic.atomgroup
from Bio.PDB import PDBIO, MMCIFParser, PDBExceptions

warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)


class DownloadProtein:
    rcsb_download_url = 'https://files.rcsb.org/download/'

    @classmethod
    def download_native_pdb(cls, pdb_id: str, chain: str, output_path: str):
        """Download the pdb file. If the pdb format is not publicly available, download the mmcif format onece,
        select only the specified chain, and convert it to pdb.

        Args:
            pdb_id (str): PDB ID.
            chain (str): Chain name.
            output_path (str): Path to the output file.
        """
        try:
            cls.download_pdb(pdb_id, output_path)
        except urllib.request.HTTPError:
            cls.download_mmcif(pdb_id, output_path)
            cls.mmcif2pdb(str(output_path), chain, str(output_path))

    @classmethod
    def download_template_pdb(cls, pdb_id: str, chain: str, output_dir: str = '.') -> str:
        """Download the pdb with the specified PDB ID and chain.
        If the pdb format is publicly available, download the pdb as (output_dir/PDBID.pdb).
        If the pdb format is not publicly available, download the mmcif format once,
        select only the specified chain, and convert it to pdb.
        The output file name is (output_dir/PDBID_chain.pdb).

        Args:
            pdb_id (str): PDB ID.
            chain (str): Chain name.
            output_dir (str, optional): Path to the output directory.

        Returns:
            str: Path to the output file.
            str: Name of the chain. If chain name is two or more characters, the first character is returned.
        """
        try:
            output_filename = pdb_id.upper() + '.pdb'
            output_path = Path(output_dir) / output_filename
            cls.download_pdb(pdb_id, output_path)
        except urllib.request.HTTPError:
            print('Download mmcif')
            output_filename = pdb_id.upper() + '_' + chain + '.pdb'
            output_path = Path(output_dir) / output_filename
            cls.download_mmcif(pdb_id, output_path)
            print('Convert mmcif to pdb')
            cls.mmcif2pdb(str(output_path), chain, str(output_path))
            print('Finish converting')
        return str(output_path)

    @staticmethod
    def _download_and_decompress_file_from_url(url: str, output_path: str):
        urllib.request.urlretrieve(url, output_path)
        with gzip.open(output_path, 'rt') as f:
            lines = f.readlines()
        with open(output_path, 'w') as f:
            f.writelines(lines)

    @classmethod
    def download_pdb(cls, pdb_id: str, output_path: str = None):
        pdb_url = cls.rcsb_download_url + pdb_id.upper() + '.pdb.gz'
        if output_path is None:
            output_path = pdb_id.upper() + '.pdb'
        cls._download_and_decompress_file_from_url(pdb_url, output_path)

    @classmethod
    def download_mmcif(cls, pdb_id: str, output_path: str = None):
        mmcif_url = cls.rcsb_download_url + pdb_id.upper() + '.cif.gz'
        if output_path is None:
            output_path = pdb_id.upper() + '.cif'
        cls._download_and_decompress_file_from_url(mmcif_url, output_path)

    @classmethod
    def mmcif2pdb(cls, input_mmcif_path: str, chain: str, output_pdb_path: str) -> str:
        parser = MMCIFParser()
        structure = parser.get_structure('', input_mmcif_path)
        try:
            sel_structure = structure[0][chain]
        except KeyError:
            print('chain list:', sorted([chain.id for chain in structure.get_chains()]))
            raise KeyError('The specified chain "{}" does not exist in template mmcif'.format(chain))

        if len(chain) > 1:
            if chain[0] in [chain.id for chain in structure.get_chains()]:
                structure[0][chain[0]].id = None
            sel_structure.id = chain[0]
        io = PDBIO()
        io.set_structure(sel_structure)
        io.save(output_pdb_path)

    @staticmethod
    def correct_resnums(mol: prody.atomic.atomgroup.AtomGroup) -> np.ndarray:
        """correct the residue number if different residues are assigned the same residue number.

        Args:
            mol (prody.Atom.atomic): mol read by prody.

        Returns:
            np.ndarray: residue numbers after modification.
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
