import io
import os
import shutil
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

from modeller import *
from modeller.automodel import *


class modeller_modeling:
    def __init__(self, target_name: str, model_num: int = 5):
        """Modeling class for modeller.

        Args:
            target_name (str): The name of the target.
            model_num (int, optional): Number of models to generate for a single pir file. Defaults to 5.
        """
        self.target_name = target_name
        self.model_num = model_num

    def modeling(self, pir_path: str, template_pdb_path: str, out_pdb_dir: str) -> None:
        """Modeling with modeller.

        Args:
            pir_path (str): Path to the pir file.
            template_pdb_path (str): Path to the template pdb directory.
            out_pdb_dir (str): Path to the output directory.
        """
        pir_path = Path(pir_path).resolve()
        template_pdb_path = Path(template_pdb_path)
        template_pdb_dir = template_pdb_path.parent.resolve()
        template_pdb_stem = template_pdb_path.stem
        out_pdb_dir = Path(out_pdb_dir).resolve()
        cwd = Path.cwd()
        with tempfile.TemporaryDirectory() as dname:
            try:
                os.chdir(dname)
                f = io.StringIO()
                with redirect_stdout(f):
                    self.run_modeller(pir_path, template_pdb_dir, template_pdb_stem)
            except ModellerError as e:
                print(e)
                print(f.getvalue())
            else:
                self.mv_pdb_files(pir_path.stem, out_pdb_dir)
            finally:
                os.chdir(cwd)

    def run_modeller(self, pir_path: Path, template_pdb_dir: Path, template_pdb_stem: str) -> None:
        """Running modeller.

        Args:
            pir_path (Path): Path to the pir file.
            template_pdb_dir (Path): Path to the template pdb directory.
            template_pdb_stem (str): File name (without suffix) of the template structure.
        """
        log.verbose()
        env = environ()
        env.io.hetatm = False
        env.io.atom_files_directory = [str(template_pdb_dir)]
        print(template_pdb_stem)

        a = automodel(env,
                      alnfile=str(pir_path),
                      knowns=template_pdb_stem,
                      sequence='TARGET')
        a.starting_model = 1
        a.ending_model = self.model_num
        a.make()

    def mv_pdb_files(self, template_name: str, out_dir: str) -> None:
        """Move generated model structures from temporary direcotry.

        Args:
            template_name (str): Name of the template structure.
            out_dir (str): Path to the output directory.
        """
        for pdb in Path('.').glob('TARGET*.pdb'):
            out_filename = self.target_name + '_' + template_name + '_' + pdb.stem[-1] + '.pdb'
            shutil.move(pdb, Path(out_dir) / out_filename)

    def modeling_dir(self, pir_dir: str, template_pdb_dir: str, out_pdb_dir: str) -> None:
        """Modeling for all pir files under a directory.

        Args:
            pir_dir (str): Path to the directory containing pir files.
            template_pdb_dir (str): Path to the template pdb directory.
            out_pdb_dir (str): Path to the output directory.
        """
        for pir_path in Path(pir_dir).glob('*.pir'):
            out_filename = self.target_name + '_' + pir_path.stem + '_' + str(self.model_num) + '.pdb'
            if not (Path(out_pdb_dir) / out_filename).exists():
                print(pir_path.stem)
                template_pdb_path = (Path(template_pdb_dir) / pir_path.stem.rsplit('_', 1)[0]).with_suffix('.pdb')
                if template_pdb_path.exists():
                    print('template exists at:', template_pdb_path)
                    self.modeling(pir_path, template_pdb_path, out_pdb_dir)
                else:
                    template_pdb_path = (Path(template_pdb_dir) / pir_path.stem[: 4]).with_suffix('.pdb')
                    print('template exists at:', template_pdb_path)
                    self.modeling(pir_path, template_pdb_path, out_pdb_dir)
