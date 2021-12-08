import argparse
import pandas as pd
from pathlib import Path

from modeller import *
from modeller.scripts import complete_pdb
from modeller import soap_protein_od


def dope_soap(pdb_dir):
    env = environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    sp = soap_protein_od.Scorer(library='./soap_protein_od.hdf5')

    model_array = []
    dope_score_array = []
    soap_score_array = []
    for pdb in Path(pdb_dir).glob('*'):
        mdl = complete_pdb(env, str(pdb))
        atmsel = selection(mdl.chains[0])

        try:
            dope_score = atmsel.assess_dope()
            soap_score = atmsel.assess(sp)
            model_array.append(pdb.stem)
            dope_score_array.append(dope_score)
            soap_score_array.append(soap_score)
        except ModellerError as e:
            print(e)

    return dope_score_array, soap_score_array, model_array


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dir', type=str,
                        help='target pdb directory of a dataset ex) ../../../pdb/target_10/1A3Q')
    args = parser.parse_args()

    target_dir = Path(args.target_dir)
    dataset_name = target_dir.parent.stem
    output_csv_path = (Path('../../../score') / dataset_name / 'dope' / target_dir.stem).with_suffix('.csv')
    output_csv_path.parent.mkdir(exist_ok=True, parents=True)
    target_sampling_dir = target_dir / 'sampling'

    dope_score, soap_score, model = dope_soap(target_sampling_dir)
    df = pd.DataFrame({'dope': dope_score, 'soap': soap_score}, index=model)
    df.to_csv(output_csv_path)
