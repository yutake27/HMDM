"""
    combining datasets with no template filtering and datasets with template filtering and sampling from combined.
"""

import argparse
import shutil
import warnings
from pathlib import Path

import pandas as pd
from prody import parsePDB, LOGGER
from prody.proteins.pdbfile import PDBParseError

LOGGER.verbosity = 'none'

warnings.simplefilter('ignore', UserWarning)


def sample_not_exceed_df_length(df: pd.DataFrame, n: int, seed: int) -> pd.Series:
    num_sample = n if len(df) > n else len(df)
    sample_index = df.sample(n=num_sample, random_state=seed).index
    return sample_index


def sample_restrict_around_best_model(df: pd.DataFrame, seed: int = 0,
                                      num_models_per_target: int = 150,
                                      gdtts_thre_from_best_model: float = 0.03,
                                      num_select_around_best_model: int = 10,
                                      gdtts_minimum_threshold: float = 0.4) -> pd.DataFrame:
    fil_gdtts_df = df.query('GDT_TS >= @gdtts_minimum_threshold')  # filtering models with GDT_TS lower than threshold
    if len(fil_gdtts_df) == 0: # if empty
        if 'SF-DOMID' in df.columns:
            target_name, sf_domid = df.iloc[0][['target', 'SF-DOMID']]
            dummy_df = fil_gdtts_df.append({'target': target_name, 'SF-DOMID': sf_domid}, ignore_index=True)
            print(dummy_df)
            return dummy_df
        else:
            return None
    selected_index_list = []
    best_gdtts, best_model = fil_gdtts_df['GDT_TS'].max(), fil_gdtts_df['GDT_TS'].idxmax()
    selected_index_list.append(best_model)
    df_around_best_model = fil_gdtts_df.query('GDT_TS >= @best_gdtts - @gdtts_thre_from_best_model and index != @best_model')
    selected_index = sample_not_exceed_df_length(df_around_best_model, num_select_around_best_model, seed)
    selected_index_list.extend(selected_index)
    df_away_from_best_model = fil_gdtts_df.query('GDT_TS < @best_gdtts - @gdtts_thre_from_best_model')
    num_models_selected = len(selected_index_list)
    num_models_should_select = num_models_per_target - num_models_selected
    selected_index = sample_not_exceed_df_length(df_away_from_best_model, num_models_should_select, seed)
    selected_index_list.extend(selected_index)
    selected_df = fil_gdtts_df.query('index in @selected_index_list')
    return selected_df


def sample(df_concat: pd.DataFrame):
    sample_df_restrict_list = []
    num_models_per_target = 150
    gdtts_thre = 0.03
    num_select_around_the_best = 10
    for name, target_df in df_concat.groupby('target'):
        sample_df_restrict\
            = sample_restrict_around_best_model(target_df, gdtts_thre_from_best_model=gdtts_thre,
                                                num_models_per_target=num_models_per_target,
                                                num_select_around_best_model=num_select_around_the_best)
        sample_df_restrict_list.append(sample_df_restrict)
    df_sample = pd.concat(sample_df_restrict_list)
    return df_sample


def copy_file(from_path: str, to_path: str):
    shutil.copy(from_path, to_path)


def check_and_copy_pdb(dataset_dir: Path, df_sample: pd.DataFrame):
    broken_pdb_list = []
    for target_name, target_df in df_sample.groupby('target'):
        print('Target:', target_name)
        target_df = target_df.dropna()
        target_dir = dataset_dir / target_name / 'sampling'
        target_dir.mkdir(parents=True, exist_ok=True)
        for index, row in target_df.iterrows():
            from_path = (dataset_dir / target_name / row['Method'] / row['model']).with_suffix('.pdb')
            to_path = (target_dir / row['model']).with_suffix('.pdb')
            if not to_path.exists():
                try:
                    parsePDB(str(from_path))
                except PDBParseError:
                    print(from_path)
                    broken_pdb_list.append(row['model'])
                except OSError:
                    shutil.rmtree(str(target_dir.parent))
                    break
                else:
                    copy_file(str(from_path), str(to_path))
    return broken_pdb_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('score_csv_path', type=str,
                        help='Path to the csv containing GDT_TS of the entire dataset.\
                              Ex) ../score/dataset/dataset.csv')
    args = parser.parse_args()
    score_csv_path = Path(args.score_csv_path)
    dataset_name = score_csv_path.stem

    df_concat = pd.read_csv(score_csv_path, index_col=0).reset_index().rename(
        columns={'index': 'model'}).drop_duplicates(subset=['model', 'GDT_TS'])
    df_sample = sample(df_concat)

    pdb_dir = Path('../pdb')
    dataset_dir = pdb_dir / dataset_name
    broken_pdb_list = check_and_copy_pdb(dataset_dir, df_sample)

    output_csv_path = (score_csv_path.parent / (dataset_name + '_sampling')).with_suffix('.csv')
    df_sample_final = df_sample.query('model not in @broken_pdb_list')
    df_sample_final.set_index('model', drop='True').to_csv(output_csv_path)


if __name__ == '__main__':
    main()

