import argparse
import pandas as pd
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tmscore_dir', type=str, help='Path to the directory of tmscore. ex) ../tmscore/')
    args = parser.parse_args()
    
    tmscore_dir = Path(args.tmscore_dir)

    tmscore_array = []
    for csv in tmscore_dir.glob('*.csv'):
        df = pd.read_csv(csv, index_col=0)
        tmscore_array.append(df)
    tmscore_df = pd.concat(tmscore_array)

    score_dir = Path('../score') / tmscore_dir.stem
    score_dir.mkdir(exist_ok=True)
    out_path = (score_dir / tmscore_dir.stem).with_suffix('.csv')
    tmscore_df.to_csv(out_path)
