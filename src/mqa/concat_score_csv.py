import argparse
import pandas as pd
from pathlib import Path


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('score_dir', type=str, help='score directory path. ex)../../score/target_10/sbrod')
    args = parser.parse_args()

    score_dir = Path(args.score_dir)
    df_array = []
    for csv_path in score_dir.glob('*.csv'):
        df_array.append(pd.read_csv(csv_path, index_col=0))

    df = pd.concat(df_array)
    out_csv_path = score_dir.with_name(f'{score_dir.stem}.csv')
    df.to_csv(out_csv_path)
