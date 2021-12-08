import argparse
from pathlib import Path
import pandas as pd


def parse_score(txt):
    with open(txt, 'r') as f:
        lines = f.readlines()
    score = lines[1].split()
    return score


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('score_dir', type=str, help='proq3 score dir for one target')
    args = parser.parse_args()

    proq3_score_dir = Path(args.score_dir)

    model_array = []
    score_array = []
    for txt in proq3_score_dir.glob('*.global'):
        model_name = str(txt.name).split('.')[0]
        model_array.append(model_name)
        score = parse_score(txt)
        score_array.append(score)
    df = pd.DataFrame(score_array, index=model_array, columns=['ProQ2D', 'ProQRosCenD', 'ProQRosFAD', 'ProQ3D'])

    df.to_csv(proq3_score_dir.with_suffix('.csv'))
