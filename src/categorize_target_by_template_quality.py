import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd


def categorize_by_label_distribution(group: pd.DataFrame,
                                     label: str,
                                     dif_threshold: float = 0.1,
                                     top_percent: float = 0.2) -> str:
    """Classify category that based on distribution of the label

    Args:
        group (pd.DataFrame): DataFrame of each target.
        label (str): The name of the label. e.g. identity(%), coverage(%).
        dif_threshold (float, optional): Threshold of top group. Defaults to 0.1.
        top_percent (float, optional): Threshold of percent of top group. Defaults to 0.2.

    Returns:
        str : The name of the category.
    """
    identity_sorted = list(group[label].sort_values())
    max_identity = identity_sorted[-1]
    second_identity = identity_sorted[-2]
    n_ten_percent = int(top_percent * len(group))
    ten_percent_identity = identity_sorted[-n_ten_percent]
    dif = max_identity - second_identity
    NAME_CATEGORY = ['top', 'multi top', 'the others']
    if dif > dif_threshold:
        return NAME_CATEGORY[0]
    elif max_identity - ten_percent_identity > dif_threshold:
        return NAME_CATEGORY[1]
    else:
        return NAME_CATEGORY[2]


def categorize_by_max_quality_template(group: pd.DataFrame,
                                       label: str = 'identity(%)',
                                       threshold_list: List = [0.4, 0.6, 0.8]) -> str:
    """Classify category that based on maximum value of the label

    Args:
        group (pd.DataFrame): DataFrame of each target.
        label (str, optional): The name of the label. Defaults to 'identity(%)'.
        thresholds_list (List, optional): Threshold of label for each category. Defaults to [0.5, 0.8].

    Returns:
        str: The name of the category.
    """
    max_quality = group[label].max()
    if len(threshold_list) == 2:
        NAME_CATEGORY = ['Low', 'Middle', 'High']
        if max_quality < threshold_list[0]:
            return NAME_CATEGORY[0]
        elif max_quality < threshold_list[1]:
            return NAME_CATEGORY[1]
        else:
            return NAME_CATEGORY[2]
    elif len(threshold_list) == 3:
        NAME_CATEGORY = ['Low', 'Mid-low', 'Mid-high', 'High']
        if max_quality < threshold_list[0]:
            return NAME_CATEGORY[0]
        elif max_quality < threshold_list[1]:
            return NAME_CATEGORY[1]
        elif max_quality < threshold_list[2]:
            return NAME_CATEGORY[2]
        else:
            return NAME_CATEGORY[3]
    else:
        raise ValueError('threshold_list should be 2 or 3')


def categorize_by_label(template_df: pd.DataFrame,
                        label: str,
                        max_value_threshold_list: List,
                        top_group_dif_threshold: float = 0.1,
                        top_group_percent: float = 0.2) -> pd.DataFrame:
    """Categorize target by the specified label.

    Args:
        template_df (pd.DataFrame): The DataFrame of the template.
        label (str): The name of the label. e.g. "identity(%)"
        top_group_dif_threshold (float):  Threshold of top group that based on distribution.
        top_group_percent (float): Threshold of percent of top group that based on distribution.
        max_value_threshold_list (List): Threshold of label for each category (Category of max Quality).

    Returns:
        pd.DataFrame: DataFrame of category.
    """
    label_extract_percent = label[: -3]
    # Categorize by distribution of the label
    dist_category = template_df.groupby('target').apply(
        lambda x: categorize_by_label_distribution(x,
                                                   label=label,
                                                   dif_threshold=top_group_dif_threshold,
                                                   top_percent=top_group_percent))
    dist_category_df = pd.DataFrame(dist_category).rename({0: label_extract_percent + '_dist_category'}, axis=1)
    # if label related to 'coverage', return only dist category
    if 'coverage' in label:
        return dist_category_df
    # Categorize by maximum value of the label
    quality_category = template_df.groupby('target').apply(
        lambda x: categorize_by_max_quality_template(x, label, threshold_list=max_value_threshold_list))
    quality_category_df = pd.DataFrame(quality_category).rename(
        {0: label_extract_percent + '_quality_category'}, axis=1)
    return pd.concat([dist_category_df, quality_category_df], axis=1)


def categorize_target(tmscore_df: pd.DataFrame) -> pd.DataFrame:
    """Categorize target by alignment quality of template (identity, positive, and coverage)

    Args:
        tmscore_df (pd.DataFrame): DataFrame of tmscore.

    Returns:
        pd.DataFrame: DataFrame of categorized target.
    """
    alignment_quality_row_columns = ['identity', 'positive', 'coverage']
    # alignment_quality_row_columns = ['identity', 'positive', 'coverage',
    #                                  'identity(-misres)', 'positive(-misres)', 'coverage(-misres)']
    alignment_quality_columns = [c + '(%)' for c in alignment_quality_row_columns]
    for c in alignment_quality_row_columns:
        tmscore_df[c + '(%)'] = tmscore_df[c] / tmscore_df['seq_len']
    template_df = tmscore_df.groupby(['target', 'template']).head(1).drop(['GDT_TS', 'GDT_HA'], axis=1)
    # Categorize by each label
    labels = alignment_quality_columns
    category_df_list = []
    for label in labels:
        if 'identity' in label:
            max_value_threshold_list = [0.4, 0.6, 0.8]
        elif 'positive' in label:
            max_value_threshold_list = [0.6, 0.8]
        else:
            max_value_threshold_list = [0.9, 0.95]
        label_category_df = categorize_by_label(template_df, label,
                                                max_value_threshold_list=max_value_threshold_list)
        category_df_list.append(label_category_df)
    category_df = pd.concat(category_df_list, axis=1)
    return category_df


def get_target_stat(target_df, label='GDT_TS'):
    d = {'Num models': len(target_df), 'mean ' + label: target_df[label].mean(), 'max ' + label: target_df[label].max(),
         'med ' + label: target_df[label].median(), 'min ' + label: target_df[label].min(),
         'var ' + label: np.var(target_df[label]), 'std ' + label: np.std(target_df[label])}
    return pd.Series(d)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('tmscore', type=str, help='path to the csv of tmscore')
    args = parser.parse_args()

    tmscore_path = Path(args.tmscore)
    tmscore_df = pd.read_csv(tmscore_path, index_col=0)
    # create category
    category_df = categorize_target(tmscore_df)
    print(category_df)
    # target info
    need_columns = ['target', 'seq_len', 'Class', 'Domain_num']
    drop_columns = set(tmscore_df.columns) - set(need_columns)
    target_df = tmscore_df.groupby('target').head(1).drop(drop_columns, axis=1)
    target_cdf = pd.merge(target_df, category_df, on='target')
    print(target_cdf)
    print(target_cdf.columns)
    # aggregate tmscore
    agg_target_df = tmscore_df.groupby('target').apply(get_target_stat)
    final_df = pd.merge(target_cdf, agg_target_df, on='target')
    int_columns = ['seq_len', 'Num models']
    for column in int_columns:
        final_df[column] = final_df[column].astype('int')
    final_df = final_df.rename({'target': 'Target', 'seq_len': 'Length'}, axis=1).set_index('Target')
    print(final_df)
    # save
    output_path = tmscore_path.parent / 'target_list.csv'
    final_df.to_csv(output_path, float_format='%.3f')


if __name__ == '__main__':
    main()
