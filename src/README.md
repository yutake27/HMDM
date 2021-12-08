# How to make dataset

## Requirement

* python 3.6 or later
  * biopython version 1.78 or later
  * prody version 2.0 or later
  * numpy
  * pandas
  * tqdm
* MODELLER version 9.24
* PSI-BLAST
* CD-HIT version 4.8.1
* TMscore

## Single-domain (scop base)

1. Select target

    ```bash
    $ cd src/scop/target_select
    $ python run_psiblast4sf.py
    $ python select_target.py
    ```

2. Modeling

    ```bash
    $ cd src
    $ python run_scop_modeling.py
    ```

    Set target list and dataset name in `run_scop_modeling.py`.

3. Concatenate

    ```bash
    $ cd src
    $ python concat_tmscore.py ../tmscore/dataset
    ```

4. Sampling models

    ```bash
    $ cd src
    $ python sampling_model.py ../score/dataset/dataset.csv
    ```

5. Replace targets that do not meet the criteria

    ```bash
    $ cd src/scop/target_select
    $ python replace_target.py --sampling_score_csv_path ../../../score/dataset/dataset_sampling.csv
    ```

    Targets to be replaced will be shown by this operation.

    You can generate new target list using below command.

    ```bash
    $ cd src/scop/target_select
    $ python replace_target.py --sampling_score_csv_path ../../../score/dataset/dataset_sampling.csv --out_csv_path ../../../scop/version_date/new_target_list_name.csv --save
    ```

    Target list will be saved in `new_target_list_name.csv`.

6. Iterate operation 2-5 until there is no target to be replaced

## Multi-domain (pisces base)

1. select target
2. modeling
3. concatenate
4. sampling
5. select targets that meet the criteria
