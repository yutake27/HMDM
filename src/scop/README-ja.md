# SCOP dataset 作成手順

1. target select

    ターゲットの合計数を100とする．

    SCOPのスーパーファミリーのエントリー数が多いものからターゲットを選択していく．

    スーパーファミリー上位100件を選択するのではなく，クラスごとに25件ずつ選択する．

    まずスーパーファミリーごとにpsiblastを行いヒット数を確認する．

    ```bash
    $ cd src/scop/target_select
    $ python run_psiblast4sf.py
    ```

    スーパーファミリーごとにヒット数が最も多いエントリーをターゲットとして選択する．

    ```bash
    $ cd src/scop/target_select
    $ python select_target.py
    ```

2. modeling

    ```bash
    $ cd src
    $ python run_scop_modeling.py
    ```

3. concat tmscore

    ```bash
    $ cd src
    $ python concat_tmscore.py ../tmscore/dataset_name
    ```
4. sampling

    ```bash
    $ cd src
    $ python sampling_model.py ../pdb/dataset_name
    ```

5. replace target

    ターゲットごとの予測結果を確認して基準を満たしていないターゲットは入れ替える．

    ```bash
    # replace
    $ cd src/scop/target_select
    $ python replace_target.py
    # modeling
    $ cd src
    $ python scop_modeling.py
    # concatenate
    $ python concat_tmscore.py ../tmscore/dataset_name
    # sampling
    $ python sampling_model.py ../pdb/dataset_name
    ```
    `replace_target.py`で入れ替えるターゲットがなくなるまでこれを繰り返す．
