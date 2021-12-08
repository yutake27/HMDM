# HMDM-multi

Created on 2021-10-15

This dataset is licensed under CC BY 4.0. (https://creativecommons.org/licenses/by/4.0/)


## Dataset structure

```
Dataset
├── fasta
│   ├── target00.fasta
│   └── target01.fasta
├── native_pdb
│   ├── target00.pdb
│   └── target01.pdb
├── pdb
│   ├── target00
│   │   ├── target00_model_structure000.pdb
│   │   └── target00_model_structure001.pdb
│   └── target01
└── data
    ├── label.csv
    ├── score.csv
    └── target.csv

```


## fasta
Target sequence in fasta format


## native_pdb
Native protein structure in pdb format


## pdb
Predicted model structure in pdb format


## data

### label.csv
A csv file with the label for each model structure such as GDT_TS

Header

* Model: Model structure name
* Target: Target name (PDB ID + chain name)
* Template: Template name (PDB ID + chain name + iteration of PSI-BLAST)
* GDT_TS: GDT_TS score
* GDT_HA: GDT_HA score
* SeqLength: Sequence length
* identity: Number of residues that match within the aligned sequence
* positive: Number of residues with a positive alignment score
* coverage: Number of residues covered by the template sequence
* identity(-misres): identity without missing residue
* positive(-misres): positive without missing residue
* coverage(-misres): coverage without missing residue
* num_misres: Number of missing residue


### mqa_score.csv
A csv file with the MQA score for each model structure

Header

* Model: Model structure name
* Target: Target name
* identity(%): Percentage of identity (identity / SeqLength)
* positive(%): Percentage of positive
* coverage(%): Percentage of coverage
* identity(-misres)(%): Percentage of identity(-misres)
* positive(-misres)(%): Percentage of positive(-misres)
* coverage(-misres)(%): Percentage of coverage(-misres)
* DOPE
* SOAP
* SBROD
* ProQ2D, ProQRosCenD, ProQRosFAD, ProQ3D
* P3CMQA
* DeepAccNet, DeepAccNet-Bert


### target.csv
A csv file with the information for each target

Header

* Target: Target name
* SeqLength: Sequence length
* IDs
* Exptl.: Experimental method
* resolution
* R-factor
* FreeRvalue
* PDB_ID
* Chain
* DomainNum: Domain number in a target protein
