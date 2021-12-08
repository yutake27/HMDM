!/bin/bash

# Download pdb sequences from PDB

wget https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gzip -d pdb_seqres.txt.gz