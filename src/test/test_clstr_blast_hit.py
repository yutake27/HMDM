import os
import sys


sys.path.append(os.path.abspath('..'))
import clstr_blast_hit

clstr_file_path = '../../blast-xml/pdbaa_20200712/1bxo_1_cdhit.fasta.clstr'
print(clstr_blast_hit.make_df_from_clstr_file(clstr_file_path))