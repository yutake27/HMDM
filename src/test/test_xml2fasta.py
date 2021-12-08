import os
import sys


sys.path.append(os.path.abspath('..'))
import xml2fasta

xml_path = '../../blast-xml/pdbaa_20200712/1bxo_1.xml'
fasta_path = '../../blast-xml/pdbaa_20200712/1bxo_1.fasta'
xml2fasta.xml2fasta(xml_path, fasta_path)

xml_path = '../../blast-xml/pdbaa_20200712/4gg1_1.xml'
fasta_path = '../../blast-xml/pdbaa_20200712/4gg1_1.fasta'
xml2fasta.xml2fasta(xml_path, fasta_path)