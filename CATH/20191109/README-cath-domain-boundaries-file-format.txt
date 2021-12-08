CATH Domall File (CDF) Format 2.0
---------------------------------
The CATH Domall file describes domain boundaries for entries in the CATH database. 
All PDB chains in CATH that contain 1 or more domains have a CathDomall entry. Whole
chain domains can be identified where the number of domains is 1 and the
number of fragments is 0.

Comment lines start with a '#' character.

Segments are continuous sequence regions of domains.
Fragments are small regions of the protein chain that are excluded from the domain definition.

Column 1: Chain name (5 characters)
Column 2: Number of domains (formatted 'D%02d')
Column 3: Number of fragments (formatted 'F%02d')

The formatting of a CATH Domall file is best explained using examples.

Example CathDomall Entries
--------------------------

KEY:
N  = Number of segments
C  = Chain character
I  = Insert character/code ('-' indicates no insert character)
S  = Start PDB number
E  = End PDB number
NR = number of residues (fragment information only)

1chmA  D02 F00  1  A    2 - A  156 -  1  A  157 - A  402 -
                N |C    S I C    E I| N |C    S I C    E I|
               |<----Domain One---->|<-----Domain Two---->|
                  |<--Segment One-->|   |<--Segment One-->|

This translates to:
1chmA01 = Chain A; 2-156
1chmA02 = Chain A; 157-402

1cnsA  D02 F00  2  A    1 - A   87 -  A  146 - A  243 -  1  A   88 - A  145 -
                N |C    S I C    E I| C    S I C    E I| N |C    S I C    E I|
               |<--------------Domain One------------->|<-----Domain Two---->|
                  |<--Segment One-->|<---Segment Two-->|   |<--Segment One-->|

This translates to:
1cnsA01 = Chain A; 1-87, 146-243
1cnsA02 = Chain A; 88-145

Fragment Information
--------------------

Fragments are small regions of the protein chain that are not included
in the domain definition. These residue ranges are tagged on the end of the
segment information. The format is different from the segment range information.

1amg0  D02 F01  1  0    1 - 0  360 -  1  0  362 - 0  417 -  0  361 - 0  361 - (1)
                N |C    S I C    E I| N |C    S I C    E I| C    S I C    E I  NR|
               |<----Domain One---->|<-----Domain Two---->|<---Fragment One----->|
                  |<--Segment One-->|   |<--Segment One-->|

This translates to:
1amg001 = No chain character; 1-360
1amg002 = No chain character; 362-417
Fragment = 361

1bcmA0 D02 F02  1  A  257 - A  487 -  1  A  492 - A  559 -  A  488 - A  491 - (4)  A  560 - A  560 - (1)
                N |C    S I C    E I| N |C    S I C    E I| C    S I C    E I  NR| C    S I C    E I  NR|
               |<----Domain One---->|<-----Domain Two---->|<---Fragment One----->|<---Fragment Two----->|
                  |<--Segment One-->|   |<--Segment One-->|

This translates to:
1bcmA01 = Chain A; 257-487
1bcmA02 = Chain A; 492-559
Fragments = 488-491, 560


Cath Chain Names
----------------
The chain names have five characters (e.g. 1oaiA).

CHARACTERS 1-4: PDB Code
The first 4 characters determine the PDB code e.g. 1oai

CHARACTER 5: Chain Character
This determines which PDB chain is represented.
Chain characters of zero ('0') indicate that the PDB file has no chain field.



New to CDF Format 2.0
=====================
Now all chains that are in CATH are represented in the CATH Domall file whether
the chain has been chopped into domains or is a whole chain domain.

Whole chain domains can be identied where the number of domains is 1 and the
number of fragments is 0. 

Chain names are now represented by 5 characters instead of 6
