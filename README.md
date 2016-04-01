# DepictStructBioinf

This script aids in the visualization of aggregated structural bioinformatics information. Like a tape player, it's a general framework that can take "project casettes", containing relevant sequence alignments, 3D structures, and functional annotations. Its main purpose is to help compare functional annotations within families of proteins, thereby enabling researchers to quickly aggregate knowledge, take project-specific notes in a relevant context, and maintain a high level of domain knowledge on multiple systems.

## Project Casettes
A "project casette" is a directory containing:

### aligned_structures
A directory containing all relevant protein structures to this project. Pymol will load structures directly from here, so be sure that any cleaning/3D alignment is already performed before the structures get here. For multi-chain proteins, I recommend splitting the chains into separate files and using a naming scheme like 3o02_A and 3o02_B

### Annotations.py
Provides the dictionary "annotations". For example

```
annotations = {}
annotations['test1'] = {'parent':'a3a',
                       'label':'Test Region',
                       'posns':[9,10,11,12,13,14,15],
                       'checkAA':'RHLMDPH',
                       'help':'An arbitrary region of the human A3A sequence.'
                       }
annotations['test2'] = {'parent':'a3b_ntd',
                       'label':'Another Test Region',
                       'posns':[0,1,2,3,4,5],
                       'checkAA':'MNPQIR',
                       'help':'Another arbitrary region, this time in the A3B N terminal domain'
                       }

...
```

The keys in this dictionary (ie. test1) are the names that will be presented to the user as annotation options in the help message.
The subdictionary key-value pairs are used as follows

- **parent**: The name of a sequence in the sequence alignment that this annotation belongs to.
- **label**: The label that will appear as floating text in the pymol visualization.
- **posns**: The sequence position(s) that this annotation refers to. Note: Indexing begins at 0, and ignores gaps.
- **checkAA**: The single-letter amino acid codes that SHOULD line up to your sequence. Used for internal checking (**TODO: implement this as a regular test**).
- **help (optional)**: A help message for other users. This currently isn't used by the program.

Note - Negative indices in "posns" will be ignored, as pymol interprets negative numbers as the end of a range (eg. -27 is interpreted as 0-27)

### Families.py
Provides the dictionary "annotations". For example
```
families = {}
families['all_a3a'] = ['4xxo_A','4xxo_B','2m65_A']
...
```


The keys in this dictionary (ie. all_a3a) are presented to the user as choices for loading and annotating in the pymol script. Structure families can be helpful if you want to quickly load groups of proteins with common characteristics.


### SequenceAlignment.fasta
The sequence alignment in fasta format. **This alignment must include all sequences referenced in Annotations.py, and sequences corresponding to structures that you will load**. For structures to load, the sequence title in this file must be the structure name, excluding the path (assumed to be "project_casette/aligned_structures/") and exluding the suffix ".pdb"

For example:
```
a3a
MEASPASGPRHLMDPHIFTSNFNNG---IG---RHKTYLCYEVERLDNGTSVKMDQHRGF
LHNQ--------AKNLLCGFYGRHAELRFLDLVPSL-QLDPAQIYRVT--WFIS--WSPC
FSWGCAGEVRAFLQENTHVRLRIFAARIY-DY-DPLYKEALQMLRD----AGAQVSIMTY
DEFKH----CWDTFVDHQGC----PFQPWDGLDEHSQALSGRLRAILQNQGN--------
------
>a3b_ntd
-MNPQIRNPMERMYRDTFYDNFENEPILYG---RSYTWLCYEVKIKRGRSNLLWDT--GV
FRGQ------------VYFKPQYHAEMCFLSWFCGN-QLPAYKCFQIT--WFVS--WTPC
PD--CVAKLAEFLSEHPNVTLTISAARLY-YYWERDYRRALCRLSQ----AGARVTIMDY
EEFAY----CWENFVYNEGQ----QFMPWYKFDENYAFLHRTLKEILRYL----------
------
>2m65_A
MEASPASGPRHLMDPHIFTSNFNNG---IG---RHKTYLCYEVERLDNGTSVKMDQHRGF
LHNQAKN--------LLCGFYGRHAELRFLDLVPSL-QLDPAQIYRVT--WFIS--WSPC
FSWGCAGEVRAFLQENTHVRLRIFAARIY--DYDPLYKEALQMLRD----AGAQVSIMTY
DEFKH----CWDTFVDHQGC----PFQPWDGLDEHSQALSGRLRAILQNQGN--------
------
>4xxo_A
-------GPRHLMDPHIFTSNFNNG---IG---RHKTYLCYEVERLDS---VKMDQHRGF
LHNQAKN--------LLCGFYGRHAALRFLDLVPSL-QLDPAQIYRVT--WFIS--WSPC
FSWGCAGEVRAFLQENTHVRLRIFAARIY--DYDPLYKEALQMLRD----AGAQVSIMTY
DEFKH----CWDTFVDHQGA----PFQPWDGLDEHSQALSGRLRAILQN-----------
------
>4xxo_B
-------GPRHLMDPHIFTSNFNNG---IG---RHKTYLCYEVERLDN-TSVKMDQHRGF
LHNQAKN---------LLGFYGRHAALRFLDLVPSL-QLDPAQIYRVT--WFIS--WSPC
FSWGCAGEVRAFLQENTHVRLRIFAARIY--DYDPLYKEALQMLRD----AGAQVSIMTY
DEFKH----CWDTFVDHQGA----PFQPWDGLDEHSQALSGRLRAILQN-----------
------
```


##General notes
- Be very careful when extracting sequences from pdb files! Different structure-to-sequence programs will interpret gaps and variable-occupany residues differently (eg. ALYS and BLYS for placements of a lysine residue may be double counted in the sequence)
- Think carefully and use your judgement when dealing with systems involving homodimers or repeating/homologous domains.
