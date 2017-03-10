# DepictStructBioinf

![alt_tag](https://github.com/j-wags/DepictStructBioinf/blob/master/test_region_image.png)

DepictStructBioinf allows researchers of a common protein domain to quickly aggregate knowledge, record findings in a relevant format, and compare hypotheses. By default, this package comes loaded with information about the APOBEC3 family of proteins. 

## Quick start (for Mac/Linux)

Ensure [Pymol](https://sourceforge.net/projects/pymol/) and [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) are installed on your system, then in a terminal, run:

```
git clone https://github.com/j-wags/DepictStructBioinf.git
pip install biopython numpy mdanalysis
cd DepictStructBioinf
python depict.py -sf a3b_ctd
pymol visualize.pml
```
Depending on your system, the "pip" command may require a "sudo" at the beginning of the line. 

Mac notes:
* The pip-affiliated python may be in ```/usr/local/bin/python``` instead of ```/usr/bin/python```
* The Pymol executable may be called "MacPyMol" on the command line (instead of just "pymol")
* I've received reports that biopython and mdanalysis installations are giving errors on some platforms. I don't yet have a machine that can reproduce these errors.

## Project Cassettes
A "project cassette" is a directory containing:

### aligned_structures
A directory containing all relevant protein structures to this project. Pymol will load structures directly from here, so be sure that any cleaning/3D alignment is already performed before the structures get here. For multi-chain proteins, I recommend splitting the chains into separate files and using a naming scheme like 3o02_A and 3o02_B

### Annotations.py
Provides the dictionary "annotations". The core of each annotation is a text label, parent sequence (eg. a3b_ctd), and list of sequence positions. For example

```
annotations = {}
annotations['asp_314_in_a3b_ctd'] = {'parent':'a3b_ctd',
                                     'label':'Asp314 in A3Bctd',
                                     'posns':[i-192 for i in [314]], 
                                     'checkAA':'D',
                                     'help':'Aspartate which has been shown to be important for -1 base specificity.'
                                  }


annotations['zn_coordinating'] = {'parent':'a3g_ctd',
                                  'label':'Zinc-coordinating amino acids',
                                  'posns':[i-197 for i in [288,257,291]],
                                  'checkAA':'HCC',
                                  'help':'Zinc-coordinating amino acids. '
                                  }
```

The keys in this dictionary (ie. test1) are the names that will be presented to the user as annotation options when they type ```python depict.py -h```.

The key-value pairs in the annotation dictionary are formatted as follows

- **parent**: The name of a sequence in the sequence alignment that this annotation belongs to.
- **label**: The label that will appear as floating text in the pymol visualization.
- **posns**: The sequence position(s) that this annotation refers to. Note: Indexing begins at 0, is for each domain (not the entire protein!) and ignores gaps.
- **checkAA**: The single-letter amino acid codes that SHOULD line up to your sequence. Used for internal checking (**TODO: implement this as a regular test**).
- **help (optional)**: A help message for other users. This currently isn't used by the program.

Note - Negative indices in "posns" will be ignored, as pymol interprets negative numbers as the end of a range (eg. -27 is interpreted as 0-27). 



### Families.py
Provides the dictionary "families". Structure families can be helpful if you want to quickly load groups of proteins with common characteristics. For example
```
families = {}
families['a3a'] = ['4xxo_A','4xxo_B','2m65_A','5sww_AE','5sww_BF','5sww_CG','5sww_DH','5td5_AC']
families['dna_bound'] = ['5sww_AE','5sww_BF','5sww_CG','5sww_DH','5td5_AC','5k83_AH','5k83_CI','5k83_EJ']

```


The keys in this dictionary (ie. a3a) are presented to the user as choices under the "-sf" argument when the user runs "python depict.py -h".


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
