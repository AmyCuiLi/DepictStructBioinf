## It's necessary to have allSequences_aydin_alignment.aln in .aln format so it can be compared easily to the review paper alignment, however clustalw can't read it so we need to convert to .fasta format.
## This requires downloading the Seq package. 

. ~/env_00/venv/bin/activate
python Seq/scripts/aln2fasta.py allSequences_aydin_alignment.aln > allSequences_aydin_alignment.fasta
deactivate

