# mask alpha locus
$ nano alphamat.bed
$ XL280alpha_Chr4    1522691     1656182
$ bedtools maskfasta \
-fi modified_sorted_XL280alpha_sequence_rename+MATa.fasta \
-bed alphamat.bed \
-fo XL280_alphamasked.fasta
# mask a locus
$ nano amat.bed
$ MATa  0   148067
$ bedtools maskfasta \
-fi modified_sorted_XL280alpha_sequence_rename+MATa.fasta \
-bed amat.bed \
-fo XL280_amasked.fasta
$ samtools faidx modified_sorted_XL280alpha_sequence_rename+MATa.fasta
$ samtools faidx XL280_alphamasked.fasta
$ samtools faidx XL280_amasked.fasta
