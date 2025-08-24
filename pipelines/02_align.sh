hisat2-build -p 16 modified_sorted_XL280alpha_sequence_rename+MATa.fasta XL280_merged &&
hisat2-build -p 16 XL280_amasked.fasta XL280_amasked &&
hisat2-build -p 16 XL280_alphamasked.fasta XL280_alphamasked &&

# MATalpha solo
for f in JH-V8_72h_XL280alpha_rep*R1_trimmed* JH-V8_72h_JHG4_rep*R1_trimmed.fastq*
do
i="${f%_*R1_trimmed.fastq*}"
hisat2 -q -x XL280_amasked \
-1 "$i"_R1_trimmed.fastq* \
-2 "$i"_R2_trimmed.fastq* \
--summary-file "$i"_HISAT2_summary.txt | \
samtools sort -o "$i"_HISAT2_aln.bam
done

# MATa solo
for f in JH-V8_72h_XL280a_rep*R1_trimmed* JH-V8_72h_JHG6_rep*R1_trimmed.fastq*
do
i="${f%_*R1_trimmed.fastq*}"
hisat2 -q -x XL280_alphamasked \
-1 "$i"_R1_trimmed.fastq* \
-2 "$i"_R2_trimmed.fastq* \
--summary-file "$i"_HISAT2_summary.txt | \
samtools sort -o "$i"_HISAT2_aln.bam
done

# crosses 
for f in JH-V8_72h_XL280alpha-a*R1_trimmed* JH-V8_72h_XL280alpha-JHG6*R1_trimmed.fastq* JH-V8_72h_JHG4-XL280a*R1_trimmed.fastq* JH-V8_72h_JHG4-JHG6*R1_trimmed.fastq*
do
i="${f%_*R1_trimmed.fastq*}"
hisat2 -q -x XL280_merged \
-1 "$i"_R1_trimmed.fastq* \
-2 "$i"_R2_trimmed.fastq* \
--summary-file "$i"_HISAT2_summary.txt | \
samtools sort -o "$i"_HISAT2_aln.bam
done
