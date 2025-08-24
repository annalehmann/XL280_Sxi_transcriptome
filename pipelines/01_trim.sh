for f in *R1*.fastq*
do
i="${f%_S*}"
n="${f%_R1*}"
cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o "$i"_R1_trimmed.fastq.gz \
    -p "$i"_R2_trimmed.fastq.gz \
    "$n"_R1_001.fastq.gz "$n"_R2_001.fastq.gz
done
