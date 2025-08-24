for f in *.bam
do
i="${f%_HISAT2_aln.bam}"
featureCounts \
-a XL280_liftmerged.gff \
-o "$i"_featurecounts.txt \
-t gene \
-g ID \
-p "$f"
done
