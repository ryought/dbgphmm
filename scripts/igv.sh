prefix="n/p01_u500_n3"
# prefix="n/p1_u500_n6"
# seqkit faidx $prefix.genome.fa g0 > $prefix.genome.uniq.fa
# cat $prefix.genome.fa > $prefix.genome.uniq.fa
seqkit faidx $prefix.genome.fa g0:301-800 > $prefix.genome.uniq.fa
# minimap2 -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.paf
minimap2 -x ava-pb --dual=yes -c --cs=long $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.paf
# maf or lastz-cigar
paftools.js view -f maf -l 1000 $prefix.genome.uniq.paf

# module load samtools
# # minimap2 -a -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.sam
# minimap2 -x ava-pb --dual=yes -a $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.sam
minimap2 -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.paf
paftools.js view $prefix.genome.uniq.paf

# module load samtools
# minimap2 -a -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.sam
# samtools view -bS $prefix.genome.uniq.sam > $prefix.genome.uniq.bam
# samtools sort $prefix.genome.uniq.bam > $prefix.genome.uniq.sorted.bam
# samtools faidx $prefix.genome.uniq.fa
# samtools index $prefix.genome.uniq.sorted.bam
# module unload samtools
# echo $prefix.genome.uniq.sorted.bam
# echo $prefix.genome.uniq.fa
