prefix="n/p1_u500_n4"

module load samtools

seqkit faidx $prefix.genome.fa g0:301-800 > $prefix.unit.fa
samtools faidx $prefix.unit.fa

minimap2 -cx asm20 --cs -a $prefix.unit.fa $prefix.reads.fa > $prefix.reads.sam
samtools view -bS $prefix.reads.sam > $prefix.reads.bam
samtools sort $prefix.reads.bam > $prefix.reads.sorted.bam
samtools index $prefix.reads.sorted.bam

minimap2 -cx asm20 --cs -a $prefix.unit.fa $prefix.genome.fa > $prefix.genome.sam
samtools view -bS $prefix.genome.sam > $prefix.genome.bam
samtools sort $prefix.genome.bam > $prefix.genome.sorted.bam
samtools index $prefix.genome.sorted.bam

module unload samtools

echo Generated:
echo $prefix.reads.sorted.bam
echo $prefix.genome.sorted.bam
echo $prefix.unit.fa

echo Next:
echo Load unit.fa as genome and sorted.bam as alignment in IGV.
echo Turn off Quick consensus mode and Hide small indels. Turn on Show mismatched bases and use collapsed view.

# prefix="n/p1_u500_n6"
# seqkit faidx $prefix.genome.fa g0 > $prefix.genome.uniq.fa
# cat $prefix.genome.fa > $prefix.genome.uniq.fa
# seqkit faidx $prefix.genome.fa g0:301-800 > $prefix.genome.uniq.fa
# minimap2 -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.paf
# minimap2 -x ava-pb --dual=yes -c --cs=long $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.paf
# maf or lastz-cigar
# paftools.js view -f maf -l 1000 $prefix.genome.uniq.paf

# module load samtools
# # minimap2 -a -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.sam
# minimap2 -x ava-pb --dual=yes -a $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.sam
# minimap2 -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.paf
# paftools.js view $prefix.genome.uniq.paf

# module load samtools
# minimap2 -a -c --cs $prefix.genome.uniq.fa $prefix.reads.fa > $prefix.genome.uniq.sam
# samtools view -bS $prefix.genome.uniq.sam > $prefix.genome.uniq.bam
# samtools sort $prefix.genome.uniq.bam > $prefix.genome.uniq.sorted.bam
# samtools faidx $prefix.genome.uniq.fa
# samtools index $prefix.genome.uniq.sorted.bam
# module unload samtools
# echo $prefix.genome.uniq.sorted.bam
# echo $prefix.genome.uniq.fa
