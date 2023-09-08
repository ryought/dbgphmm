# (1) map whole
minimap2 -a -c --cs -t 16 -DP -x asm20 alts.fa alts.fa > alts.sam
python ../scripts/kir/ref_trim.py alts.sam --start 680000 --end 860000 --name LRC_KIR_Primary
# worked, but some end breakpoint is wrong (by checking dotplot)

# (2) map prefix/suffix only
seqkit subseq -r 679000:681000 LRC_KIR.fa > prefix.fa
seqkit subseq -r 859000:861000 LRC_KIR.fa > suffix.fa
minimap2 -c --cs -t 16 -DP -x asm20 alts.fa prefix.fa > prefix.paf
minimap2 -c --cs -t 16 -DP -x asm20 alts.fa suffix.fa > suffix.paf
python ../scripts/kir/lift.py prefix.sam --pos 1001
python ../scripts/kir/lift.py suffix.sam --pos 1001
# this gives alts.txt

# generate haps.fa
seqkit faidx alts.fa --region-file ../scripts/kir/alts.txt > haps.fa
# check mutations by visualizing alignment
minimap2 -a -c --cs -t 16 -x asm20 chr19.fa haps.fa > haps.sam
