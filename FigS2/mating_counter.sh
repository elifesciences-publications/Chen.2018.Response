echo "Reads mapped against scaffold511:" > A4_mating_types.txt
samtools view    A4/A4.combined.bam scaffold511:7128-10167 | sed "s/^.*RG:Z://m" | sort | uniq -c | sort -nr >> A4_mating_types.txt
samtools view -h A4/A4.combined.bam scaffold511:7128-10167 > A4_511.sam
echo "Reads mapped against scaffold5308:" >> A4_mating_types.txt
samtools view    A4/A4.combined.bam scaffold5308:1-2739 | sed "s/^.*RG:Z://m" | sort | uniq -c | sort -nr >> A4_mating_types.txt
samtools view -h A4/A4.combined.bam scaffold5308:1-2739 > A4_5307.sam


echo "Reads mapped against scaffold4735:" > A5_mating_types.txt
samtools view    A5/A5.combined.bam scaffold4735:1-2739 | sed "s/^.*RG:Z://m" | sort | uniq -c | sort -nr >> A5_mating_types.txt
samtools view -h A5/A5.combined.bam scaffold4735:1-2739 > A5_4735.sam
echo "Reads mapped against scaffold2225:" >> A5_mating_types.txt
samtools view    A5/A5.combined.bam scaffold2225:11810-14510 | sed "s/^.*RG:Z://m" | sort | uniq -c | sort -nr >> A5_mating_types.txt
samtools view -h A5/A5.combined.bam scaffold2225:11810-14510 > A5_2225.sam


echo "Reads mapped against scaffold11994:" > SL1_mating_types.txt
samtools view    SL1_spades/SL1.combined.bam scaffold_11994:1-2690 | sed "s/^.*RG:Z://m" | sort | uniq -c | sort -nr >> SL1_mating_types.txt
samtools view -h SL1_spades/SL1.combined.bam scaffold_11994:1-2690 > SL1_11994.sam
echo "Reads mapped against scaffold11899:" >> SL1_mating_types.txt
samtools view    SL1_spades/SL1.combined.bam scaffold_11899:1-2719 | sed "s/^.*RG:Z://m" | sort | uniq -c | sort -nr >> SL1_mating_types.txt
samtools view -h SL1_spades/SL1.combined.bam scaffold_11899:1-2719 > SL1_11899.sam
