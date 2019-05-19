echo "starting"
cd /mnt/scratch/auxie001/rhizo/genomes/A4

for i in {51..64}
do /mnt/scratch/auxie001/programs/fastp/fastp -w 8 --stdout \
-i SRR74164$i\_1.fastq.gz -I SRR74164$i\_2.fastq.gz | \
/mnt/scratch/auxie001/programs/bwa/bwa mem -t 8 -M A4_spades.fa - > SRR74164$i.sam &&
samtools view -Sbhu SRR74164$i.sam | samtools sort -@ 8 - SRR74164$i.sorted &&
samtools rmdup SRR74164$i.sorted.bam SRR74164$i.nodup.bam &&
rm SRR74164$i.sorted.bam SRR74164$i.sam
done



for i in 37 39 40
do /mnt/scratch/auxie001/programs/fastp/fastp -w 8 --stdout \
-i SRR27261$i\_1.fastq.gz -I SRR27261$i\_2.fastq.gz | \
/mnt/scratch/auxie001/programs/bwa/bwa mem -t 8 -M -R '@RG\tID:SRR27261'"$i"'\tSM:SRR27261'"$i"'\tLB:A4_singles' \
A4_spades.fa - > SRR27261$i.sam &&
samtools view -Sbhu SRR27261$i.sam | samtools sort -@ 8 - SRR27261$i.sorted &&
samtools rmdup SRR27261$i.sorted.bam SRR27261$i.nodup.bam &&
rm SRR27261$i.sorted.bam SRR27261$i.sam
done




