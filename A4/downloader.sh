cd /mnt/scratch/auxie001/rhizo/genomes/

cd /mnt/scratch/auxie001/rhizo/genomes/A4_singles


#this downloads the 14 single-nuclei datasets for A4
echo "starting"
for i in {51..65}
do wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR741/SRR74164$i/SRR74164$i.sra &&
/mnt/scratch/auxie001/programs/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump --split-3 --gzip SRR74164$i.sra
done

#This downloads the 3 libraries of bulk illumina data for this isolate
echo "starting"
for i in 37 39 40
do wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR272/SRR27261$i/SRR27261$i.sra &&
/mnt/scratch/auxie001/programs/sratoolkit.2.9.2-ubuntu64/bin/fastq-dump --split-3 --gzip SRR27261$i.sra
done
