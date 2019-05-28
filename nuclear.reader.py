from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import pysam
import csv

#these are the A4 single nuclear reads, they can be reordered to give the desired output format ordering
samples = ["SN03","SN09","SN19","SN02","SN11","SN12","SN20","SN21","SN22","SN24","SN26","SN13","SN14","SN16","SN27"]

print "there are " + str(len(samples)) + " samples to be analyzed"

#now we need to read in the reference to check the alleles from supp6, with the A4 page
reference = pysam.Fastafile("/mnt/scratch/auxie001/rhizo/genomes/A4/A4_spades.fa")
print "now have reference sequence"

#now loading the .bam alignment
bamfile = pysam.AlignmentFile("A4.combined.bam","rb")

with open("supp6.A4.raw.csv") as A4csv:
	with open("A4_recomb_depths.csv","wb") as csv_writer:
		depth_writer = csv.writer(csv_writer,delimiter=",")
		allele_reader = csv.reader(A4csv,delimiter=",")
		depth_writer.writerow(["scaffold","position","individual","A","C","T","G","N","blast1","blast2","blast3","blast4","blast5"])
		for row in allele_reader:
			if row[0][0] != "#":
				contig = row[0]
				blast1_count, blast2_count, blast3_count, blast4_count, blast5_count = 0,0,0,0,0
				try:
					region = reference.fetch(contig,int(row[2])-100,int(row[2])+100)
					temp_fasta = open("temp.fasta", "wb")
					temp_fasta.write(region)
					temp_fasta.close()
					blastn_cline = NcbiblastnCommandline(cmd="blastn", out="temp.xml", outfmt=5, query = "temp.fasta", db="A4_spades.fa", num_threads = 4, evalue=0.001)
					stdout, stderr = blastn_cline()
					result_handle = open("temp.xml","rb")
					blast_records = NCBIXML.parse(result_handle)
					for blast in blast_records:
						for hit in blast.alignments:
							blast1_count += 1
				except ValueError:
					blast1_count = "NA"

				try:
                                        region = reference.fetch(contig,int(row[2])-250,int(row[2])+250)
                                        temp_fasta = open("temp.fasta", "wb")
                                        temp_fasta.write(region)
                                        temp_fasta.close()
                                        blastn_cline = NcbiblastnCommandline(cmd="blastn", out="temp.xml", outfmt=5, query = "temp.fasta", db="A4_spades.fa", num_threads = 4, evalue=0.001)
                                        stdout, stderr = blastn_cline()
                                        result_handle = open("temp.xml","rb")
                                        blast_records = NCBIXML.parse(result_handle)
                                        for blast in blast_records:
                                                for hit in blast.alignments:
                                                        blast2_count += 1
				except ValueError:
					blast2_count = "NA"
                                try:
                                        region = reference.fetch(contig,int(row[2])-375,int(row[2])+375)
                                        temp_fasta = open("temp.fasta", "wb")
                                        temp_fasta.write(region)
                                        temp_fasta.close()
                                        blastn_cline = NcbiblastnCommandline(cmd="blastn", out="temp.xml", outfmt=5, query = "temp.fasta", db="A4_spades.fa", num_threads = 4, evalue=0.001)
                                        stdout, stderr = blastn_cline()
                                        result_handle = open("temp.xml","rb")
                                        blast_records = NCBIXML.parse(result_handle)
                                        for blast in blast_records:
                                                for hit in blast.alignments:
                                                        blast3_count += 1
                                except ValueError:
                                        blast3_count = "NA"
                                try:
                                        region = reference.fetch(contig,int(row[2])-500,int(row[2])+500)
                                        temp_fasta = open("temp.fasta", "wb")
                                        temp_fasta.write(region)
                                        temp_fasta.close()
                                        blastn_cline = NcbiblastnCommandline(cmd="blastn", out="temp.xml", outfmt=5, query = "temp.fasta", db="A4_spades.fa", num_threads = 4, evalue=0.001)
                                        stdout, stderr = blastn_cline()
                                        result_handle = open("temp.xml","rb")
                                        blast_records = NCBIXML.parse(result_handle)
                                        for blast in blast_records:
                                                for hit in blast.alignments:
                                                        blast4_count += 1
                                except ValueError:
                                        blast4_count = "NA"
                                try:
                                        region = reference.fetch(contig,int(row[2])-1000,int(row[2])+1000)
                                        temp_fasta = open("temp.fasta", "wb")
                                        temp_fasta.write(region)
                                        temp_fasta.close()
                                        blastn_cline = NcbiblastnCommandline(cmd="blastn", out="temp.xml", outfmt=5, query = "temp.fasta", db="A4_spades.fa", num_threads = 4, evalue=0.001)
                                        stdout, stderr = blastn_cline()
                                        result_handle = open("temp.xml","rb")
                                        blast_records = NCBIXML.parse(result_handle)
                                        for blast in blast_records:
                                                for hit in blast.alignments:
                                                        blast5_count += 1
                                except ValueError:
                                        blast5_count = "NA"
				print row[0], contig, row[2], int(row[2]) + 1, "blast count", blast1_count,blast2_count,blast3_count,blast4_count,blast5_count
				print "gives: ", row[1], reference.fetch(contig,int(row[2])-1,int(row[2]))
				for pileupcolumn in bamfile.pileup(contig,int(row[2]),int(row[2])+1):
					if pileupcolumn.pos == int(row[2])-1:
						d1 = dict()
						for i in samples: #ignore last one since nucleus 13 stripped from sample
							#counts ACGT in that order
							d1[i] = [0,0,0,0,0]
						for pileupread in pileupcolumn.pileups:
							#print pileupread.query_position
							id = pileupread.alignment.get_tag("RG")
							try:
								base = pileupread.alignment.query_sequence[pileupread.query_position]
							except TypeError:
								print "missing base"
								base = "N"
								d1[id][4] += 1
							#print id,base,
							if base.upper() == "A":
								d1[id][0] += 1
							if base.upper() == "C":
								d1[id][1] += 1
							if base.upper() == "G":
								d1[id][2] += 1
							if base.upper() == "T":
								d1[id][3] += 1
				for id in samples:
					#print id, d1[id]
					string = [contig,row[2],id] + d1[id] + [blast1_count,blast2_count,blast3_count,blast4_count,blast5_count]
					#print(string)
					depth_writer.writerow(string)
				print ""
				#raw_input("Press Enter")
reference.close()
bamfile.close()
