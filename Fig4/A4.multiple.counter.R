setwd("/Users/user/Desktop/Chen.analysis")
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(ape)
library(gridExtra)

#The A5 csv has been modified to add scaffold to the scaffold number
A4 <- data.frame(read.csv("supp6.A4.raw.csv",skip=4,na.strings="",stringsAsFactors = FALSE))

#remove SN13 because that is what Chen et al did, then make all character class
A4[,4:18] <- lapply(A4[,4:18], as.character)
A4 <- A4[,colnames(A4) != "SN13"]

#I need this order for comparison with the two csv files. This is the order
#from Chen et al figure 2
A4_name_order <- c("SN19","SN03","SN09","SN27","SN14","SN16","SN26",
                   "SN20","SN02","SN12","SN11","SN24","SN21","SN22")
#make missing values NA so they are ignore by the recombination algorithm
A4[A4 == ""]<- NA
#now make long form of A4
data_long <- data.frame(gather(A4, individual, genotype, SN03:SN27,factor_key = FALSE),stringsAsFactors = FALSE)
data_long <- data_long[data_long$individual != "SN13",]
data_long$newname <- paste0(data_long$scaffold,data_long$position)
data_long[] <- lapply(data_long,as.character)
data_long <- data_long[order(data_long$scaffold,data_long$position,match(data_long$individual,A4_name_order)),]

A4_recomb_depths <- data.frame(read.csv("A4_recomb_depths.csv"),stringsAsFactors = FALSE)
A4_recomb_depths <- A4_recomb_depths[A4_recomb_depths$individual != "SN13",]
A4_recomb_depths$total <- rowSums(A4_recomb_depths[,5:9])
#calculate max of reads that are not N
A4_recomb_depths$max <- apply(A4_recomb_depths[,5:8],1,max)
#maf is the major allele fraction. I will keep things above 0.90 if there are above 10 reads
A4_recomb_depths$maf <- (A4_recomb_depths$max)/(A4_recomb_depths$total)
A4_recomb_depths$scaffold <- as.character(A4_recomb_depths$scaffold)
A4_recomb_depths$position <- as.character(A4_recomb_depths$position)
A4_recomb_depths <- A4_recomb_depths[order(A4_recomb_depths$scaffold,A4_recomb_depths$position,match(A4_recomb_depths$individual,A4_name_order)),]
A4_recomb_depths$pair_ratio <- A4_recomb_depths$bad_reads/A4_recomb_depths$total_reads
#this makes an identifier to be used later
A4_recomb_depths$newname <- paste0(A4_recomb_depths$scaffold,A4_recomb_depths$position)

#this checks that the newname column is the same, there should be 33068 values
sum(data_long$newname == A4_recomb_depths$newname)
sum(data_long$individual == A4_recomb_depths$individual)

#removing low read depths, with less than 5 coverage
data_long$no_low <- data_long$genotype
data_long$no_low[A4_recomb_depths$total < 5] <- NA

#now remove positions where at least one individual is heterozygous.
#the justification is that due to incomplete coverage, not all heterozygotes will
#be caught. First we get het_sites which has all the heterozygous sites
het_sites <- data_long[which(A4_recomb_depths$maf  < 0.9 & A4_recomb_depths$total > 9),]
#then we get a new column, and match if both the the position and the scaffold are in het sites
#based on the newname column, which is the scaffold and position pasted together
data_long$no_het <- data_long$genotype
data_long$no_het[data_long$newname %in% het_sites$newname]<- NA

data_long$mega200<- data_long$genotype
data_long$mega200[A4_recomb_depths$mega200 > 1] <- NA

data_long$mega500 <- data_long$genotype
data_long$mega500[A4_recomb_depths$mega500 > 1] <- NA

data_long$mega750 <- data_long$genotype
data_long$mega750[A4_recomb_depths$mega750 > 1] <- NA

data_long$mega1000 <- data_long$genotype
data_long$mega1000[A4_recomb_depths$mega1000 > 1] <- NA

data_long$mega1500 <- data_long$genotype
data_long$mega1500[A4_recomb_depths$meag1500 > 1] <- NA

data_long$no_three <- data_long$genotype
data_long$no_three[is.na(data_long$no_low)] <- NA
data_long$no_three[is.na(data_long$no_het)] <- NA
data_long$no_three[is.na(data_long$mega200)] <- NA

colnames(data_long)

#Now lets make new version of the matrices and count differences
A4_mock      <- spread(data_long[,c(1,3,6,7)],individual,genotype)
A4_mock      <- A4_mock[,c("scaffold","position",A4_name_order)]
A4_no_low    <- spread(data_long[,c(1,3,6,9)],individual,no_low)
A4_no_low    <- A4_no_low[,c("scaffold","position",A4_name_order)]
A4_no_het   <- spread(data_long[,c(1,3,6,10)],individual,no_het)
A4_no_het   <- A4_no_het[,c("scaffold","position",A4_name_order)]
A4_mega200 <- spread(data_long[,c(1,3,6,11)],individual,mega200)
A4_mega200 <- A4_mega200[,c("scaffold","position",A4_name_order)]
A4_mega500 <- spread(data_long[,c(1,3,6,12)],individual,mega500)
A4_mega500 <- A4_mega500[,c("scaffold","position",A4_name_order)]
A4_mega750 <- spread(data_long[,c(1,3,6,13)],individual,mega750)
A4_mega750 <- A4_mega750[,c("scaffold","position",A4_name_order)]
A4_mega1000 <- spread(data_long[,c(1,3,6,14)],individual,mega1000)
A4_mega1000 <- A4_mega1000[,c("scaffold","position",A4_name_order)]
A4_mega1500 <- spread(data_long[,c(1,3,6,15)],individual,mega1500)
A4_mega1500 <- A4_mega1500[,c("scaffold","position",A4_name_order)]
A4_no_three <- spread(data_long[,c(1,3,6,16)],individual,no_three)
A4_no_three <- A4_no_three[,c("scaffold","position",A4_name_order)]




A4_mock$mat_a_unique <- apply(A4_mock[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_mock$mat_b_unique <- apply(A4_mock[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_mock$mat_a_unique[A4_mock$mat_a_unique == 2])
length(A4_mock$mat_b_unique[A4_mock$mat_b_unique == 2])

A4_no_low$mat_a_unique <- apply(A4_no_low[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_no_low$mat_b_unique <- apply(A4_no_low[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_no_low$mat_a_unique[A4_no_low$mat_a_unique == 2])
length(A4_no_low$mat_b_unique[A4_no_low$mat_b_unique==2])

A4_no_het$mat_a_unique <- apply(A4_no_het[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_no_het$mat_b_unique <- apply(A4_no_het[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_no_het$mat_a_unique[A4_no_het$mat_a_unique == 2])
length(A4_no_het$mat_b_unique[A4_no_het$mat_b_unique==2])

A4_mega200$mat_a_unique <- apply(A4_mega200[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_mega200$mat_b_unique <- apply(A4_mega200[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_mega200$mat_a_unique[A4_mega200$mat_a_unique == 2])
length(A4_mega200$mat_b_unique[A4_mega200$mat_b_unique==2])

A4_two_blast$mat_a_unique <- apply(A4_two_blast[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_two_blast$mat_b_unique <- apply(A4_two_blast[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_two_blast$mat_a_unique[A4_two_blast$mat_a_unique == 2])
length(A4_two_blast$mat_b_unique[A4_two_blast$mat_b_unique==2])

A4_three_blast$mat_a_unique <- apply(A4_three_blast[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_three_blast$mat_b_unique <- apply(A4_three_blast[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_three_blast$mat_a_unique[A4_three_blast$mat_a_unique == 2])
length(A4_three_blast$mat_b_unique[A4_three_blast$mat_b_unique==2])

A4_four_blast$mat_a_unique <- apply(A4_four_blast[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_four_blast$mat_b_unique <- apply(A4_four_blast[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_four_blast$mat_a_unique[A4_four_blast$mat_a_unique == 2])
length(A4_four_blast$mat_b_unique[A4_four_blast$mat_b_unique==2])

A4_five_blast$mat_a_unique <- apply(A4_five_blast[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_five_blast$mat_b_unique <- apply(A4_five_blast[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_five_blast$mat_a_unique[A4_five_blast$mat_a_unique == 2])
length(A4_five_blast$mat_b_unique[A4_five_blast$mat_b_unique==2])

A4_no_three$mat_a_unique <- apply(A4_no_three[,3:5]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A4_no_three$mat_b_unique <- apply(A4_no_three[,6:16] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A4_no_three$mat_a_unique[A4_no_three$mat_a_unique == 2])
length(A4_no_three$mat_b_unique[A4_no_three$mat_b_unique==2])

A4_mock[A4_mock$mat_a_unique == 2 | A4_mock$mat_b_unique == 2,]
A4_no_three[A4_no_three$mat_a_unique == 2 | A4_no_three$mat_b_unique == 2,]




A4_m_mock    <- as.DNAbin(matrix(A4_mock[,3:16]))
A4_m_no_three  <- as.DNAbin(matrix(A4_no_three[,3:16]))

A4_mock_melted   <- melt(1- dist.dna(A4_m_mock,   model="raw",as.matrix=TRUE,pairwise.deletion = TRUE))
A4_no_three_melted <- melt(1- dist.dna(A4_m_no_three, model="raw",as.matrix=TRUE,pairwise.deletion = TRUE))

A4_name_order <- sub("SN","",A4_name_order)
png("/Users/user/Desktop/Chen.analysis/A4.matrices.png", width=600 ,height=300)
p1 <- ggplot(data = A4_mock_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(show.legend=FALSE) + theme_classic() +
  scale_x_continuous(labels=A4_name_order, breaks=seq(1,length(A4_name_order)),expand=c(0,0))+
  scale_y_continuous(labels=A4_name_order, breaks=seq(1,length(A4_name_order)),expand=c(0,0))+
  scale_fill_gradient(low = "#ecf5fc", high = "#668dc8", space = "Lab", 
                      na.value = "black", guide = "colourbar",aesthetics = "fill")+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust =0.5),
        axis.text.y = element_text(size=16),
        axis.title =element_blank(),
        axis.line = element_blank())
p2 <- ggplot(data = A4_no_three_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(show.legend=FALSE) + theme_classic() +
  scale_x_continuous(labels=A4_name_order, breaks=seq(1,length(A4_name_order)),expand=c(0,0))+
  scale_y_continuous(labels=A4_name_order, breaks=seq(1,length(A4_name_order)),expand=c(0,0))+
  scale_fill_gradient(low = "#ecf5fc", high = "#668dc8", space = "Lab", 
                      na.value = "black", guide = "colourbar",aesthetics = "fill")+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust =0.5),
        axis.text.y = element_text(size=16),
        axis.title =element_blank(),
        axis.line = element_blank())
grid.arrange(p1,p2,nrow=1)
dev.off()
