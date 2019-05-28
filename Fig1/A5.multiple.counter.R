setwd("/Users/user/Desktop/Chen.analysis")
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(ape)

#The A5 csv has been modified to add scaffold to the scaffold number
A5 <- data.frame(read.csv("supp6.A5.raw.csv",skip=4,na.strings="",stringsAsFactors = FALSE))
#make genotypes character vector
A5[,4:13] <- lapply(A5[,4:13], as.character)


#I need this order for comparison with the two csv files
A5_name_order <- c("SN03E","SN05E","SN01H","SN02A","SN05B","SN06C","SN03C","SN05D")

A5_recomb_depths <- data.frame(read.csv("A5_recomb_depths.csv"),stringsAsFactors = FALSE)
A5_recomb_depths$total <- rowSums(A5_recomb_depths[,5:9])
A5_recomb_depths$max <- apply(A5_recomb_depths[,5:8],1,max)
#maf is the major allele fraction. I will keep things above 0.90
A5_recomb_depths$maf <- (A5_recomb_depths$max)/(A5_recomb_depths$total)
A5_recomb_depths$scaffold <- as.character(A5_recomb_depths$scaffold)
A5_recomb_depths$position <- as.numeric(A5_recomb_depths$position)
A5_recomb_depths$newname <- paste0(A5_recomb_depths$scaffold,A5_recomb_depths$position)

A5_recomb_depths <- A5_recomb_depths[order(A5_recomb_depths$scaffold,A5_recomb_depths$position,match(A5_recomb_depths$individual,A5_name_order)),]


data_long <- gather(A5, individual, genotype, SN03E:SN02A, factor_key=FALSE)
data_long <- data_long[order(data_long$scaffold,data_long$position,match(data_long$individual,A5_name_order)),]
data_long$newname <- paste0(data_long$scaffold,data_long$position)

#this checks the newname of the two, should sum to 19864
sum(A5_recomb_depths$newname == data_long$newname)
sum(A5_recomb_depths$individual == data_long$individual)
#removing low read depth
#now we remove the sites with less than 4 coverage
data_long$no_low <- data_long$genotype
data_long$no_low[A5_recomb_depths$total < 5] <- NA

#now remove positions where at least one individual is heterozygous.
#the justification is that due to incomplete coverage, not all heterozygotes will
#be caught. First we get het_sites which has all the heterozygous sites
het_sites <- data_long[which(A5_recomb_depths$maf  < 0.9 & A5_recomb_depths$total > 9),]
#then we get a new column, and match if both the the position and the scaffold are in het sites
#while this may potentially cause some false matches since a position could be heterozyous in one position on one scaffold
#and not on another, since the data is so sparse, I do not think this is a problem. 

data_long$no_het <- data_long$genotype
data_long$no_het[data_long$newname %in% het_sites$newname] <- NA

data_long$mega200 <- data_long$genotype
data_long$mega200[which(A5_recomb_depths$mega200 > 1)] <- NA

data_long$mega500 <- data_long$genotype
data_long$mega500[which(A5_recomb_depths$mega500 > 1)] <- NA

data_long$mega750 <- data_long$genotype
data_long$mega750[which(A5_recomb_depths$mega750 > 1)] <- NA

data_long$mega1000 <- data_long$genotype
data_long$mega1000[which(A5_recomb_depths$mega1000 > 1)] <- NA

data_long$mega1500 <- data_long$genotype
data_long$mega1500[which(A5_recomb_depths$mega1500 > 1)] <- NA

data_long$no_three <- data_long$genotype
data_long$no_three[A5_recomb_depths$total < 5] <- NA
data_long$no_three[data_long$newname %in% het_sites$newname] <- NA
data_long$no_three[which(A5_recomb_depths$mega200 > 1)] <- NA
colnames(data_long)

#Now lets make new version of the matrices and count differences
A5_mock <- spread(data_long[,c(1,3,6,7)],individual,genotype)
A5_mock <- A5_mock[,c("scaffold","position",A5_name_order)]
A5_no_low <- spread(data_long[,c(1,3,6,9)],individual,no_low)
A5_no_low <- A5_no_low[,c("scaffold","position",A5_name_order)]
A5_no_het <- spread(data_long[,c(1,3,6,10)],individual,no_het)
A5_no_het <- A5_no_het[,c("scaffold","position",A5_name_order)]
A5_mega200 <- spread(data_long[,c(1,3,6,11)],individual,mega200)
A5_mega200 <- A5_mega200[,c("scaffold","position",A5_name_order)]
A5_two_blast <- spread(data_long[,c(1,3,6,12)],individual,two_blast)
A5_two_blast <- A5_two_blast[,c("scaffold","position",A5_name_order)]
A5_three_blast <- spread(data_long[,c(1,3,6,13)],individual,three_blast)
A5_three_blast <- A5_three_blast[,c("scaffold","position",A5_name_order)]
A5_four_blast <- spread(data_long[,c(1,3,6,14)],individual,four_blast)
A5_four_blast <- A5_four_blast[,c("scaffold","position",A5_name_order)]
A5_five_blast <- spread(data_long[,c(1,3,6,15)],individual,five_blast)
A5_five_blast <- A5_five_blast[,c("scaffold","position",A5_name_order)]
A5_no_three <- spread(data_long[,c(1,3,6,16)],individual,no_three)
A5_no_three <- A5_no_three[,c("scaffold","position",A5_name_order)]




A5_mock$mat_a_unique <- apply(A5_mock[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_mock$mat_b_unique <- apply(A5_mock[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_mock$mat_a_unique[A5_mock$mat_a_unique == 2])
length(A5_mock$mat_b_unique[A5_mock$mat_b_unique == 2])

A5_no_low$mat_a_unique <- apply(A5_no_low[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_no_low$mat_b_unique <- apply(A5_no_low[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_no_low$mat_a_unique[A5_no_low$mat_a_unique == 2])
length(A5_no_low$mat_b_unique[A5_no_low$mat_b_unique==2])

A5_no_het$mat_a_unique <- apply(A5_no_het[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_no_het$mat_b_unique <- apply(A5_no_het[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_no_het$mat_a_unique[A5_no_het$mat_a_unique == 2])
length(A5_no_het$mat_b_unique[A5_no_het$mat_b_unique==2])

A5_mega200$mat_a_unique <- apply(A5_mega200[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_mega200$mat_b_unique <- apply(A5_mega200[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_mega200$mat_a_unique[A5_mega200$mat_a_unique == 2])
length(A5_mega200$mat_b_unique[A5_mega200$mat_b_unique==2])

A5_two_blast$mat_a_unique <- apply(A5_two_blast[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_two_blast$mat_b_unique <- apply(A5_two_blast[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_two_blast$mat_a_unique[A5_two_blast$mat_a_unique == 2])
length(A5_two_blast$mat_b_unique[A5_two_blast$mat_b_unique==2])

A5_three_blast$mat_a_unique <- apply(A5_three_blast[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_three_blast$mat_b_unique <- apply(A5_three_blast[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_three_blast$mat_a_unique[A5_three_blast$mat_a_unique == 2])
length(A5_three_blast$mat_b_unique[A5_three_blast$mat_b_unique==2])

A5_four_blast$mat_a_unique <- apply(A5_four_blast[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_four_blast$mat_b_unique <- apply(A5_four_blast[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_four_blast$mat_a_unique[A5_four_blast$mat_a_unique == 2])
length(A5_four_blast$mat_b_unique[A5_four_blast$mat_b_unique==2])

A5_five_blast$mat_a_unique <- apply(A5_five_blast[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_five_blast$mat_b_unique <- apply(A5_five_blast[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_five_blast$mat_a_unique[A5_five_blast$mat_a_unique == 2])
length(A5_five_blast$mat_b_unique[A5_five_blast$mat_b_unique==2])

A5_no_three$mat_a_unique <- apply(A5_no_three[,3:4]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
A5_no_three$mat_b_unique <- apply(A5_no_three[,5:10] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(A5_no_three$mat_a_unique[A5_no_three$mat_a_unique == 2])
length(A5_no_three$mat_b_unique[A5_no_three$mat_b_unique==2])

A5_no_three[A5_no_three$mat_a_unique > 1 | A5_no_three$mat_b_unique > 1,]
doubled_five_blast <- A5_five_blast[A5_five_blast$mat_a_unique > 1 | A5_five_blast$mat_b_unique > 1,]


A5_m_mock    <- as.DNAbin(matrix(A5_mock[,3:10]))
A5_m_no_three  <- as.DNAbin(matrix(A5_no_three[,3:10]))

A5_mock_melted   <- melt(1-dist.dna(A5_m_mock, model="raw",as.matrix=TRUE,pairwise.deletion = TRUE))
A5_no_three_melted <- melt(1-dist.dna(A5_m_no_three, model="raw",as.matrix=TRUE,pairwise.deletion = TRUE))
A5_name_order <- sub("SN","",A5_name_order)

png("/Users/user/Desktop/Chen.analysis/A5.matrices.png",width=600, height=300)
p1 <- ggplot(data = A5_mock_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(show.legend=FALSE) + theme_classic() +
  scale_x_continuous(labels=A5_name_order, breaks=seq(1,length(A5_name_order)),expand=c(0,0))+
  scale_y_continuous(labels=A5_name_order, breaks=seq(1,length(A5_name_order)),expand=c(0,0))+
  scale_fill_gradient(low = "#ecf5fc", high = "#668dc8", space = "Lab", 
                      na.value = "black", guide = "colourbar",aesthetics = "fill")+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust =0.5),
        axis.text.y = element_text(size=16),
        axis.title =element_blank(),
        axis.line = element_blank())
p2 <- ggplot(data = A5_no_three_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(show.legend=FALSE) + theme_classic() +
  scale_x_continuous(labels=A5_name_order, breaks=seq(1,length(A5_name_order)),expand=c(0,0))+
  scale_y_continuous(labels=A5_name_order, breaks=seq(1,length(A5_name_order)),expand=c(0,0))+
  scale_fill_gradient(low = "#ecf5fc", high = "#668dc8", space = "Lab", 
                      na.value = "black", guide = "colourbar",aesthetics = "fill")+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust =0.5),
        axis.text.y = element_text(size=16),
        axis.title =element_blank(),
        axis.line = element_blank())
grid.arrange(p1,p2,nrow=1)


dev.off()
