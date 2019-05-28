setwd("/Users/user/Desktop/Chen.analysis")
library(ggplot2)
library(tidyr)
library(dplyr)
library(ape)
library(reshape2)
library(gridExtra)
library(grid)

### Note that supp6.SL1.raw.csv has multiple duplicate rows. On scaffold 66 and 98. This should not be possible for a .vcf file
### no idea if this means anything or not.

#The SL1 csv has been modified to add scaffold to the scaffold number
SL1 <- data.frame(read.csv("supp6.SL1.raw.csv",skip=4,na.strings="",stringsAsFactors = FALSE))
#make genotypes character vector
SL1[,4:18] <- lapply(SL1[,4:18], as.character)
SL1 <- SL1[,colnames(SL1) != "SN13"]
#I need to reorder the rows to match the Figure 3 order

SL1_name_order <- c("SN21","SN18","SN16","SN25","SN15","SN06","SN14",
                    "SN10","SN04","SN17","SN05","SN12","SN23","SN07","SN08")

SL1_recomb_depths <- data.frame(read.csv("SL1_recomb_depths.csv"),stringsAsFactors = FALSE)
SL1_recomb_depths <- SL1_recomb_depths[SL1_recomb_depths$individual != "SN13",]
SL1_recomb_depths$total <- rowSums(SL1_recomb_depths[,5:9])
SL1_recomb_depths$max <- apply(SL1_recomb_depths[,5:8],1,max)
#maf is the major allele fraction. I will keep things above 0.90
SL1_recomb_depths$maf <- (SL1_recomb_depths$max)/(SL1_recomb_depths$total)
SL1_recomb_depths$scaffold <- as.character(SL1_recomb_depths$scaffold)
SL1_recomb_depths$newname <- paste0(SL1_recomb_depths$scaffold,SL1_recomb_depths$position)

SL1_recomb_depths <- SL1_recomb_depths[order(SL1_recomb_depths$scaffold,SL1_recomb_depths$position,match(SL1_recomb_depths$individual,SL1_name_order)),]


data_long <- gather(SL1, individual, genotype, SN06:SN25, factor_key=FALSE)
data_long <- data_long[order(data_long$scaffold,data_long$position,match(data_long$individual,SL1_name_order)),]
data_long$newname <- paste0(data_long$scaffold,data_long$position)

#should sum to 18465
sum(data_long$newname == SL1_recomb_depths$newname)
sum(data_long$individual == SL1_recomb_depths$individual)

#removing low read depth
#now we remove the sites with less than 4 coverage
data_long$no_low <- data_long$genotype
data_long$no_low[SL1_recomb_depths$total < 5] <- NA
#reomve sites that are heterozygous

#now remove positions where at least one individual is heterozygous.
#the justification is that due to incomplete coverage, not all heterozygotes will
#be caught. First we get het_sites which has all the heterozygous sites
het_sites <- data_long[which(SL1_recomb_depths$maf  < 0.9 & SL1_recomb_depths$total > 9),]
#then we get a new column, and match if both the the position and the scaffold are in het sites
#while this may potentially cause some false matches since a position could be heterozyous in one position on one scaffold
#and not on another, since the data is so sparse, I do not think this is a problem. 
data_long$no_het <- data_long$genotype
data_long$no_het[data_long$newname %in% het_sites$newname]<- NA

data_long$mega200 <- data_long$genotype
data_long$mega200[SL1_recomb_depths$mega200 > 1] <- NA

data_long$no_three <- data_long$genotype
data_long$no_three[is.na(data_long$no_low)] <- NA
data_long$no_three[is.na(data_long$no_het)] <- NA
data_long$no_three[SL1_recomb_depths$mega200 > 1] <- NA

data_long$mega500 <- data_long$genotype
data_long$mega500[is.na(data_long$no_low)] <- NA
data_long$mega500[is.na(data_long$no_het)] <- NA
data_long$meag500[SL1_recomb_depths$mega500 > 1] <- NA

data_long$meag750 <- data_long$genotype
data_long$mega750[is.na(data_long$no_low)] <- NA
data_long$mega750[is.na(data_long$no_het)] <- NA
data_long$mega750[SL1_recomb_depths$mega750 > 1] <- NA

data_long$mega1000 <- data_long$genotype
data_long$mega1000[is.na(data_long$no_low)] <- NA
data_long$mega1000[is.na(data_long$no_het)] <- NA
data_long$mega1000[SL1_recomb_depths$mega1500 > 1] <- NA

data_long$mega1500 <- data_long$genotype
data_long$mega1500[is.na(data_long$no_low)] <- NA
data_long$mega1500[is.na(data_long$no_het)] <- NA
data_long$mega1500[SL1_recomb_depths$mega1500 > 1] <- NA


colnames(data_long)

#Now lets make new version of the matrices and count differences
SL1_mock <- spread(data_long[,c(1,3,6,7)],individual,genotype)
SL1_mock <- SL1_mock[,c("scaffold","position",SL1_name_order)]
SL1_no_low <- spread(data_long[,c(1,3,6,9)],individual,no_low)
SL1_no_low <- SL1_no_low[,c("scaffold","position",SL1_name_order)]
SL1_no_het <- spread(data_long[,c(1,3,6,10)],individual,no_het)
SL1_no_het <- SL1_no_het[,c("scaffold","position",SL1_name_order)]
SL1_mega200 <- spread(data_long[,c(1,3,6,11)],individual,mega200)
SL1_mega200 <- SL1_mega200[,c("scaffold","position",SL1_name_order)]
SL1_two_blast <- spread(data_long[,c(1,2,3,4,6,12)],individual,two_blast)
SL1_two_blast <- SL1_two_blast[,c("scaffold","position",SL1_name_order)]
SL1_three_blast <- spread(data_long[,c(1,2,3,4,6,13)],individual,three_blast)
SL1_three_blast <- SL1_three_blast[,c("scaffold","position",SL1_name_order)]
SL1_four_blast <- spread(data_long[,c(1,2,3,4,6,14)],individual,four_blast)
SL1_four_blast <- SL1_four_blast[,c("scaffold","position",SL1_name_order)]
SL1_five_blast <- spread(data_long[,c(1,2,3,4,6,15)],individual,five_blast)
SL1_five_blast <- SL1_five_blast[,c("scaffold","position",SL1_name_order)]
SL1_no_three <- spread(data_long[,c(1,3,6,16)],individual,no_three)
SL1_no_three <- SL1_no_three[,c("scaffold","position",SL1_name_order)]


SL1_mock$mat_a_unique <- apply(SL1_mock[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_mock$mat_b_unique <- apply(SL1_mock[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_mock$mat_a_unique[SL1_mock$mat_a_unique == 2])
length(SL1_mock$mat_b_unique[SL1_mock$mat_b_unique == 2])

SL1_no_low$mat_a_unique <- apply(SL1_no_low[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_no_low$mat_b_unique <- apply(SL1_no_low[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_no_low$mat_a_unique[SL1_no_low$mat_a_unique == 2])
length(SL1_no_low$mat_b_unique[SL1_no_low$mat_b_unique==2])

SL1_no_het$mat_a_unique <- apply(SL1_no_het[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_no_het$mat_b_unique <- apply(SL1_no_het[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_no_het$mat_a_unique[SL1_no_het$mat_a_unique == 2])
length(SL1_no_het$mat_b_unique[SL1_no_het$mat_b_unique==2])

SL1_mega200$mat_a_unique <- apply(SL1_mega200[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_mega200$mat_b_unique <- apply(SL1_mega200[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_mega200$mat_a_unique[SL1_mega200$mat_a_unique == 2])
length(SL1_mega200$mat_b_unique[SL1_mega200$mat_b_unique==2])

SL1_mega500$mat_a_unique <- apply(SL1_two_blast[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_mega500$mat_b_unique <- apply(SL1_two_blast[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_two_blast$mat_a_unique[SL1_two_blast$mat_a_unique == 2])
length(SL1_two_blast$mat_b_unique[SL1_two_blast$mat_b_unique==2])

SL1_three_blast$mat_a_unique <- apply(SL1_three_blast[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_three_blast$mat_b_unique <- apply(SL1_three_blast[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_three_blast$mat_a_unique[SL1_three_blast$mat_a_unique == 2])
length(SL1_three_blast$mat_b_unique[SL1_three_blast$mat_b_unique==2])

SL1_four_blast$mat_a_unique <- apply(SL1_four_blast[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_four_blast$mat_b_unique <- apply(SL1_four_blast[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_four_blast$mat_a_unique[SL1_four_blast$mat_a_unique == 2])
length(SL1_four_blast$mat_b_unique[SL1_four_blast$mat_b_unique==2])

SL1_five_blast$mat_a_unique <- apply(SL1_five_blast[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_five_blast$mat_b_unique <- apply(SL1_five_blast[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_five_blast$mat_a_unique[SL1_five_blast$mat_a_unique == 2])
length(SL1_five_blast$mat_b_unique[SL1_five_blast$mat_b_unique==2])

SL1_no_three$mat_a_unique <- apply(SL1_no_three[,3:9]  ,1,function(x) sum(!is.na(unique(toupper(x)))))
SL1_no_three$mat_b_unique <- apply(SL1_no_three[,10:17] ,1,function(x) sum(!is.na(unique(toupper(x)))))
length(SL1_no_three$mat_a_unique[SL1_no_three$mat_a_unique == 2])
length(SL1_no_three$mat_b_unique[SL1_no_three$mat_b_unique==2])

SL1_no_three[SL1_no_three$mat_a_unique == 2 | SL1_no_three$mat_b_unique == 2,]

SL1_five_blast[which(is.na(SL1_five_blast$SN21) == FALSE),]

SL1_m_mock    <- as.DNAbin(matrix(SL1_mock[,3:17]))
SL1_m_no_three <- as.DNAbin(matrix(SL1_no_three[,3:17]))

SL1_mock_melted   <- melt(1-dist.dna(SL1_m_mock, model="raw",as.matrix=TRUE,pairwise.deletion = TRUE))
SL1_no_three_melted<- melt(1-dist.dna(SL1_m_no_three,model="raw",as.matrix=TRUE,pairwise.deletion = TRUE))

png("/Users/user/Desktop/Chen.analysis/SL1.matrices.png",width=600,height=300)
SL1_name_order <- sub("SN","",SL1_name_order)
p1 <- ggplot(data = SL1_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(show.legend=FALSE) + theme_classic() + 
  scale_x_continuous(labels=SL1_name_order, breaks=seq(1,length(SL1_name_order)),expand=c(0,0))+
  scale_y_continuous(labels=SL1_name_order, breaks=seq(1,length(SL1_name_order)),expand=c(0,0))+
  scale_fill_gradient(low = "#ecf5fc", high = "#668dc8", space = "Lab", 
                      na.value = "black", guide = "colourbar",aesthetics = "fill")+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust =0.5),
        axis.text.y = element_text(size=16),
        axis.title =element_blank(),
        axis.line = element_blank())
p2 <- ggplot(data = SL1_no_three_melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(show.legend=FALSE) + theme_classic() +
  scale_x_continuous(labels=SL1_name_order, breaks=seq(1,length(SL1_name_order)),expand=c(0,0))+
  scale_y_continuous(labels=SL1_name_order, breaks=seq(1,length(SL1_name_order)),expand=c(0,0))+
  scale_fill_gradient(low = "#ecf5fc", high = "#668dc8", space = "Lab", 
                      na.value = "black", guide = "colourbar",aesthetics = "fill")+
  theme(axis.text.x = element_text(size=16,angle = 90, vjust =0.5),
      axis.text.y = element_text(size=16),
      axis.title =element_blank(),
      axis.line = element_blank())
grid.arrange(p1,p2,nrow=1)


dev.off()

png("/Users/user/Desktop/Chen.analysis/matrix.scale.png")
p2 + geom_tile(show.legend=TRUE) + theme(legend.text=element_text(size=20),
                                         legend.title=element_text(size=20),
                                         legend.key.size = unit(1.5,"cm"))
dev.off()


