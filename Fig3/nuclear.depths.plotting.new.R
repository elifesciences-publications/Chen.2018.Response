library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)

set.seed(3634567)

A4_r <- read.csv("/Users/user/Desktop/Chen.analysis/A4_recomb_depths.csv",header=TRUE,stringsAsFactors=FALSE)
A4_r$total <- rowSums(A4_r[,5:9])
A4_r$max <- apply(A4_r[,5:8],1,max)
A4_recomb_sites <- read.delim("/Users/user/Desktop/Chen.analysis/supp6.A4.recombine.marked.csv",sep="\t",header=TRUE,skip=2,stringsAsFactors = FALSE)
A4_recomb_melt <- melt(A4_recomb_sites, id=c("scaffold","position","reference"))
#A4_recomb_melt <- A4_recomb_melt[A4_recomb_melt$value=="R",]
colnames(A4_recomb_melt) <- c("scaffold","position","reference","individual","value")
A4_recomb_melt$scaffold <- as.character(A4_recomb_melt$scaffold)
A4_recomb_melt$individual <- as.character(A4_recomb_melt$individual)
A4_recomb <- inner_join(A4_r,A4_recomb_melt)
head(A4_recomb)
head(A4_r)


A5_r <- read.csv("/Users/user/Desktop/Chen.analysis/A5_recomb_depths.csv",header=TRUE)
A5_r$total <- rowSums(A5_r[,5:9])
A5_r$max <- apply(A5_r[,5:8],1,max)
A5_recomb_sites <- read.delim("/Users/user/Desktop/Chen.analysis/supp6.A5.recombine.marked.csv",sep="\t",header=TRUE,skip=4)
A5_recomb_melt <- melt(A5_recomb_sites, id=c("scaffold","reference","position"))
colnames(A5_recomb_melt) <- c("scaffold","reference", "position","individual","value")
A5_recomb <- inner_join(A5_r,A5_recomb_melt)
head(A5_recomb)

SL1_r <- read.csv("/Users/user/Desktop/Chen.analysis/SL1_recomb_depths.csv",header=TRUE)
SL1_r$total <- rowSums(SL1_r[,5:9])
SL1_r$max <- apply(SL1_r[,5:8],1,max)
SL1_recomb_sites <- read.delim("/Users/user/Desktop/Chen.analysis/supp6.SL1.recombine.marked.csv",header=TRUE,skip=4,sep="\t")
SL1_recomb_melt <- melt(SL1_recomb_sites, id=c("scaffold","position"))
#SL1_recomb_melt <- SL1_recomb_melt[SL1_recomb_melt$value=="R",]
colnames(SL1_recomb_melt) <- c("scaffold","position","individual","value")
SL1_recomb <- inner_join(SL1_r,SL1_recomb_melt)
head(SL1_recomb)


theme_update(text = element_text(size=26))
g1 <- ggplot(A4_recomb) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,2500)
g2 <- ggplot(A5_recomb) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,2500)
g3 <- ggplot(SL1_recomb) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,2500)+
  labs(x=expression("Log"[10]*" of read depth per site")) 
g4 <- ggplot(A4_recomb[A4_recomb$value == "R",]) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,165)
g5 <- ggplot(A5_recomb[A5_recomb$value == "R",]) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,165)
g6 <- ggplot(SL1_recomb[SL1_recomb$value == "R",]) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits=c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,165) +
    labs(x=expression("Log"[10]*" of read depth per site"))
png("/Users/user/Desktop/Chen.analysis/Fig2.recombined.sites.test.png",width=1500,height=1000)
grid.arrange(g1,g2,g3,g4,g5,g6,nrow=2)
dev.off()


nrow(A4_recomb[A4_recomb$value == "R",])
A4_sub <- A4_recomb[A4_recomb$total>0,]
A4_sub <- A4_sub[A4_sub$value != "R",]
A4_sub <- A4_sub[sample(nrow(A4_sub),284),]

nrow(A5_recomb[A5_recomb$value == "R",])
A5_sub <- A5_recomb[A5_recomb$total>0,]
A5_sub <- A5_sub[A5_sub$value != "R",]
A5_sub <- A5_sub[sample(nrow(A5_sub),146),]

nrow(SL1_recomb[SL1_recomb$value == "R",])
SL1_sub <- SL1_recomb[SL1_recomb$total>0,]
SL1_sub <- SL1_sub[SL1_sub$value != "R",]
SL1_sub <- SL1_sub[sample(nrow(SL1_sub),171),]



theme_update(axis.text = element_text(size=26))
g1 <- ggplot(A4_sub) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,150)
g2 <- ggplot(A5_sub) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,150)
g3 <- ggplot(SL1_sub) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,150)
g4 <- ggplot(A4_recomb[A4_recomb$value == "R",]) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,150)
g5 <- ggplot(A5_recomb[A5_recomb$value == "R",]) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits = c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,150)
g6 <- ggplot(SL1_recomb[SL1_recomb$value == "R",]) + geom_histogram(aes(x=total),bins=100) + scale_x_log10(limits=c(0.75,15000),breaks=c(1,10,100,1000)) + ylim(0,150)
png("/Users/user/Desktop/Chen.analysis/Fig2.recombined.sites.subsampled.Sep1.png",width=800,height=550)
grid.arrange(g1,g2,g3,g4,g5,g6,nrow=2)
dev.off()

wilcox.test(A4_sub$total,A4_recomb$total[A4_recomb$value == "R"],paired=FALSE)
wilcox.test(A5_sub$total,A5_recomb$total[A5_recomb$value == "R"],paired=FALSE)
wilcox.test(SL1_sub$total,SL1_recomb$total[SL1_recomb$value == "R"],paired=FALSE)
