library(ggplot2)
library(gridExtra)
library(reshape2)


A4_r <- read.csv("/Users/user/Desktop/Chen.analysis/A4_recomb_depths.csv",header=TRUE,stringsAsFactors=FALSE)
A4_r$total <- rowSums(A4_r[,4:8])
A4_r$max <- apply(A4_r[,4:7],1,max)
A4_recomb_sites <- read.csv("/Users/user/Desktop/Chen.analysis/supp6.A4.recombine.marked.csv",header=TRUE,skip=4,stringsAsFactors = FALSE)
A4_recomb_melt <- melt(A4_recomb_sites, id=c("scaffold","position","reference"))
#A4_recomb_melt <- A4_recomb_melt[A4_recomb_melt$value=="R",]
colnames(A4_recomb_melt) <- c("scaffold","position","reference","individual","value")
A4_recomb_melt$scaffold <- as.character(A4_recomb_melt$scaffold)
A4_recomb_melt$individual <- as.character(A4_recomb_melt$individual)
A4_recomb <- inner_join(A4_r,A4_recomb_melt)
head(A4_recomb)
head(A4_r)


A5_r <- read.csv("/Users/user/Desktop/Chen.analysis/A5_recomb_depths.csv",header=TRUE)
A5_r$total <- rowSums(A5_r[,4:8])
A5_r$max <- apply(A5_r[,4:7],1,max)
A5_recomb_sites <- read.csv("/Users/user/Desktop/Chen.analysis/supp6.A5.recombine.marked.csv",header=TRUE,skip=4,sep="\t")
A5_recomb_melt <- melt(A5_recomb_sites, id=c("scaffold","reference","position"))
colnames(A5_recomb_melt) <- c("scaffold","reference", "position","individual","value")
A5_recomb <- inner_join(A5_r,A5_recomb_melt)
head(A5_recomb)

SL1_r <- read.csv("/Users/user/Desktop/Chen.analysis/SL1_recomb_depths.csv",header=TRUE)
SL1_r$total <- rowSums(SL1_r[,4:8])
SL1_r$max <- apply(SL1_r[,4:7],1,max)
SL1_recomb_sites <- read.csv("/Users/user/Desktop/Chen.analysis/supp6.SL1.recombine.marked.csv",header=TRUE,skip=4)
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
png("/Users/user/Desktop/Chen.analysis/recombined.sites.png",width=900,height=600)
grid.arrange(g1,g2,g3,g4,g5,g6,nrow=2)
dev.off()

# max.len<-max(length(A4_r$total),length(A5_r$total),length(SL1_r$total))
# 
# excel_sites <- data.frame(A4 = c(A4_r$total, rep(NA, max.len - length(A4_r$total))),
#                           A5 = c(A5_r$total, rep(NA, max.len - length(A5_r$total))),
#                           SL1 = c(SL1_r$total, rep(NA, max.len - length(SL1_r$total))))
# 
# excel_melt <- melt(excel_sites)
# 
# max.rec.len <- max(length(A4_recomb$total),length(A5_recomb$total),length(SL1_recomb$total))
# 
# recomb_sites <- data.frame(A4 = c(A4_recomb$total, rep(NA, max.len - length(A4_recomb$total))),
#                            A5 = c(A5_recomb$total, rep(NA, max.len - length(A5_recomb$total))),
#                            SL1 = c(SL1_recomb$total, rep(NA, max.len - length(SL1_recomb$total))))
# 
# recomb_melt <- melt(recomb_sites)
# 
# 
# p4 <- ggplot(recomb_melt, aes(x=value,y=..ndensity..)) + 
#   #scale_y_continuous(expand=c(0,0.8)) +
#   geom_histogram(bins=50)# +
# #geom_vline(data=summed_means,aes(xintercept=means),linetype=1)+
# #geom_vline(data=summed_means,aes(xintercept=medians),linetype=2)+
# #geom_text(data=summed_means,aes(label=means_text ,x=38000,y=14000),inherit.aes = FALSE,hjust=1)+
# #geom_text(data=summed_means,aes(label=medians_text,x=38000,y=12000),inherit.aes = FALSE,hjust=1)
# #png("/Users/user/Downloads/recombined.depth.png")
# p4 + facet_grid(rows=vars(variable)) + theme_classic() + 
#   ggtitle("Sequencing Depth for Sites With 'Wrong'  \n Color Indicating Recombination") +
#   labs(y ="Number of loci",x="") + scale_x_log10(breaks=c(0.5,1,10,100,1000,10000),label=c("0","1","10","100","1000","10000"))
# 
# print(A4_nozero,vp = viewport(width=0.5,height=0.22,x=0.65,y=0.91,just="top"))
# print(A5_nozero,vp = viewport(width=0.5,height=0.22,x=0.65,y=0.64,just="top"))
# print(SL1_nozero,vp = viewport(width=0.5,height=0.22,x=0.65,y=0.355,just="top"))
# dev.off()
# 
# 
# nrow(recomb_melt[recomb_melt$variable == "A4" & recomb_melt$value > 9,])
# nrow(recomb_melt[recomb_melt$variable == "A4" & recomb_melt$value > 9,])
# nrow(recomb_melt[recomb_melt$variable == "A4" & recomb_melt$value > 9,])
# 
# 
# hist(log10(A4_r$total[A4_r$total > 0]),xlim=c(0,5),breaks=20,xlab="log10 of depth",main="")
# title("Sequencing Depth of all non-zero \n sites in excel file for A4")
# 
# 
# hist(log10(A5_r$total[A5_r$total > 0]),xlim=c(0,5),breaks=20,xlab="log10 of depth",main="")
# title("Sequencing Depth of all non-zero \n sites in excel file for A5")
# hist(log10(recomb$total[A4_recomb$total]),xlim=c(0,5),breaks=20,xlab="log10 of depth",main="")
# title("Sequencing Depth of  \n 'interesting' sites for A5")
# 
# 
# hist(log10(SL1_r$total[SL1_r$total > 0]),xlim=c(0,5),breaks=20,xlab="log10 of depth",main="")
# title("Sequencing Depth of all non-zero \n sites in excel file for SL1")
# hist(log10(recomb$total[recomb$total]),xlim=c(0,5),breaks=20,xlab="log10 of depth",main="")
# title("Sequencing Depth of  \n 'interesting' sites for SL1")

