library(ggplot2)
library(gridExtra)
SL1_mat5 = c("17","05","12","23")
SL1_mat1 = c("14","15","06","07")
A4_mat2 = c("03","19","09")
A4_mat1 = c("26","20","02","12","11","24","21","22")
A5_mat3 = c("3E","5E")
A5_mat6 = c("1H","5B","6C","3C","5D")
setwd("/Users/user/Desktop/Chen.analysis")

A4 <- data.frame(
labels = c("19","03","09","27","14","16","26","20","02","12","11","24","21","22",""),
m511 =  c(  15,   2,  10,   0,   0,   0,  60, 318,   0,1273, 828, 459, 167, 285,0),
m5308 = c( 125,   0, 191,   0,   0,   0,   0,   0,   0,   2,   0,   1,   0,   0,0))

A4$log511 <- log10(A4$m511)
A4$log5308 <- -log10(A4$m5308)

A5 <- data.frame(
  labels=c("3E","5E","1H","2A","5B","6C","3C","5D","","","","","","",""),
  m4735=c(    8, 268,   5,   1,   8,  15,   2,   0,NA,NA,NA,NA,NA,NA,NA),
  m2225=c(    3,   5,1947,   1, 270, 407,  58, 386,NA,NA,NA,NA,NA,NA,NA))

A5$log4735 <- log10(A5$m4735)
A5$log2225 <- -log10(A5$m2225)

SL1 <- data.frame(
  labels=c("21","18","16","25","15","06","14","10","04","17","05","12","23","07","08", ""),
  m11994=c(   0,   0,   0,   0,   0,   0,   0,   9,   0, 964,   1,  71,   0,   0,   0,  0),
  m11899=c(   0,   0,   3,   0,  28,   8,  54,   0,   0,   4,   2,   0,   0,   0,   0,  0)
)

SL1$log11994 <- log10(SL1$m11994)
SL1$log11899 <- -log10(SL1$m11899)

theme_set(theme_classic())
theme_update(axis.title.y=element_blank(),
             axis.text.y=element_text(size=16),
             plot.title = element_text(size = 18),
             axis.text.x=element_text(size=12))

p1 <- ggplot() + 
  geom_bar(data = A4[which(A4$m511>0),], aes(x=labels, y=log511),stat = "identity",fill="#ed7facff",color="black") +
  geom_bar(data = A4[which(A4$m5308>0),], aes(x=labels, y=log5308),stat = "identity",fill="#9eddbcff",color="black") +
  coord_flip() + scale_x_discrete(limits=A4$labels) + 
  theme(axis.text.y=element_text(size=16),plot.title = element_text(size = 18)) + 
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),labels=c("1000","100","10","1","10","100","1000"),limits=c(-4.5,4))+
  ylab('mat2(scaffold_5308)              mat-1(scaffold_511) \n Number of reads mapped to:') + ggtitle("Mapped reads to alternate mating \nloci for A4 nuclei")+
  geom_rect(aes(ymin=-4.5,ymax=-3.8,
                xmin=which(A4$labels %in% A4_mat2)-0.4,
                xmax=which(A4$labels %in% A4_mat2)+0.4),fill="#9eddbcff")+
  geom_rect(aes(ymin=-4.5,ymax=-3.8,
                xmin=which(A4$labels %in% A4_mat1)-0.4,
                xmax=which(A4$labels %in% A4_mat1)+0.4),fill="#ed7facff")
p2 <- ggplot() + 
  geom_bar(data = A5[which(A5$m4735>0),], aes(x=labels, y=log4735),stat = "identity",fill="#fab364ff",color="black") +
  geom_bar(data = A5[which(A5$m2225>0),], aes(x=labels, y=log2225),stat = "identity",fill="#d0b5e0ff",color="black") +
  coord_flip() + scale_x_discrete(limits=A5$labels) +
  theme(axis.text.y=element_text(size=16),plot.title = element_text(size = 18))+ 
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),labels=c("1000","100","10","1","10","100","1000"),limits=c(-4.5,4))+
  ylab('mat-6(scaffold_2225)           mat-3(scaffold_4735) \n Number of reads mapped to:') + ggtitle("Mapped reads to alternate mating\nloci for A5 nuclei")+
  geom_rect(aes(ymin=-4.5,ymax=-3.8,
                xmin=which(A5$labels %in% A5_mat3)-0.4,
                xmax=which(A5$labels %in% A5_mat3)+0.4),fill="#fab364ff")+
  geom_rect(aes(ymin=-4.5,ymax=-3.8,
                xmin=which(A5$labels %in% A5_mat6)-0.4,
                xmax=which(A5$labels %in% A5_mat6)+0.4),fill="#d0b5e0ff")
p3 <- ggplot() + 
  geom_bar(data =SL1[which(SL1$m11899>0),], aes(x=labels, y=log11899),stat = "identity",fill="#ed7facff",color="black") +
  geom_bar(data =SL1[which(SL1$m11994>0),], aes(x=labels, y=log11994),stat = "identity",fill="#fae041ff",color="black") +
  scale_x_discrete(limits=SL1$labels) + coord_flip() +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3),labels=c("1000","100","10","1","10","100","1000"),limits=c(-4.5,4))+
  ylab('mat-1(scaffold_11994)           mat-5(scaffold_11899) \n Number of reads mapped to:') + ggtitle("Mapped reads to alternate mating\nloci for SL1 nuclei")+
  geom_rect(aes(ymin=-4.5,ymax=-3.8,
                xmin=which(SL1$labels %in% SL1_mat1)-0.4,
                xmax=which(SL1$labels %in% SL1_mat1)+0.4),fill="#ed7facff")+
  geom_rect(aes(ymin=-4.5,ymax=-3.8,
                xmin=which(SL1$labels %in% SL1_mat5)-0.4,
                xmax=which(SL1$labels %in% SL1_mat5)+0.4),fill="#fae041ff")
png("/Users/user/Desktop/Chen.analysis/mating_read.Sep1.png",width=1500,height=500)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()

