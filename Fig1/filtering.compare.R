library(ggplot2)
library(RColorBrewer)

col <- brewer.pal(5,"Set1")
col <- col[3:5]

sites <- data.frame(read.csv("/Users/user/Desktop/Chen.analysis/doubled.sites.filtered.csv",header=TRUE,sep="\t"))
sites$pos <- c(1,2,3,4,5)
p1 <- ggplot(sites) + 
  geom_point(aes(x=pos-0.02,y=A4_doubled),color=col[1],shape=16,size=5.5) + #geom_line(aes(x=pos,y=A4),color="red",lwd=0.35) +
  geom_point(aes(x=pos,y=A5_doubled),color=col[2],shape=17,size=5.5) + #geom_line(aes(x=pos,y=A5),color="blue",lwd=0.35) +
  geom_point(aes(x=pos+0.02,y=SL1_doubled),color=col[3],shape=18,size=5.5) + #geom_line(aes(x=pos,y=SL1),color="green",lwd=0.35) +
  geom_text(aes(x=4.5,y=330,hjust=0),label="A4",size=6)+ geom_point(aes(x=4.35,y=330),color=col[1],shape=16,size=4.5)+
  geom_text(aes(x=4.5,y=300,hjust=0),label="A5",size=6)+ geom_point(aes(x=4.35,y=300),color=col[2],shape=17,size=4.5)+
  geom_text(aes(x=4.5,y=270,hjust=0),label="SL1",size=6)+ geom_point(aes(x=4.35,y=270),color=col[3],shape=18,size=4.5)+
  geom_rect(aes(xmin=4.1,ymin=240,xmax=5,ymax=350),color="black",fill=NA)+
  scale_y_continuous(expand = c(0, 1),limits=c(0,400))+
  theme_classic() + labs(y="Number of recombined positions",x=NULL)+
  ggtitle("Recombined Sites")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        text = element_text(size=24)) + scale_x_continuous(breaks=c(1,2,3,4,5),
                                                                                labels=c("Original",
                                                                                         "No Heterozygous Loci", 
                                                                                         "No Repeat Regions",
                                                                                         "No Low Coverage Sites",
                                                                                         "All Three Filters"))

p2 <- ggplot(sites) + 
  geom_point(aes(x=pos-0.02,y=A4_all),color=col[1],shape=16,size=5.5) + #geom_line(aes(x=pos,y=A4),color="red",lwd=0.35) +
  geom_point(aes(x=pos,y=A5_all),color=col[2],shape=17,size=5.5) + #geom_line(aes(x=pos,y=A5),color="blue",lwd=0.35) +
  geom_point(aes(x=pos+0.02,y=SL1_all),color=col[3],shape=18,size=5.5) + #geom_line(aes(x=pos,y=SL1),color="green",lwd=0.35) +
  scale_y_continuous(expand = c(0, 1),limits=c(0,25000))+
  theme_classic() + labs(y="Total number of basepair calls",x=NULL)+
  ggtitle("All Sites")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        text = element_text(size=24)) + scale_x_continuous(breaks=c(1,2,3,4,5),
                                                           labels=c("Original",
                                                                    "No Heterozygous Loci", 
                                                                    "No Repeat Regions",
                                                                    "No Low Coverage Sites",
                                                                    "All Three Filters"))
png("/Users/user/Desktop/Chen.analysis/filtering_both.png",width=750,height=500)
grid.arrange(p2,p1,nrow=1)
dev.off()



