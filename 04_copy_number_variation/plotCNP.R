#@author:aniafijarczyk
library(ggplot2)
library(dplyr)


### Plotting copy number variants in relatives of Jean-Talon

#chrom_conv = read.table("chromLengthConvertion.out", sep = "\t", header = TRUE)
chrom_conv = data.frame("chrom"=c("ref|NC_001133|","ref|NC_001134|","ref|NC_001135|","ref|NC_001136|","ref|NC_001137|","ref|NC_001138|","ref|NC_001139|","ref|NC_001140|","ref|NC_001141|","ref|NC_001142|","ref|NC_001143|","ref|NC_001144|","ref|NC_001145|","ref|NC_001146|","ref|NC_001147|","ref|NC_001148|","ref|NC_001224|"),
                        "length"=c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066,85779),
                        "addition"=c(0,230218,1043402,1360022,2891955,3468829,3738990,4829930,5392573,5832461,6578212,7245028,8323205,9247636,10031969,11123260,12071326))

chrf <- chrom_conv %>% filter(chrom != "ref|NC_001224|")
chrf["end"] = chrf["addition"] + chrf["length"]
chrf["midpoint"] = (chrf["addition"] + chrf["end"])/2

cnp <- read.table("./input_files/getCNVIntervals_all.tab", sep="\t",header=FALSE,col.names=c("chrom","start","stop","CN","strain"))

cnm <- merge(cnp,chrf,by = "chrom",sort=FALSE)
cnm["nstart"] <- cnm["start"] + cnm["addition"]
cnm["nstop"] <- cnm["stop"] + cnm["addition"]

strains <- as.vector(unique(cnm$strain))
ds <- data.frame("strain" = strains, order = c(9,6,1,7,3,8,4,5,2))
do <- ds[order(ds$order),]
cns <- merge(cnm,ds,on="strain",sort=FALSE)


big_chr <- chrf[c(1,2,3,4)]
big_chr["strain"] = strains[1]
for (i in 2:length(strains)) {
  newchrf <- big_chr[c(1,2,3,4)]
  newchrf["strain"] <- strains[i]
  big_chr = rbind(big_chr,newchrf)
}


big_chrm <- merge(big_chr,ds,on="strain",sort=FALSE)
big_chrm["midpoint"] <- (big_chrm["addition"]+big_chrm["end"])/2
chr15 <- chrf %>% filter(chrom == "ref|NC_001147|")

cobalt = data.frame(strain = c('Jean-Talon'), chrom = c('ref|NC_001147|'),
                    start = c(906236+chr15$addition-20000), stop = c(907555+chr15$addition+20000))


p <- ggplot(cns) + 
  geom_rect(data = big_chrm, aes(xmin=addition, xmax=end, ymin=order-0.25, ymax=order+0.25, fill = chrom)) +
  geom_rect(data = cns, aes(xmin=nstart, xmax=nstop, ymin=order-0.25, ymax=order+0.25, fill = as.factor(CN))) +
  geom_rect(data = cobalt, aes(xmin=start, xmax=stop, ymin=0.75,ymax=1.25,fill=chrom),fill="black") +
  geom_text(data = cobalt, aes(x = start, y = 0.4, label="COT1"), colour="black") +
  #scale_fill_manual(values = c("#0072B2","#56B4E9","#E69F00","#D55E00",rep(c("grey80","grey60"),8))) +
  scale_fill_manual(breaks = c("2","3","5","6"), values = c("#0072B2","#56B4E9","#E69F00","#D55E00",rep(c("grey80","grey60"),8)))+
  scale_x_reverse(breaks=as.vector(chrf$midpoint),
                  labels=c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"))+
  scale_y_continuous(breaks = c(1:length(strains)),
                   labels = as.vector(do$strain)) +
  coord_flip() +
  labs(x = "chromosomes") +
  guides(fill=guide_legend(title="copy number")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle=45, hjust=1),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
p

pdf("./output/plotCNP.pdf",w=6,h=5)
p
dev.off()
