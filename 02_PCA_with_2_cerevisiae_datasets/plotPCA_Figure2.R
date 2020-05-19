#author@aniafijarczyk
library("SNPRelate")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)


### Running PCA on 1011 yeast collection + Jean-Talon


# Reading vcf files
#vcf.fn <- "genotypes_1011_JT_biall.recode_AF_PL.vcf.gz"
#vcf.fn <- "FileS3.vcf.gz"
vcf.fn <- "./input_files/sample_1011_JT_biall.recode_AF_PL.vcf.gz"
snpgdsVCF2GDS_R(vcf.fn, "genotypes_1011_JT_biall.recode.gds", method="biallelic.only",option=snpgdsOption(chromosome1=1, chromosome2=2, chromosome3=3, chromosome4=4, chromosome5=5, chromosome6=6, chromosome7=7, chromosome8=8, chromosome9=9, chromosome10=10, chromosome11=11, chromosome12=12, chromosome13=13, chromosome14=14, chromosome15=15, chromosome16=16, chromosomeMT=17))
#snpgdsClose(genofile)
genofile <- snpgdsOpen("genotypes_1011_JT_biall.recode.gds")

# Running PCA
pca <- snpgdsPCA(genofile, num.thread=3)
pc.percent <- pca$varprop*100

# Getting sample ids and population info

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.csv("./input_files/pops_PCA.csv",sep=",",header=TRUE)

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  EV7 = pca$eigenvect[,7],
                  EV8 = pca$eigenvect[,8],
                  stringsAsFactors = FALSE)

tabp <- merge(tab, pop_code, by.x = 'sample.id', by.y = 'VcfName', sort = FALSE)


# Setting colors & shapes

#colors <- read.table('./input_files/pop_colors_1011.csv', sep='\t', header=TRUE)
colors <- data.frame("group" = c("1.Wine/European","1.Wine/European(subclade1)","1.Wine/European(subclade2)","1.Wine/European(subclade3)","1.Wine/European(subclade4)","2.Alpechin","3.Brazilianbioethanol","4.Mediterraneanoak","5.Frenchdairy","6.Africanbeer","7.Mosaicbeer","8.Mixedorigin","9.Mexicanagave","10.FrenchGuianahuman","11.Alebeer","12.WestAfricancocoa","13.Africanpalmwine","14.CHNIII","15.CHNII","16.CHNI","17.Taiwanese","18.FarEastAsia","19.Malaysian","20.CHNV","21.Ecuadorean","22.FarEastRussian","23.NorthAmericanoak","24.Asianislands","25.Sake","26.Asianfermentation","Jean_Talon","M1.Mosaicregion1","M2.Mosaicregion2","M3.Mosaicregion3","Unclustered"),
                     "color" = c("#D55E00","#D55E00","#D55E00","#D55E00","#D55E00","#009E73","#F0E442","#009E73","#009E73","#0072B2","#0072B2","#E69F00","#F0E442","#F0E442","#56B4E9","black","black","black","black","black","black","black","black","black","#F0E442","black","#009E73","black","#009E73","#009E73","black","#CC79A7","#CC79A7","#CC79A7","grey"),
                     "fill" = c("#D55E00","#D55E00","#D55E00","#D55E00","#D55E00","#009E73","#F0E442","#009E73","#009E73","#0072B2","#0072B2","#E69F00","#F0E442","#F0E442","#56B4E9","black","black","black","black","black","black","black","black","black","#F0E442","black","#009E73","black","#009E73","#009E73","black","#CC79A7","#CC79A7","#CC79A7","grey"),
                     "shape" = c(1,0,2,4,3,2,17,1,3,2,4,19,3,19,1,1,6,4,13,8,2,3,0,12,1,5,19,9,0,15,19,1,2,0,1),
                     "pop" = c("01.W","01.W1","01.W2","01.W3","01.W4","02.A","03.B","04.M","05.F","06.A","07.M","08.M","09.M","10.F","11.A","12.W","13.A","14.C","15.C","16.C","17.T","18.F","19.M","20.C","21.E","22.R","23.N","24.A","25.S","26.A","Jean-Talon","M1.M","M2.M","M3.M","unk"))

tabpc <- merge(tabp, colors, by.x = 'Pop', by.y = 'pop', sort = FALSE)
tab_cols <- tabpc %>% select("Pop","color","fill","shape") %>% distinct() %>% arrange(Pop)


# Plotting PC7 vs. PC8 by pop

p <- ggplot(tabp) + aes(x = EV7, y = EV8, color = Pop, fill = Pop, shape = Pop) +
  geom_point(size=2) +
  scale_shape_manual(values = as.vector(tab_cols$shape)) +
  scale_fill_manual(values = as.vector(tab_cols$fill)) +
  scale_color_manual(values = as.vector(tab_cols$color)) +
  scale_x_continuous(limits = c(-0.08,0.19)) +
  labs(x = paste0("PC7 = ",round(pc.percent, 2)[7],"%"), y = paste0("PC8 = ",round(pc.percent, 2)[8],"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        legend.background = element_blank(),
        legend.key =  element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.85,0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=16))
p

# Saving plot

pdf("./output/sample_plotPCA_1011.pdf",w=7,h=5)
p
dev.off()






### Running PCA on Beer collection + Jean-Talon


# Reading vcf files
#vcf.fn <- "genotypes_Beer_JT_biall.recode_AF_PL.vcf.gz"
#vcf.fn <- "FileS2.vcf.gz"
vcf.fn <- "./input_files/sample_Beer_JT_biall.recode_AF_PL.vcf.gz"
#snpgdsClose(genofile)
snpgdsVCF2GDS_R(vcf.fn, "genotypes_Beer_JT_biall.recode.gds", method="biallelic.only",option=snpgdsOption(chr1=1, chr2=2, chr3=3, chr4=4, chr5=5, chr6=6, chr7=7, chr8=8, chr9=9, chr10=10, chr11=11, chr12=12, chr13=13, chr14=14, chr15=15, chr16=16, chrmt=17))
genofile <- snpgdsOpen("genotypes_Beer_JT_biall.recode.gds")


# Running PCA
pca <- snpgdsPCA(genofile, num.thread=3)
pc.percent <- pca$varprop*100


# Getting sample ids and population info

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.table("./input_files/pops_PCA.csv",sep=",",header=TRUE)

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
tabp <- merge(tab, pop_code, by.x = 'sample.id', by.y = "VcfName", sort = FALSE)


# Setting colors & shapes

popcols <- c("black","#56B4E9","#0072B2","#009E73","#E69F00",
             "black","#F0E442","#D55E00","black","black",
             "black","#CC79A7","#0072B2","#009E73","#D55E00",
             "#D55E00","grey")
popfils <- popcols
popshapes <- c(1,1,3,0,19,
               4,1,1,2,3,
               19,1,0,1,15,
               17,1)

# Plotting PC2 vs. PC1 by population

p <- ggplot(tabp) + aes(x = EV1, y = EV2, color = Pop, fill = Pop, shape = Pop) +
  geom_point(size=2) +
  scale_color_manual(values = popcols) +
  scale_shape_manual(values = popshapes) +
  scale_fill_manual(values = popfils) +
  scale_x_continuous(limits = c(-0.08,0.16)) +
  labs(x = paste0("PC1 = ",format(round(pc.percent[1],2)),"%"), y = paste0("PC2 = ",format(round(pc.percent[2],2)),"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        legend.background = element_blank(),
        legend.key =  element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9,0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=16))

p


# Plotting PC3 vs. PC2 by population

p2 <- ggplot(tabp) + aes(x = EV2, y = EV3, color = Pop, fill = Pop, shape = Pop) +
  geom_point(size=2) +
  scale_color_manual(values = popcols) +
  scale_shape_manual(values = popshapes) +
  scale_fill_manual(values = popfils) +
  scale_x_continuous(limits = c(-0.13,0.16)) +
  labs(x = paste0("PC2 = ",format(round(pc.percent[2],2)),"%"), y = paste0("PC3 = ",format(round(pc.percent[3],2)),"%")) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = FALSE),
        legend.background = element_blank(),
        legend.key =  element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9,0.5),
        axis.title = element_text(size=20),
        axis.text = element_text(size=16))

p2


# Saving plot

pdf("./output/sample_plotPCA_beer.pdf",w=7,h=5)
p2
dev.off()
