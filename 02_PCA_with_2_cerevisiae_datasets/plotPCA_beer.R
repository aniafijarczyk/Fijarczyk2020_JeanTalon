setwd("/media/anna/Volume/Yeast_JeanTalon/2020_02_figshare/02_Figure2_PCA_with_2_cerevisiae_datasets")


# author:aniafijarczyk
library("SNPRelate")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)



### Running PCA on Beer collection + Jean-Talon



# Reading vcf files
vcf.fn <- "./input_files/genotypes_Beer_JT_biall.recode_AF_PL.vcf.gz"
snpgdsVCF2GDS_R(vcf.fn, "genotypes_Beer_JT_biall.recode.gds", method="biallelic.only",option=snpgdsOption(chr1=1, chr2=2, chr3=3, chr4=4, chr5=5, chr6=6, chr7=7, chr8=8, chr9=9, chr10=10, chr11=11, chr12=12, chr13=13, chr14=14, chr15=15, chr16=16, chrmt=17))
genofile <- snpgdsOpen("genotypes_Beer_JT_biall.recode.gds")


# Running PCA
pca <- snpgdsPCA(genofile, num.thread=3)
pc.percent <- pca$varprop*100


# Getting sample ids and population info

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.table("./input_files/pop_info_beer.txt",sep="\t",header=TRUE)
colnames(pop_code) <- c("Strain","sample.id","Location","Country","Continent","Type","Source","pop_code")
pop_codev <- as.vector(pop_code$pop_code)
country <- as.vector(pop_code$Country)
head(cbind(sample.id, pop_codev))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)
tabp <- merge(tab, pop_code, by = 'sample.id', sort = FALSE)


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

p <- ggplot(tabp) + aes(x = EV1, y = EV2, color = pop_code, fill = pop_code, shape = pop_code) +
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

p2 <- ggplot(tabp) + aes(x = EV2, y = EV3, color = pop_code, fill = pop_code, shape = pop_code) +
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

pdf("plotPCA_beer.pdf",w=7,h=5)
p2
dev.off()
