setwd("/media/anna/Volume/Yeast_JeanTalon/2020_02_figshare/02_Figure2_PCA_with_2_cerevisiae_datasets")


# author:aniafijarczyk
library("SNPRelate")
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)


### Running PCA on 1011 yeast collection + Jean-Talon


# Reading vcf files
vcf.fn <- "./input_files/genotypes_1011_JT_biall.recode_AF_PL.vcf.gz"
snpgdsVCF2GDS_R(vcf.fn, "genotypes_1011_JT_biall.recode.gds", method="biallelic.only",option=snpgdsOption(chromosome1=1, chromosome2=2, chromosome3=3, chromosome4=4, chromosome5=5, chromosome6=6, chromosome7=7, chromosome8=8, chromosome9=9, chromosome10=10, chromosome11=11, chromosome12=12, chromosome13=13, chromosome14=14, chromosome15=15, chromosome16=16, chromosomeMT=17))
#snpgdsClose(genofile)
genofile <- snpgdsOpen("genotypes_1011_JT_biall.recode.gds")


# Running PCA
pca <- snpgdsPCA(genofile, num.thread=3)
pc.percent <- pca$varprop*100


# Getting sample ids and population info

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.csv("./input_file/pop_info_1011.csv",sep=",",header=TRUE)
pop_codev <- as.vector(pop_code$Abbreviation)
country <- as.vector(pop_code$Geographical.origins)
head(cbind(sample.id, pop_codev))

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

colors <- read.table('pop_colors_1011.csv', sep='\t', header=TRUE)
tabpc <- merge(tabp, colors, by.x = 'Abbreviation', by.y = 'pop', sort = FALSE)
tab_cols <- tabpc[c(1,20,21,22)] %>% distinct() %>% arrange(Abbreviation)


# Plotting PC7 vs. PC8 by pop

p <- ggplot(tabp) + aes(x = EV7, y = EV8, color = Abbreviation, fill = Abbreviation, shape = Abbreviation) +
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

pdf("plotPCA_1011.pdf",w=7,h=5)
p
dev.off()

