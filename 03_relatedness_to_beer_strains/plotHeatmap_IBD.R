setwd("/media/anna/Volume/Yeast_JeanTalon/2020_02_figshare/03_relatedness_to_beer_strains")

#author:aniafijarczyk
library("SNPRelate")
library(dplyr)
library(RColorBrewer)



### Reading vcf file with beer strains

vcf.fn <- "./input_files/sample_VarFiltr.vcf.gz"
snpgdsClose(genofile)
snpgdsVCF2GDS_R(vcf.fn, "genotypes_allsamples_VarFiltr.gds", method="biallelic.only",
                option=snpgdsOption("ref|NC_001133|"=1, "ref|NC_001134|"=2, 
                                    "ref|NC_001135|"=3, "ref|NC_001136|"=4, 
                                    "ref|NC_001137|"=5, "ref|NC_001138|"=6, 
                                    "ref|NC_001139|"=7, "ref|NC_001140|"=8, 
                                    "ref|NC_001141|"=9, "ref|NC_001142|"=10, 
                                    "ref|NC_001143|"=11, "ref|NC_001144|"=12, 
                                    "ref|NC_001145|"=13, "ref|NC_001146|"=14, 
                                    "ref|NC_001147|"=15, "ref|NC_001148|"=16, "ref|NC_001224|"=17))
genofile <- snpgdsOpen("genotypes_allsamples_VarFiltr.gds")


### Getting table with pop info
samples <- read.table("./input_files/samples_pop.txt", sep = "\t", header=TRUE)





### Calculating IBD KINGMAN

# Selecting samples

bf <- samples %>% filter(Groups == "Beer/baking" | Groups == "Beer/Bread" | Groups == "08.M" | Groups == "Mixed" | Groups == "Jean-Talon")

# Calculating matrix

ibd.robust <- snpgdsIBDKING(genofile, sample.id=bf$sampleID, useMatrix=TRUE)
plot(ibd.robust$IBS0, ibd.robust$kinship, xlab="Proportion of Zero IBS",ylab="Estimated Kinship Coefficient (KING-robust)")
samps <- bf$sampleID
mat <- ibd.robust$kinship  # dspMatrix
colnames(mat) <- samps
rownames(mat) <- samps
mat[mat < 0] <- as.numeric(0)
kin <- 2*mat


# Setting group colors

pal5 <- c("grey20","grey40","grey60","grey80","grey90")
pal3 <- colorRampPalette(brewer.pal(9, "PuBu"))(256)
sg <- read.table("./input_files/beerbakery_groups.txt",sep="\t",header=TRUE)
samps_gr <- merge(bf,sg, by.x = "sampleID", by.y = "sample",all.x = TRUE, sort = FALSE)
levels(samps_gr$gr) <- c(levels(samps_gr$gr), "GroupX") 
samps_gr$gr[is.na(samps_gr$gr)] <- "GroupX"
colcode <- data.frame("gr" = c("Group1","Group2","Group3","Group4","GroupX"),"color" = c(pal5))
samps_gr_col <- merge(samps_gr, colcode, by = "gr", sort=FALSE)
sgc <- merge(bf, samps_gr_col, by = "sampleID", sort=FALSE)

types2 <- read.table("./input_files/beer_types2.txt",sep="\t",header=TRUE)
sgt <- merge(sgc, types2, by.x = "sampleID",by.y="sample", sort=TRUE)
typecode <- data.frame("Type2" = c("Baking","Beer","Wild","Clinical","Spirits","Wine","unk","Bio-ethanol","Sake","Lab"),
                       "color" = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","grey90","grey80","grey50","grey30"))
stt <- merge(sgt,typecode,by.x = "Type2.y", by.y="Type2",sort=FALSE)
so <- merge(bf,stt[c(1,2,5,14,15)],by.x="sampleID",by.y="sampleID",sort=FALSE)


# Plotting kinship matrix

pdf("plotHeatmap_IBD_sample.pdf",h=5,w=5)
heatmap(as.matrix(kin), scale="none", col = pal3,
        RowSideColors = as.vector(so$color.y),
        ColSideColors = as.vector(so$color.x))
dev.off()

# Table of kinship coefficients of each strain with Jean-Talon

df <- data.frame("kinship"=as.vector(kin[,'Jean-Talon']),"sample"=samps)
write.table(df, file="IBDKINGmethod_kinship_with_JT_sample.out",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE)





### Calculating IBD KINGMAN for all samples

# Calculating matrix

bf <- samples
ibd.robust <- snpgdsIBDKING(genofile, sample.id=bf$sampleID, useMatrix=TRUE)
samps <- bf$sampleID
mat <- ibd.robust$kinship  # dspMatrix
colnames(mat) <- samps
rownames(mat) <- samps
mat[mat < 0] <- as.numeric(0)
kin <- 2*mat

# Matrix with kinship coefficients

df_kin <- as.data.frame(as.matrix(kin))
write.table(df_kin,file = "IBDKINGmethod_kinship_matrix_sample.out",quote=FALSE, sep=",")




