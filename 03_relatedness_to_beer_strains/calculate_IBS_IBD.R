#author@aniafijarczyk
library("SNPRelate")
library(ape)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)




### Reading vcf file with beer strains

#vcf.fn <- "allsamples_VarFiltr.vcf.gz"
#vcf.fn <- "FileS4.vcf.gz"
vcf.fn <- "./input_files/sample_VarFiltr.vcf.gz"
snpgdsVCF2GDS_R(vcf.fn, "genotypes_allsamples_VarFiltr.gds", method="biallelic.only",
                option=snpgdsOption("ref|NC_001133|"=1, "ref|NC_001134|"=2, 
                                    "ref|NC_001135|"=3, "ref|NC_001136|"=4, 
                                    "ref|NC_001137|"=5, "ref|NC_001138|"=6, 
                                    "ref|NC_001139|"=7, "ref|NC_001140|"=8, 
                                    "ref|NC_001141|"=9, "ref|NC_001142|"=10, 
                                    "ref|NC_001143|"=11, "ref|NC_001144|"=12, 
                                    "ref|NC_001145|"=13, "ref|NC_001146|"=14, 
                                    "ref|NC_001147|"=15, "ref|NC_001148|"=16, "ref|NC_001224|"=17))
#snpgdsClose(genofile)
genofile <- snpgdsOpen("genotypes_allsamples_VarFiltr.gds")




### Getting table with pop info
samples <- read.table("./input_files/table_S2.csv", sep = ",", header=TRUE) %>% 
  select(Strain.name,Genetic.group) %>% 
  unite(group, c(Strain.name, Genetic.group), sep = "_", remove = FALSE)




### Calculating identity by state matrix and writing down tree

ibs <- snpgdsIBS(genofile)
distSamples <- data.frame("Samples" = ibs$sample.id)
newDF <- merge(distSamples, samples, by.x = "Samples", by.y = "Strain.name", sort = FALSE)
names <- newDF$group
mat_ibs <- as.matrix(ibs$ibs)
dimnames(mat_ibs) <- list(names, names)
newmat <- 1 - mat_ibs
tree <- bionj(newmat)
plot(tree, type = "unrooted")
write.tree(tree, file = "./output/sample_VarFiltr_ibs.nwk")



### Calculating IBD KINGMAN for all samples

ibd.robust <- snpgdsIBDKING(genofile, sample.id=samples$Strain.name, useMatrix=TRUE)
samps <- samples$Strain.name
mat <- ibd.robust$kinship  # dspMatrix
colnames(mat) <- ibd.robust$sample.id
rownames(mat) <- ibd.robust$sample.id
mat[mat < 0] <- as.numeric(0)
kin <- 2*mat

# Matrix with kinship coefficients

df_kin <- as.data.frame(as.matrix(kin))
write.table(df_kin,file = "./output/sample_IBDKINGmethod_kinship_matrix.out",quote=FALSE, sep=",")




### Plotting heatmap of kinship coefficients only for the Beer/baking group

# Selecting samples

bf <- read.table("./input_files/table_S2.csv", sep = ",", header=TRUE,na.strings = c("","NA")) %>% 
  select(Strain.name,Genetic.group,Source,BB_clade) %>% 
  filter(Genetic.group == "Beer/baking" | Genetic.group == "Beer/Bread" | Genetic.group == "08.M" | Genetic.group == "Mixed" | Genetic.group == "Jean-Talon")

# Calculating matrix

ibd.robust <- snpgdsIBDKING(genofile, sample.id=bf$Strain.name, useMatrix=TRUE)
plot(ibd.robust$IBS0, ibd.robust$kinship, xlab="Proportion of Zero IBS",ylab="Estimated Kinship Coefficient (KING-robust)")
mat <- ibd.robust$kinship  # dspMatrix
colnames(mat) <- ibd.robust$sample.id
rownames(mat) <- ibd.robust$sample.id
mat[mat < 0] <- as.numeric(0)
kin <- 2*mat
pal3 <- colorRampPalette(brewer.pal(9, "PuBu"))(256)
heatmap(as.matrix(kin), scale="none", col = pal3)

# Setting group colors

levels(bf$BB_clade) <- c(levels(bf$BB_clade), "cladeX") 
bf$BB_clade[is.na(bf$BB_clade)] <- "cladeX"
levels(bf$Source) <- c(levels(bf$Source), "unk") 
bf$Source[is.na(bf$Source)] <- "unk"

colcode <- data.frame("BB_clade" = c("clade1A","clade1B","clade1C","clade1D","cladeX"),
                      "color_clade" = c("grey20","grey40","grey60","grey80","white"))
bf_col <- merge(bf, colcode, by = "BB_clade", sort=FALSE)

typecode <- data.frame("Source" = c("Baking","Beer","Wild","Clinical","Spirits","Wine","unk","Bio-ethanol","Sake","Lab"),
                       "color_source" = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","grey90","grey80","grey50","grey30"))
bf_col2 <- merge(bf_col,typecode,by = "Source",sort=FALSE)

# Putting sample IDs in order

ordered_names <- data.frame("Strain.name" = ibd.robust$sample.id)
so <- merge(ordered_names, bf_col2, by = "Strain.name",sort=FALSE)


# Plotting heatmap

pal3 <- colorRampPalette(brewer.pal(9, "PuBu"))(256)
pdf("./output/sample_plotHeatmap_IBD_sample.pdf",h=5,w=5)
heatmap(as.matrix(kin), scale="none", col = pal3,
        RowSideColors = as.vector(so$color_source),
        ColSideColors = as.vector(so$color_clade))
dev.off()





### Plotting estimates of diversity and divergence

# Data frame with diversity estimates for 5713 single exon non-overlapping genes

dt <- data.frame("Strain"=c("A.2565","A.Muntons","A.S-33","A.T-58","BE005","CFI","CFN","CFP","CHK","Jean-Talon"),
                 "Group_order"=c(8,1,0,7,2,3,4,6,5,9),
                 "PiW/nt"=c(0.000968,0.003076,0.002023,0.001833,0.008096,0.008516,0.008481,0.008551,0.008417,0.008108),
                 "PiA/nt"=c(0.005904,0.002859,0.002264,0.003272,0.005718,0.005865,0.005867,0.007652,0.005817,NA))

df <- dt %>% filter(Strain != "Jean-Talon")


# Plotting estimates

ds <- df[order(df$Group_order),]
ordered <- as.vector(ds$Strain)
df$Strain <- factor(df$Strain, levels = ordered)

p <- ggplot(df) + aes(x = Strain, y = PiW.nt) +
  geom_point(size=5, colour="black") +
  geom_hline(yintercept = as.numeric(dt$PiW.nt[10]), linetype="dashed") +
  geom_point(data=df, aes(x = Strain, y = PiA.nt), size=5, shape=2) +
  labs(y = "synonymus divergence/ \n nucleotide diversity") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_text(size = 18, angle = 45, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top")
p

pdf("./output/plotDiversity.pdf",w=6,h=5)
p
dev.off()






