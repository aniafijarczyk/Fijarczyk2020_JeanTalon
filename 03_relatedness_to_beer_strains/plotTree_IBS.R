setwd("/media/anna/Volume/Yeast_JeanTalon/2020_02_figshare/03_relatedness_to_beer_strains")

#author:aniafijarczyk
library("SNPRelate")
library(ape)




### Reading vcf file with beer strains

vcf.fn <- "./input_files/sample_VarFiltr.vcf.gz"
#snpgdsClose(genofile)
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


### Calculating identity by state matrix and writing down tree

ibs <- snpgdsIBS(genofile)
distSamples <- data.frame("Samples" = ibs$sample.id)
newDF <- merge(distSamples, samples, by.x = "Samples", by.y = "sampleID", sort = FALSE)
names <- newDF$group
mat_ibs <- as.matrix(ibs$ibs)
dimnames(mat_ibs) <- list(names, names)
newmat <- 1 - mat_ibs
tree <- bionj(newmat)
plot(tree, type = "unrooted")
write.tree(tree, file = "sample_VarFiltr_ibs.nwk")




