# Analyses related to estimating relatedness of the Jean-Talon strain to other beer strains (Figure 3)


## Calculating IBS, IBD and plotting figures in Figure 3 
### Complete vcf file for these analyses can be downloaded from figshare as File S4

calculate_IBS_IBD.R


### NJ tree

./input_files/sample_VarFiltr.vcf.gz	# complete vcf file can be downloaded from figshare as File S4
./input_files/table_S2.csv		# metadata for all strains
./output/sample_VarFiltr_ibs.nwk	# tree based on IBS matrix
genotypes_allsamples_VarFiltr.gds	# temporary file


### Matrix of kinship coefficients & heatmap

./input_files/sample_VarFiltr.vcf.gz
./input_files/table_S2.csv

./output/sample_IBDKINGmethod_kinship_matrix.out	# matrix of kinship coefficients between all samples
./output/sample_plotHeatmap_IBD_sample.pdf		# heatmap of kinship coefficients in Beer/baking strains
genotypes_allsamples_VarFiltr.gds			# temporary file


### Estimates of synonymous diversity and divergence in 5713 genes

./output/plotDiversity.pdf		# plot with synonymous summary statistics




# Manipulating fasta files
### Selecting non-overlapping and single exon genes; 
### if vcf with genotypes of relatives is present (can be downloaded from figshare - File S5), 
### after running several bcftools and seqtk commands one can generate msa of concatenated genes

manipulateFasta.ipynb

./input_files/saccharomyces_cerevisiae.gff.gz	# script downloads the file automatically from sgd
./input_files/relatives_annot_Filtered2.samples	# names of haplotypes from the .tab file with genotypes

./output/manipulateFasta_nonoverlappingCDS.bed	# bed file with single exon non-overlapping list of genes




# Calculating divergence time between Jean-Talon and relatives using S. cerevisiae S288C as an outgroup

calcRelativeTime_Scer.ipynb

./input_files/relatives_annot_Filtered2.samples		# names of haplotypes from the .tab file with genotypes
./input_files/relatives_annot_Filtered2_01.tab.gz	# .tab file with genotypes
./input_files/relatives_annot_synonymous_snpEff.tab.gz	# list of synonymous variants

./output/calcRelativeTime_Scer_sample.out
