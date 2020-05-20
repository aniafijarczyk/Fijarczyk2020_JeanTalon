# Code for analysing structural variation with long reads

S288c.all_feature.gff	# GFF file of all genomic annotations produced by Yue et al. Nat. Genet. 2017. 
			# Downloaded from https://yjx1217.github.io/Yeast_PacBio_2016/data/
S288c.genome.fa		# Genome of the S.cerevisiae strain S288C produced by Yue et al. Nat. Genet. 2017. 
			# Downloaded from https://yjx1217.github.io/Yeast_PacBio_2016/data/
SAM.pkl			# python pickle serialized file containing the SAM records of supplementary alignments of long-reads for the detection of translocations
SAM_pivot.pkl		# Reformatted version of SAM.pkl
SVs.ipynb		# Main script
rl.pkl			# python pickle serialized file containing the distributions of read lengths for long read datasets
svim/A2565.pb/final_results.vcf			# results of SVIM for strain A2565
svim/T58.pb/final_results.vcf			# results of SVIM for strain T58
svim/barcode11.ont/final_results.vcf		# results of SVIM for strain Jean-Talon
svim/barcode11_sub.ont/final_results.vcf	# results of SVIM for strain Jean-Talon (subset of reads)
svim/s288c.pb/final_results.vcf			# results of SVIM for strain S288C
