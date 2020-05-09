setwd("/media/anna/Volume/Yeast_JeanTalon/2020_02_figshare/figure1")

# author:aniafijarczyk
library(dplyr)
library(ggplot2)


df <- read.table("./input_files/SNP_frequencies.tab", sep = "\t", header=TRUE)


p <- ggplot(df) + aes(x = FreqAlt) +
  geom_histogram(colour = "white", fill = "black") +
  labs(x = "Variant frequency", y = "No. variants") +
  #ggtitle(label = " Frequency of mapped variants") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20))
p


pdf("plotFreq.pdf",w=7,h=7)
p
dev.off()

