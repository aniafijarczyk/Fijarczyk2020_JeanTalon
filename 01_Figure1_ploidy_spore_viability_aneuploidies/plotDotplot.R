# author:aniafijarczyk
library(ggplot2)
library(dplyr)
library(tidyr)


df <- read.table("./input_files/sporulation_data.tab", sep = "\t", header = TRUE)
dd <- df %>% gather(key = "Replicate", value = "Spores", JT1:JT3)


p <-ggplot(dd) + aes(x = Replicate, y = Spores) +
  geom_boxplot(fill = "grey80") +
  geom_dotplot(binaxis='y', stackdir='center',dotsize=1.2, stackratio = 1) +
  labs(x = "Replicates", y = "No. viable spores per tetrad") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18))
p

pdf("plotDotplot.pdf",w=7,h=4)
p
dev.off()
