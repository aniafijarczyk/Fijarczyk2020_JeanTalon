#author:aniafijarczyk
library(ggplot2)
library(dplyr)
library(tidyr)



# Reading files with estimates

dt <- read.table("./input_files/diversity_6041genes.tab", sep = "\t", header = TRUE)
df <- dt %>% filter(Strain != "Jean-Talon")


# Plotting estimates

ds <- df[order(df$Pop),]
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

pdf("plotDiversity.pdf",w=6,h=5)
p
dev.off()

