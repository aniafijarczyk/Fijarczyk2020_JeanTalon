
### Plotting ploidy from DNA staining of the Jean-Talon colonies

# author:Souhir Marsit
library(ggplot2)
library(dplyr)


tab_total = read.table(gzfile('./input_files/Ploidy_JT_tot.txt.gz'), header = T, sep = '\t')
tab_total= filter (tab_total, GRN.B.HLog> 2000)
tab_ech= filter (tab_total, ech != 90)
tab_ech= filter (tab_ech, ech != 93)
tab_ech= filter (tab_ech, ech !=91)
tab_ech= filter (tab_ech, ech != 92)
tab_ech= filter (tab_ech, ech != 89)
tab_ech= filter (tab_ech, ech != 96)
tab_ech= filter (tab_ech, ech != 95)
tab_ech= filter (tab_ech, ech != 94)
tab_ech= filter (tab_ech, GRN.B.HLog> 5000)
tab_ech= filter (tab_ech, ech != "ech")
tab_ech = mutate (tab_ech, "rep"= "sampl")
tab_ctl1= filter(tab_total, ech==90 )
tab_ctl2= filter(tab_total, ech==93 )
tab_ctl1 = mutate (tab_ctl1, "rep"= "ctl1")
tab_ctl2 = mutate (tab_ctl2, "rep"= "ctl2")
tab_tot = rbind (tab_ctl1, tab_ctl2, tab_ech)


p <- ggplot(tab_tot, aes(x = GRN.B.HLog, colour = rep, group = ech, linetype = rep)) +
  geom_density() + 
  scale_colour_manual(values = c(ctl1 = "black", ctl2 = "black", sampl= "black", name="")) + 
  scale_linetype_manual(values = c(ctl1 = "dotted", ctl2 = "longdash", sampl = "solid", name="")) +
  xlim(0, 35000) +
  labs(x = "Fluorescence (a.u.)", y = "Density") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "None",
        legend.title = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))
p


pdf("./output/plotPloidy.pdf", onefile = T, width = 7, height = 7)
p
dev.off()



### Plotting spore viability

# author:aniafijarczyk
library(ggplot2)
library(dplyr)
library(tidyr)

#df <- read.table("./input_files/sporulation_data.tab", sep = "\t", header = TRUE)

df <- data.frame("JT1" = c(1,2,0,2,0,2,1,2,1,2,0,1,3,1,2,1,3,0,0,2,1,2,1,0),
                 "JT2" = c(2,3,1,1,2,2,1,2,2,0,1,1,2,2,1,1,0,1,2,1,0,2,1,0),
                 "JT3" = c(2,1,0,2,1,4,1,2,2,2,0,0,0,2,2,1,2,1,1,2,0,2,1,1))
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

pdf("./output/plotDotplot.pdf",w=7,h=4)
p
dev.off()



### Plotting frequencies of the Jean-Talon variants


df <- read.table(gzfile("./input_files/SNP_frequencies.tab.gz"), sep = "\t", header=TRUE)

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

pdf("./output/plotFreq.pdf",w=7,h=7)
p
dev.off()



### Plotting coverage profile across Jean-Talon genome


# Reading & filtering input

df <- read.table(gzfile("./input_files/formatInput_Jean-Talon.tab.gz"), sep = "\t", header = TRUE)
dff <- df %>% dplyr::filter(row_number() %% 3 != 1) ## Delete every 3rd row starting from 1


# Plotting coverage
p <- ggplot(dff)  +
  geom_rect(data=df,aes(xmin = win_stop-250, xmax = win_stop, ymin = -Inf, ymax = Inf, fill = chrom)) +
  scale_fill_manual(values = rep(c("white","grey80"),8)) +
  geom_point(aes(x = win_stop-125, y = cn, colour = type), size = 0.1) +
  scale_colour_manual(values = c("#D55E00", "#0072B2","grey20")) +
  #facet_wrap(~chrom, ncol = 4) +
  scale_x_continuous(expand = c(0, 0), breaks = c(115132.5,637066.5,1202495,2127210,3182043,3606023,4286949.5,5114076,5615542,6208748.5,6915625.5,7788439.5,8789909.5,9644448.5,10582321.5,11602066.5),
                     labels = c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")) +
  scale_y_continuous(limits = c(0,8)) +
  #geom_text(data = dtext, aes(x = pos.x, y = pos.y, label = label)) +
  labs(x = "Chromosomes", y = "Norm. copy number") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank())
p


# Saving plot
pdf("./output/plotCNVsInGenome.pdf",w=12,h=3)
p
dev.off()

