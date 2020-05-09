setwd("/media/anna/Volume/Yeast_JeanTalon/2020_02_figshare/figure1")

# author:aniafijarczyk
library(ggplot2)
library(dplyr)


# Reading & filtering input
fname = "./input_files/formatInput_Jean-Talon.tab"
df <- read.table(fname, sep = "\t", header = TRUE)
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
pdf("plotCNVsInGenome.pdf",w=12,h=3)
p
dev.off()

