# author:Souhir Marsit
library(ggplot2)
library(dplyr)


tab_total = read.table('./input_files/Ploidy_JT_tot.txt', header = T, sep = '\t')
#tab_total= filter (tab_total, GRN.B.HLog> 2000)
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


pdf("plotPloidy.pdf", onefile = T, width = 7, height = 7)
p
dev.off()
