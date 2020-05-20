
##R version 3.6.1  was used to run the following scripts
## author: Souhir Marsit
#####Libraries used
library(lemon)
library(ggplot2)
library(cowplot)
library(lattice)
library(dplyr)
library(plyr)
library(lemon)
library(tidyr)
library(data.table)
library(magrittr)
library(stringr)
#library(UpSetR)
#library(ggpubr)
library(gdata)
#library(filesstrings)
library(gplots)
rm(list=ls())


#calculate max slope
data_all = read.table("./input_files/Cobalt_resistance_JT_raw_15_05_2020.txt", header = T, sep = '\t')
data_all = filter (data_all, strain != "bord")

cal_slope <- function(size_vect, time_vect, window.regression = 5){
  #I create a temporary vector for the slopes and the intercepts
  slopes = NA
  #it corresponds to the mean of the window
  fact_nas = sum(!is.na(size_vect))/length(size_vect)
  if (fact_nas > 0.8){ 
    
    for(j in 1:(length(size_vect) - window.regression + 1)){
      #I determine the start and end positions of the current window
      start = j
      end = j + window.regression - 1
      #I calculate the regression
      regression = lm(size_vect[start:end] ~ time_vect[start:end])
      slope = regression$coefficient[2]
      #I fill the temporary vectors
      slopes = c(slopes, slope)
    }
    #assemble regression parameters and coose the defined percentile for each
    parameters = data.frame(slopes)
    parameters = parameters[order(parameters$slopes),]
    slope = parameters[length(parameters)-1]
  }
  
  else { slope = NA }
  return(slope)
}


#first transform time from second to hour
data_all$time=as.numeric(data_all$time)
datah = mutate (data_all, "timeh"= (time/3600))
fdata = select (datah, strain, sample, concentration, timeh, OD, plate)

#calculate slopes for all data
#slope per strain, per plate through time

df_slopes = fdata %>% group_by(plate, sample, strain, concentration) %>%
  dplyr::summarise(max_slope = cal_slope(OD, timeh),
                   max_size = nth(OD, -3,order_by=OD))%>% arrange(plate, sample, strain, concentration)

write.table(df_slopes, file="./output/phenotyp_JT_slope_tot.txt", sep = '\t') 


#calculate median of all replicates

df_slopes = read.table("./output/phenotyp_JT_slope_tot.txt", header = T, sep = '\t')

df_sum <- df_slopes %>% group_by(strain, concentration) %>%
  dplyr::summarize(median_max_slope = median(max_slope),
                   median_max_size = median(max_size),
                   sd_slope= sd(max_slope),
                   sd_size= sd(max_size)) 


write.table(df_sum, file="./output/phenotyp_JT_slope_sd.txt", sep = '\t') 

##Figure 4 B####Growth rate##############################################################

data = read.table("./output/phenotyp_JT_slope_sd.txt", header = T, sep = '\t')
data $concentration_f = factor(data $concentration, levels=c("0", "2", "4", "6", "8", "10"))
data $Strain = factor(data $strain, levels=c("BY_1n", "cot1", "BY_2n", "LL13_40", "LL13_54", "Jean_Talon", "London", "Windsor"))


pdf(("./output/Phenotypage_JT_fig_11_05_2020.pdf"), height = 7, width =11 )

ggplot(data, aes(x = concentration_f, y = as.numeric(median_max_slope), ymin=median_max_slope-sd_slope, ymax=median_max_slope+sd_slope, group=Strain, color=Strain)) + 
  geom_line(alpha=5/10, size=1) +
  geom_pointrange(fill = "black", alpha=7/10, size=1) + 
  scale_color_manual(values=c("royalblue2", "grey20", "navyblue", "seagreen", "yellowgreen", "Orangered2", "goldenrod3", "orange2"))+
  ylab("Growth rate (OD/h)")+ theme(axis.text=element_text(size=18))+
  xlab("Cobalt concentration (mM)")+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) +
  theme(panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(legend.title=element_text(size=18), legend.text=element_text(size=16))
dev.off()

##Figure_4C##############Growth curve####################################################################

data_all = read.table("./input_files/Cobalt_resistance_JT_raw_15_05_2020.txt", header = T, sep = '\t')
data_all = filter (data_all, strain != "bord")
data_all = mutate (data_all, "timeh"= (time/3600))
data_6mM = filter (data_all, concentration == "6")

data_6mM $Strain = factor(data_6mM $strain, levels=c("BY_1n", "cot1", "BY_2n", "LL13_40", "LL13_54", "Jean_Talon", "London", "Windsor"))

COLS = c("royalblue2", "grey20", "navyblue", "seagreen", "yellowgreen", "Orangered2", "goldenrod3", "orange2")
names(COLS) = c("BY_1n","cot1","BY_2n", "LL13_40", "LL13_54", "Jean_Talon", "London", "Windsor")

pdf("./output/Phenotypage_JT_fig_G_curve_11_05_2020.pdf", height = 7, width =11 )
ggplot(data_6mM)  + 
  geom_smooth(aes(x = timeh, y = OD, colour= Strain, fill=Strain)) + 
  #scale_color_manual(values=c("royalblue2", "grey20", "navyblue", "seagreen", "yellowgreen", "Orangered2", "goldenrod3", "orange2"))+
  scale_fill_manual(name="Strain",values=COLS)+
  scale_color_manual(name="Strain",values=COLS)+
  ylab("OD")+ theme(axis.text=element_text(size=11))+
  xlab("Time (h)")+
  theme(axis.title.y=element_text(size=18))+
  theme(axis.title.x=element_text(size=18))+
  theme(axis.text.x=element_text(size=18, colour = "black"))+ theme(axis.text.y=element_text(size=18, colour = "black")) +
  theme(panel.background = element_blank(), axis.line=element_line(colour="black"))+
  theme(legend.title=element_text(size=18), 
        legend.text=element_text(size=16))
dev.off()
