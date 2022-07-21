#WP4 qPCR data (R-script, R version 4.1.1)
setwd("S:/Mibi/Produkte/Sainur/WP4/qPCR_data_WP")
library(car)
library(dplyr)
library(agricolae)

data = read.csv("S:/Mibi/Produkte/Sainur/WP4/qPCR_data_WP/16S_Bacterial_qPCR.CSV")
head(data)


#Day as factor
n_Day  <- factor(data$Day)
n_Temp <- factor(data$Temp.)
library(Rmisc)
tg <- summarySE(data, measurevar="B_no_gc_per_g_soil", groupvars=c("Sample","n_Day","n_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)
pd <- position_dodge(0.001) # move them .05 to the left and right
p=ggplot(tg, aes(x=n_Day, y=B_no_gc_per_g_soil , group=n_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "free_y") +
  
  geom_errorbar(aes(ymin= B_no_gc_per_g_soil -se, ymax= B_no_gc_per_g_soil +se), colour = "black", position = pd)+
  geom_line(aes(linetype= n_Temp, colour= n_Temp),position = pd )+
  geom_point(aes(shape=n_Temp, colour= n_Temp)) +
  
  xlab ("Day")+
  scale_x_discrete()+
  
  
  ylab ("Bacterial abundance per g of Dry soil")+
  expand_limits(y=0) +
  #  scale_y_continuous()+
  scale_y_continuous()+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 7, angle =0))

pp <- p+ scale_colour_manual(values = c("#EF6C00", "#0000FF"))
  
pp  
  


#reorder Bacteria

#Day as factor
b_Day  <- factor(data$Day)
b_Temp <- factor(data$Temp.)
library(Rmisc)
tg <- summarySE(data, measurevar="X1ng_log10", groupvars=c("Sample","b_Day","b_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)


data_new <- tg
data_new$Sample <- factor(data_new$Sample, levels= c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20"))

pd <- position_dodge(0.001) # move them .05 to the left and right
p=ggplot(data_new, aes(x=b_Day, y=X1ng_log10, group=b_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "fixed") +
  
  geom_errorbar(aes(ymin= X1ng_log10-se, ymax= X1ng_log10+se), colour = "black", position = pd)+
  geom_line(aes(linetype= b_Temp, colour= b_Temp),position = pd )+
  geom_point(aes(shape=b_Temp, colour= b_Temp)) +
  scale_linetype_manual(values=c("solid", "solid"))+
  
  xlab ("Day")+
  scale_x_discrete()+
  
  
  ylab ("Bacterial abundance per ng of DNA (log10)")+
  expand_limits(y=0) +
  #  scale_y_continuous()+
  scale_y_continuous(limits = c(2, 6.75))+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 8, angle =0))

pp <- p+ scale_colour_manual(values = c("#0000FF", "#FF0000"))

pp  

#print size 10x8

#per g soil
#Day as factor
b_Day  <- factor(data$Day)
b_Temp <- factor(data$Temp.)
library(Rmisc)
tg <- summarySE(data, measurevar="B_no_gc_per_g_soil_log10", groupvars=c("Sample","b_Day","b_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)


data_new <- tg
data_new$Sample <- factor(data_new$Sample, levels= c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20"))

pd <- position_dodge(0.001) # move them .05 to the left and right
p=ggplot(data_new, aes(x=b_Day, y=B_no_gc_per_g_soil_log10, group=b_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "fixed") +
  
  geom_errorbar(aes(ymin= B_no_gc_per_g_soil_log10-se, ymax= B_no_gc_per_g_soil_log10+se), colour = "black", position = pd)+
  geom_line(aes(linetype= b_Temp, colour= b_Temp),position = pd )+
  geom_point(aes(shape=b_Temp, colour= b_Temp)) +
  scale_linetype_manual(values=c("solid", "solid"))+
  
  xlab ("Day")+
  scale_x_discrete()+
  
  
  ylab ("Bacterial abundance per g of dry soil (log10)")+
  expand_limits(y=0) +
  #  scale_y_continuous()+
  scale_y_continuous(limits = c(8, 11))+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 8, angle =0))

pp <- p+ scale_colour_manual(values = c("#0000FF", "#FF0000"))

pp  


  
  
#Archeal

data2 = read.csv("S:/Mibi/Produkte/Sainur/WP4/qPCR_data_WP/16S_Archaeal_qPCR.csv")
head(data2)


#figure size: 10x6 8()


#reorder (Archaea)
#Day as factor
head(data2)
a_Day  <- factor(data2$Day)
a_Temp <- factor(data2$Temp.)
library(Rmisc)
tg <- summarySE(data2, measurevar="X1ng_log10", groupvars=c("Sample","a_Day","a_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)

data_new <- tg
data_new$Sample <- factor(data_new$Sample, levels= c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20"))


pd <- position_dodge(0.001) # move them .05 to the left and right
p=ggplot(data_new, aes(x=a_Day, y=X1ng_log10, group=a_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "fixed") +
  
  geom_errorbar(aes(ymin= X1ng_log10-se, ymax= X1ng_log10+se), colour = "black", position = pd)+
  geom_line(aes(linetype= a_Temp, colour= a_Temp),position = pd )+
  geom_point(aes(shape=a_Temp, colour= a_Temp)) +
  scale_linetype_manual(values=c("solid", "solid"))+
  
  xlab ("Day")+
  scale_x_discrete()+
  
  
  ylab ("Archeal abundance per ng of DNA (log10)")+
  expand_limits(y=0) +
  #  scale_y_continuous()+
  scale_y_continuous(limits = c(2, 5))+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 8, angle =0))

pp <- p+ scale_colour_manual(values = c("#0000FF", "#FF0000"))
  
pp  

#print 10x8
#per g of osil
head(data2)
a_Day  <- factor(data2$Day)
a_Temp <- factor(data2$Temp.)
library(Rmisc)
tg <- summarySE(data2, measurevar="A_no_gc_per_g_soil_log10", groupvars=c("Sample","a_Day","a_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)

data_new <- tg
data_new$Sample <- factor(data_new$Sample, levels= c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20"))


pd <- position_dodge(0.001) # move them .05 to the left and right
p=ggplot(data_new, aes(x=a_Day, y=A_no_gc_per_g_soil_log10, group=a_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "fixed") +
  
  geom_errorbar(aes(ymin= A_no_gc_per_g_soil_log10-se, ymax= A_no_gc_per_g_soil_log10+se), colour = "black", position = pd)+
  geom_line(aes(linetype= a_Temp, colour= a_Temp),position = pd )+
  geom_point(aes(shape=a_Temp, colour= a_Temp)) +
  scale_linetype_manual(values=c("solid", "solid"))+
  
  xlab ("Day")+
  scale_x_discrete()+
  
  
  ylab ("Archeal abundance per of dry soil (log10)")+
  expand_limits(y=0) +
  #  scale_y_continuous()+
#  scale_y_continuous(limits = c(0, 10))+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 8, angle =0))

pp <- p+ scale_colour_manual(values = c("#0000FF", "#FF0000"))

pp 




#Fungal abundance (qPCR)

data3 = read.csv("S:/Mibi/Produkte/Sainur/WP4/qPCR_data_WP/ITS_Fungal_qPCR.csv")
head(data3)


#figure size: 10x6 8()


#Day as factor
head(data3)
f_Day  <- factor(data3$Day)
f_Temp <- factor(data3$Temp.)
library(Rmisc)
tg <- summarySE(data3, measurevar="X1ng_log10", groupvars=c("Sample","f_Day","f_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)
#neworder <- c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20")
pd <- position_dodge(0.001) # move them .05 to the left and right

data_new <- tg
data_new$Sample <- factor(data_new$Sample, levels= c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20"))



p=ggplot(data_new, aes(x=f_Day, y=X1ng_log10, group=f_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "fixed") +
  
  geom_errorbar(aes(ymin= X1ng_log10-se, ymax= X1ng_log10+se), colour = "black", position = pd)+
  geom_line(aes(linetype= f_Temp, colour= f_Temp),position = pd )+
  geom_point(aes(shape=f_Temp, colour= f_Temp)) +
  scale_linetype_manual(values=c("solid", "solid"))+
  
  xlab ("Day")+
  scale_x_discrete()+
  
  
  ylab ("Fungal abundance per ng of DNA (log10)")+
  expand_limits(y=0) +
#  scale_y_continuous()+
  scale_y_continuous(limits = c(2, 5))+
  
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 8, angle =0))

pp <- p+ scale_colour_manual(values = c("#0000FF", "#FF0000"))

pp  

#print size (10x8)


#per g of dry soil

#Day as factor
head(data3)
f_Day  <- factor(data3$Day)
f_Temp <- factor(data3$Temp.)
library(Rmisc)
tg <- summarySE(data3, measurevar="F_no_gc_per_g_soil_log10", groupvars=c("Sample","f_Day","f_Temp"))
tg
require(grid)
require(plyr)
library(ggplot2)
#neworder <- c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20")
pd <- position_dodge(0.001) # move them .05 to the left and right

data_new <- tg
data_new$Sample <- factor(data_new$Sample, levels= c("S1", "S2", "S3", "S4","S5", "S6", "S7", "S8","S9", "S10", "S11", "S12","S13", "S14", "S15", "S16","S17", "S18", "S19", "S20"))



p=ggplot(data_new, aes(x=f_Day, y=F_no_gc_per_g_soil_log10, group=f_Temp)) + 
  
  facet_wrap(~Sample, ncol = 4, scales = "fixed") +
  
  geom_errorbar(aes(ymin= F_no_gc_per_g_soil_log10-se, ymax= F_no_gc_per_g_soil_log10+se), colour = "black", position = pd)+
  geom_line(aes(linetype= f_Temp, colour= f_Temp),position = pd )+
  geom_point(aes(shape=f_Temp, colour= f_Temp)) +
  scale_linetype_manual(values=c("solid", "solid"))+
  
  xlab ("Day")+ 
  scale_x_discrete()+
  
  
  ylab ("Fungal abundance per g of dry soil (log10)")+
  expand_limits(y=0) +
  #  scale_y_continuous()+
  scale_y_continuous(limits = c(6, 9))+
  
  theme_bw()+
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 8, angle =0))

pp <- p+ scale_colour_manual(values = c("#0000FF", "#FF0000"))

pp  
