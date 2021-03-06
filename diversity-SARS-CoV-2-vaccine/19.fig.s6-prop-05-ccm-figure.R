
library(showtext)
showtext_auto()
font_add("ArialN","arial.ttf")
font.families()
rm(list = ls())
# par(mar = c(4,4,2,4))
# par(mfrow = c(3,1))
library(ggpubr) # ???ذ?
library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(readxl)
library(lubridate)
library(gtable)
library(ggplot2)
library(reshape)
display.brewer.all(type="seq")
library(scales)
library(showtext)
showtext_auto()
font_add("ArialN","arial.ttf")
font.families()
c1="#f79c1d"  #PHSM 
c3="#1426a4"  #SI 
c4="#3a9679"  #GR 
c2= (brewer.pal(8, "Greys"))[7]
c5= (brewer.pal(8, "Blues"))[7]
windowsFonts(A=windowsFont("Times New Roman"),B=windowsFont("Arial"))

DHS1<- read.csv("ccm_results2315/GR_xmap_NPI_rho_statistics.csv",header=FALSE)
DHS2<- read.csv("ccm_results2315/NPI_xmap_GR_rho_statistics.csv",header=FALSE)


meanP1=DHS1$V3
lowP1=DHS1$V4
highP1=DHS1$V2
meanP2=DHS2$V3
lowP2=DHS2$V4
highP2=DHS2$V2

x=DHS1$V1

dfR=data.frame(x,meanP1,meanP2,lowP1,highP1,lowP2,highP2)
dfR[dfR<0]<-0
#A
gorder<-c("PHSM forces Growth Rate","Growth Rate forces PHSM")


TH1<-ggplot(dfR,aes(x)) +
  geom_ribbon(aes(ymin = lowP1, ymax = highP1),fill=c1,alpha=0.4) +
  geom_ribbon(aes(ymin = lowP2, ymax = highP2),fill=c2,alpha=0.4)+
  geom_line(aes(x=x,y=meanP1,color="PHSM forces Growth Rate"),size=0.8)+
  geom_line(aes(x=x,y=meanP2,color="Growth Rate forces PHSM"),size=0.8)+
  theme_classic()+
  xlab(NULL) + 
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(name="Cross Mapping Skill(p)",limits = c(0, 0.8),
                     breaks = seq(0,1,0.2))+
  theme(panel.background = element_rect( colour = "black", size = 0.8)) +
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(name="Library Size",   limits = c(0, 1180),
                     breaks = seq(0,1180,400),expand = c(0, 0))+
  scale_colour_manual(values=c(c1,c2),labels=gorder,breaks=gorder)+
  theme(axis.text.x = element_text(size = 14, color="black",vjust = 0.5, angle = 0))+
  theme(axis.text.y = element_text(size = 14, color="black",vjust = 0.5))+
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position = c(0.295,0.96),legend.direction = "vertical",legend.key.size=unit(0.4,'cm'))+
  labs(color="")+
  theme(plot.title = element_text(hjust = 0.5,size=12)) +
  theme(plot.margin=unit(c(1,0.3,1,0.3),'lines'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(text=element_text(size=14,  family="serif"))
TH1


DHS1<- read.csv("ccm_results2315/GR_xmap_Vero_rho_statistics.csv",header=FALSE)
DHS2<- read.csv("ccm_results2315/Vero_xmap_GR_rho_statistics.csv",header=FALSE)


meanP1=DHS1$V3
lowP1=DHS1$V4
highP1=DHS1$V2
meanP2=DHS2$V3
lowP2=DHS2$V4
highP2=DHS2$V2

x=DHS1$V1

dfR=data.frame(x,meanP1,meanP2,lowP1,highP1,lowP2,highP2)
dfR[dfR<0]<-0
#B
gorder<-c("Vaccine Coverage forces Growth Rate","Growth Rate forces Vaccine Coverage")


TH2<-ggplot(dfR,aes(x)) +
  geom_ribbon(aes(ymin = lowP1, ymax = highP1),fill=c4,alpha=0.4) +
  geom_ribbon(aes(ymin = lowP2, ymax = highP2),fill=c2,alpha=0.4)+
  geom_line(aes(x=x,y=meanP1,color="Vaccine Coverage forces Growth Rate"),size=0.8)+
  geom_line(aes(x=x,y=meanP2,color="Growth Rate forces Vaccine Coverage"),size=0.8)+
  theme_classic()+
  xlab(NULL) + 
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(name="Cross Mapping Skill(p)",   limits = c(0, 0.8),
                     breaks = seq(0,1,0.2))+
  theme(panel.background = element_rect( colour = "black", size = 0.8)) +
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(name="Library Size",   limits = c(0, 1180),
                     breaks = seq(0,1180,400),expand = c(0, 0))+
  scale_colour_manual(values=c(c4,c2),labels=gorder,breaks=gorder)+
  theme(axis.text.x = element_text(size = 14, color="black",vjust = 0.5, angle = 0))+
  theme(axis.text.y = element_text(size = 14, color="black",vjust = 0.5))+
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position = c(0.4,0.96),legend.direction = "vertical",legend.key.size=unit(0.4,'cm'))+
  labs(color="")+
  theme(plot.title = element_text(hjust = 0.5,size=12)) +
  theme(plot.margin=unit(c(1,0.3,1,0.3),'lines'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(text=element_text(size=14,  family="serif"))
TH2

DHS1<- read.csv("ccm_results2315/GR_xmap_TR_rho_statistics.csv",header=FALSE)
DHS2<- read.csv("ccm_results2315/TR_xmap_GR_rho_statistics.csv",header=FALSE)


meanP1=DHS1$V3
lowP1=DHS1$V4
highP1=DHS1$V2
meanP2=DHS2$V3
lowP2=DHS2$V4
highP2=DHS2$V2

x=DHS1$V1

dfR=data.frame(x,meanP1,meanP2,lowP1,highP1,lowP2,highP2)
dfR[dfR<0]<-0
#C
gorder<-c("International forces Travel Growth Rate","Growth Rate forces International Travel")


TH3<-ggplot(dfR,aes(x)) +
  geom_ribbon(aes(ymin = lowP1, ymax = highP1),fill=c5,alpha=0.4) +
  geom_ribbon(aes(ymin = lowP2, ymax = highP2),fill=c2,alpha=0.4)+
  geom_line(aes(x=x,y=meanP1,color="International forces Travel Growth Rate"),size=0.8)+
  geom_line(aes(x=x,y=meanP2,color="Growth Rate forces International Travel"),size=0.8)+
  theme_classic()+
  xlab(NULL) + 
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(name="Cross Mapping Skill(p)",   limits = c(0, 0.8),
                     breaks = seq(0,1,0.2))+
  theme(panel.background = element_rect( colour = "black", size = 0.8)) +
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(name="Library Size",   limits = c(0, 1180),
                     breaks = seq(0,1180,400),expand = c(0, 0))+
  scale_colour_manual(values=c(c5,c2),labels=gorder,breaks=gorder)+
  theme(axis.text.x = element_text(size = 14, color="black",vjust = 0.5, angle = 0))+
  theme(axis.text.y = element_text(size = 14, color="black",vjust = 0.5))+
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position = c(0.415,0.96),legend.direction = "vertical",legend.key.size=unit(0.4,'cm'))+
  labs(color="")+
  theme(plot.title = element_text(hjust = 0.5,size=12)) +
  theme(plot.margin=unit(c(1,0.3,1,0.3),'lines'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(text=element_text(size=14,  family="serif"))
TH3


c1="#f79c1d"
c3="#1426a4"
c4="#3a9679"
c5= (brewer.pal(8, "Blues"))[7]
c2= (brewer.pal(8, "Greys"))[7]


DHS1<- read.csv("ccm_results2315/SI_xmap_NPI_rho_statistics.csv",header=FALSE)
DHS2<- read.csv("ccm_results2315/NPI_xmap_SI_rho_statistics.csv",header=FALSE)


meanP1=DHS1$V3
lowP1=DHS1$V4
highP1=DHS1$V2
meanP2=DHS2$V3
lowP2=DHS2$V4
highP2=DHS2$V2

x=DHS1$V1

dfR=data.frame(x,meanP1,meanP2,lowP1,highP1,lowP2,highP2)
dfR[dfR<0]<-0
#D
gorder<-c("PHSM forces Shannon's Index","Shannon's Index forces PHSM")

TH4<-ggplot(dfR,aes(x)) +
  geom_ribbon(aes(ymin = lowP2, ymax = highP2),fill=c2,alpha=0.4)+
  geom_ribbon(aes(ymin = lowP1, ymax = highP1),fill=c1,alpha=0.4) +
  geom_line(aes(x=x,y=meanP2,color="Shannon's Index forces PHSM"),size=0.8)+
  geom_line(aes(x=x,y=meanP1,color="PHSM forces Shannon's Index"),size=0.8)+
  
  theme_classic()+
  xlab(NULL) + 
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(name="Cross Mapping Skill(p)",limits = c(0,0.8),
                     breaks = seq(0,1,0.2))+
  theme(panel.background = element_rect( colour = "black", size = 0.8)) +
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(name="Library Size",   limits = c(0, 1180),
                     breaks = seq(0,1180,400),expand = c(0, 0))+
  scale_colour_manual(values=c(c1,c2),labels=gorder,breaks=gorder)+
  theme(axis.text.x = element_text(size = 14, color="black",vjust = 0.5, angle = 0))+
  theme(axis.text.y = element_text(size = 14, color="black",vjust = 0.5))+
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position = c(0.335,0.96),legend.direction = "vertical",legend.key.size=unit(0.4,'cm'))+
  labs(color="")+
  theme(plot.title = element_text(hjust = 0.5,size=12)) +
  theme(plot.margin=unit(c(1,0.3,1,0.3),'lines'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(text=element_text(size=14,  family="serif"))
TH4


DHS1<- read.csv("ccm_results2315/SI_xmap_Vero_rho_statistics.csv",header=FALSE)
DHS2<- read.csv("ccm_results2315/Vero_xmap_SI_rho_statistics.csv",header=FALSE)



meanP1=DHS1$V3
lowP1=DHS1$V4
highP1=DHS1$V2
meanP2=DHS2$V3
lowP2=DHS2$V4
highP2=DHS2$V2


x=DHS1$V1

dfR=data.frame(x,meanP1,meanP2,lowP1,highP1,lowP2,highP2)
dfR[dfR<0]<-0
#E

gorder<-c("Vaccine Coverage forces Shannon's Index","Shannon's Index forces Vaccine Coverage")


TH5<-ggplot(dfR,aes(x=x)) +
  geom_ribbon(aes(ymin = lowP1, ymax = highP1),fill=c4,alpha=0.4) +
  geom_ribbon(aes(ymin = lowP2, ymax = highP2),fill=c2,alpha=0.4)+
  geom_line(aes(y=meanP1,color="Vaccine Coverage forces Shannon's Index"),size=0.8)+
  geom_line(aes(y=meanP2,color="Shannon's Index forces Vaccine Coverage"),size=0.8)+
  theme_classic()+
  xlab(NULL) +
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(name="Cross Mapping Skill(p)",limits = c(0,0.8),
                     breaks = seq(0,1,0.2))+
  theme(panel.background = element_rect( colour = "black", size = 0.8)) +
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(name="Library Size",   limits = c(0, 1180),
                     breaks = seq(0,1180,400),expand = c(0, 0))+
  scale_fill_manual(values=c(c4,c2))+
  scale_colour_manual(values=c(c4,c2),labels=gorder,breaks=gorder)+
  theme(axis.text.x = element_text(size = 14, color="black",vjust = 0.5, angle = 0))+
  theme(axis.text.y = element_text(size = 14, color="black",vjust = 0.5))+
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position = c(0.44,0.96),legend.direction = "vertical",legend.key.size=unit(0.4,'cm'))+
  labs(color="")+
  theme(plot.title = element_text(hjust = 0.5,size=12)) +
  theme(plot.margin=unit(c(1,0.3,1,0.3),'lines'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(text=element_text(size=14,  family="serif"))



TH5


DHS1<- read.csv("ccm_results2315/SI_xmap_TR_rho_statistics.csv",header=FALSE)
DHS2<- read.csv("ccm_results2315/TR_xmap_SI_rho_statistics.csv",header=FALSE)


meanP1=DHS1$V3
lowP1=DHS1$V4
highP1=DHS1$V2
meanP2=DHS2$V3
lowP2=DHS2$V4
highP2=DHS2$V2

x=DHS1$V1

dfR=data.frame(x,meanP1,meanP2,lowP1,highP1,lowP2,highP2)
dfR[dfR<0]<-0
#F
gorder<-c("International Travel forces Shannon's Index","Shannon's Index forces International Travel")


TH6<-ggplot(dfR,aes(x)) +
  geom_ribbon(aes(ymin = lowP1, ymax = highP1),fill=c5,alpha=0.4) +
  geom_ribbon(aes(ymin = lowP2, ymax = highP2),fill=c2,alpha=0.4)+
  geom_line(aes(x=x,y=meanP1,color="International Travel forces Shannon's Index"),size=0.8)+
  geom_line(aes(x=x,y=meanP2,color="Shannon's Index forces International Travel"),size=0.8)+
  theme_classic()+
  xlab(NULL) + 
  theme(axis.text=element_text(size=12))+
  scale_y_continuous(name="Cross Mapping Skill(p)",limits = c(0,0.8),
                     breaks = seq(0,1,0.2))+
  theme(panel.background = element_rect( colour = "black", size = 0.8)) +
  guides(fill=guide_legend(title=NULL))+
  scale_x_continuous(name="Library Size",   limits = c(0, 1180),
                     breaks = seq(0,1180,400),expand = c(0, 0))+
  scale_colour_manual(values=c(c5,c2),labels=gorder,breaks=gorder)+
  theme(axis.text.x = element_text(size = 14, color="black",vjust = 0.5, angle = 0))+
  theme(axis.text.y = element_text(size = 14, color="black",vjust = 0.5))+
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.position = c(0.45,0.96),legend.direction = "vertical",legend.key.size=unit(0.4,'cm'))+
  labs(color="")+
  theme(plot.title = element_text(hjust = 0.5,size=12)) +
  theme(plot.margin=unit(c(1,0.3,1,0.3),'lines'),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  theme(text=element_text(size=14,  family="serif"))
TH6



library(ggpubr) 
library(ggplot2) 

e<-ggarrange(TH1,TH2,TH3,TH4,TH5,TH6,ncol =3, nrow = 2,
             labels = c("A","B","C","D","E","F"), 
             font.label = list(size = 16, face = "bold",family="ArialN"))
ggsave(file="ccm_results2315/fig_ccm_sub_shannon_prop_05_mean.pdf", plot=e, width=12, height=8)