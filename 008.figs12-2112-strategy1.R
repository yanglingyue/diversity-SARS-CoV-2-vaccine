setwd("D:/code_upload/vaccine-SARS-CoV-2-diversity")
rm(list = ls())
library(ggplot2);library(ggsci);library(RColorBrewer);library(tsModel);
library(tidyverse);library(tidyr);library(lubridate);library(mgcv);
library(MASS);library(splines);library(dlnm);library(showtext);
library(car);library(R330);library(boot);library(dplyr)
showtext_auto()
font_add("ArialN","arial.ttf")
font_families()
dat1 <-read.csv("01_data/subsample_diversity_vaccine_data.csv")# Load data

#Rename the Shannon's index of lineage diversity data, calculated based on subsampling strategy 1, as "sub_lineage"
colnames(dat1)[13] <- "sub_lineage" 

#Rename the richness of lineage diversity data, calculated based on subsampling strategy 1, as "sub_richness"
colnames(dat1)[14] <- "sub_richness" 

#Rename the evenness of lineage diversity data, calculated based on subsampling strategy 1, as "sub_evenness"
colnames(dat1)[15] <- "sub_evenness" 

#Rename the Simpson' indec of lineage diversity data, calculated based on subsampling strategy 1, as "sub_evenness"
colnames(dat1)[19] <- "simpson_lineage" 

dat1 <- dat1 %>% 
  group_by(month) %>% 
  mutate(country_index = row_number())

dat1$country_index <- as.factor(dat1$country_index)
dat1$travel2 <- log(dat1$travel_route)
data <- dat1
y1 <- which(data$month=="2019-12") 
y2 <- which(data$month=="2020-01") 
y3 <- which(data$month=="2020-02") 
y4 <- which(data$month=="2022-09"|data$month=="2022-08"|data$month=="2022-07"|data$month=="2022-06"
            |data$month=="2022-05"|data$month=="2022-04"|data$month=="2022-03"|data$month=="2022-02"
            |data$month=="2022-01") 

data <-data[-c(y1,y2,y3,y4),]

dir.create("02_fig/figs12_2112_strategy1")
outpath <-  "02_fig/figs12_2112_strategy1/figs12-2112"


source('002.load-data-cb-function.R')

# set data for models of growth rate of cases and lineage diversity
Y1  <- data$GR_NCM_log_diff # response variable (growth rate of cases)
Y2  <- data$sub_lineage # response variable (Shannon's index of lineage diversity based on the subsampling strategy 1)
Y3  <- data$sub_richness # response variable (the richness of lineage diversity based on the subsampling strategy 1)
Y4  <- data$sub_evenness # response variable (the evenness of lineage diversity based on the subsampling strategy 1)
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
M <- data$month # The month index of the dataset 

df <- data.frame(Y1, Y2, Y3, Y4, P, V, N, C, M)
A <- round(max(df$V,na.rm = T),2)*10



model_gr <-  gam(Y1 ~ npi_gr_cb+FBW_gr_cb+lineage_gr_cb+travel_gr_cb+N+s(C, bs = "re")+P:V
                 ,method = "REML",data=df)
summary(model_gr)
AIC(model_gr)
cp_npi_gr <- crosspred(npi_gr_cb, model_gr, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr <- crosspred(travel_gr_cb, model_gr, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
cp_lineage_gr <- crosspred(lineage_gr_cb, model_gr,cen=0,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_FI_gr <- crosspred(FBW_gr_cb, model_gr, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


model_si <-gam(Y2 ~ Lag(df$Y2,1)+travel_si_cb+gr_cb+npi_si_cb+FBW_si_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_si)
AIC(model_si)

cp_npi_si <- crosspred(npi_si_cb, model_si,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
#cp_travel_si <- crosspred(travel_si_cb, model_si, cen=0,at=0:15000/1000,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si <- crosspred(travel_si_cb, model_si, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_gr_si <- crosspred(gr_cb, model_si, cen=0,by=0.1,at=c(-35:0/10,0:95/10),bylag=0.1,cumul=TRUE) #slices
cp_FI_si <- crosspred(FBW_si_cb, model_si, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


model_gr_rich <- gam(Y1 ~ FBW_gr_cb+richness_gr_cb+npi_gr_cb+travel_gr_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_gr_rich)
cp_npi_gr_rich <- crosspred(npi_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr_rich <- crosspred(travel_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
cp_FI_gr_rich <- crosspred(FBW_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


model_si_rich <- gam(Y3 ~ Lag(df$Y3,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_si_rich)
cp_npi_si_rich <- crosspred(npi_si_cb, model_si_rich,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si_rich <- crosspred(travel_si_cb, model_si_rich, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_FI_si_rich <- crosspred(FBW_si_cb, model_si_rich, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


model_gr_even <- gam(Y1 ~ FBW_gr_cb+evenness_gr_cb+npi_gr_cb+travel_gr_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_gr_even)
cp_npi_gr_even <- crosspred(npi_gr_cb, model_gr_even, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr_even <- crosspred(travel_gr_cb, model_gr_even, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
cp_FI_gr_even <- crosspred(FBW_gr_cb, model_gr_even, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


model_si_even <- gam(Y4 ~ Lag(df$Y4,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+N+s(C, bs = "re")+P:V,method = "REML",data=data)
summary(model_si_even)
cp_npi_si_even <- crosspred(npi_si_cb, model_si_even,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si_even <- crosspred(travel_si_cb, model_si_even, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_FI_si_even <- crosspred(FBW_si_cb, model_si_even, cen=0,by=0.1,at=0:A/10,bylag=0.1,cumul=TRUE) #slices


#si_m <- mean(data$sub_lineage,na.rm=T)
#it_m <- mean(data$travel2,na.rm=T)

#r_l <- cp_travel_si$cumlow["1.593",]/si_m 
#r_f <- cp_travel_si$cumfit["1.593",]/si_m 
#r_h <- cp_travel_si$cumhigh["1.593",]/si_m 

#0.09753737/si_m 

pdf(file=paste0(outpath,"_a-f.pdf"), width = 21, height = 21/3*2)

par(mfrow=c(2,3),mar = c(8, 8, 5, 2),oma = c(0, 2, 0, 0))

# get exposures values
vars <- cp_travel_gr$predvar

# obtain growth rate fit and upper and lower confidence limits for all exposure variables
gr <- cp_travel_gr$cumfit
gr.lci <- cp_travel_gr$cumlow
gr.uci <- cp_travel_gr$cumhigh

# set growth rate range 
g1 <- min(range(gr, gr.lci, gr.uci))
g2 <- max(range(gr, gr.lci, gr.uci))

si <- cp_travel_si$cumfit
si.lci <- cp_travel_si$cumlow
si.uci <- cp_travel_si$cumhigh

# set shannon index range 
s1 <- min(range(si, si.lci, si.uci))
s2 <- max(range(si, si.lci, si.uci))

# get selected lag variable positions
lag_gr <- 4
lag_si <- 4

# define colours
col1 <- "#00468BFF"
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/3, max = 255))

col2 <- "#ED0000FF"
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))


# define x values (var, by var)
varby <- seq(0, 15, 0.1)

# growth rate
#par(mar = c(5, 6.5, 3, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "International travel", ylab = "", main = "", 
     ylim = range(-2, 2), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1,at=0:5*3,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-2:0*1,2:0*1),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
xx <- c(varby, rev(varby))
yy <- c(gr.lci[,lag_gr], rev(gr.uci[,lag_gr]))
polygon(xx, yy, col = tcol1, border = tcol1)
abline(h = 0, lty = 3)
# diversity
lines(varby, si[,lag_si], col = col2, lwd = 3)
xx <- c(varby, rev(varby))
yy <- c(si.lci[,lag_si],rev(si.uci[,lag_si]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 0, lty = 3)
mtext("Effect on growth rate of cases and \nShannon's index of lineage diversity",side = 2,at=-0, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Shannon's index of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#      y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 2.32, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


# get exposures values
vars <- cp_travel_gr_even$predvar

# obtain growth rate fit and upper and lower confidence limits for all exposure variables
gr <- cp_travel_gr_even$cumfit
gr.lci <- cp_travel_gr_even$cumlow
gr.uci <- cp_travel_gr_even$cumhigh

# set growth rate range 
g1 <- min(range(gr, gr.lci, gr.uci))
g2 <- max(range(gr, gr.lci, gr.uci))

si <- cp_travel_si_even$cumfit
si.lci <- cp_travel_si_even$cumlow
si.uci <- cp_travel_si_even$cumhigh

# set shannon index range 
s1 <- min(range(si, si.lci, si.uci))
s2 <- max(range(si, si.lci, si.uci))

# get selected lag variable positions
lag_gr <- 4
lag_si <- 4
#mx2 <- 1

# define colours
col1 <- "#00468BFF"
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/3, max = 255))

col2 <-  "seagreen" 
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))


# define x values (var, by var)
varby <- seq(0, 15, 0.1)

# growth rate
#par(mar = c(5, 6.5, 3, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "International travel", ylab = "", main = "", 
     ylim = range(-2, 2), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1,at=0:5*3,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-2:0*1,2:0*1),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
xx <- c(varby, rev(varby))
yy <- c(gr.lci[,lag_gr], rev(gr.uci[,lag_gr]))
polygon(xx, yy, col = tcol1, border = tcol1)
abline(h = 0, lty = 3)
# diversity
lines(varby, si[,lag_si], col = col2, lwd = 3)
xx <- c(varby, rev(varby))
yy <- c(si.lci[,lag_si],rev(si.uci[,lag_si]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 0, lty = 3)
mtext("Effect on growth rate of cases and \nevenness of lineage diversity",side = 2,at=0, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Evenness of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 2.32, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 5.5)

# get exposures values
vars <- cp_travel_gr_rich$predvar

# obtain growth rate fit and upper and lower confidence limits for all exposure variables
gr <- cp_travel_gr_rich$cumfit
gr.lci <- cp_travel_gr_rich$cumlow
gr.uci <- cp_travel_gr_rich$cumhigh

# set growth rate range 
g1 <- min(range(gr, gr.lci, gr.uci))
g2 <- max(range(gr, gr.lci, gr.uci))

si <- cp_travel_si_rich$cumfit
si.lci <- cp_travel_si_rich$cumlow
si.uci <- cp_travel_si_rich$cumhigh

# set shannon index range 
s1 <- min(range(si, si.lci, si.uci))
s2 <- max(range(si, si.lci, si.uci))

# get selected lag variable positions
lag_gr <- 4
lag_si <- 4
#mx2 <- 1

# define colours
col1 <- "#00468BFF"
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/3, max = 255))

col2 <- "#E18727B2"

tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))


# define x values (var, by var)
varby <- seq(0, 15, 0.1)

# growth rate
#par(mar = c(5, 6.5, 3, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "International travel", ylab = "", main = "", 
     ylim = range(0,20), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1,at=0:5*3,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(4:0*5),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
xx <- c(varby, rev(varby))
yy <- c(gr.lci[,lag_gr], rev(gr.uci[,lag_gr]))
polygon(xx, yy, col = tcol1, border = tcol1)
abline(h = 0, lty = 3)
# diversity
lines(varby, si[,lag_si], col = col2, lwd = 3)
xx <- c(varby, rev(varby))
yy <- c(si.lci[,lag_si],rev(si.uci[,lag_si]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 0, lty = 3)
mtext("Effect on growth rate of cases and \nrichness of lineage diversity",side = 2,at=10, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Richness of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 22, text = "C", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


# get exposures values
vars <- cp_npi_gr$predvar

# obtain growth rate fit and upper and lower confidence limits for all exposure variables
gr <- cp_npi_gr$cumfit
gr.lci <- cp_npi_gr$cumlow
gr.uci <- cp_npi_gr$cumhigh

# set growth rate range 
g1 <- min(range(gr, gr.lci, gr.uci))
g2 <- max(range(gr, gr.lci, gr.uci))

si <- cp_npi_si$cumfit
si.lci <- cp_npi_si$cumlow
si.uci <- cp_npi_si$cumhigh

# set shannon index range 
s1 <- min(range(si, si.lci, si.uci))
s2 <- max(range(si, si.lci, si.uci))


# get selected lag variable positions
lag_gr <- 3
lag_si <- 3
#mx2 <- 1

# define colours
col1 <- "#00468BFF"
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/3, max = 255))

col2 <-  "#ED0000FF"
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))


# define x values (var, by var)
varby <- seq(0, 100, 0.1)

# growth rate
#par(mar = c(5, 6.5, 5, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "Stringency index", ylab = "", main = "", 
     ylim = range(-3, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*1),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
xx <- c(varby, rev(varby))
yy <- c(gr.lci[,lag_gr], rev(gr.uci[,lag_gr]))
polygon(xx, yy, col = tcol1, border = tcol1)
abline(h = 0, lty = 3)
# diversity
lines(varby, si[,lag_si], col = col2, lwd = 3)
xx <- c(varby, rev(varby))
yy <- c(si.lci[,lag_si],rev(si.uci[,lag_si]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 0, lty = 3)
mtext("Effect on growth rate of cases and\nShannon's index of lineage diversity",side = 2,at=-1.5, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Shannon's index of lineage diversity"),
#                                   paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#      y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 0.42, text = "D", las = 2, cex = 2.5, family = "ArialN",line = 5.5)

# get exposures values
vars <- cp_npi_gr_even$predvar

# obtain growth rate fit and upper and lower confidence limits for all exposure variables
gr <- cp_npi_gr_even$cumfit
gr.lci <- cp_npi_gr_even$cumlow
gr.uci <- cp_npi_gr_even$cumhigh

# set growth rate range 
g1 <- min(range(gr, gr.lci, gr.uci))
g2 <- max(range(gr, gr.lci, gr.uci))

si <- cp_npi_si_even$cumfit
si.lci <- cp_npi_si_even$cumlow
si.uci <- cp_npi_si_even$cumhigh

# set shannon index range 
s1 <- min(range(si, si.lci, si.uci))
s2 <- max(range(si, si.lci, si.uci))


# get selected lag variable positions
lag_gr <- 3
lag_si <- 3
#mx2 <- 1

# define colours
col1 <- "#00468BFF"
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/3, max = 255))

col2 <- "seagreen" 
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))


# define x values (var, by var)
varby <- seq(0, 100, 0.1)

# growth rate
#par(mar = c(5, 6.5, 5, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "Stringency index", ylab = "", main = "", 
     ylim = range(-3, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*1),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
xx <- c(varby, rev(varby))
yy <- c(gr.lci[,lag_gr], rev(gr.uci[,lag_gr]))
polygon(xx, yy, col = tcol1, border = tcol1)
abline(h = 0, lty = 3)
# diversity
lines(varby, si[,lag_si], col = col2, lwd = 3)
xx <- c(varby, rev(varby))
yy <- c(si.lci[,lag_si],rev(si.uci[,lag_si]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 0, lty = 3)
mtext("Effect on growth rate of cases and\nevenness of lineage diversity",side = 2,at=-1.5, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Evenness of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#      y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 0.42, text = "E", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


# get exposures values
vars <- cp_npi_gr_rich$predvar

# obtain growth rate fit and upper and lower confidence limits for all exposure variables
gr <- cp_npi_gr_rich$cumfit
gr.lci <- cp_npi_gr_rich$cumlow
gr.uci <- cp_npi_gr_rich$cumhigh

# set growth rate range 
g1 <- min(range(gr, gr.lci, gr.uci))
g2 <- max(range(gr, gr.lci, gr.uci))

si <- cp_npi_si_rich$cumfit
si.lci <- cp_npi_si_rich$cumlow
si.uci <- cp_npi_si_rich$cumhigh

# set shannon index range 
s1 <- min(range(si, si.lci, si.uci))
s2 <- max(range(si, si.lci, si.uci))


# get selected lag variable positions
lag_gr <- 3
lag_si <- 3
#mx2 <- 1

# define colours
col1 <- "#00468BFF"
tcol1 <- do.call(rgb, c(as.list(col2rgb(col1)), alpha = 255/3, max = 255))

col2 <-"#E18727B2"
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))


# define x values (var, by var)
varby <- seq(0, 100, 0.1)

# growth rate
#par(mar = c(5, 6.5, 5, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "Stringency index", ylab = "", main = "", 
     ylim = range(-18, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*6),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
xx <- c(varby, rev(varby))
yy <- c(gr.lci[,lag_gr], rev(gr.uci[,lag_gr]))
polygon(xx, yy, col = tcol1, border = tcol1)
abline(h = 0, lty = 3)
# diversity
lines(varby, si[,lag_si], col = col2, lwd = 3)
xx <- c(varby, rev(varby))
yy <- c(si.lci[,lag_si],rev(si.uci[,lag_si]))
polygon(xx, yy, col = tcol2, border = tcol2)
abline(h = 0, lty = 3)
mtext("Effect on growth rate of cases and \nrichness of lineage diversity",side = 2,at=-9, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Richness of lineage diversity"),
#                                    paste0("Growth rate of cases")),
#      col = c(col2, col1), 
#    lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 0.42, text = "F", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


dev.off()





library(data.table);library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(ggpubr)
mytheme <- theme_bw()+
  theme(
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.y=element_blank(),
    panel.grid.minor.y=element_blank(),
    # plot.title=element_text(size=12, vjust = 0.5, hjust = 0.5),
    legend.position = c("none"),
    # legend.key.size = unit(0.25,"cm"),
    #legend.key.width = unit(0.2,"cm"),
    #legend.key.height  = unit(0.5,"cm"),
    legend.background = element_rect(fill=NA, size=0,color=NA),
    legend.text=element_text(size=20),
    legend.title=element_text(face="bold",size=20),
    axis.line.y=element_line(linetype=1,color='black',size=0.5),
    axis.line.x=element_line(linetype=1,color='black',size=0.5),
    axis.ticks = element_line(linetype=2,color='black'),
    panel.grid=element_line(linetype=2,color='grey'),
    plot.title = element_text(hjust = 0.5,size=22, face = "bold"),
    axis.title.y.left = element_text(size = 20,color="black",vjust=2,face = "bold"),
    axis.title.y.right =element_text(size = 20,color="black",vjust=0,angle=90),
    axis.title.x = element_text(size = 20, color="black",vjust = -0.5,face = "bold"),
    axis.text.y.left = element_text(size = 20,color="black",vjust = 0,angle = 0),
    axis.text.y.right = element_text(size = 20,color="black",vjust = 0,angle = 0),
    axis.text.x = element_text(size = 20, color="black",vjust = 0.5, angle = 0),
    axis.ticks.length=unit(0.15,"cm"),
    axis.ticks.y=element_line(color="black",size=.5),
    axis.ticks.x=element_line(color="black",size=.5),
    plot.margin=unit(c(2.5,2.5,2.5,2.5),'lines'),
    panel.border = element_rect(colour = "white", fill=NA, size=1)
  )



vero_low <- cp_FI_si$cumlow["75","lag1"]
vero_fit <- cp_FI_si$cumfit["75","lag1"]
vero_high <- cp_FI_si$cumhigh["75","lag1"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$Type <- "SI"
result$coverage1 <- "75"
result$coverage <- "High (75)"
result_vero_75 <- result 

vero_low <- cp_FI_si$cumlow["50","lag1"]
vero_fit <- cp_FI_si$cumfit["50","lag1"]
vero_high <- cp_FI_si$cumhigh["50","lag1"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$Type <-  "SI"
result$coverage1 <- "50"
result$coverage <- "Medium (50)"
result_vero_50 <- result

vero_low <- cp_FI_si$cumlow["25","lag1"]
vero_fit <- cp_FI_si$cumfit["25","lag1"]
vero_high <- cp_FI_si$cumhigh["25","lag1"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$Type <- "SI"
result$coverage1 <- "25"
result$coverage <- "Low (25)"
result_vero_25 <- result


result_vero <- rbind(result_vero_25,result_vero_50,result_vero_75)
result_vero$Lag <- "lag1"
result_vero_shannon <- result_vero


result_vero_shannon$coverage <-factor(result_vero_shannon$coverage,levels= c('Low (25)','Medium (50)','High (75)',
                                                                             ordered=TRUE))

p1  <- ggplot(result_vero_shannon,aes(x=coverage, y=vero_fit,group=coverage))+ 
  geom_point(aes(x=coverage,y=vero_fit,color=coverage),position=position_dodge(0.5),size=8,alpha=1)+
  geom_errorbar(aes(x=coverage, y=vero_fit,ymin=vero_low,ymax=vero_high,color=coverage),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6)+
  #facet_wrap(~Type,scales = "free",nrow = 1)+
  labs(x="Adjusted vaccine coverage (%)", y= "Effect on Shannon's index of lineage diversty" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+ 
  scale_y_continuous(limits = c(-1.5,0.5),breaks = c(-1.5,-1,-0.5,0,0.5))+
  mytheme+
  theme(axis.line.x = element_line(color="black", size = 0.8),axis.line.y = element_line(color="black", size = 0.8))

p1

ggsave(file=paste0(outpath,"-h",".pdf"),
       p1,width = 7.2,height =7.2)

vero_low <- cp_FI_gr$cumlow["75","lag1"]
vero_fit <- cp_FI_gr$cumfit["75","lag1"]
vero_high <- cp_FI_gr$cumhigh["75","lag1"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$Type <- "GR_SI"
result$coverage1 <- "75"
result$coverage <- "High (75)"
result_vero_75 <- result 

vero_low <- cp_FI_gr$cumlow["50","lag1"]
vero_fit <- cp_FI_gr$cumfit["50","lag1"]
vero_high <- cp_FI_gr$cumhigh["50","lag1"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$Type <-  "GR_SI"
result$coverage1 <- "50"
result$coverage <- "Medium (50)"
result_vero_50 <- result

vero_low <- cp_FI_gr$cumlow["25","lag1"]
vero_fit <- cp_FI_gr$cumfit["25","lag1"]
vero_high <- cp_FI_gr$cumhigh["25","lag1"]
result <- cbind(vero_low,vero_fit,vero_high)
result <- as.data.frame(result)
result$Lag <- rownames(result)
result$Type <- "GR_SI"
result$coverage1 <- "25"
result$coverage <- "Low (25)"
result_vero_25 <- result

result_vero <- rbind(result_vero_25,result_vero_50,result_vero_75)
result_vero$Lag <- "lag1"
result_vero_shannon_gr <- result_vero
result_vero_shannon_gr$coverage <-factor(result_vero_shannon_gr$coverage,levels= c('Low (25)','Medium (50)','High (75)',
                                                                                   ordered=TRUE))

p2  <- ggplot(result_vero_shannon_gr,aes(x=coverage, y=vero_fit,group=coverage))+ 
  geom_point(aes(x=coverage,y=vero_fit,color=coverage),position=position_dodge(0.5),size=8,alpha=1)+
  geom_errorbar(aes(x=coverage, y=vero_fit,ymin=vero_low,ymax=vero_high,color=coverage),
                width=0.2,stat="identity",position = position_dodge(0.5),lwd=0.6)+
  #facet_wrap(~Type,scales = "free",nrow = 1)+
  labs(x="Adjusted vaccine coverage (%)", y= "Effect on growth rate of cases" )+  
  geom_hline(aes(yintercept=0),lty='dashed',colour ="black",size=0.5)+
  scale_color_manual(name = 'Adjusted vaccine coverage',values = c("#92C5DE","#4393C3","#2166AC"))+ 
  scale_y_continuous(limits = c(-4,0),breaks = c(-4,-3,-2,-1,0))+
  mytheme+
  theme(axis.line.x = element_line(color="black", size = 0.8),axis.line.y = element_line(color="black", size = 0.8))



p2

ggsave(file=paste0(outpath,"-g",".pdf"),
       p2,width = 7.2,height =7.2)




p<-ggarrange(p1, p2,ncol =2, nrow = 1,widths = c(1,1),heights = c(1,1),
             labels = c("A","B"),  vjust=1.3,align = "hv",
             #common.legend=TRUE,legend = "bottom",
             font.label = list(size = 24))
p


ggsave(file=paste0(outpath,"_growth_rate_shannon",".pdf"),
       p,width = 16,height =8)

