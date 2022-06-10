rm(list = ls())
library(ggplot2);library(ggsci);library(tidyr);library(car)
library(corrplot);library(ggcorrplot);library(RColorBrewer);library(tsModel)
library(tidyverse);library(naniar);library(mgcv);library(MASS);library(Cairo)
library(lubridate);library(dlnm);library(splines);library(showtext)
showtext_auto()
font_add("ArialN","arial.ttf")
font.families()


data <-read.csv("01_data/subsample_total_2000.csv")
#data <-read.csv("01_data/subsample_prop_05.csv")

data<- data[order(data$country,data$month,decreasing=F),]


data$travel2 <- log(data$travel_route)
y1 <- which(data$month=="2019-12") 
y2 <- which(data$month=="2020-01") 
y3 <- which(data$month=="2020-02") 
data <-data[-c(y1,y2,y3),]
# Set maximum lag
nlag = 3
rownames(data) <- data$record

colnames(data)
lag_npi <- tsModel::Lag(data$NPI, group = data$country, k = 0:nlag)
lag_lineage <- tsModel::Lag(data$sub_lineage, group = data$country, k = 0:nlag)
lag_travel <- tsModel::Lag(data$travel2, group = data$country, k = 0:nlag)
lag_gr <- tsModel::Lag(data$GR_NCM_log_diff, group = data$country, k = 0:nlag)
lag_FBW_vero <- tsModel::Lag(data$F_I_BWP, group = data$country, k = 0:nlag)




# Define cross basis matrix (combining nonlinear exposure and hysteresis function)
# Set the lag section
lagknot1 = equalknots(0:nlag, 1)
lagknot2 = equalknots(0:nlag, 2)
var <- lag_npi
npi_gr_cb <- crossbasis(var,
                        argvar = list(fun = "poly",degree=2),
                        arglag = list(fun = "ns", knots = lagknot2))


var <- lag_FBW_vero
FBW_gr_cb <- crossbasis(var,
                        argvar = list(fun = "ns",knots = equalknots(data$F_I_BWP, 2)),
                        arglag = list(fun = "ns", knots =lagknot1))




var <- lag_lineage
lineage_gr_cb <- crossbasis(var,
                            argvar = list(fun = "poly",degree=3),
                            arglag = list(fun = "ns", knots = lagknot1))


var <- lag_travel
travel_gr_cb <- crossbasis(var,
                           argvar = list(fun = "lin"),
                           arglag = list(fun = "ns", knots = lagknot2))








var <- lag_gr
gr_cb <- crossbasis(var,
                    argvar = list(fun = "poly",degree =3),
                    arglag = list(fun = "ns",knots =lagknot1))

var <- lag_npi
npi_si_cb <- crossbasis(var,
                        argvar = list(fun = "ns",Boundary.knots = range(data$NPI,na.rm=T)),
                        arglag = list(fun = "ns", knots = lagknot2))



var <- lag_FBW_vero
FBW_si_cb <- crossbasis(var,
                        argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns", knots =lagknot1))


var <- lag_travel
travel_si_cb <- crossbasis(var,
                           argvar = list(fun = "ns",Boundary.knots = range(data$travel2,na.rm=T)),
                           arglag = list(fun = "ns",knots = lagknot2))

# Specifies a unique column name for the inla() model
# Note: GLM (), GAM () or GLM.nb () models are not required
#colnames(npi_cb) = paste0("npi_cb.", colnames(npi_cb))
#colnames(vero_cb) = paste0("vero_cb.", colnames(vero_cb))
#colnames(T_vero_cb) = paste0("T_vero_cb.", colnames(T_vero_cb))
#colnames(F_vero_cb) = paste0("F_vero_cb.", colnames(F_vero_cb))
#colnames(shannon_cb) = paste0("shannon_cb.", colnames(npi_cb))



