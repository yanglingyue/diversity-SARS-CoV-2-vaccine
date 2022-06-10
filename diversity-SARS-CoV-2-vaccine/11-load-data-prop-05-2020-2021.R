rm(list = ls())
library(ggplot2);library(ggsci);library(tidyr);library(car)
library(corrplot);library(ggcorrplot);library(RColorBrewer);library(tsModel)
library(tidyverse);library(naniar);library(mgcv);library(MASS);library(Cairo)
library(lubridate);library(dlnm);library(splines);library(showtext)
showtext_auto()
font_add("ArialN","arial.ttf")
font.families()


data <-read.csv("01_data/subsample_prop_05.csv")
data<- data[order(data$country,data$month,decreasing=F),]
data$date <- parse_date_time(data$month, "ym")
data$date <-as.Date(as.POSIXct(data$date,tz="Hongkong"))
data$month1 <- month(data$date)
data$year <- year(data$date)
#table <- as.data.frame(table(data$country))
#colnames(table) <- c("country","month")
data$travel2 <- log(data$travel_route)
y1 <- which(data$month=="2019-12") 
y2 <- which(data$month=="2020-01") 
y3 <- which(data$month=="2020-02") 
data <-data[-c(y1,y2,y3),]

y4 <- which(data$month=="2020-09") 
y5 <- which(data$month=="2020-10")
y6 <- which(data$month=="2020-11") 
y7 <- which(data$month=="2020-12")
#y8 <- which(data$month=="2020-01") 
#y3 <- which(data$month=="2020-02") 
data_delta <- data[c(y4,y5,y6,y7),]

data_2021_delta = subset(data,data$year== "2021")

data_2021 <- rbind(data_delta,data_2021_delta)

# Set maximum lag
nlag = 3
rownames(data_2021) <- data_2021$record

colnames(data_2021)
lag_npi <- tsModel::Lag(data_2021$NPI, group = data_2021$country, k = 0:nlag)
lag_lineage <- tsModel::Lag(data_2021$sub_lineage, group = data_2021$country, k = 0:nlag)
lag_travel <- tsModel::Lag(data_2021$travel2, group = data_2021$country, k = 0:nlag)
lag_gr <- tsModel::Lag(data_2021$GR_NCM_log_diff, group = data_2021$country, k = 0:nlag)
lag_FBW_vero <- tsModel::Lag(data_2021$F_I_BWP, group = data_2021$country, k = 0:nlag)




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



