rm(list = ls())
#Sys.setenv(LANGUAGE = "en")
#Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")
library(ggplot2);library(ggsci);library(RColorBrewer);library(tsModel);
library(tidyverse);library(tidyr);library(lubridate);library(mgcv);
library(MASS);library(splines);library(dlnm);library(showtext);
library(car);library(R330);library(boot)
showtext_auto()
font_add("ArialN","arial.ttf")
font_families()

data <-read.csv("01_data/subsample_prop_05_lineage_diversity.csv")# Load data
data$travel2 <- log(data$travel_route)
data$country <- as.factor(data$country)
data<- data[order(data$country,data$month,decreasing=F),]

y1 <- which(data$month=="2019-12") 
y2 <- which(data$month=="2020-01") 
y3 <- which(data$month=="2020-02") 
data <-data[-c(y1,y2,y3),]

rownames(data) <- data$record

# Set maximum lag
nlag = 3

# Creating lagged variables
# Define the monthly mean covariate lag term matrix
lag_npi <- tsModel::Lag(data$NPI, group = data$country, k = 0:nlag)
lag_lineage <- tsModel::Lag(data$sub_lineage, group = data$country, k = 0:nlag)
lag_travel <- tsModel::Lag(data$travel2, group = data$country, k = 0:nlag)
lag_gr <- tsModel::Lag(data$GR_NCM_log_diff, group = data$country, k = 0:nlag)
lag_FW_vero <- tsModel::Lag(data$F_I_BWP, group = data$country, k = 0:nlag)

# Define cross basis matrix (combining nonlinear exposure and hysteresis function)

# Set the lag section

var <- lag_npi
npi_cb <- crossbasis(var,
                    argvar = list(fun = "lin"),
                    arglag = list(fun = "ns",df=2))
summary(npi_cb)

###########################international travel######################################

var <- lag_travel
travel_cb <- crossbasis(var,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "ns", df=2))
summary(travel_cb)



###########################lineage diversity######################################
var <- lag_lineage
lineage_cb1 <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns",df=2))
summary(lineage_cb1)


var <- lag_lineage
lineage_cb2 <- crossbasis(var,
                        argvar = list(fun = "poly",degree=2),
                        arglag = list(fun = "ns", df=2))
summary(lineage_cb2)


var <- lag_lineage
lineage_cb3 <- crossbasis(var,
                        argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns", df=2))
summary(lineage_cb3)


var <- lag_lineage
lineage_cb4 <- crossbasis(var,
                        argvar = list(fun = "ns",
                                      knots = equalknots(data$sub_lineage, 1)),
                        arglag = list(fun = "ns", df=2))
summary(lineage_cb4)


var <- lag_lineage
lineage_cb5 <- crossbasis(var,
                        argvar = list(fun = "ns",
                                      knots = equalknots(data$sub_lineage, 2)),
                        arglag = list(fun = "ns", df=2))
summary(lineage_cb5)



############adjusted vaccine coverage consider immunity waning and booster#######################################

var <- lag_FW_vero
FW_cb1 <- crossbasis(var,
                      argvar = list(fun = "lin"),
                      arglag = list(fun = "ns", df=2))
summary(FW_cb1)


var <- lag_FW_vero
FW_cb2 <- crossbasis(var,
                      argvar = list(fun = "poly",degree=2),
                      arglag = list(fun = "ns", df=2))
summary(FW_cb2)


var <- lag_FW_vero
FW_cb3 <- crossbasis(var,
                      argvar = list(fun = "poly",degree=3),
                      arglag = list(fun = "ns", df=2))
summary(FW_cb3)


var <- lag_FW_vero
FW_cb4 <- crossbasis(var,
                      argvar = list(fun = "ns",
                                    knots = equalknots(data$F_I_BWP, 1)),
                      arglag = list(fun = "ns", df=2))
summary(FW_cb4)


var <- lag_FW_vero
FW_cb5 <- crossbasis(var,
                      argvar = list(fun = "ns",
                                    knots = equalknots(data$F_I_BWP, 2)),
                      arglag = list(fun = "ns", df=2))
summary(FW_cb5)


############growth rate of cases#######################################

var <- lag_gr
gr_cb1 <- crossbasis(var,
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns", df=2))
summary(gr_cb1)


var <- lag_gr
gr_cb2 <- crossbasis(var,
                     argvar = list(fun = "poly",degree=2),
                     arglag = list(fun = "ns", df=2))
summary(gr_cb2)


var <- lag_gr
gr_cb3 <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                     arglag = list(fun = "ns", df=2))
summary(gr_cb3)


var <- lag_gr
gr_cb4 <- crossbasis(var,
                     argvar = list(fun = "ns",
                                   knots = equalknots(data$GR_NCM_log_diff, 1)),
                     arglag = list(fun = "ns", df=2))
summary(gr_cb4)


var <- lag_gr
gr_cb5 <- crossbasis(var,
                     argvar = list(fun = "ns",
                                   knots = equalknots(data$GR_NCM_log_diff, 2)),
                     arglag = list(fun = "ns", df=2))
summary(gr_cb5)


# assign unique column names to cross-basis matrix for gam() model
# note: not necessary for glm(), gam() or glm.nb() models
colnames(npi_cb) = paste0("npi_cb.", colnames(npi_cb))

colnames(travel_cb) = paste0("travel_cb.", colnames(travel_cb))

colnames(lineage_cb1) = paste0("lineage_cb1.", colnames(lineage_cb1))
colnames(lineage_cb2) = paste0("lineage_cb2.", colnames(lineage_cb2))
colnames(lineage_cb3) = paste0("lineage_cb3.", colnames(lineage_cb3))
colnames(lineage_cb4) = paste0("lineage_cb4.", colnames(lineage_cb4))
colnames(lineage_cb5) = paste0("lineage_cb5.", colnames(lineage_cb5))

colnames(FW_cb1) = paste0("FW_cb1.", colnames(FW_cb1))
colnames(FW_cb2) = paste0("FW_cb2.", colnames(FW_cb2))
colnames(FW_cb3) = paste0("FW_cb3.", colnames(FW_cb3))
colnames(FW_cb4) = paste0("FW_cb4.", colnames(FW_cb4))
colnames(FW_cb5) = paste0("FW_cb5.", colnames(FW_cb5))

colnames(gr_cb1) = paste0("gr_cb1.", colnames(gr_cb1))
colnames(gr_cb2) = paste0("gr_cb2.", colnames(gr_cb2))
colnames(gr_cb3) = paste0("gr_cb3.", colnames(gr_cb3))
colnames(gr_cb4) = paste0("gr_cb4.", colnames(gr_cb4))
colnames(gr_cb5) = paste0("gr_cb5.", colnames(gr_cb5))




            
