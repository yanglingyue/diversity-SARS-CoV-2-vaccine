rm(list = ls())
#Sys.setenv(LANGUAGE = "en")
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")
library(ggplot2);library(ggsci);library(RColorBrewer);library(tsModel);
library(tidyverse);library(tidyr);library(lubridate);
library(mgcv);library(MASS);library(splines);library(dlnm);library(showtext);
library(car);library(R330);
library(boot)
showtext_auto()
font_add("ArialN","arial.ttf")
font_families()

data <-read.csv("01_data/subsample_prop_05.csv")# Load data
data<- data[order(data$country,data$month,decreasing=F),]
#data$date <- parse_date_time(data$month, "ym")
#data$year <- year(data$date)
data$travel2 <- log(data$travel_route)

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
lagknot1 = equalknots(0:nlag, 1)
lagknot2 = equalknots(0:nlag, 2)

var <- lag_npi
npi_cb1 <- crossbasis(var,
                    argvar = list(fun = "lin"),
                    arglag = list(fun = "ns",knots = lagknot2))
summary(npi_cb1)

var <- lag_npi
npi_cb2 <- crossbasis(var,
                        argvar = list(fun = "poly",degree=2),
                        arglag = list(fun = "ns", knots = lagknot2))
summary(npi_cb2)

var <- lag_npi
npi_cb3 <- crossbasis(var,
                        argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns", knots = lagknot2))
summary(npi_cb3)

var <- lag_npi
npi_cb4 <- crossbasis(var,
                        argvar = list(fun = "ns",knots = equalknots(data$NPI, 1)),
                        arglag = list(fun = "ns", knots = lagknot2))
summary(npi_cb4)

var <- lag_npi
npi_cb5 <- crossbasis(var,
                        argvar = list(fun = "ns",knots = equalknots(data$NPI, 2)),
                        arglag = list(fun = "ns", knots = lagknot2))
summary(npi_cb5)



###########################S_I######################################


var <- lag_lineage
lineage_cb1 <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns",knots = lagknot1))
summary(lineage_cb1)



var <- lag_lineage
lineage_cb2 <- crossbasis(var,
                        argvar = list(fun = "poly",degree=2),
                        arglag = list(fun = "ns", knots = lagknot1))
summary(lineage_cb2)

var <- lag_lineage
lineage_cb3 <- crossbasis(var,
                        argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns", knots = lagknot1))
summary(lineage_cb3)

var <- lag_lineage
lineage_cb4 <- crossbasis(var,
                        argvar = list(fun = "ns",
                                      knots = equalknots(data$sub_lineage, 1)),
                        arglag = list(fun = "ns", knots = lagknot1))
summary(lineage_cb4)

var <- lag_lineage
lineage_cb5 <- crossbasis(var,
                        argvar = list(fun = "ns",
                                      knots = equalknots(data$sub_lineage, 2)),
                        arglag = list(fun = "ns", knots = lagknot1))
summary(lineage_cb5)

###########################I_T######################################


var <- lag_travel
travel_cb1 <- crossbasis(var,
                          argvar = list(fun = "lin"),
                          arglag = list(fun = "ns", knots = lagknot2))
summary(travel_cb1)


var <- lag_travel
travel_cb2 <- crossbasis(var,
                          argvar = list(fun = "poly",degree=2),
                          arglag = list(fun = "ns", knots = lagknot2))
summary(travel_cb2)

var <- lag_travel
travel_cb3 <- crossbasis(var,
                          argvar = list(fun = "poly",degree=3),
                          arglag = list(fun = "ns", knots = lagknot2))
summary(travel_cb3)

var <- lag_travel
travel_cb4 <- crossbasis(var,
                          argvar = list(fun = "ns",
                                        knots = equalknots(data$travel2, 1)),
                          arglag = list(fun = "ns", knots = lagknot2))
summary(travel_cb4)

var <- lag_travel
travel_cb5 <- crossbasis(var,
                          argvar = list(fun = "ns",
                                        knots = equalknots(data$travel2, 2)),
                          arglag = list(fun = "ns", knots = lagknot2))
summary(travel_cb5)


############vero_immunity_waning#######################################

var <- lag_FW_vero
FW_cb1 <- crossbasis(var,
                      argvar = list(fun = "lin"),
                      arglag = list(fun = "ns", knots = lagknot1))
summary(FW_cb1)

var <- lag_FW_vero
FW_cb2 <- crossbasis(var,
                      argvar = list(fun = "poly",degree=2),
                      arglag = list(fun = "ns", knots = lagknot1))
summary(FW_cb2)


var <- lag_FW_vero
FW_cb3 <- crossbasis(var,
                      argvar = list(fun = "poly",degree=3),
                      arglag = list(fun = "ns", knots =lagknot1))
summary(FW_cb3)

var <- lag_FW_vero
FW_cb4 <- crossbasis(var,
                      argvar = list(fun = "ns",
                                    knots = equalknots(data$F_I_BWP, 1)),
                      arglag = list(fun = "ns", knots = lagknot1))
summary(FW_cb4)

var <- lag_FW_vero
FW_cb5 <- crossbasis(var,
                      argvar = list(fun = "ns",
                                    knots = equalknots(data$F_I_BWP, 2)),
                      arglag = list(fun = "ns", knots=lagknot1))
summary(FW_cb5)



############growth_rate#######################################

var <- lag_gr
gr_cb1 <- crossbasis(var,
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "ns", knots=lagknot1))
summary(gr_cb1)

var <- lag_gr
gr_cb2 <- crossbasis(var,
                     argvar = list(fun = "poly",degree=2),
                     arglag = list(fun = "ns", knots=lagknot1))
summary(gr_cb2)

var <- lag_gr
gr_cb3 <- crossbasis(var,
                     argvar = list(fun = "poly",degree=3),
                     arglag = list(fun = "ns", knots=lagknot1))
summary(gr_cb3)

var <- lag_gr
gr_cb4 <- crossbasis(var,
                     argvar = list(fun = "ns",
                                   knots = equalknots(data$GR_NCM_log_diff, 1)),
                     arglag = list(fun = "ns", knots=lagknot1))
summary(gr_cb4)

var <- lag_gr
gr_cb5 <- crossbasis(var,
                     argvar = list(fun = "ns",
                                   knots = equalknots(data$GR_NCM_log_diff, 2)),
                     arglag = list(fun = "ns", knots=lagknot1))
summary(gr_cb5)



#urban_ind1 <- data$NPI - quantile(data$NPI, p = 0.75,na.rm = T) # highly urbanised 
#NPI_ind <- data$NPI
#FI_ind <- data$F_I

#FI_basis6_npi <- FI_cb6*npi_cb6
#npi_basis6_FI <- npi_cb6*FI_ind

# Specifies a unique column name for the inla() model
# Note: GLM (), GAM () or GLM.nb () models are not required
#colnames(npi_cb) = paste0("npi_cb.", colnames(npi_cb))
#colnames(vero_cb) = paste0("vero_cb.", colnames(vero_cb))
#colnames(T_vero_cb) = paste0("T_vero_cb.", colnames(T_vero_cb))
#colnames(F_vero_cb) = paste0("F_vero_cb.", colnames(F_vero_cb))
#colnames(shannon_cb) = paste0("shannon_cb.", colnames(npi_cb))



