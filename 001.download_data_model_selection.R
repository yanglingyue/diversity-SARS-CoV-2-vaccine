setwd("D:/code_upload/vaccine-SARS-CoV-2-diversity")
rm(list = ls())
#Sys.setenv(LANGUAGE = "en")
#Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")
library(ggplot2);library(ggsci);library(RColorBrewer);library(tsModel);
library(tidyverse);library(tidyr);library(lubridate);library(mgcv);
library(MASS);library(splines);library(dlnm);library(showtext);
library(car);library(R330);library(boot);library(dplyr)
showtext_auto()
font_add("ArialN","arial.ttf")
font_families()

data <-read.csv("01_data/subsample_diversity_vaccine_data.csv")# Load data
colnames(data)[13] <- "sub_lineage" #Select the data for Shannon's index of lineage diversity calculated based on sampling strategy 1
data <- data %>% 
  group_by(month) %>% 
  mutate(country_index = row_number())


data$travel2 <- log(data$travel_route)
data$country <- as.factor(data$country)
data$country_index <- as.factor(data$country_index)

data<- data[order(data$country,data$month,decreasing=F),]

y1 <- which(data$month=="2019-12") 
y2 <- which(data$month=="2020-01") 
y3 <- which(data$month=="2020-02") 
y4 <- which(data$month=="2022-09")

data <-data[-c(y1,y2,y3,y4),]

rownames(data) <- data$record

# Set maximum lag
nlag = 3

# Creating lagged variables
# Define the monthly mean covariate lag term matrix
lag_npi <- tsModel::Lag(data$NPI_NV, group = data$country_index, k = 0:nlag)
lag_lineage <- tsModel::Lag(data$sub_lineage, group = data$country_index, k = 0:nlag)
lag_travel <- tsModel::Lag(data$travel2, group = data$country_index, k = 0:nlag)
lag_gr <- tsModel::Lag(data$GR_NCM_log_diff, group = data$country_index, k = 0:nlag)
lag_FW_vero <- tsModel::Lag(data$FW_BW_P, group = data$country_index, k = 0:nlag)

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
                                    knots = equalknots(data$FW_BW_P, 1)),
                      arglag = list(fun = "ns", df=2))
summary(FW_cb4)


var <- lag_FW_vero
FW_cb5 <- crossbasis(var,
                      argvar = list(fun = "ns",
                                    knots = equalknots(data$FW_BW_P, 2)),
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




# set data for models of growth rate of cases and lineage diversity
Y1  <- data$GR_NCM_log_diff # response variable (growth rate of cases)
Y2  <- data$sub_lineage # response variable (Shannon's index of lineage diversity based on the subsampling strategy 1)
P  <- data$NPI_NV # forpublic health and social measure (PHSM)
V <- data$FW_BW_P # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country_index # for country interaction with month random effect
M <- data$month # The month index of the dataset 

df <- data.frame(Y1, Y2, P, V, N, C, M)


# Stepwise methods were used to screen growth rate model variables using AIC criterion---lineage diversity and adjusted vaccine coverage
mylist <- list(lineage_cb1,lineage_cb2,lineage_cb3,lineage_cb4,lineage_cb5,
               FW_cb1,FW_cb2,FW_cb3,FW_cb4,FW_cb5)

names(mylist)<- c( "lineage_cb1","lineage_cb2","lineage_cb3","lineage_cb4","lineage_cb5",
                   "FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")

gam_list <- lapply(mylist, function(x){                                                              
  mgcv::gam(Y1 ~ + x +travel_cb+npi_cb
            +N+s(C, bs = "re")+P:V
            , data = df,   method = "REML") 
})   

result <- sapply(X = gam_list , FUN = AIC)
result <- as.data.frame(result) 
result
min(result[,1])
rownames(result)[which.min(result$result)]



mylist <- list(FW_cb1,FW_cb2,FW_cb3,FW_cb4,FW_cb5)

names(mylist)<- c("FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")

gam_list <- lapply(mylist, function(x){                                                              
  mgcv::gam(Y1 ~ + x +travel_cb+npi_cb+lineage_cb5
            +N+s(C, bs = "re")+P:V
            , data = df, method = "REML") 
})   

result <- sapply(X = gam_list , FUN = AIC)
result <- as.data.frame(result) 
result
min(result[,1])
rownames(result)[which.min(result$result)]



#Select the final model for growth rate with the minimum AIC 
model.gr <- gam(Y1  ~ travel_cb + npi_cb + lineage_cb5 + FW_cb3 +
                  N+s(C, bs = "re") + P:V,  data = df, method = "REML")
summary(model.gr)
AIC(model.gr)



#pdf("report/plot/report_growth_rate_fig5.pdf",height=6,width=10)
par(oma = c(0, 0, 1, 0),cex.axis=1,cex.lab=1,family = "ArialN")
plot(model.gr,family = "ArialN",cex.lab = 1,cex.main=1,cex.axis=1)
#dev.off()

#pdf("report/plot/report_growth_rate_fig6.pdf",height=5,width=10)

model.gr.res <- residuals(model.gr,type="deviance")
par(mfrow = c(1,2), mar = c(5,4,1,2)+0.2)
hist(model.gr.res,
     xlab = "Residuals", main = "",
     nclass = 20, prob = T,xlim=c(-2.2,4.4))
model.gr.res =
  model.gr.res[which(model.gr.res < abs(min(model.gr.res)))]
normFitDens = dnorm(-22:44/10, mean =  8.497202e-15, sd = 0.8218212)
lines(-22:44/10, normFitDens)
qqnorm(model.gr.res, main ="", col = "#00000077", pch = 16,ylim=c(-2,2))
qqline(model.gr.res, col="red")
#dev.off()



###################################################################################################
# Stepwise methods were used to screen lineage diversity model variables using AIC criterion---the case growth rate and adjusted vaccine coverage

mylist <- list(gr_cb1,gr_cb2,gr_cb3,gr_cb4,gr_cb5,
               FW_cb1,FW_cb2,FW_cb3,FW_cb4,FW_cb5)

names(mylist)<- c(  "gr_cb1","gr_cb2","gr_cb3","gr_cb4","gr_cb5",
                    "FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")



gam_list <- lapply(mylist, function(x){                                                              
  mgcv::gam(Y2~ + x + travel_cb + npi_cb +
              +Lag(df$Y2,1)+N+s(C, bs = "re")+P:V
            , data = df,   method = "REML") 
})   

result <- sapply(X = gam_list , FUN = AIC)
result <- as.data.frame(result) 
result
min(result[,1])
rownames(result)[which.min(result$result)]


mylist <- list(gr_cb1,gr_cb2,gr_cb3,gr_cb4,gr_cb5)

names(mylist)<- c("gr_cb1","gr_cb2","gr_cb3","gr_cb4","gr_cb5")

gam_list <- lapply(mylist, function(x){                                                              
  mgcv::gam(Y2~ + x + travel_cb + npi_cb + FW_cb3
            +Lag(df$Y2,1)+N+s(C, bs = "re")+P:V
            , data = df,   method = "REML") 
})   
result <- sapply(X = gam_list , FUN = AIC)
result <- as.data.frame(result) 
result
min(result[,1])
rownames(result)[which.min(result$result)]




#Select the final model for lineage diversity with the minimum AIC 
model.si <- gam(Y2  ~ travel_cb + npi_cb + gr_cb2 + FW_cb3 +
                  + Lag(df$Y2,1) + N + s(C, bs = "re") + P:V, data = df, method = "REML")
summary(model.si)
AIC(model.si)





#pdf("report/plot/report_diversity_fig6.pdf",height=5,width=10)
model.si.res <- residuals(model.si,type="deviance")
par(mfrow = c(1,2), mar = c(5,4,1,2)+0.2)
hist(model.si.res,
     xlab = "Residuals", main = "",
     nclass = 20, prob = T,xlim=c(-2,2),ylim=c(0,1))
model.si.res =
  model.si.res[which(model.si.res < abs(min(model.si.res)))]
normFitDens = dnorm(-20:20/10, mean = -3.493873e-16, sd = 0.4989637)
lines(-20:20/10, normFitDens)
qqnorm(model.si.res, main ="", col = "#00000077", pch = 16,ylim=c(-2,2))
qqline(model.si.res, col="red")
#dev.off()

