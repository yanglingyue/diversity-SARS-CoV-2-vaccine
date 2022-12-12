source('001.download_data_parameter.R')
# set data for models of growth rate of cases and lineage diversity
Y1  <- data$GR_NCM_log_diff # response variable (growth rate of cases)
Y2  <- data$sub_lineage # response variable (Shannon's index of lineage diversity based on the subsampling strategy 1)
P  <- data$NPI # forpublic health and social measure (PHSM)
V <- data$F_I_BWP # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country # for country interaction with month random effect
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
  mgcv::gam(Y1 ~ + x +travel_cb+npi_cb+lineage_cb3
            +N+s(C, bs = "re")+P:V
            , data = df, method = "REML") 
})   

result <- sapply(X = gam_list , FUN = AIC)
result <- as.data.frame(result) 
result
min(result[,1])
rownames(result)[which.min(result$result)]



#Select the final model for growth rate with the minimum AIC 
model.gr <- gam(Y1  ~ travel_cb + npi_cb + lineage_cb3 + FW_cb3 +
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
model.si <- gam(Y2  ~ travel_cb + npi_cb + gr_cb3 + FW_cb3 +
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
