
rm(list = ls())
library(ggplot2);library(ggsci);library(tidyr);library(car)
library(corrplot);library(ggcorrplot);library(RColorBrewer)
library(tidyverse);library(naniar);library(mgcv);library(MASS)
library(lubridate);library(mvabund);library(dlnm);library(rEDM)

data<-read.csv("subsample_prop_05.csv")
data1 <- data
#view(data)
data$date <- parse_date_time(data$month, "ym")
data$date <-as.Date(as.POSIXct(data$date,tz="Hongkong"))
data$month1 <- month(data$date)
data$year <- year(data$date)
colnames(data)[29]<-"Growth_Rate"
colnames(data)[55]<-"Shannon"
colnames(data)[23]<-"vaccine"
y1 <- which(data$month=="2019-12")
y2 <- which(data$month=="2020-01")
y3 <- which(data$month=="2020-02")
data <-data[-c(y1,y2,y3),]
data$travel2 <- log(data$travel_route)
data$travel2<-as.numeric(data$travel2)

normalize <- function(x,na.rm=TRUE, ...) {
  (x - mean(x,na.rm=TRUE, ...))/sd(x,na.rm=TRUE, ...)
}
# separate time column from data
vars <- c("Shannon", "Growth_Rate", "NPI", "vaccine","travel2")
composite_ts <- data[, vars]

# normalize each time series within a country
data_by_country <- split(composite_ts, data$country)
normalized_data <- lapply(data_by_country, function(df) sapply(df, normalize))
w=(1:24)

data$month1=rep(w,60)
composite_ts <- cbind(month1 = data$month1, data.frame(do.call(rbind, normalized_data)))


segments_end <- cumsum(sapply(data_by_country, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)


# Choose random segments for prediction
set.seed(2315)

rndlib <- sample(1:NROW(segments), floor(NROW(segments) * 0.75))
composite_lib <- segments[rndlib, ]
composite_pred <- segments[-rndlib, ]

#Quantifying predictability 

simplex_out <- lapply(vars, function(var) {
  simplex(composite_ts[, c("month1", var)], E = 2:10, lib = composite_lib, pred = composite_pred)
})

names(simplex_out) <- vars

par(mfrow = c(3, 2))
for (var in names(simplex_out)) {
  plot(simplex_out[[var]]$E, simplex_out[[var]]$rho, type = "l", xlab = "Embedding Dimension (E)", 
       ylab = "Forecast Skill (rho)", main = var)
}

best_E <- sapply(simplex_out, function(df) {
  df$E[which.max(df$rho)]
})
best_E
# dev.off()
#quantifying nonlinearity
smap_out <- lapply(vars, function(var) {
  s_map(composite_ts[, c("month1", var)], E = best_E[var], lib = composite_lib, 
        pred = composite_pred)
})
names(smap_out) <- names(simplex_out)

par(mfrow = c(3, 2))
for (var in names(smap_out)) {
  plot(smap_out[[var]]$theta, smap_out[[var]]$rho, type = "l", xlab = "Nonlinearity (theta)", 
       ylab = "Forecast Skill (rho)", main = var)
}




#Convergent Cross Mapping 


lib_sizes <- c(seq(20, 300, by = 20), seq(320, 1440, by = 80))
GR_xmap_SI <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Growth_Rate", 
                  target_column = "Shannon", E = best_E["Growth_Rate"], lib_sizes = lib_sizes, 
                  silent = TRUE)
SI_xmap_GR <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Shannon", 
                  target_column = "Growth_Rate", E = best_E["Shannon"], lib_sizes = lib_sizes, 
                  silent = TRUE)
#1
write.table(GR_xmap_SI[,c(7,9)], file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/GR_xmap_SI.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table(SI_xmap_GR[,c(7,9)], file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/SI_xmap_GR.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
GR_xmap_SI_means <- ccm_means(GR_xmap_SI)
SI_xmap_GR_means <- ccm_means(SI_xmap_GR)

plot(GR_xmap_SI_means$lib_size, pmax(0, GR_xmap_SI_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0,1), lwd = 2)
lines(SI_xmap_GR_means$lib_size, pmax(0, SI_xmap_GR_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("GR xmap SI", 
                                                                  "SI xmap GR"), inset = 0.02, bty = "n", cex = 0.8)
abline(h = 0, lty = 3)

####
GR_xmap_NPI <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Growth_Rate", 
                   target_column = "NPI", E = best_E["Growth_Rate"], lib_sizes = lib_sizes, 
                   silent = TRUE)
NPI_xmap_GR <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "NPI", 
                   target_column = "Growth_Rate", E = best_E["NPI"], lib_sizes = lib_sizes, 
                   silent = TRUE)
write.table( GR_xmap_NPI[,c(7,9)],file="C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/GR_xmap_NPI.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table( NPI_xmap_GR[,c(7,9)],file="C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/NPI_xmap_GR.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
GR_xmap_NPI_means <- ccm_means(GR_xmap_NPI)
NPI_xmap_GR_means <- ccm_means(NPI_xmap_GR)

plot(GR_xmap_NPI_means$lib_size, pmax(0, GR_xmap_NPI_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(NPI_xmap_GR_means$lib_size, pmax(0, NPI_xmap_GR_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("GR xmap NPI", 
                                                                  "NPI xmap GR"), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)
####
GR_xmap_Vero <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Growth_Rate", 
                    target_column = "vaccine", E = best_E["Growth_Rate"], lib_sizes = lib_sizes, 
                    silent = TRUE)
Vero_xmap_GR <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "vaccine", 
                    target_column = "Growth_Rate", E = best_E["vaccine"], lib_sizes = lib_sizes, 
                    silent = TRUE)
write.table( GR_xmap_Vero[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/GR_xmap_Vero.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table( Vero_xmap_GR[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/Vero_xmap_GR.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
GR_xmap_Vero_means <- ccm_means(GR_xmap_Vero)
Vero_xmap_GR_means <- ccm_means(Vero_xmap_GR)

plot(GR_xmap_Vero_means$lib_size, pmax(0, GR_xmap_Vero_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(Vero_xmap_GR_means$lib_size, pmax(0, Vero_xmap_GR_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("GR xmap VP", 
                                                                  "VP xmap GR"), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)
####
SI_xmap_Vero <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Shannon", 
                    target_column = "vaccine", E = best_E["Shannon"], lib_sizes = lib_sizes, 
                    silent = TRUE)
Vero_xmap_SI <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "vaccine", 
                    target_column = "Shannon", E = best_E["vaccine"], lib_sizes = lib_sizes, 
                    silent = TRUE)
write.table( SI_xmap_Vero[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/SI_xmap_Vero.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table( Vero_xmap_SI[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/Vero_xmap_SI.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
SI_xmap_Vero_means <- ccm_means(SI_xmap_Vero)
Vero_xmap_SI_means <- ccm_means(Vero_xmap_SI)

plot(SI_xmap_Vero_means$lib_size, pmax(0, SI_xmap_Vero_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(Vero_xmap_SI_means$lib_size, pmax(0, Vero_xmap_SI_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("SI xmap VP", 
                                                                  "VP xmap SI"), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)
####
SI_xmap_NPI <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Shannon", 
                   target_column = "NPI", E = best_E["Shannon"], lib_sizes = lib_sizes, 
                   silent = TRUE)
NPI_xmap_SI <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "NPI", 
                   target_column = "Shannon", E = best_E["NPI"], lib_sizes = lib_sizes, 
                   silent = TRUE)
write.table( SI_xmap_NPI[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/SI_xmap_NPI.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table( NPI_xmap_SI[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/NPI_xmap_SI.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
SI_xmap_NPI_means <- ccm_means(SI_xmap_NPI)
NPI_xmap_SI_means <- ccm_means(NPI_xmap_SI)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) #
plot(SI_xmap_NPI_means$lib_size, pmax(0, SI_xmap_NPI_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(NPI_xmap_SI_means$lib_size, pmax(0, NPI_xmap_SI_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("SI xmap NPI", 
                                                                  "NPI xmap SI "), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)
####
GR_xmap_TR <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "travel2", 
                  target_column = "Growth_Rate", E = best_E["travel2"], lib_sizes = lib_sizes, 
                  silent = TRUE)
TR_xmap_GR <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Growth_Rate", 
                  target_column = "travel2", E = best_E["Growth_Rate"], lib_sizes = lib_sizes, 
                  silent = TRUE)
write.table(  GR_xmap_TR[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/GR_xmap_TR.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table(  TR_xmap_GR[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/TR_xmap_GR.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
GR_xmap_TR_means <- ccm_means(GR_xmap_TR)
TR_xmap_GR_means <- ccm_means(TR_xmap_GR)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) #
plot(GR_xmap_TR_means$lib_size, pmax(0, GR_xmap_TR_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(TR_xmap_GR_means$lib_size, pmax(0, TR_xmap_GR_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("GR_xmap_TR", 
                                                                  "TR_xmap_GR "), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)
####
SI_xmap_TR <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "travel2", 
                  target_column = "Shannon", E = best_E["travel2"], lib_sizes = lib_sizes, 
                  silent = TRUE)
TR_xmap_SI <- ccm(composite_ts, lib = segments, pred = segments, lib_column = "Shannon", 
                  target_column = "travel2", E = best_E["Shannon"], lib_sizes = lib_sizes, 
                  silent = TRUE)
write.table(  SI_xmap_TR[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/SI_xmap_TR.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
write.table(  TR_xmap_SI[,c(7,9)],file = "C:/Users/nodiet/Desktop/Global_ccm/prop_05_country_0517/ccm_results2315/TR_xmap_SI.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")
SI_xmap_TR_means <- ccm_means(SI_xmap_TR)
TR_xmap_SI_means <- ccm_means(TR_xmap_SI)
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0)) #
plot(SI_xmap_TR_means$lib_size, pmax(0, SI_xmap_TR_means$rho), type = "l", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", col = "red", ylim = c(0, 
                                                                                  1), lwd = 2)
lines(TR_xmap_SI_means$lib_size, pmax(0, TR_xmap_SI_means$rho), col = "blue", 
      lwd = 2)
legend(x = "topleft", col = c("red", "blue"), lwd = 2, legend = c("SI_xmap_TR", 
                                                                  "TR_xmap_SI"), inset = 0.02, bty = "n", cex = 0.8)

abline(h = 0, lty = 3)