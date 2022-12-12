rm(list = ls())

data <-read.csv("01_data/subsample_total_2000_lineage_diversity.csv")

y1 <- which(data$month=="2019-12") 
y2 <- which(data$month=="2020-01") 
y3 <- which(data$month=="2020-02") 
data <-data[-c(y1,y2,y3),]

source('003.load-data-cb-function.R')

# set data for models of growth rate of cases and lineage diversity
Y1  <- data$GR_NCM_log_diff # response variable (growth rate of cases)
Y2  <- data$sub_lineage # response variable (Shannon's index of lineage diversity based on the subsampling strategy 1)
Y3  <- data$sub_richness # response variable (the richness of lineage diversity based on the subsampling strategy 1)
Y4  <- data$sub_evenness # response variable (the evenness of lineage diversity based on the subsampling strategy 1)
P  <- data$NPI # forpublic health and social measure (PHSM)
V <- data$F_I_BWP # for adjusted vaccine coverage
N <- log(data$new_cases_month) # Log transformation of the number of new cases per month 
C <- data$country # for country interaction with month random effect
M <- data$month # The month index of the dataset 

df <- data.frame(Y1, Y2, Y3, Y4, P, V, N, C, M)



model_gr <-  gam(Y1 ~ npi_gr_cb+FBW_gr_cb+lineage_gr_cb+travel_gr_cb+N+s(C, bs = "re")+P:V
                 ,method = "REML",data=df)
summary(model_gr)
#AIC(model_gr)
cp_npi_gr <- crosspred(npi_gr_cb, model_gr, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr <- crosspred(travel_gr_cb, model_gr, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
#cp_lineage_gr <- crosspred(lineage_gr_cb, model_gr,at=0:40/10,cen=0,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_FI_gr <- crosspred(FBW_gr_cb, model_gr, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


model_si <-gam(Y2 ~ Lag(df$Y2,1)+travel_si_cb+gr_cb+npi_si_cb+FBW_si_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_si)
cp_npi_si <- crosspred(npi_si_cb, model_si,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si <- crosspred(travel_si_cb, model_si, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_gr_si <- crosspred(gr_cb, model_si, cen=0,by=0.1,at=c(-35:0/10,0:95/10),bylag=0.1,cumul=TRUE) #slices
cp_FI_si <- crosspred(FBW_si_cb, model_si, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


model_gr_rich <- gam(Y1 ~ FBW_gr_cb+richness_gr_cb+npi_gr_cb+travel_gr_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_gr_rich)
cp_npi_gr_rich <- crosspred(npi_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr_rich <- crosspred(travel_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
cp_FI_gr_rich <- crosspred(FBW_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


model_si_rich <- gam(Y3 ~ Lag(df$Y3,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_si_rich)
cp_npi_si_rich <- crosspred(npi_si_cb, model_si_rich,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si_rich <- crosspred(travel_si_cb, model_si_rich, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_FI_si_rich <- crosspred(FBW_si_cb, model_si_rich, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


model_gr_even <- gam(Y1 ~ FBW_gr_cb+evenness_gr_cb+npi_gr_cb+travel_gr_cb+N+s(C, bs = "re")+P:V,method = "REML",data=df)
summary(model_gr_even)
cp_npi_gr_even <- crosspred(npi_gr_cb, model_gr_even, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr_even <- crosspred(travel_gr_cb, model_gr_even, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
cp_FI_gr_even <- crosspred(FBW_gr_cb, model_gr_even, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


model_si_even <- gam(Y4 ~ Lag(df$Y4,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+N+s(C, bs = "re")+P:V,method = "REML",data=data)
summary(model_si_even)
cp_npi_si_even <- crosspred(npi_si_cb, model_si_even,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si_even <- crosspred(travel_si_cb, model_si_even, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_FI_si_even <- crosspred(FBW_si_cb, model_si_even, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


#si_m <- mean(data$lineage,na.rm=T)
#it_m <- mean(data$travel2,na.rm=T)

#r_l <- cp_travel_si $cumlow["0.332",]/si_m 
#r_f <- cp_travel_si $cumfit["0.332",]/si_m 
#r_h <- cp_travel_si $cumhigh["0.332",]/si_m 


#0.09753737/si_m 

pdf("02_fig/figs7-2000-a-f.pdf", width = 21, height = 21/3*2)

par(mfrow=c(2,3),mar = c(8, 8, 5, 4),oma = c(0, 2, 0, 0))

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
#mx2 <- 1

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
     ylim = range(-3, 2), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.8,line =5)
axis(1,at=0:5*3,cex.axis=2.8,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*1,2:0*1),cex.axis=2.8,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and \nShannon's index of lineage diversity",side = 2,at=-0.5, line =4,cex=1.9,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Shannon's index of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.4)
mtext(side = 2, at = 2.4, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 5.5)






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
     ylim = range(-3, 1), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.8,line =5)
axis(1,at=0:5*3,cex.axis=2.8,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*1,1:0*1),cex.axis=2.8,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and \nevenness of lineage diversity",side = 2,at=-1, line =4,cex=1.9,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Evenness of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.4)
mtext(side = 2, at = 1.35, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 5.5)

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
     ylim = range(-4, 12), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.8,line =5)
axis(1,at=0:5*3,cex.axis=2.8,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-1:0*4,3:0*4),cex.axis=2.8,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and \nrichness of lineage diversity",side = 2,at=4, line =4,cex=1.9,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Richness of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.4)
mtext(side = 2, at = 13, text = "C", las = 2, cex = 2.5, family = "ArialN",line = 5.5)




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
     ylim = range(-5, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.8,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.8,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-5:0*1),cex.axis=2.8,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and\nShannon's index of lineage diversity",side = 2,at=-2.5, line =4,cex=1.9,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Shannon's index of lineage diversity"),
                                   paste0("Growth rate of cases")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
      y.intersp = 1.5, horiz = F,cex=2.4)
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
     ylim = range(-5, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.8,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.8,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-5:0*1),cex.axis=2.8,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and\nevenness of lineage diversity",side = 2,at=-2.5, line =4,cex=1.9,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Evenness of lineage diversity"),
                  paste0("Growth rate of cases")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.4)
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
     ylim = range(-6, 4), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.8,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.8,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*2,2:0*2),cex.axis=2.8,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and \nrichness of lineage diversity",side = 2,at=-1, line =4,cex=1.9,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Richness of lineage diversity"),
                                    paste0("Growth rate of cases")),
      col = c(col2, col1), 
    lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.4)
mtext(side = 2, at = 5, text = "F", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


dev.off()






pdf(file="02_fig/figs7-2000-g.pdf",width=7.5,height=7.5)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_gr$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_gr$cumfit)
max(z)
min(z)

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-8,0), 50)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0,max(y), length.out = ncol(z)),
               z,col=cols,level=levels,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect on growth rate of cases",
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=1.8,
               key.axes=axis(4,family = "ArialN",cex.axis=1.8,las=1,seq(-8,0,by=2)),
               plot.axes = { axis(1,cex.axis=1.8, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=1.8,family = "ArialN",las=3)})
abline(h = 70, lty = 5,lwd=2.5)
mtext(side = 2, at = 110, text = "G", las = 2, cex = 2.5, family = "ArialN",line = 2)

dev.off()







pdf(file="02_fig/figs7-2000-h.pdf",width=7.5,height=7.5)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_si$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_si$cumfit)
#z <- round(t(cp_FI_si$cumfit),3)
#z1 <- t(z)
max(z)
min(z)
max(y)
pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-1.9,0.6), 60)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0, max(y), length.out = ncol(z)),
               z,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect of Shannon's index\n of lineage diversity",
               col = cols,level=levels,
               key.title=title("",cex.main=2,las=3,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=1.8,
               key.axes=axis(4,family = "ArialN",cex.axis=1.8,las=1,seq(-1.9,0.6,by=0.5)),
               plot.axes = { axis(1,cex.axis=1.8, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=1.8,family = "ArialN",las=3) })
abline(h = 70, lty = 5,lwd=2.5)
mtext(side = 2, at = 110, text = "H", las = 2, cex = 2.5, family = "ArialN",line = 2)

dev.off()



pdf(file="02_fig/figs7-2000-i.pdf",width=7.65,height=8.4)
par(mar = c(8, 8, 5, 2),oma = c(0, 1, 0, 1),family="ArialN")

#par(mfrow=c(1,2),mar = c(5,5, 1, 5),oma = c(1, 1, 1, 1))

cen_gr1 <- 0
min_gr1 <-  min(data$GR_NCM_log_diff, na.rm = T)
max_gr1 <- max(data$GR_NCM_log_diff, na.rm = T)
y <- cp_gr_si$predvar

# Forecast overall accumulation
gr_cp <- crosspred(gr_cb, model_si, cen=0,by=0.1,bylag=0.1,cumul=TRUE)
gr_low_cp <- crosspred(gr_cb,model_si, bylag = 0.1, at=seq(c(min_gr1-1),cen_gr1,by=0.1),cen=cen_gr1,cumul=TRUE) 
gr_high_cp <- crosspred(gr_cb,model_si, bylag = 0.1, at=seq(cen_gr1,max_gr1,by=0.1),cen=cen_gr1,cumul=TRUE) 


col2 <- brewer.pal(11, "RdBu")[3]
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/4, max = 255))

# Plot
#par(mar = c(5, 4, 4, 4))
plot(gr_low_cp,"slices",lag=3,cumul=T,col="black",ylim=c(-14,7),xlim=c(-5,10),axes=F,ann=F,lwd=3,cex.axis=1.9,
     ci.arg=list(col = grey(0.85)))
lines(gr_high_cp,"slices",lag=3,cumul=T,ci="area",col="black",lwd=3, ci.arg=list(col = grey(0.85)))
abline(h=0,col="black",lty = 3)
axis(1,at=c(seq(-5,10,3)),cex.axis=1.8,family = "ArialN",mgp=c(3,1,0))
axis(2,at=c(-2:0*7,1:0*7),cex.axis=1.8,family = "ArialN",las=3,mgp=c(3,1,0))
title(xlab= "Growth rate of cases",cex.lab=1.8,family = "ArialN",line =3)
#mtext("Time lag: 0 months",side = 3,line =1,cex=1.8,family = "ArialN",font=2)
mtext("Effect on Shannon's index of lineage diversity",side = 2,at=-3.5, line =3,cex=1.8,family = "ArialN")
mtext(side = 2, at = 15, text = "I", las = 2, cex = 2.5, family = "ArialN",line = 5.5)
dev.off()






