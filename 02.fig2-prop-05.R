model_gr <-  glm(GR_NCM_log_diff ~ npi_gr_cb+FBW_gr_cb+lineage_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)
summary(model_gr)

model_si <-glm(sub_lineage ~ Lag(data$sub_lineage,1)+travel_si_cb+gr_cb+npi_si_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)
summary(model_si)


cp_npi_gr <- crosspred(npi_gr_cb, model_gr, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr <- crosspred(travel_gr_cb, model_gr, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices
cp_lineage_gr <- crosspred(lineage_gr_cb, model_gr,at=0:40/10,cen=0,by=0.1,bylag=0.1,cumul=TRUE) #slices

cp_npi_si <- crosspred(npi_si_cb, model_si,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si <- crosspred(travel_si_cb, model_si, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_gr_si <- crosspred(gr_cb, model_si, cen=0,by=0.1,at=c(-35:0/10,0:95/10),bylag=0.1,cumul=TRUE) #slices

model_gr_rich <- glm(GR_NCM_log_diff ~ FBW_gr_cb+richness_gr_cb+npi_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)
summary(model_gr_rich)

model_si_rich <- glm(sub_richness ~ Lag(data$sub_richness,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)
summary(model_si_rich)

cp_npi_gr_rich <- crosspred(npi_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr_rich <- crosspred(travel_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices

cp_npi_si_rich <- crosspred(npi_si_cb, model_si_rich,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si_rich <- crosspred(travel_si_cb, model_si_rich, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices



model_gr_even <- glm(GR_NCM_log_diff ~ FBW_gr_cb+evenness_gr_cb+npi_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)
summary(model_gr_even)

model_si_even <- glm(sub_evenness ~ Lag(data$sub_evenness,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)
summary(model_si_even)

cp_npi_gr_even <- crosspred(npi_gr_cb, model_gr_even, cen=0,by=0.1,at=0:1000/10,bylag=0.1,cumul=TRUE) #slices
cp_travel_gr_even <- crosspred(travel_gr_cb, model_gr_even, cen=0,by=0.1,at=0:150/10,bylag=0.1,cumul=TRUE) #slices

cp_npi_si_even <- crosspred(npi_si_cb, model_si_even,cen=0,at=0:1000/10,by=0.1,bylag=0.1,cumul=TRUE) #slices
cp_travel_si_even <- crosspred(travel_si_cb, model_si_even, cen=0,at=0:150/10,by=0.1,bylag=0.1,cumul=TRUE) #slices


#si_m <- mean(data$lineage,na.rm=T)
#it_m <- mean(data$travel2,na.rm=T)

#r_l <- cp_travel_si $cumlow["0.332",]/si_m 
#r_f <- cp_travel_si $cumfit["0.332",]/si_m 
#r_h <- cp_travel_si $cumhigh["0.332",]/si_m 

#0.09753737/si_m 

pdf("02_fig/prop-05/fig3-3.pdf", width = 21, height = 21/3*2)

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
mtext("Effect on growth rate of cases and \nShannon's index of lineage diversity",side = 2,at=0, line =4,cex=1.7,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Shannon's index of lineage diversity"),
                  
                  paste0("Growth rate of cases")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 2.52, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


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

legend("bottomleft",
       legend = c(paste0("Evenness of lineage diversity"),
                  
                  paste0("Growth rate of cases")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 2.45, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 5.5)

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
     ylim = range(-8, 24), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1,at=0:5*3,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-1:0*8,3:0*8),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and \nrichness of lineage diversity",side = 2,at=8, line =4,cex=1.7,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Richness of lineage diversity"),
                  
                  paste0("Growth rate of cases")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 28, text = "C", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


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
     ylim = range(-4, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-4:0*1),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and\nShannon's index of lineage diversity",side = 2,at=-2, line =4,cex=1.7,family = "ArialN")

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
     ylim = range(-4, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-4:0*1),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and\nevenness of lineage diversity",side = 2,at=-2, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Evenness of lineage diversity"),
#                  paste0("Growth rate of cases")),
#       col = c(col2, col1), 
#       lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.5)
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
     ylim = range(-24, 8), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*8,1:0*8),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate of cases and \nrichness of lineage diversity",side = 2,at=-8, line =4,cex=1.7,family = "ArialN")

#legend("bottomleft",
#       legend = c(paste0("Richness of lineage diversity"),
#                                    paste0("Growth rate of cases")),
#      col = c(col2, col1), 
#    lwd = 3, lty = 1, bty = "n", 
#       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 11, text = "F", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


dev.off()





