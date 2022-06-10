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


si_m <- mean(data$lineage,na.rm=T)
it_m <- mean(data$travel2,na.rm=T)

r_l <- cp_travel_si $cumlow["1",]/si_m 
r_f <- cp_travel_si $cumfit["1",]/si_m 
r_h <- cp_travel_si $cumhigh["1",]/si_m 


0.09753737/si_m 

pdf("02_fig/prop-05/fig3.pdf", width = 18, height = 6)

par(mfrow=c(1,3),mar = c(6.5,6.5, 4, 5),oma = c(2, 2, 3, 1),family="ArialN")

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
plot(gr_low_cp,"slices",lag=3,cumul=T,col="black",ylim=c(-30,10),xlim=c(-5,10),axes=F,ann=F,lwd=3,cex.axis=5,
     ci.arg=list(col = grey(0.85)))
lines(gr_high_cp,"slices",lag=3,cumul=T,ci="area",col="black",lwd=3, ci.arg=list(col = grey(0.85)))
abline(h=0,col="black",lty = 3)
axis(1,at=c(seq(-5,10,3)),cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*10,1:0*10),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
title(xlab= "Growth rate",cex.lab=2.5,family = "ArialN",line =5)
#mtext("Time lag: 0 months",side = 3,line =1,cex=1.8,family = "ArialN",font=2)
mtext("Effect on lineage diversity",side = 2,at=-10, line =4,cex=1.7,family = "ArialN")

mtext(side = 2, at = 17.5, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 5.5)



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
     xlab = "International Travel", ylab = "", main = "", 
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
mtext("Effect on growth rate and \n lineage diversity",side = 2,at=0, line =4,cex=1.7,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Lineage diversity"),
                  
                  paste0("Growth rate")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 2.72, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 5.5)





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

col2 <- "#ED0000FF"
tcol2 <- do.call(rgb, c(as.list(col2rgb(col2)), alpha = 255/3, max = 255))

# define x values (var, by var)
varby <- seq(0, 100, 0.1)

# growth rate
#par(mar = c(5, 6.5, 5, 2),family="ArialN")
plot(varby, gr[,lag_gr], col = col1, type = "l", lwd = 3, 
     xlab = "Stringency index", ylab = "", main = "", 
     ylim = range(-8, 0), frame.plot = F, axes = F,family = "ArialN",cex.lab=2.5,line =5)
axis(1, at = 0:5*20, labels = 0:5*20,cex.axis=2.5,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-4:0*2),cex.axis=2.5,family = "ArialN",mgp=c(3,1,0))
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
mtext("Effect on growth rate and \n lineage diversity",side = 2,at=-4, line =4,cex=1.7,family = "ArialN")

legend("bottomleft",
       legend = c(paste0("Lineage diversity"),
                  
                  paste0("Growth rate")),
       col = c(col2, col1), 
       lwd = 3, lty = 1, bty = "n", 
       y.intersp = 1.5, horiz = F,cex=2.5)
mtext(side = 2, at = 1.48, text = "C", las = 2, cex = 2.5, family = "ArialN",line = 5.5)



dev.off()





