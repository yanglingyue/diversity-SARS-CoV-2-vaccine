#########################################################################################################################################
#multi var model
#########################################################################################################################################
model_gr <- glm(GR_NCM_log_diff ~ FBW_gr_cb+lineage_gr_cb+npi_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)

summary(model_gr)

cp_FI_gr <- crosspred(FBW_gr_cb, model_gr, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


pdf(file="02_fig/fig3_a.pdf",width=8,height=8)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_gr$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_gr$cumfit)
max(z)
min(z)

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-3.5,0.5), 50)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0,max(y), length.out = ncol(z)),
               z,col=cols,level=levels,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect on growth rate of cases",
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=2,
               key.axes=axis(4,family = "ArialN",cex.axis=2,seq(-3.5,0.5,by=1)),
               plot.axes = { axis(1,cex.axis=2, at = 0:nlag, c(0:nlag),family = "ArialN") 
                       axis(2,cex.axis=2,family = "ArialN")})
abline(h = 70, lty = 5,lwd=2.5)
dev.off()





model_si <- glm(sub_lineage ~ Lag(data$sub_lineage,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)

summary(model_si)

cp_FI_si <- crosspred(FBW_si_cb, model_si, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices

pdf(file="02_fig/fig3_b.pdf",width=8,height=8)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_si$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_si$cumfit)
max(z)
min(z)
max(y)
pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-3.5,0.5), 60)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0, max(y), length.out = ncol(z)),
               z,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect on Shannon diversity",
               col = cols,level=levels,
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=2,
               key.axes=axis(4,family = "ArialN",cex.axis=2,seq(-3.5,0.5,by=1)),
               plot.axes = { axis(1,cex.axis=2, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=2,family = "ArialN") })
abline(h = 70, lty = 5,lwd=2.5)
dev.off()



pdf(file="02_fig/fig3_c.pdf",width=8,height=8)
par(mar = c(5, 4.5, 4, 3),oma = c(0, 1, 0, 1),family="ArialN")

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
plot(gr_low_cp,"slices",lag=3,cumul=T,col="black",ylim=c(-30,10),xlim=c(-5,10),axes=F,ann=F,lwd=3,cex.axis=2,
     ci.arg=list(col = grey(0.85)))
lines(gr_high_cp,"slices",lag=3,cumul=T,ci="area",col="black",lwd=3, ci.arg=list(col = grey(0.85)))
abline(h=0,col="black",lty = 3)
axis(1,at=c(seq(-5,10,3)),cex.axis=2,family = "ArialN",mgp=c(3,2,0))
axis(2,at=c(-3:0*10,1:0*10),cex.axis=2,family = "ArialN",las=1,mgp=c(3,1,0))
title(xlab= "Growth rate of cases",cex.lab=2,family = "ArialN",line =4)
#mtext("Time lag: 0 months",side = 3,line =1,cex=1.8,family = "ArialN",font=2)
mtext("Effect on Shannon's index of lineage diversity",side = 2,at=-10, line =4,cex=2,family = "ArialN")
#par(new=T)
#hist(data$GR_NCM_log_diff,xlim=c(0,max(y)),ylim=c(0,850),axes=F,ann=F,col=tcol2,breaks=20)
#axis(2,at=c(0:2*160),cex.axis=2,family = "ArialN",line = 1)
#mtext("Frequency",side =2,at=160, line =4,cex=1.8,family = "ArialN")
mtext(side = 2, at = 17.5, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 5.5)
dev.off()



