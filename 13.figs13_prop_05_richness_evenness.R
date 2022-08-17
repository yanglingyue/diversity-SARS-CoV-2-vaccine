model_gr_rich <- glm(GR_NCM_log_diff ~ FBW_gr_cb+richness_gr_cb+npi_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)
model_si_rich <- glm(sub_richness ~ Lag(data$sub_richness,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)

cp_FI_gr_rich <- crosspred(FBW_gr_cb, model_gr_rich, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices
cp_FI_si_rich <- crosspred(FBW_si_cb, model_si_rich, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


pdf(file="02_fig/figs13_a-richness.pdf",width=7.5,height=7.5)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_gr_rich$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_gr_rich$cumfit)
max(z)
min(z)

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-3,1), 50)
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
               key.axes=axis(4,family = "ArialN",cex.axis=1.8,seq(-3,1,by=1)),
               plot.axes = { axis(1,cex.axis=1.8, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=1.8,family = "ArialN")})
abline(h = 70, lty = 5,lwd=2.5)
mtext(side = 2, at = 110, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 2)
dev.off()




pdf(file="02_fig/figS13_richness.pdf",width=7.5,height=7.57)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_si_rich$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_si_rich$cumfit)
#z <- round(t(cp_FI_si$cumfit),3)
#z1 <- t(z)
max(z)
min(z)
max(y)
pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-60.2,30), 60)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0, max(y), length.out = ncol(z)),
               z,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect of richness\n of lineage diversity",
               col = cols,level=levels,
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=2,
               key.axes=axis(4,family = "ArialN",cex.axis=1.8,seq(-60,30,by=30)),
               plot.axes = { axis(1,cex.axis=1.8, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=1.8,family = "ArialN") })
abline(h = 70, lty = 5,lwd=2.5)
#mtext(side = 2, at = 110, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 2)

dev.off()








model_gr_even <- glm(GR_NCM_log_diff ~ FBW_gr_cb+evenness_gr_cb+npi_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)

model_si_even <- glm(sub_evenness ~ Lag(data$sub_evenness,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)


cp_FI_gr_even <- crosspred(FBW_gr_cb, model_gr_even, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices
cp_FI_si_even <- crosspred(FBW_si_cb, model_si_even, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


pdf(file="02_fig/figs13_a-evenness.pdf",width=7.5,height=7.5)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_gr_even$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_gr_even$cumfit)
max(z)
min(z)

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-3,1), 50)
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
               key.axes=axis(4,family = "ArialN",cex.axis=1.8,seq(-3,1,by=1)),
               plot.axes = { axis(1,cex.axis=1.8, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=1.8,family = "ArialN")})
abline(h = 70, lty = 5,lwd=2.5)
mtext(side = 2, at = 110, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 2)
dev.off()




pdf(file="02_fig/figS13_evenness.pdf",width=7.5,height=7.57)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_si_even$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_si_even$cumfit)
#z <- round(t(cp_FI_si$cumfit),3)
#z1 <- t(z)
max(z)
min(z)
max(y)
pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-1,0.5), 60)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0, max(y), length.out = ncol(z)),
               z,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect of evenness\n of lineage diversity",
               col = cols,level=levels,
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=2,
               key.axes=axis(4,family = "ArialN",cex.axis=1.8,seq(-1,0.5,by=0.5)),
               plot.axes = { axis(1,cex.axis=1.8, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=1.8,family = "ArialN") })
abline(h = 70, lty = 5,lwd=2.5)
#mtext(side = 2, at = 110, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 2)

dev.off()



