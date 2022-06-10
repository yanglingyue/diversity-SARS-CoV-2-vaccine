#########################################################################################################################################
#multi var model
#########################################################################################################################################
model_gr <- glm(GR_NCM_log_diff ~ FBW_gr_cb+lineage_gr_cb+npi_gr_cb+travel_gr_cb+NPI:F_I_BWP,family=gaussian,data=data)

summary(model_gr)

cp_FI_gr <- crosspred(FBW_gr_cb, model_gr, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices


pdf(file="02_fig/prop-05/fig4_a.pdf",width=8,height=8)
par(mar = c(5, 4.5, 4, 3))

nlag = 3
y <- cp_FI_gr$predvar
x <- seq(0, nlag, 0.1)
z <- t(cp_FI_gr$cumfit)
max(z)
min(z)

pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-3,0), 50)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0,max(y), length.out = ncol(z)),
               z,col=cols,level=levels,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect on growth rate",
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=2,
               key.axes=axis(4,family = "ArialN",cex.axis=2,seq(-3,0,by=1)),
               plot.axes = { axis(1,cex.axis=2, at = 0:nlag, c(0:nlag),family = "ArialN") 
                       axis(2,cex.axis=2,family = "ArialN")})
abline(h = 70, lty = 5,lwd=1)
dev.off()





model_si <- glm(sub_lineage ~ Lag(data$sub_lineage,1)+travel_si_cb+npi_si_cb+gr_cb+FBW_si_cb+NPI:F_I_BWP,family=gaussian,data=data)
#model1 <- glm(lineage ~ travel_cb+ npi_cb+FI_vero_cb+NPI:F_I,family=gaussian,data=data)

summary(model_si)

cp_FI_si <- crosspred(FBW_si_cb, model_si, cen=0,by=0.1,at=0:100,bylag=0.1,cumul=TRUE) #slices

#pdf("02.fig/Fig_3/Fig_3B-c1-contour.pdf", width=5, height=5)
pdf(file="02_fig/prop-05/fig4_b.pdf",width=8,height=8)
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
levels <- pretty(c(-3.5,0.5), 60)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
filled.contour(x= seq(0, 3, length.out = nrow(z)),
               y=seq(0, max(y), length.out = ncol(z)),
               z,
               ylim = range(y, finite = TRUE),
               xlab = "Lag, months", ylab = "Adjusted vaccine coverage", main = "Effect on lineage diversity",
               col = cols,level=levels,
               key.title=title("",cex.main=2,las=1,family = "ArialN"),
               family = "ArialN",cex.main=2,cex.lab=2,
               key.axes=axis(4,family = "ArialN",cex.axis=2,seq(-3.5,0.5,by=1)),
               plot.axes = { axis(1,cex.axis=2, at = 0:nlag, c(0:nlag),family = "ArialN") 
                 axis(2,cex.axis=2,family = "ArialN") })
abline(h = 70, lty = 5,lwd=1)
dev.off()



r_l <- cp_FI_gr $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag0"]
r_f <- cp_FI_gr $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag0"]
r_h <- cp_FI_gr $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag0"]

r_f_a <- cp_FI_gr $allfit[c("10","20","30","40","50","60","70","80","90","100")]
r_l_a <- cp_FI_gr $alllow[c("10","20","30","40","50","60","70","80","90","100")]

result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag0","r_l_lag0","r_h_lag0","CI_lag0") 
result_lag0 <- result



r_l <- cp_FI_gr $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag1"]
r_f <- cp_FI_gr $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag1"]
r_h <- cp_FI_gr $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag1"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag1","r_l_lag1","r_h_lag1","CI_lag1") 
result_lag1 <- result


r_l <- cp_FI_gr $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag2"]
r_f <- cp_FI_gr $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag2"]
r_h <- cp_FI_gr $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag2"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag2","r_l_lag2","r_h_lag2","CI_lag2") 
result_lag2 <- result



r_l <- cp_FI_gr $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag3"]
r_f <- cp_FI_gr $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag3"]
r_h <- cp_FI_gr $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag3"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag3","r_l_lag3","r_h_lag3","CI_lag3") 
result_lag3 <- result


result_gr <- cbind(result_lag0,result_lag1,result_lag2,result_lag3)
result_gr <- result_gr[,c("CI_lag0","CI_lag1","CI_lag2","CI_lag3")]
result_gr$Coverage <- rownames(result_gr)
result_gr <- result_gr[,c(5,1,2,3,4)]
colnames(result_gr) <- c("Coverages","Lag0","lag1","lag2","lag3")

write.csv(result_gr,"02_fig/prop-05/fig4_a.csv",row.names = F)



r_l <- cp_FI_si $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag0"]
r_f <- cp_FI_si $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag0"]
r_h <- cp_FI_si $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag0"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag0","r_l_lag0","r_h_lag0","CI_lag0") 
result_lag0 <- result



r_l <- cp_FI_si $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag1"]
r_f <- cp_FI_si $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag1"]
r_h <- cp_FI_si $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag1"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag1","r_l_lag1","r_h_lag1","CI_lag1") 
result_lag1 <- result


r_l <- cp_FI_si $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag2"]
r_f <- cp_FI_si $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag2"]
r_h <- cp_FI_si $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag2"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag2","r_l_lag2","r_h_lag2","CI_lag2") 
result_lag2 <- result



r_l <- cp_FI_si $cumlow[c("10","20","30","40","50","60","70","80","90","100"),"lag3"]
r_f <- cp_FI_si $cumfit[c("10","20","30","40","50","60","70","80","90","100"),"lag3"]
r_h <- cp_FI_si $cumhigh[c("10","20","30","40","50","60","70","80","90","100"),"lag3"]
result <- cbind(r_l,r_f,r_h)
result<-data.frame(result,CI=NA)

for (i in 1:10){
  result[i,4]<- paste(format(round(result[i,2],2))," (",format(round(result[i,1],2)),", ",format(round(result[i,3],2)),")",sep = "")}
colnames(result) <-c("r_l_lag3","r_l_lag3","r_h_lag3","CI_lag3") 
result_lag3 <- result


result_si <- cbind(result_lag0,result_lag1,result_lag2,result_lag3)
result_si <- result_si[,c("CI_lag0","CI_lag1","CI_lag2","CI_lag3")]
result_si$Coverage <- rownames(result_si)
result_si <- result_si[,c(5,1,2,3,4)]
colnames(result_si) <- c("Coverages","Lag0","lag1","lag2","lag3")

write.csv(result_si,"02_fig/prop-05/fig4_b.csv",row.names = F)






# set value
cen_vero_gr <- 0
min_vero_gr <- min(data$F_I_BWP, na.rm = T)
max_vero_gr <- max(data$F_I_BWP, na.rm = T)

# Forecast overall accumulation
FI_gr_cp <- crosspred(FBW_gr_cb, model_gr, cen=0,by=0.1,bylag=0.1,cumul=TRUE)
FI_low_gr_cp <- crosspred(FBW_gr_cb,model_gr, bylag = 0.1, at=seq(cen_vero_gr,min_vero_gr),by=0.1,cen=cen_vero_gr,cumul=T) 
FI_high_gr_cp <- crosspred(FBW_gr_cb,model_gr, bylag = 0.1, at=seq(min_vero_gr,100,by=0.1),cen=cen_vero_gr,cumul=T) 


# Plot
pdf("02_fig/prop-05/FIG-s1.pdf", width=15, height=7.5)

par(mfrow=c(1,2),mar = c(6.5,6.5, 4, 5),oma = c(2, 2, 3, 1),family="ArialN")
plot(FI_low_gr_cp,"slices",lag=1,cumul=T,col="black",ylim=c(-6,0),xlim=c(0,100),axes=F,ann=F,lwd=2.5,cex.axis=5)
lines(FI_high_gr_cp,"slices",lag=1,cumul=T,ci="area",col="black",lwd=3, ci.arg=list(col = grey(0.85)))
abline(h=0,col="black")
axis(1,at=0:4*25,cex.axis=2,family = "ArialN")
axis(2,at=c(-3:0*2),cex.axis=2,family = "ArialN")
title(xlab= "Adjusted vaccine coverage",cex.lab=2,family = "ArialN")
mtext("Effect on growth Rate",side = 2,at=-3, line =2.5,cex=2,family = "ArialN")
mtext(side = 2, at = 0.7, text = "A", las = 2, cex = 2.5, family = "ArialN",line = 5.5)


# set value
cen_vero_si <- 0
min_vero_si <- min(data$F_I_BWP, na.rm = T)
max_vero_si <- max(data$F_I_BWP, na.rm = T)

# Forecast overall accumulation
FI_si_cp <- crosspred(FBW_si_cb, model_si, cen=0,by=0.1,bylag=0.1,cumul=TRUE)
FI_low_si_cp <- crosspred(FBW_si_cb,model_si, bylag = 0.1, at=seq(cen_vero_si,min_vero_si),by=0.1,cen=cen_vero_si,cumul=TRUE) 
FI_high_si_cp <- crosspred(FBW_si_cb,model_si, bylag = 0.1, at=seq(min_vero_si,100,by=0.1),cen=cen_vero_si,cumul=TRUE) 


plot(FI_low_si_cp,"slices",lag=1,cumul=T,col="black",ylim=c(-6,2),xlim=c(0,100),axes=F,ann=F,lwd=2.5,cex.axis=5)
lines(FI_high_si_cp,"slices",lag=1,cumul=T,ci="area",col="black",lwd=3, ci.arg=list(col = grey(0.85)))
abline(h=0,col="black")
axis(1,at=0:4*25,cex.axis=2,family = "ArialN")
axis(2,at=c(-3:0*2,1:0*2),cex.axis=2,family = "ArialN")
title(xlab= "Adjusted vaccine coverage",cex.lab=2,family = "ArialN")
mtext("Effect on lineage diversity",side = 2,at=-2, line =2.5,cex=2,family = "ArialN",mgp=c(3,0.8,0))
mtext(side = 2, at = 3, text = "B", las = 2, cex = 2.5, family = "ArialN",line = 5.5)
dev.off()



