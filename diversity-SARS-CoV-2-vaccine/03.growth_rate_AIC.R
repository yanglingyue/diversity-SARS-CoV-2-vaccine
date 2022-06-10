vars <- c("npi_cb1","npi_cb2","npi_cb3","npi_cb4","npi_cb5",
          "lineage_cb1","lineage_cb2","lineage_cb3","lineage_cb4","lineage_cb5",
          "FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")
result <- c()
for(i in 1:length(vars)){
  model <- glm(substitute(GR_NCM_log_diff  ~ x +travel_cb1+NPI:F_I_BWP, 
                          list(x=as.name(vars[i]))), family="gaussian",
                   data=data)
  result <- rbind(result,c(vars[i],model[["aic"]]))
}
result                  
min(result[,2])


vars <- c("npi_cb1","npi_cb2","npi_cb3","npi_cb4","npi_cb5",
          "FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")
result <- c()
for(i in 1:length(vars)){
  model <- glm(substitute(GR_NCM_log_diff  ~ x +travel_cb1+lineage_cb3+NPI:F_I_BWP, 
                          list(x=as.name(vars[i]))),family="gaussian",
               data=data)
  result <- rbind(result,c(vars[i],model[["aic"]]))
}
result                  
min(result[,2])


vars <- c("FW_cb1","FW_cb2","FW_cb3","FW_cb4","FW_cb5")
result <- c()
for(i in 1:length(vars)){
  model <- glm(substitute(GR_NCM_log_diff  ~ x+travel_cb1+lineage_cb3+npi_cb2+NPI:F_I_BWP,
                          list(x=as.name(vars[i]))),family="gaussian",
               data=data)
  result <- rbind(result,c(vars[i],model[["aic"]]))
}
result                  
min(result[,2])





model.gr <- glm(GR_NCM_log_diff  ~ travel_cb1 + lineage_cb3+npi_cb2+FW_cb5+
               NPI:F_I_BWP, family="gaussian",data=data)
summary(model.gr)


pdf("report/plot/report_growth_rate_fig5.pdf",height=6,width=10)
par(mfrow=c(2,2),oma = c(0, 0, 1, 0),cex.axis=1,cex.lab=1,family = "ArialN")
plot(model.gr,family = "ArialN",cex.lab = 1,cex.main=1,cex.axis=1)
dev.off()

pdf("report/plot/report_growth_rate_fig6.pdf",height=5,width=10)

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
dev.off()

shapiro.test(model.gr.res)



data$number <- rownames(data)
data1 <- data[,c("month","country","number")]
model.gr.res <- residuals(model.gr,type="deviance")
residuals <- as.data.frame(model.gr.res)
residuals$number <- rownames(residuals)
result.res <- merge(data1,residuals,by=c("number"),all.x = T)



pal <- rev(brewer.pal(11, "RdBu"))
levels <- pretty(c(-5, 5), 40)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))


p <- ggplot(result.res, aes(month,country)) +
  geom_tile(aes(fill=model.gr.res),color = "white",size=0.6) +
  scale_fill_gradientn(colors=cols, limits = c(-5,5),na.value="grey88") +
  theme_bw()+theme(panel.grid=element_blank(),
                   axis.line=element_line(size=0.3,colour="black"))+
  theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 0),
        axis.text.y = element_text(size = 9,color="black",vjust = 0.5))+
  scale_x_discrete(breaks=c("2020-04","2020-07","2020-10","2021-01","2021-04","2021-07","2021-10","2022-01"),
                   labels=c("Apr\n2020","Jul\n2020","Oct\n2020","Jan\n2020","Apr\n2021","Jul\n2021","Oct\n2021","Jan\n2022"))+
  xlab("") + 
  ylab("") + 
  guides(fill=guide_colorbar("Residuals",title.position = "right",
                             title.theme = element_text(size = 20, angle = 90,hjust=0.5),
                             label.theme= element_text(size = 20),
                             
                             title.vjust = 0.9)) +
  theme(legend.background = element_rect(fill = NA))+
  theme(legend.key.size=unit(2,'cm'),legend.key.width=unit(0.5,'cm'))+
  
  theme(plot.margin=unit(c(0.5,3.5,0.5,0.5),'lines'),##?????????????????????
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  
  theme(axis.title.y =element_text(family = "ArialN",size=20,colour = "black"),axis.text.y=element_text(family = "ArialN",size=20,colour = "black"),
        axis.title.x =element_text(family = "ArialN",size=20,colour = "black"),axis.text.x=element_text(family = "ArialN",size=20,color = 'black'))


pdf("report/plot/report_growth_rate_fig7.pdf",height=15,width=12)
p
dev.off()
