#data$country <- as.factor(data$country)
#data$travel2 <- log(data$travel_route)
data<- data[order(data$country,data$month,decreasing=F),]

# Set maximum lag
nlag = 3
#rownames(data) <- data$record

colnames(data)
# Setting the lag for public health and social measure (PHSM) #
lag_npi <- tsModel::Lag(data$NPI_NV, group = data$country_index, k = 0:nlag)  

# Setting the lag for Shannon's index of lineage diversity based on PANGOLIN#      
lag_lineage <- tsModel::Lag(data$sub_lineage, group = data$country_index, k = 0:nlag) 

# Setting the lag for international travel#  
lag_travel <- tsModel::Lag(data$travel2, group = data$country_index, k = 0:nlag)     

# Setting the lag for the growth rate of cases#  
lag_gr <- tsModel::Lag(data$GR_NCM_log_diff, group = data$country_index, k = 0:nlag) 

# Setting the lag for adjusted vaccine coverage#  
lag_FBW_vero <- tsModel::Lag(data$FW_BW_P, group = data$country_index, k = 0:nlag)

# Setting the lag for the richness of lineage diversity# 
lag_richness <- tsModel::Lag(data$sub_richness, group = data$country_index, k = 0:nlag) 

# Setting the lag for the evenness of lineage diversity# 
lag_evenness <- tsModel::Lag(data$sub_evenness, group = data$country_index, k = 0:nlag) 

# Setting the lag for the simpson' index of lineage diversity# 
lag_simpson <- tsModel::Lag(data$simpson_lineage, group = data$country, k = 0:nlag) 



# Define cross basis matrix (combining nonlinear exposure and hysteresis function)
# Set the lag section

#Construct a cross-basis of the growth rate of cases and PHSM
var <- lag_npi
npi_gr_cb <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns",df=2))
summary(npi_gr_cb)


#Construct a cross-basis of lineage diversity and PHSM
var <- lag_npi
npi_si_cb <- crossbasis(var,
                        argvar = list(fun = "lin"),
                        arglag = list(fun = "ns",df=2))
summary(npi_si_cb)


#Construct a cross-basis of the growth rate of cases and adjusted vaccine coverage
var <- lag_FBW_vero
FBW_gr_cb <- crossbasis(var,
                        argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns",df=2))
summary(FBW_gr_cb)


#Construct a cross-basis of the lineage diversity and adjusted vaccine coverage
var <- lag_FBW_vero
FBW_si_cb <- crossbasis(var,
                        argvar = list(fun = "poly",degree=3),
                        arglag = list(fun = "ns",df=2))
summary(FBW_si_cb)


#Construct a cross-basis of the growth rate of cases and adjusted vaccine coverage
var <- lag_travel
travel_gr_cb <- crossbasis(var,
                         argvar = list(fun = "lin"),
                         arglag = list(fun = "ns",df=2))
summary(travel_gr_cb)

#Construct a cross-basis of the lineage diversity of cases and adjusted vaccine coverage

var <- lag_travel
travel_si_cb <- crossbasis(var,
                           argvar = list(fun = "lin"),
                           arglag = list(fun = "ns",df=2))
summary(travel_si_cb)


#Construct a cross-basis of the Shannon's index of lineage diversity and the growth rate
var <- lag_lineage
lineage_gr_cb <- crossbasis(var,
                            argvar = list(fun = "ns",
                                          knots = equalknots(data$sub_lineage, 2)),
                            arglag = list(fun = "ns",df=2))
summary(lineage_gr_cb)


#Construct a cross-basis of the richness of lineage diversity and the growth rate
var <- lag_richness
richness_gr_cb <- crossbasis(var,
                             argvar = list(fun = "ns",
                                           knots = equalknots(data$sub_richness, 2)),
                             arglag = list(fun = "ns",df=2))
summary(richness_gr_cb)


#Construct a cross-basis of the evenness of lineage diversity and the growth rate
var <- lag_evenness
evenness_gr_cb <- crossbasis(var,
                             argvar = list(fun = "ns",
                                           knots = equalknots(data$sub_evenness, 2)),
                             arglag = list(fun = "ns",df=2))
summary(evenness_gr_cb)


#Construct a cross-basis of the Simpson' index of lineage diversity and the growth rate
var <- lag_simpson
simpson_gr_cb <- crossbasis(var,
                            argvar = list(fun = "ns",
                                          knots = equalknots(data$simpson_lineage, 2)),
                            arglag = list(fun = "ns",df=2))
summary(simpson_gr_cb)


#Construct a cross-basis of the growth rate based on comfirmed cases and lineage diversity
var <- lag_gr
gr_cb <- crossbasis(var,
                    argvar = list(fun = "poly",degree =2),
                    arglag = list(fun = "ns",df=2))
summary(gr_cb)






# Specifies a unique column name for the gam() model
# Note: GLM (), GAM () or GLM.nb () models are not required
#colnames(npi_cb) = paste0("npi_cb.", colnames(npi_cb))
#colnames(vero_cb) = paste0("vero_cb.", colnames(vero_cb))
#colnames(T_vero_cb) = paste0("T_vero_cb.", colnames(T_vero_cb))
#colnames(F_vero_cb) = paste0("F_vero_cb.", colnames(F_vero_cb))
#colnames(shannon_cb) = paste0("shannon_cb.", colnames(npi_cb))



