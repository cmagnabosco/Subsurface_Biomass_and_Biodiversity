## Depth and Temperature Fits

gridCells = read.csv("metadata_by_grid.csv",stringsAsFactors = FALSE)
GreenlandFID = c(3774:3776,3759:3762,3737:3742,3711:3715,3681:3685,3639:3643,3584:3588,3522:3525,3452:3454,3378:3380)
AntarcFID = 3791:4163 

df.trimmed = read.csv("cores_with_PCR.csv")
# Select direct measurements only
all = df.trimmed[which(df.trimmed$MethodCM=="direct"),]
all$Depth = as.numeric(all$Depth) # in meters
all$cellsPer = as.numeric(all$cellsPer) # in cell cm-3


# load indices
myIndices = as.matrix(read.csv("1000_indices_for_bootstrap.csv",header=FALSE))
myIndices=myIndices
bootstraps = nrow(myIndices)
depthsToIterate = gridCells$Z122_Med_HF_km*1000 # in meters


AntarcFID = 3791:4163 # These are FID_1 on TC's documents already corrected in CM
GreenlandFID = c(3774:3776,3759:3762,3737:3742,3711:3715,3681:3685,3639:3643,3584:3588,3522:3525,3452:3454,3378:3380)





####################################
## 4) DEFINE OUTPUTS
####################################
temperature.error = vector(length=bootstraps)
cv.error = vector(length=bootstraps)
lm.error = vector(length=bootstraps)

temperature.biomass = vector(length=bootstraps)
lm.biomass = vector(length=bootstraps)
cv.biomass = vector(length=bootstraps)

temperature.byGridResult = data.frame(matrix(NA, nrow =nrow(gridCells), ncol = bootstraps))
cv.byGridResult = data.frame(matrix(NA, nrow =nrow(gridCells), ncol = bootstraps))
lm.byGridResult = data.frame(matrix(NA, nrow =nrow(gridCells), ncol = bootstraps))

temperature.parameters = data.frame(matrix(NA,nrow=bootstraps,ncol=2))
lm.parameters = data.frame(matrix(NA,nrow=bootstraps,ncol=2))

####################################
### 5) BEGIN BOOTSTRAP LOOP
####################################

ptm <- proc.time()
for (n in 1:nrow(myIndices)) {
  if(n%%5==0){ 
    print(n)
  }
  ####################################
  ## 5a) Setup training and test datasets
  ####################################
  train_ind <- myIndices[n,]
  
  
  
  ## MAKE temperature vector based on maps. This model is good. when used to predict temperatures in real dataset,Adjusted R-squared = 0.91
  all$temperature = all$mast + all$Depth*(all$medianHF/1000)/all$tcon
  

  
  ## Training set for Temperature
  train.temp = all[train_ind,]
  test.temp = all[-train_ind,]
  
  ## Training set for Depth
  pf.train = all[train_ind,]
  pf.test = all[train_ind,]
  
  ## Temperature Fit
  train.temp$cellsPer = log10(train.temp$cellsPer) 
  test.temp = test.temp[which(complete.cases(test.temp$T..oC.)),]
  test.temp$cellsPer = log10(test.temp$cellsPer) 
  test.temp$T..oC. = test.temp$temperature
  
  powerFit= lm(cellsPer~T..oC.,train.temp)
  
  powerResult = predict(powerFit,test.temp)
  
  
  
  a =  powerFit$coefficients[1]
  slope = powerFit$coefficients[2]
  
  temperature.parameters[n,1] =a
  temperature.parameters[n,2] =slope
  
  biomass = 0
  temp = matrix(NA,nrow=nrow(gridCells),ncol=1)
  for (i in 1:nrow(gridCells)){
    mast = as.numeric(gridCells$MEAN_Annual_Temp[i])
    medianHF = as.numeric(gridCells$Median_HF[i])
    tcon = as.numeric(gridCells$Therm_Cond__W_per_m_K[i])
    if (gridCells$FID[i]%in%GreenlandFID){      
      integralFun = function(x) {10^(slope*(0 + (x)*(medianHF/1000)/tcon)+a)}
      temp[i]=gridCells$grid_area_m2[i]*100*100*(15000000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
      biomass = biomass + gridCells$grid_area_m2[i]*100*100*(15000000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
      
    }
    else if (gridCells$FID[i]%in%AntarcFID) {
      integralFun = function(x) {10^(slope*(0 + (x)*(medianHF/1000)/tcon)+a)}
      temp[i]=gridCells$grid_area_m2[i]*100*100*(2150000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
      biomass = biomass + gridCells$grid_area_m2[i]*100*100*(2150000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
    }
    else{
      integralFun = function(x) {10^(slope*(mast + (x)*(medianHF/1000)/tcon)+a)}
      temp[i]=gridCells$grid_area_m2[i]*100*100*integrate(integralFun,1,depthsToIterate[i])$value*100
      biomass = biomass + gridCells$grid_area_m2[i]*100*100*integrate(integralFun,1,depthsToIterate[i])$value*100
    }
  }
  
  temperature.byGridResult[,n] =temp 
  temperature.biomass[n]=biomass
  temperature.error[n] = mean((powerResult-test.temp$cellsPer)^2)
  #print("temperature error")
  #print(mean((powerResult-test.temp$cellsPer)^2))
  print("temperature biomass")
  print(biomass)
  
  ## Depth PowerFit 
  pf.train$Depth = log10(as.numeric(pf.train$Depth))
  pf.train$cellsPer = log10(as.numeric(pf.train$cellsPer))
  pf.test$Depth = log10(as.numeric(pf.test$Depth))
  pf.test$cellsPer = log10(as.numeric(pf.test$cellsPer))
  powerFit = lm(cellsPer~Depth,data=pf.train)
  #  pf.x = data.frame(test$Depth)
  powerResult = predict(powerFit,pf.test)
  a =  powerFit$coefficients[1]
  b = powerFit$coefficients[2]
  
  
  lm.parameters[n,1] =a
  lm.parameters[n,2] =b
  #### Integration step
  biomass = 0
  dep = matrix(NA,nrow=nrow(gridCells),ncol=1)
 for (i in 1:nrow(gridCells)){
    if (gridCells$FID[i]%in%GreenlandFID){
      integralFun = function(x) {(10^a)*x^b}
      biomass = biomass + gridCells$grid_area_m2[i]*100*100*(15000000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
      dep[i]=gridCells$grid_area_m2[i]*100*100*(15000000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
    }
    else if (gridCells$FID[i]%in%AntarcFID) {
      integralFun = function(x) {(10^6)*x^b}
      biomass = biomass + gridCells$grid_area_m2[i]*100*100*(2150000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
      dep[i]=gridCells$grid_area_m2[i]*100*100*(2150000000+integrate(integralFun,1,depthsToIterate[i])$value*100)
    }
    else{
      integralFun = function(x) {(10^a)*x^b}
      biomass = biomass + gridCells$grid_area_m2[i]*100*100*integrate(integralFun,1,depthsToIterate[i])$value*100
      dep[i] =  gridCells$grid_area_m2[i]*100*100*integrate(integralFun,1,depthsToIterate[i])$value*100
    }
   }
  
  lm.byGridResult[,n] = dep 
  lm.biomass[n]=biomass
  lm.error[n] = mean((powerResult-pf.test$cellsPer)^2)
  print("depth biomass")
  print(biomass)
}
proc.time() - ptm

