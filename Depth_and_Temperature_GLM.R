library(doParallel)
registerDoParallel(cores = 16)
library(foreach)

library(glmnet)
library(fields)
library(nlstools)

####################################
### 1) OPEN AND FORMAT FILES
####################################
gridCells = read.csv("metadata_by_grid.csv",stringsAsFactors = FALSE)
GreenlandFID = c(3774:3776,3759:3762,3737:3742,3711:3715,3681:3685,3639:3643,3584:3588,3522:3525,3452:3454,3378:3380)
AntarcFID = 3791:4163 
gridCells$rechargeType = gridCells$rechargeDepth
gridCells$combined = factor(paste(gridCells$rechargeVolume,gridCells$crustScheme2))
gridCells$rechargeVolume = factor(gridCells$rechargeVolume)
gridCells$rechargeType = factor(gridCells$rechargeDepth)
gridCells$crustScheme2 = factor(gridCells$crustScheme2)
gridCells$rechargeFull = factor(gridCells$Descriptio)
gridCells$combined_cv = factor(gridCells$combined_cv)
gridCells$combined_cr = factor(gridCells$combined_cr)
gridCells$rechargeShort = factor(gridCells$rechargeShort)


df.trimmed = read.csv("cores_with_PCR.csv")
# Select direct measurements only
all = df.trimmed[which(df.trimmed$MethodCM=="direct"),]
# makes sure all parameters are in correct format
all$Depth = as.numeric(all$Depth) # in meters
all$cellsPer = as.numeric(all$cellsPer) # in cell cm-3

####################################
#### 3) SELECT TRAINING DATASET and BOOTSTRAPS
####################################
# load indices selected for bootstraps
myIndices = as.matrix(read.csv("1000_indices_for_bootstrap.csv",header=FALSE))
bootstraps = nrow(myIndices)
depthsToIterate = gridCells$Z122_Med_HF_km*1000 # in meters

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

  ## Set recharge and crust as factors
  all$combined = factor(paste(all$rechargeVolume,all$crustScheme2))
  all$rechargeType = factor(all$rechargeType)
  all$rechargeVolume = factor(all$rechargeVolume)
  all$crustScheme2 = factor(all$crustScheme2)
  all$rechargeFull = factor(all$rechargeFull)
  all$combined_cv = factor(all$combined_cv)
  all$combined_cr = factor(all$combined_cr)
  all$rechargeShort = factor(all$rechargeShort)
  
  
  ## Training set for GLM 
    keep = c("cellsPer","Depth","temperature")
  
    train <- all[train_ind,which(colnames(all)%in%keep)]
    test <- all[-train_ind,which(colnames(all)%in%keep)]
    train$Depth = log10(train$Depth)
    test$Depth = log10(test$Depth)
    train$cellsPer = log10(train$cellsPer)
    test$cellsPer = log10(test$cellsPer)
  
    train = train[which(complete.cases(train)),]
    train.x =  model.matrix(cellsPer~.,train)[,-1]
    train.y = train$cellsPer
    test = test[which(complete.cases(test)),]
    test.x =  model.matrix(cellsPer~.,test)[,-1]
    
  ## Training set for Temperature
    train.temp = all[train_ind,]
    test.temp = all[-train_ind,]
    
  ## Training set for Depth
    pf.train = all[train_ind,]
    pf.test = all[train_ind,]
####################################
## 6) RUN GLMnet FIT
####################################
    glmfit = glmnet(train.x,train.y)
    cv.fit = cv.glmnet(train.x,train.y)
    print(coef(cv.fit,cv.fit$lambda.min))
    cv.result = predict(cv.fit, newx = test.x, s = "lambda.min")
 ###################################
 ### 6b) integrate GLMnet FIT
 ################################
  
    
    results = foreach(i = 1:length(depthsToIterate), .combine = "c") %dopar% {

      mySlices = seq(1,depthsToIterate[i],0.01) # 1 cm slices from 1 m to isotherm
      patch = data.frame(gridCells[i,])
      patch = patch[rep(seq_len(nrow(patch)), each=length(mySlices)),]
      patch$Depth = as.numeric(mySlices)
      patch$Depth = as.numeric(log10(mySlices))
      if(patch$FID%in%c(GreenlandFID,AntarcFID)){patch$temperature = 0 + mySlices*(gridCells[i,]$Median_HF/1000)/gridCells[i,]$Therm_Cond__W_per_m_K}
      else{
      patch$temperature = gridCells[i,]$MEAN_Annual_Temp + mySlices*(gridCells[i,]$Median_HF/1000)/gridCells[i,]$Therm_Cond__W_per_m_K}
      patch$dummy = rep(1,nrow(patch))
      patch = patch[,which(colnames(patch)%in%append(keep,"dummy"))]
      patch = model.matrix(dummy~.,patch)
      patch = patch[,which(colnames(patch)%in%colnames(train.x))]
      
    if (length(setdiff(colnames(train.x),colnames(patch)))>0){ # control in case of descepancy of categorical factor in training/test sets
      tmp = data.frame(rep(0,nrow(patch)))
      colnames(tmp) = setdiff(colnames(train.x),colnames(patch))
      patch = cbind(patch,tmp)
    }
   sum((10^predict(cv.fit, newx = data.matrix(patch), s = "lambda.min")))
   }

  #results is a vector of of length "depthsToIterate"
  print (sum(results*gridCells$grid_area_m2)*100*100) # Total Biomass
  cv.biomass[n]=sum(results*gridCells$grid_area_m2)*100*100 
  cv.byGridResult[,n]=results*gridCells$grid_area_m2*100*100
  cv.error[n]=mean((cv.result-test$cellsPer)^2) 
  print("cv error")
  print(mean((cv.result-test$cellsPer)^2) )
}

proc.time() - ptm

write.table(cv.biomass,file = "glm_depthtempZ122_Med_HF_cv.biomass.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(cv.error,file = "glm_depthtempZ122_Med_HF_cv.error.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(cv.byGridResult, file = 'glm_depthtempZ122_Med_HF_cvGridResult.csv',sep=',',row.names=FALSE,col.names=FALSE)





