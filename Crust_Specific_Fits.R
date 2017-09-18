## Depth power fit that is regionally specific 

library(ggplot2)
gridCells = read.csv("metadata_by_grid.csv",stringsAsFactors = FALSE)
GreenlandFID = c(3774:3776,3759:3762,3737:3742,3711:3715,3681:3685,3639:3643,3584:3588,3522:3525,3452:3454,3378:3380)
AntarcFID = 3791:4163 # These are FID_1 on TC's documents already corrected in CM
#gridCells$rechargeFull = factor(gridCells$rechargeVolumeFull)
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

# load indices
myIndices = as.matrix(read.csv("1000_indices_for_bootstrap.csv",header=FALSE))[1:2,]
bootstraps = nrow(myIndices)
depthsToIterate = gridCells$Z122_Med_HF_km*1000 # in meters

plot(log10(all$Depth),log10(all$cellsPer))



AntarcFID = 3791:4163 # These are FID_1 on TC's documents already corrected in CM
GreenlandFID = c(3774:3776,3759:3762,3737:3742,3711:3715,3681:3685,3639:3643,3584:3588,3522:3525,3452:3454,3378:3380)



linearModel = function(var,map,train,test){
  biomass = 0
  error = matrix(nrow=0,ncol=2)
  univar = sort(unique(var))
  loopAB=matrix(nrow=1,ncol=3*length(univar))
  gridVec = vector(length=nrow(gridCells))
  for (i in 1:length(univar)){
    #print(univar[i])
    #print(lm(c~d,all[which(var==univar[i]),]))
    powerFit = lm(c~d,train[which(var==univar[i]),])
    pf.test = test[which(var==univar[i]),]
    a =  powerFit$coefficients[1]
    b = powerFit$coefficients[2]
    loopAB[,(i-1)*3+1]=a
    loopAB[,(i-1)*3+2]=b
    loopAB[,(i-1)*3+3]=summary(powerFit)$r.squared
    powerResult = predict(powerFit,pf.test)
    tempError = cbind(pf.test$c,powerResult)
    error = rbind(error,tempError)
    #      if(mean((tempError[1]-tempError[2])^2)>1){print(univar[i])}
    ###      if (is.na(b)){b=0}
    for (g in which(map==univar[i])){
      if (gridCells$FID[g]%in%GreenlandFID){      
        integralFun = function(x) {(10^7.73)*x^-0.66}
        biomass = biomass + gridCells$grid_area_m2[g]*100*100*(15000000000+integrate(integralFun,1,depthsToIterate[g])$value*100)
        gridVec[g]=gridCells$grid_area_m2[g]*100*100*(15000000000+integrate(integralFun,1,depthsToIterate[g])$value*100)     
      }
      else if (gridCells$FID[g]%in%AntarcFID) {
        integralFun = function(x) {(10^6)*x^-0.66}
        biomass = biomass + gridCells$grid_area_m2[g]*100*100*(2150000000+integrate(integralFun,1,depthsToIterate[g])$value*100)
        gridVec[g]=gridCells$grid_area_m2[g]*100*100*(2150000000+integrate(integralFun,1,depthsToIterate[g])$value*100)        
      }
      
      else{        
        integralFun = function(x) {(10^a)*x^b}
        biomass= biomass + gridCells$grid_area_m2[g]*100*100*integrate(integralFun,1,depthsToIterate[g])$value*100
        gridVec[g] = gridCells$grid_area_m2[g]*100*100*integrate(integralFun,1,depthsToIterate[g])$value*100
      }
    }
  }
  #  print(dim(error))
  mse = mean((error[1]-error[2])^2)
  biomassAndError = t(matrix(c(biomass,mse)))[1,]
  #print(biomass)
  if(mse>2){myList=list(estimate=rep(NA,2),grid=rep(NA,length=nrow(gridCells)),params=rep(NA,length(univar)*3))}
  else if (mse==0){myList=list(estimate=rep(NA,2),grid=rep(NA,length=nrow(gridCells)),params=rep(NA,length(univar)*3))}
  else{myList = list(estimate=biomassAndError,grid = gridVec,params=loopAB)}
  return(myList)
}

all$d = log10(all$Depth)
all$c = log10(all$cellsPer)

index=myIndices
newEstimate = matrix(nrow=nrow(index), ncol=2)
gridValues = data.frame(matrix(nrow=nrow(gridCells),ncol=nrow(index)))
parameters = matrix(nrow=nrow(index), ncol=3*5)
for (n in 1:nrow(index)){
  if(n%%10==0){print(n)}
  trainSet = all[index[n,],]
  testSet = all[-index[n,],]
  if(sum(table(trainSet$crustScheme2)<4)>0){print(n); newEstimate[n]=NA}
  else{
        output = linearModel(trainSet$crustScheme2,gridCells$crustScheme2,trainSet,testSet)
        newEstimate[n,]=output$estimate
        print("GLM Estimate and Error")
        print(output$estimate)
        gridValues[,n]=output$grid
        parameters[n,]=output$params
      }
  }


crustTypes = sort(unique(trainSet$crustScheme2))
median_params = apply(parameters,2,median,na.rm=TRUE)

colors = c("blue","green","red","purple","orange")
par(mfrow=c(3,2))
for (i in 1:length(crustTypes)){
  tmp = all[which(all$crustScheme2==crustTypes[i]),]
  plot(tmp$d,tmp$c,col=colors[i],xlim=c(-1,4),ylim=c(3,10),
       main=paste(paste(crustTypes[i],"; R-squared = ",sep=""),toString(round(median_params[(i-1)*3+3],digits = 2))))
  #print(summary(lm(tmp$c~tmp$d))$r.squared)
  #print(lm(tmp$c~tmp$d)$coefficients)
  tmp.lm = lm(tmp$c~tmp$d)
  lwr = predict(tmp.lm,data.frame(tmp$d), interval="predict")[,2]
  upr =  predict(tmp.lm,data.frame(tmp$d), interval="predict")[,3]
  abline(a=median_params[(i-1)*3+1],b=median_params[(i-1)*3+2],col=colors[i])
  lines(tmp$d,lwr)
  lines(tmp$d,upr)
  print(mean(upr-lwr))
}

totdat.lm = lm(all$c~all$d)
lwr = predict(totdat.lm,data.frame(all$d), interval="predict")[,2]
upr =  predict(totdat.lm,data.frame(all$d), interval="predict")[,3]
plot(all$d,all$c,main="Total Dataset; R-squared = 0.19")
abline(a=coef(totdat.lm)[1],b=coef(totdat.lm)[2],col="gray" )
lines(all$d,lwr,col="gray")
lines(all$d,upr,col="gray")
print(mean(upr-lwr))

