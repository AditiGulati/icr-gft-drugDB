##
## Functions for generating dose response curves
##

#set maximum dose
maxDose<-6

require(drc)

#### drc for a given drug
getdrc<-NULL
getdrc<-function(mydata){
  drc<- NULL
  drc<-tryCatch(
    drm(mydata$response~log10(mydata$dosage), fct=LL.3u()), logDose=10, error=function(x) NA)
  return(drc)
}

#### get drc Slope
getSlope<-NULL
getSlope<-function(drc){
  b<-NA
  b<-drc$coefficients[1]
  return(b)
}

#### get drc Einf
getEinf<-NULL
getEinf<-function(drc){
  c<-NA
  c<-drc$coefficients[2]
  return(c)
}

#### get drc EC50
getEC50<-NULL
getEC50<-function(drc){
  e<-NA
  e<-drc$coefficients[3]
  return(e)
}

#### get pvalue for the drc slope
getpval<-NULL
getpval<-function(drc){
  pval<-NA
  pval<-summary(drc)$coefficients[1,4]
  return(pval)
}

####
getSF<-NULL
getSF<-function(drc, mydata, survivalLevel){  
  sf<-NA
  if(is.na(sf)){
    sf<-maxDose
  }
  sf<-survivalFraction(drc, response=survivalLevel)
  if(min(mydata$response)>survivalLevel){
    sf<-maxDose
  }
  return(sf)  
}

####
ll3<-NULL
ll3<-function(x, drc){
  b<-drc$coefficients[1]
  c<-drc$coefficients[2]
  e<-drc$coefficients[3]
  
  lx<-log(x)
  le<-log(e)
  top<-1-c
  bott<-(1+exp(b*(lx-le)))
  c+top/bott  
}

#### calculate survival fraction
survivalFraction<-NULL
survivalFraction<-function(drc, response=0.5){  
  b<-drc$coefficients[1]
  c<-drc$coefficients[2]
  e<-drc$coefficients[3]
  val<-NA
  if(response<c){
    val<-maxDose
  }else{
  target<-NA
  target<-(1-c)/(response-c)
  val<-log(e)+1/b*log(target^(1/1)-1)
  val<-exp(val)
  }
  return(as.numeric(val))
}  

#### get auc
auc<-NULL
auc<-function(minC, maxC, drc, funct=ll3){
  ii<- NA
  ii<-tryCatch(
    integrate(funct, lower=minC, upper=maxC, drc=drc)$value, error=function(x) NA)
  
  return(ii)
}

#### get pvalue for lack of fit test
getLOF_pval<-NULL
getLOF_pval<-function(drc){
  lof<- NA
  lof<-modelFit(drc)
  lof_pval<-lof[2,5]
  return(lof_pval)
}
