##
## drc based on cellHTS2 pocscores
## dose concentration in picoMoles now 
## the script's been modified to  
## log10 transform dose concentration
## for d~r analysis
##

source("/Rnaidb_git/drugdb_devel/icr-gft-drugdb/DRC_functions_ll3u.R")

library(Hmisc)

drcResult<- function(
  path="/Users/agulati/Desktop/drug_screen_batch_processing/drug_screen_output/Prostate_2014_03_31_CTG_drugscreen_22RV1_singleArm/",
  cell_line_name="22RV1",
  drcFile="Prostate_2014_03_31_CTG_drugscreen_22RV1_singleArm_pocscores_drugresponse.txt",
  summaryName="Prostate_2014_03_31_CTG_drugscreen_22RV1_singleArm_pocscores_summary.txt",
  drcFigures="Prostate_2014_03_31_CTG_drugscreen_22RV1_singleArm_pocscore_figures.pdf"
){

  maxDose=6
  
  # do we need a closing /?
  last<- substr(path, nchar(path), nchar(path))
  if(last!="/")
  {
    path<-paste(path,"/",sep="")
  }
  require(drc)
  reportName<-paste(path, drcFile, sep="")
  scores<-read.delim(paste(path, summaryName, sep=""))
  scores$GeneID<-toupper(scores$GeneID)
  mydrugs<-NULL
  mydrugs<-unique(scores$GeneID)
  mydrugs<-na.omit(mydrugs)
  mydrugs<-sort(mydrugs)
  mydrugs <- mydrugs[ mydrugs != "EMPTY" & mydrugs != "POS" & mydrugs != "NEG" & mydrugs != "empty" & mydrugs != "pos" & mydrugs != "neg" & mydrugs != "DMSO" & mydrugs != "dmso" ]
  
  report=data.frame(
    drug=character(), 
    layout=character(), 
    sf80=numeric(), 
    sf50=numeric(), 
    sf20=numeric(),
    auc100=numeric(),
    slope=numeric(), 
    eInf=numeric(), 
    slope_pval=numeric(), 
    ec50=numeric(),
    LOF_pval=numeric(),
    auc100_actual=numeric()
  )
  scores_sem <- NULL
  for(i in 1:nrow(scores)){
    scores_sem[i] <- sd(scores[i,c("normalized_r1_ch1","normalized_r2_ch1","normalized_r3_ch1")])/(3^0.5)
  }
  scores_with_sem <- data.frame(
    scores,
    sem=scores_sem
  )
  scores_with_sem$GeneID<-toupper(scores_with_sem$GeneID)
    
  figs<-paste(path, drcFigures, sep="")
  pdf(figs, width=5, height=5)
  
  for(drug in mydrugs){
    drug_rows <- which(scores_with_sem$GeneID == drug)
    
    mean_score <- apply(scores_with_sem[drug_rows,c("normalized_r1_ch1","normalized_r2_ch1","normalized_r3_ch1")], 1, mean)
    score_conc <- log10(scores_with_sem[drug_rows,"Concentration_pM"]) 
  
    yplus <- mean_score + scores_with_sem[drug_rows,"sem"]
    yminus <- mean_score - scores_with_sem[drug_rows,"sem"]
    pl<-which(scores_with_sem$GeneID==drug)[1]
    playout<-scores_with_sem$Layout[pl]
    # scored data 
    mydata2<-subset(scores_with_sem, subset=(GeneID==drug))
    mydata<-data.frame(response=c(mydata2$normalized_r1_ch1, mydata2$normalized_r2_ch1, mydata2$normalized_r3_ch1), dosage=rep(mydata2$Concentration_pM, 3))
    mydata$dosage<-as.numeric(as.character(mydata$dosage))
    drc<-getdrc(mydata)
  
  #
  # NB - plot.drc used the mean of the replicate responses per dose...
  #
  print(mydata$response[1])
  
  # may want to calculate the mean/median value of the lowest dose value
  # mean(mydata$response[which(mydata$dosage == min(mydata$dosage))])
  
  # check if we have a model - if not fill in NAs
  if(is.na(drc))
  {
    temp<-data.frame(
      drug=drug,                      
      layout=playout, 
      sf80=NA, 
      sf50=NA, 
      sf20=NA, 
      auc100=NA, 
      slope=NA, 
      slope_pval=NA, 
      eInf=NA, 
      ec50=NA,
      LOF_pval=NA,
      auc100_actual=NA
    )    
    report=rbind(report, temp)
	plot(
	  mydata$response~log10(mydata$dosage),
	  main=paste(drug, " response in ", cell_line_name),
	  xlab="Dose (log10 pM)",
	  ylab="Surviving fraction",
	  ylim=c(0.0, 1.50),
	  lwd=2
    )
    print("skipping setting SF and AUC values - no model?")
  }else{  
  
  # set sf variables to NA before assignment
  # otherwise might get a stale value
  sf80 <- NA
  sf50 <- NA
  sf20 <- NA
  auc100 <- NA
  auc100_actual <- NA
  
  slope <- NA
  einf <- NA
  ec50 <- NA
  ps <- NA
  
  slope<-getSlope(drc)
  slope <- round(slope, digits=3)
  einf<-getEinf(drc)
  einf <- round(einf, digits=3)
  ec50<-getEC50(drc)
  ec50 <- round(ec50, digits=3) 
  ps<-getpval(drc)
  ps <- round(ps, digits=3)
  
  # try setting ps to 1 if ps is.na...
  if(is.na(ps)){
    ps <- 1
  }
  #if(ps > 0.1 | slope <= 0){
  if(slope <= 0){
    sf80 <- maxDose
    sf50 <- maxDose
    sf20 <- maxDose
   }else{
     sf80<-getSF(drc, mydata=mydata, survivalLevel=0.8)
     sf80 <- round(sf80, digits=3)
     sf50<-getSF(drc, mydata=mydata, survivalLevel=0.5)
     sf50 <- round(sf50, digits=3)
     sf20<-getSF(drc, mydata=mydata, survivalLevel=0.2)
     sf20 <- round(sf20, digits=3)  
    
    # need to set this to a big value?
    # if min response is 0.5 then we never saw 
    # a concentration that killed 50% of the cells
    # Set sf50/20 to max dose (10^3?)
    # values were 100e-6 or 100e-8
    # changed to 10^3
    
    if(min(mydata$response)>0.5){
      sf50<-maxDose
      sf20<-maxDose
    }
    
    if(min(mydata$response)>0.8){
      sf80<-maxDose
    }  
  }
  
  # get per drug per screen auc 
  auc100_actual=auc(min(log10(mydata$dosage)), max(log10(mydata$dosage)), drc)
  auc100_actual <- round(auc100_actual, digits=3)
  #auc100=auc(min(mydata$dosage), max(mydata$dosage), drc)
  
   # change auc scale to 0-1 and cap auc values > 1 - max log10 auc value for no response (0 slope) = 3.31
  if(is.na(auc100_actual)){
    auc100 = NA
  }else{
    auc100=auc100_actual/3.31
  }
  if(auc100 > 1){
    auc100 = 1    
  }
  auc100 <- round(auc100, digits=3)
  
  # get per drug per screen lack of fit test p-value
  LOF_pval<-getLOF_pval(drc)
  LOF_pval <- round(LOF_pval, digits=3)
  
  plot(
    drc,
    type="all",
    main=paste(drug, " response in ", cell_line_name),
    xlab="Dose (log10 pM)",
    #xlab="Dose (nM)",
    ylab="Surviving fraction",
    ylim=c(0.0, 1.50),
    #cex.lab=1.4,
    #cex.axis=1.4,
    #cex.main=1.4,
    lwd=2
  )
  errbar(score_conc, mean_score, yplus, yminus, col="black", errbar.col="black", pch=19, add=TRUE)
  
    if((sf80 < maxDose) | (!is.null(drc))) {
      lines(c(sf80,sf80), c(-1, 0.8), col="red")
      text(sf80, 0.86, labels="SF80", col="red")
    }   
    if((sf50 < maxDose) | (!is.null(drc))) {
      lines(c(sf50,sf50), c(-1, 0.5), col="red")
      text(sf50, 0.56, labels="SF50", col="red")
    }
    if((sf20 < maxDose) | (!is.null(drc))) {
      lines(c(sf20,sf20), c(-1, 0.2), col="red")
      text(sf20, 0.26, labels="SF20", col="red")
    }   
    temp<-data.frame(
      drug=drug, 
      layout=playout, 
      sf80=sf80, 
      sf50=sf50, 
      sf20=sf20, 
      auc100=auc100, 
      slope=slope, 
      slope_pval=ps,
      eInf=einf, 
      ec50=ec50,
      LOF_pval=LOF_pval,
      auc100_actual=auc100_actual
    )   
    report=rbind(report, temp)
  }
}  
  dev.off()

  write.table(report, file = reportName, sep = "\t", row.names=FALSE)
}