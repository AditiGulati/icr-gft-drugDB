##
## cellHTS2 analysis function 
##

zScreen <- function(
  name,
  path, 
  posControls = "pos",
  negControls = "neg",
  description = "Description.txt",
  compoundlibrary = "Compoundlibrary_p11_p12.txt",
  reportHTML = TRUE,
  plateconf = "Plateconf.txt",
  platelist = "Platelist.txt",
  zscoreName="zscore.txt",
  summaryName="zscore_summary.txt",
  reportdirName="reportdir_z"
){
  
  require(cellHTS2)
  x<-readPlateList(platelist, name=name, path=path)
  x<-configure(x, description, plateconf, path=path)

  # for poc, switch controls
  # since CELLHTS2 calculates percent of positive control
  #posControls<-tolower(posControls)
  #negControls<-tolower(negControls)
  xn<-normalizePlates(x, scale="multiplicative", log=TRUE, posControls= posControls, negControls= negControls, method="median", varianceAdjust="none")

  # summarize reps
  xz<-scoreReplicates(xn, method="zscore", sign="+")
  xsc<-summarizeReplicates(xz, summary="median")
  
  # annotate
  xa<-cellHTS2::annotate(xz, compoundlibrary, path=path)
  xa_summary<-cellHTS2::annotate(xsc, compoundlibrary, path=path)
  
  # top table
  getTopTable(list("raw"=x, "normalized"=xn, "scored"=xa_summary), file=paste(path, summaryName, sep=""))
  setSettings(list(
    plateList=list(reproducibility=list(include=TRUE, map=TRUE),
                   intensities=list(include=TRUE, map=TRUE)),
    screenSummary=list(scores=list(range=c(-20, 10), map=TRUE)))) 

  ###########################################################
  # write scores for each rep for this screen in the same folder
  scorefile<-paste(path, zscoreName, sep="")
  plates<-plate(xa)
  wells<-well(xa)
  scores_rep1<-Data(xa[,1])
  scores_rep2<-Data(xa[,2])
  scores_rep3<-Data(xa[,3])
  # prepare a simple text report
  combinedz<-data.frame(compound=geneAnno(xa), plate=plates, well=wells, zscore_rep1=scores_rep1, zscore_rep2=scores_rep2, zscore_rep3=scores_rep3)
  names(combinedz)<-c("Compound", "Plate", "Well", "Zscore_rep1", "Zscore_rep2", "Zscore_rep3")
  write.table(combinedz, scorefile, sep="\t", quote=FALSE, row.names=FALSE)

    ###########################################################	
    # write a QC report for zscores
    if(reportHTML){
  
      reportdir<-paste(path, reportdirName, sep="")
      try(writeReport(raw=x, normalized=xn, scored=xsc, outdir=reportdir, 
                force=TRUE, posControls= posControls, negControls= negControls, 
                mainScriptFile="/Rnaidb_git/icr-gft-drugdb/batch_drugScreen_cellhts2_ZSCORES.R"))
    }
}