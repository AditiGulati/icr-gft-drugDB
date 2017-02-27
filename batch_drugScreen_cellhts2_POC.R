##
## cellHTS2 analysis based on POC scores
##

pocScreen <- function(
  name,
  path, 
  posControls = "pos",
  negControls = "neg",
  descripFile = "Description.txt",
  compoundlibrary = "CompoundLibrary_p11_p12.txt",
  reportHTML = TRUE,
  plateconf = "Plateconf.txt",
  platelist = "Platelist.txt",
  pocsName="pocs.txt",
  summaryName="pocs_summary.txt",
  reportdirName="reportdir_poc",
  screenlog="Screenlog.txt" 
){

  require(cellHTS2)
  x<-readPlateList(platelist, name=name, path=path)
  #	configure
  cat("-----",name,"------ Z score\n")
  # check if we have a screenlog 
  if (file.exists( paste(path,screenlog,sep="") )){
	cat("Screenlog found.\n")
	x<-configure(x, descripFile=descripFile, confFile=plateconf, logFile=screenlog, path=path);
  } else{
	cat("No Screenlog found. Proceeding without.\n")
	x<-configure(x, descripFile=descripFile, confFile=plateconf, path=path);
  }
  
  # for poc, switch controls
  # since CELLHTS2 calculates percent of positive control
  posControls<-tolower(posControls)
  negControls<-tolower(negControls)
  xn<-normalizePlates(x, scale="multiplicate", posControls=posControls, negControls=negControls, method="negatives")

  # summarize reps
  xsc<-summarizeReplicates(xn, summary="median")

  # annotate
  xa<-cellHTS2::annotate(xsc, compoundlibrary, path=path)

  # top table
  getTopTable(list("raw"=x, "normalized"=xn, "scored"=xa), file=paste(path, summaryName, sep=""))
  setSettings(list(plateList=list(reproducibility=list(include=TRUE, map=TRUE), intensities=list(include=TRUE, map=TRUE)), screenSummary=list(scores=list(range=c(-20, 10), map=TRUE)))) 

  ###########################################################
  # write scores for this screen in the same folder
  scorefile<-paste(path, pocsName, sep="")
  plates<-plate(xa)
  wells<-well(xa)
  scores<-Data(xa)
  # prepare a simple text report
  combinedz<-data.frame(compound=geneAnno(xa), plate=plates, well=wells, pocscore=scores)
  names(combinedz)<-c("Compound", "Plate", "Well", "POCscore")
  write.table(combinedz, scorefile, sep="\t", quote=FALSE, row.names=FALSE)

  ###########################################################	
  # write a QC report for poc scores
  if(reportHTML){
    reportdir<-paste(path, reportdirName, sep="")
    try(writeReport(raw=x, normalized=xn, scored=xsc, outdir=reportdir, 
                force=TRUE, posControls= posControls, negControls= negControls, 
                mainScriptFile="/Rnaidb_git/icr-gft-drugdb/batch_drugScreen_cellhts2_POC.R"))
  }
}