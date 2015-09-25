# This script reads the guide file and 
# analyses the drug screens by calling 
# drug screen analysis pipeline

source("/Rnaidb_git/icr-gft-drugdb/xls2platefiles.R")
source("/Rnaidb_git/icr-gft-drugdb/batch_drugScreen_cellhts2_ZSCORES.R")
source("/Rnaidb_git/icr-gft-drugdb/reps_qc.R")
source("/Rnaidb_git/icr-gft-drugdb/batch_drugScreen_cellhts2_POC.R")
source("/Rnaidb_git/icr-gft-drugdb/batch_drugResponseCurve_POC_ll3u.R")

guideFile <- NULL
guideFile <- list.files(pattern="^.*guide_file.*.txt$")
print (guideFile)

guide <- read.table(
  guideFile, 
  header=TRUE,
  sep = "\t", 
  stringsAsFactors=FALSE
)

Xls=guide$xls_file
Datapath=guide$Datapath
Cell_line=guide$Cell_line
#Negcontrols = "(?i)^neg$|^neg2$|^siCON1$|^siCON2$|^allstar$"
#Poscontrols = "(?i)^pos$|^siPLK1$|plk1"
Annotationfile=guide$Compound_library_file
DescripFile=guide$Descrip_file		
Platelist=guide$Platelist_file
Plateconf=guide$Plateconf_file
ReportHTML1=guide$report_html1
ReportHTML2=guide$report_html2
ZscoreName=guide$zscore_file
PocscoreName=guide$pocscore_file
SummaryName1=guide$summary_file1
SummaryName2=guide$summary_file2
ZprimeName=guide$zprime_file
ReportdirName1=guide$reportdir_file1
ReportdirName2=guide$reportdir_file2
Controls_qc=guide$qc_file
Qc_plot_1=guide$plot_1_file
Qc_plot_2=guide$plot_2_file
Qc_plot_3=guide$plot_3_file
Corr_coeff=guide$corr_file
separateZprimeFile=guide$separate_zprime_file
Drc_file=guide$drc_file
Drc_figures=guide$drc_figures

#
## run xls2platelist.R script
#
  
script_1 <- xls2platelist(
  datapath=Datapath,
  xls=Xls
)   

#
## run cellHTS2 analysis for zscores
#

  script_2 <- zScreen(  
    name=Cell_line,
    path=paste(Datapath,"/",sep=""),
    #path=Datapath,
    posControls="pos",
    negControls="neg",
    description=DescripFile,
    compoundlibrary=Annotationfile,
    reportHTML=ReportHTML1,
    plateconf=Plateconf,
    platelist=Platelist,
    zscoreName=ZscoreName, 
    summaryName=SummaryName1,
    reportdirName=ReportdirName1   
  )

#
## run cellHTS2 analysis for POC scores
#

  script_3 <- pocScreen(
    name=Cell_line,
    path=paste(Datapath,"/",sep=""),
    posControls="pos",
    negControls="neg",
    description=DescripFile,
    compoundlibrary=Annotationfile,
    reportHTML=ReportHTML2,
    plateconf=Plateconf,
    platelist=Platelist,
    pocsName=PocscoreName, 
    summaryName=SummaryName2,
    reportdirName=ReportdirName2   
  )
  
#
## function for qc - creating boxplots and calculating correlation coefficient
#
  
  script_4 <- Zscreen_qc(
    datapath=Datapath,
    summaryName=SummaryName1,
    qc_plot_1=Qc_plot_1,
    qc_plot_2=Qc_plot_2,
    qc_plot_3=Qc_plot_3,
    corr_coeff=Corr_coeff
 )

#
## run drug response curve analysis
#

  script_5 <- drcResult(
    path=Datapath,
    cell_line_name=Cell_line,
    drcFile=Drc_file,
    summaryName=SummaryName2,
    drcFigures=Drc_figures 
  )