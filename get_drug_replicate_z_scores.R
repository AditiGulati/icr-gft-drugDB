##
## export all p11-12-13 cellHTS2 zscores 
##

z_summaryFileList <- "z_summary_file_list.txt"

z_summaryFiles <- read.table(
  z_summaryFileList,
  header=FALSE,
  sep="\t",
  stringsAsFactors=FALSE
)

for(z_summaryFile in 1:nrow(z_summaryFiles)) {
  z_summary <- read.table(
    z_summaryFiles[z_summaryFile,],
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE
  )
  z_summary <- z_summary[ order(z_summary[,"plate"],z_summary[,"well"]), ]
  
  colnames(z_summary)[22] <- "Concentration_pM"
  compound <- paste(z_summary$GeneID, "_", z_summary$Concentration_pM, sep="")
  compound <- toupper(compound)
  z_summary <- cbind(z_summary, compound)
  drugs <- as.character(z_summary$compound) 
  summary_with_z <- NULL
  summary_with_z <- cbind(
    z_summary, 
    zscore_rep1=rep(NA, times=nrow(z_summary)),
    zscore_rep2=rep(NA, times=nrow(z_summary)),
    zscore_rep3=rep(NA, times=nrow(z_summary))
  )

  drug_rows <- which(! drugs %in% c("EMPTY_NA", "NA_NA", "POS_NA", "NEG_NA", "NEG2_NA"))

  z_drug_mad_rep1 <- mad(z_summary$normalized_r1_ch1[drug_rows], na.rm =TRUE)
  z_drug_mad_rep2 <- mad(z_summary$normalized_r2_ch1[drug_rows], na.rm =TRUE)
  z_drug_mad_rep3 <- mad(z_summary$normalized_r3_ch1[drug_rows], na.rm =TRUE)
  
  for (i in 1:nrow(summary_with_z)) {
    summary_with_z$zscore_rep1[i] <- z_summary$normalized_r1_ch1[i]/z_drug_mad_rep1
    summary_with_z$zscore_rep2[i] <- z_summary$normalized_r2_ch1[i]/z_drug_mad_rep2
    summary_with_z$zscore_rep3[i] <- z_summary$normalized_r3_ch1[i]/z_drug_mad_rep3
  }
  z_summaryName <- gsub(".txt","_with_rep_zscores.txt", z_summaryFiles[z_summaryFile,], fixed=TRUE)
  write.table(
    summary_with_z,
    file=z_summaryName,
    col.names=TRUE,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )
}