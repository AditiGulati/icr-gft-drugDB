##
## export all p11-12-13 cellHTS2 zscores 
##

z_summaryFileList <- "z_summary_file_list.txt"
poc_summaryFileList <- "poc_summary_file_list.txt"

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

  summary_with_z <- NULL
  summary_with_z <- cbind(
    z_summary, 
    zscore_rep1=rep(NA, times=nrow(z_summary)),
    zscore_rep2=rep(NA, times=nrow(z_summary)),
    zscore_rep3=rep(NA, times=nrow(z_summary))
  )

  drugs <- na.omit(unique(z_summary$GeneID))
  for (drug in drugs) {
    lib_rows <- which(z_summary$GeneID == drug)
      z_drug_mad_rep1 <- mad(z_summary$normalized_r1_ch1[lib_rows], na.rm =TRUE)
      z_drug_mad_rep2 <- mad(z_summary$normalized_r2_ch1[lib_rows], na.rm =TRUE)
      z_drug_mad_rep3 <- mad(z_summary$normalized_r3_ch1[lib_rows], na.rm =TRUE)
      summary_with_z$zscore_rep1[lib_rows] <- z_summary$normalized_r1_ch1[lib_rows]/z_drug_mad_rep1
      summary_with_z$zscore_rep2[lib_rows] <- z_summary$normalized_r2_ch1[lib_rows]/z_drug_mad_rep2
      summary_with_z$zscore_rep3[lib_rows] <- z_summary$normalized_r3_ch1[lib_rows]/z_drug_mad_rep3
  }
  z_summaryName <- gsub("_zscores_summary.txt","_zscores_summary_with_rep_zscores.txt", z_summaryFiles[z_summaryFile,], fixed=TRUE)
  write.table(
    summary_with_z,
    file=z_summaryName,
    col.names=TRUE,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )
}

poc_summaryFiles <- read.table(
  poc_summaryFileList,
  header=FALSE,
  sep="\t",
  stringsAsFactors=FALSE
)

for(poc_summaryFile in 1:nrow(poc_summaryFiles)) {
  poc_summary <- read.table(
    poc_summaryFiles[poc_summaryFile,],
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE
  )

  summary_with_poc <- NULL
  summary_with_poc <- cbind(
    poc_summary, 
    pocscore_rep1=rep(NA, times=nrow(poc_summary)),
    pocscore_rep2=rep(NA, times=nrow(poc_summary)),
    pocscore_rep3=rep(NA, times=nrow(poc_summary))
  )

  drugs <- na.omit(unique(poc_summary$GeneID))
  for (drug in drugs) {
    lib_rows <- which(poc_summary$GeneID == drug)
      poc_drug_mad_rep1 <- mad(poc_summary$normalized_r1_ch1[lib_rows], na.rm =TRUE)
      poc_drug_mad_rep2 <- mad(poc_summary$normalized_r2_ch1[lib_rows], na.rm =TRUE)
      poc_drug_mad_rep3 <- mad(poc_summary$normalized_r3_ch1[lib_rows], na.rm =TRUE)
      summary_with_poc$pocscore_rep1[lib_rows] <- poc_summary$normalized_r1_ch1[lib_rows]/poc_drug_mad_rep1
      summary_with_poc$pocscore_rep2[lib_rows] <- poc_summary$normalized_r2_ch1[lib_rows]/poc_drug_mad_rep2
      summary_with_poc$pocscore_rep3[lib_rows] <- poc_summary$normalized_r3_ch1[lib_rows]/poc_drug_mad_rep3
  }
  poc_summaryName <- gsub("_pocscores_summary.txt","_pocscores_summary_with_rep_pocscores.txt", poc_summaryFiles[poc_summaryFile,], fixed=TRUE)
  write.table(
    summary_with_poc,
    file=poc_summaryName,
    col.names=TRUE,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )
}