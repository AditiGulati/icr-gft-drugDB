##
## Drug effect analysis script
##

inputFiles = list.files()

file <- grep("summaryFile1_", inputFiles, ignore.case=T)

for (i in 1:length(inputFiles)) {
  if (i == file) {
    dmso_File_Path = inputFiles[i]
  }else{
    drug_File_Path = inputFiles[i]
  }
}

drug_effect <- NULL

# read in the dmso file
dmso <- read.table(
  dmso_File_Path,
  head=T,
  sep="\t",
  check.names=FALSE,
  stringsAsFactors=FALSE
)
dmso <- dmso[ order(dmso[,1],dmso[,4]), ]

# read in the dmso file
drug <- read.table(
  drug_File_Path,
  head=T,
  sep="\t",
  check.names=FALSE,
  stringsAsFactors=FALSE
)
drug <- drug[ order(drug[,1],drug[,4]), ]

dmso$wellAnno <-  toupper(dmso$wellAnno)

# populate the output file
value <- NULL
Sample <- NULL
for (value in 1:length(dmso$wellAnno)) {
  sample <- NULL
  if (
    dmso$wellAnno[value] == "SAMPLE"
    #dmso$wellAnno[value] == "sample" | dmso$wellAnno[value] == "SAMPLE"
  ) {
    sample[value] = "SAMPLE"
  }
  if (
    dmso$wellAnno[value] == "NEG" |
    dmso$wellAnno[value] == "POS" 
  ) {    
    sample[value] = "CONTROL" 
  }
  if (
    dmso$wellAnno[value] == "EMPTY" |
    dmso$wellAnno[value] == "COMPOUND" |
    dmso$wellAnno[value] == "STAURASPORINEDRUG" |
    dmso$wellAnno[value] == "VEHICLEDRUG" |
    dmso$wellAnno[value] == "DRUG.DMSO" |
    dmso$wellAnno[value] == "DRUG.STAURO"
  ) {    
    sample[value] = "EMPTY"
  }
  Sample <- rbind(Sample, sample[value])
}

value <- NULL
for (value in 1:length(dmso$wellAnno)) {
  if (
    dmso$wellAnno[value] == "NEG" |
    dmso$wellAnno[value] == "POS" 
  ) {  
    dmso$GeneID[value] = dmso$wellAnno[value]
  }
}
value <- NULL
for (value in 1:length(dmso$wellAnno)) {  
  if (
    dmso$wellAnno[value] == "NEG" 
  ) {    
    dmso$Entrez_gene_ID[value] = "NEG"
  }  
  if (
    dmso$wellAnno[value] == "POS" 
  ) {    
    dmso$Entrez_gene_ID[value] = "POS"
  }
  if (
    dmso$wellAnno[value] == "EMPTY"
  ) {    
    dmso$Entrez_gene_ID[value] = "EMPTY"
  }
}

# add columns 1 to 11 
drug_effect <- data.frame(
  dmso[ ,"plate"], 
  dmso[ ,"well"], 
  dmso$GeneID,
  dmso$Function,
  Sample, 
  dmso[ ,"raw_r1_ch1"], 
  dmso[ ,"raw_r2_ch1"], 
  dmso[ ,"raw_r3_ch1"], 
  drug[ ,"raw_r1_ch1"], 
  drug[ ,"raw_r2_ch1"], 
  drug[ ,"raw_r3_ch1"],
  stringsAsFactors=FALSE
)

colnames(drug_effect) <- cbind(
  "Plate", 
  "Well",
  "Mature_Sanger_ID",
  "Function",
  "Sample", 
  "dmso_r1_ch1", 
  "dmso_r2_ch1", 
  "dmso_r3_ch1", 
  "drug_r1_ch1", 
  "drug_r2_ch1", 
  "drug_r3_ch1"
)

##
## log2 transform raw scores 
## for all reps in DMSO and 
## DRUG screens
##

i <- NULL 
j <- NULL
dmso_r1_log2 <- NULL
dmso_r2_log2 <- NULL
dmso_r3_log2 <- NULL
drug_r1_log2 <- NULL
drug_r2_log2 <- NULL
drug_r3_log2 <- NULL
for (j in 6:11) {
  for (i in 1:nrow(drug_effect)) {  
    dmso_r1_log2[i] <- log2(drug_effect[i,6])
    dmso_r2_log2[i] <- log2(drug_effect[i,7])
    dmso_r3_log2[i] <- log2(drug_effect[i,8])
    drug_r1_log2[i] <- log2(drug_effect[i,9])
    drug_r2_log2[i] <- log2(drug_effect[i,10])
    drug_r3_log2[i] <- log2(drug_effect[i,11])
  }
}

## add columns 12-17 to the drug_effect file ##
drug_effect <- cbind(
  drug_effect,                      
  dmso_r1_log2, 
  dmso_r2_log2, 
  dmso_r3_log2, 
  drug_r1_log2, 
  drug_r2_log2, 
  drug_r3_log2
)
drug_effect <- data.frame(drug_effect)

##
## 1. calculate per plate median for each log2 rep in
##    DMSO and DRUG screens excluding controls so that 
##    the median is not skewed by extreme effects
##
## 2. Add rows with plate centered values
##    for each of DMSO and DRUG reps
##    to the drug_effect file
##    

plates <- NULL
plates <- levels(as.factor(drug_effect$Plate))
plate <- NULL
plate_sample_rows <- NULL
plate_control_rows <- NULL
plate_sample_rows_for_controls <- NULL
pp_sample_median <- NULL

# calculate plate centered log2 values for samples only

plate_centered_log_intensities <- matrix(data=rep(NA, times=6*length(drug_effect$Sample)), ncol=6, nrow=length(drug_effect$Sample)) # assume we have 2*3 reps
column <- NULL
log2_columns <- NULL
log2_columns <- cbind(drug_effect[,12:17])
for(column in 1:ncol(log2_columns)) {
  for (plate in plates) {
    plate_control_rows <- which(
      drug_effect$Plate == plate & 
        drug_effect$Sample != "SAMPLE" &
        drug_effect$Sample != "EMPTY"
    )
    plate_sample_rows_for_controls <- which(
      drug_effect$Plate == plate & 
      drug_effect$Sample == "SAMPLE"
    )
    pp_control_median <- median(log2_columns[plate_sample_rows_for_controls, column],na.rm=TRUE) 
    row <- NULL
    for (row in plate_control_rows) {
      plate_centered_log_intensities[row,column] <- (log2_columns[row,column] - pp_control_median)
    }
      plate_sample_rows <- which(
        drug_effect$Plate == plate & 
          drug_effect$Sample == "SAMPLE"
      )
      pp_sample_median <- median(log2_columns[plate_sample_rows, column],na.rm=TRUE) 
      row <- NULL
      for (row in plate_sample_rows) {
        plate_centered_log_intensities[row,column] <- (log2_columns[row,column] - pp_sample_median)
      }
   }
}
#add columns 18-23 to the drug_effect file
drug_effect <- cbind(drug_effect, plate_centered_log_intensities)

names(drug_effect)[18] <- paste("dmso_r1_log2_plate_centered") 
names(drug_effect)[19] <- paste("dmso_r2_log2_plate_centered") 
names(drug_effect)[20] <- paste("dmso_r3_log2_plate_centered") 
names(drug_effect)[21] <- paste("drug_r1_log2_plate_centered") 
names(drug_effect)[22] <- paste("drug_r2_log2_plate_centered") 
names(drug_effect)[23] <- paste("drug_r3_log2_plate_centered") 


##
## median of plate-centered data
##

dmso_pcentered_median <- NULL
drug_pcentered_median <- NULL
de_pcentered_median <- NULL

columns18_20 <- cbind(drug_effect[,18:20])
columns18_20 <- as.matrix(columns18_20)
row <- NULL
for (row in 1:nrow(columns18_20)) {
  drug_pcentered_med <- NULL
  dmso_pcentered_med <- median(columns18_20[row,], na.rm=TRUE)
  dmso_pcentered_median <- rbind(dmso_pcentered_median, dmso_pcentered_med)
  colnames(dmso_pcentered_median) <- "dmso_median_log2_plate_centered"
}

columns21_23 <- cbind(drug_effect[,21:23])
columns21_23 <- as.matrix(columns21_23)
row <- NULL
for (row in 1:nrow(columns21_23)) {
  drug_pcentered_med <- NULL
  drug_pcentered_med <- median(columns21_23[row,], na.rm=TRUE)
  drug_pcentered_median <- rbind(drug_pcentered_median, drug_pcentered_med)
  colnames(drug_pcentered_median) <- "drug_median_log2_plate_centered"
}
#add column 24 & 25 to the drug_effect file
drug_effect <- cbind(drug_effect, dmso_pcentered_median, drug_pcentered_median)



##
## Calculate DRUG EFFECT (DE) 
## (drug_log2_plate_centered_score - dmso_log2_plate_centered_score)
##

median_pcentered_scores <- NULL
median_pcentered_scores <- cbind(drug_effect[24], drug_effect[25])
row <- NULL
DE_score <- NULL
drug_value <- NULL
dmso_value <- NULL
for (row in 1:nrow(median_pcentered_scores)){
  de_value <- NULL
  drug_value <- median_pcentered_scores[row,2]
  dmso_value <- median_pcentered_scores[row,1]
  de_value <- (drug_value) - (dmso_value)
  DE_score <- rbind(DE_score, de_value)
  nrow(DE_score)
  colnames(DE_score) <- "DE" 
}
## add column 26 to the drug_effect file ##
drug_effect <- cbind(drug_effect, DE_score)

##
## Calculate zscores  
## for each rep
##

column <- NULL
zscore_columns <- matrix(data=rep(NA, times=3*length(drug_effect$Sample)), ncol=3, nrow=length(drug_effect$Sample))
plate_centered_columns <- cbind(drug_effect[,24:26])

#control_MEDIAN <- NULL
control_MAD <- NULL

#sample_MEDIAN <- NULL
sample_MAD <- NULL

sample_rows_for_control <- NULL

for (column in 1:ncol(plate_centered_columns)) {
  control_rows <- which(
    drug_effect$Sample != "SAMPLE" & 
      drug_effect$Sample != "EMPTY"
  )
  sample_rows_for_control <- which(
    drug_effect$Sample == "SAMPLE"
  )
  #control_MEDIAN <- median(plate_centered_columns[sample_rows_for_control,column], na.rm = TRUE)
  #print(control_MEDIAN)
  control_MAD <- mad(plate_centered_columns[sample_rows_for_control,column], na.rm = TRUE)
  i <- NULL
  for (i in control_rows) {
    zscore_columns[i,column] <- (plate_centered_columns[i,column])/control_MAD
  }
    sample_rows <- NULL
    zscore_sublib_column <- NULL
    sample_rows <- which(
      drug_effect$Sample == "SAMPLE"
      )
    #sample_MEDIAN <- median(plate_centered_columns[sample_rows,column], na.rm = TRUE)
    #print(sample_MEDIAN)
    sample_MAD <- mad(plate_centered_columns[sample_rows,column], na.rm = TRUE)
    i <- NULL
    for (i in sample_rows) {
      zscore_columns[i,column] <- (plate_centered_columns[i,column])/sample_MAD
    }
}
#add columns 27_29 to the drug_effect file
drug_effect <- cbind(drug_effect, zscore_columns)

names(drug_effect)[27] <- paste("dmso_Zscore") 
names(drug_effect)[28] <- paste("drug_Zscore")
names(drug_effect)[29] <- paste("DE_Zscore")  

##
## round off the values in drug effect file to 3 decimal places
##

#drug_effect$dmso_median_zscore <- round(drug_effect$dmso_median_Zscore, digits=3)
#drug_effect$drug_median_zscore <- round(drug_effect$drug_median_Zscore, digits=3)
#drug_effect$DE_median_zscore <- round(drug_effect$DE_median_Zscore, digits=3)

##
## write drug_effect results to a text file
##

short_drug = gsub("_zscores_summary_with_rep_zscores.txt","",drug_File_Path)
shorter_drug = gsub("summaryFile2_","",short_drug)
short_dmso = gsub("_zscores_summary_with_rep_zscores.txt","",dmso_File_Path)
shorter_dmso = gsub("summaryFile1_","",short_dmso)

screen_DE_dir_name = paste("drugEffect_", shorter_dmso,"_and_",shorter_drug, sep="")

write.table(
  drug_effect, 
  file=paste(screen_DE_dir_name,".txt", sep=""), 
  row.names=FALSE,
  sep="\t"
)
