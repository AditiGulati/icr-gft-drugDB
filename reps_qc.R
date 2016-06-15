Zscreen_qc <- function(
  datapath, 
  summaryName, 
 # controls_qc, 
  qc_plot_1, 
  qc_plot_2, 
  qc_plot_3, 
  corr_coeff
) {

  datapath=paste(datapath,"/", sep="")

  summary <- paste(datapath, summaryName, sep="")

  summary_file <- read.table(file= summary, sep="\t", header=TRUE) 				

  #sample_summary_data<-subset(summary_file, subset=(Function!=""))
  sample_summary_data<-summary_file[which(summary_file$wellAnno == "sample"),]
  rep1<-sample_summary_data$normalized_r1_ch1
  rep2<-sample_summary_data$normalized_r2_ch1
  rep3<-sample_summary_data$normalized_r3_ch1


  # debug
  write.table(rep1, file="rep1_test.txt")

  r12.cor <- cor(rep1,rep2,method="pearson",use="pairwise.complete.obs")
 #r12.cor <- cor(sample_summary_data$normalized_r1_ch1,sample_summary_data$normalized_r2_ch1,method="pearson",use="pairwise.complete.obs")

  r13.cor <- cor(rep1,rep3,method="pearson",use="pairwise.complete.obs")
 #r13.cor <- cor(sample_summary_data$normalized_r1_ch1,sample_summary_data$normalized_r3_ch1,method="pearson",use="pairwise.complete.obs")

  r23.cor <- cor(rep2,rep3,method="pearson",use="pairwise.complete.obs")
 #r23.cor <- cor(sample_summary_data$normalized_r2_ch1,sample_summary_data$normalized_r3_ch1,method="pearson",use="pairwise.complete.obs")

  min.cor <- min(r12.cor, r13.cor, r23.cor)

  max.cor <- max(r12.cor, r13.cor, r23.cor)


  ## write correlation coefficient values in a text file

  combined.cor<-data.frame(min.cor, max.cor)

  names(combined.cor)<-c("Minimum_correlation", "Maximum_correlation")

  corr_coeff_file <- paste(datapath, corr_coeff, sep="")

  write.table(combined.cor, corr_coeff, sep="\t", quote=FALSE, row.names=FALSE)
	
  ##
  ## scatter plots
  ##
  	
  numberOfPlates = levels(as.factor(sample_summary_data$plate))

  plate_cols <- NULL

  for (plate in numberOfPlates) {
    plate_cols = rainbow(length(numberOfPlates), s=1, v=1, alpha=0.5)
    names(plate_cols) <- numberOfPlates
  }
 	
  ## scatter plot showing correlation between rep1 and rep2

  qc_plot_file_1 <- paste(datapath, qc_plot_1 , sep="")

  png(file=qc_plot_file_1, width=400, height=400)
  
  plot(
    NULL,
    NULL,
    xlim=c(min(rep1,na.rm=TRUE),max(rep1,na.rm=TRUE)),
    ylim=c(min(rep2,na.rm=TRUE),max(rep2,na.rm=TRUE)),
    main=paste("Rep1 vs Rep2 (r=",round(r12.cor, 2),")", sep=""),
    xlab="Normalised scores for Rep1",
    ylab="Normalised scores for Rep2",
    cex=1.5
  )
  
  for (plate in numberOfPlates) {
    points(
      rep1[which(sample_summary_data$plate == plate)], 
      rep2[which(sample_summary_data$plate == plate)],
      cex=1.5,
      pch=19,
      col=plate_cols[which(names(plate_cols) == plate)]
    )	
  }
  legend( 
    x="bottomright", 
    legend=paste("Plate",names(plate_cols), sep=""),
    col=plate_cols, lwd=2, lty=c(1,2), 
    pch=19
  )
  dev.off()

  ## scatter plot showing correlation between rep2 and rep3

  qc_plot_file_2 <- paste(datapath, qc_plot_2, sep="")

  png(file=qc_plot_file_2, width=400, height=400)

  plot(
    NULL,
    NULL,
    xlim=c(min(rep2,na.rm=TRUE),max(rep2,na.rm=TRUE)),
    ylim=c(min(rep3,na.rm=TRUE),max(rep3,na.rm=TRUE)),
    main=paste("Rep2 vs Rep3 (r=",round(r23.cor, 2),")", sep=""),
    xlab="Normalised scores for Rep2",
    ylab="Normalised scores for Rep3",
    cex=1.5
  )
  
  for (plate in numberOfPlates) {
    points(
      rep2[which(sample_summary_data$plate == plate)], 
      rep3[which(sample_summary_data$plate == plate)],
      cex=1.5,
      pch=19,
      col=plate_cols[which(names(plate_cols) == plate)]
    )	
  }
  legend( 
    x="bottomright", 
    legend=paste("Plate",names(plate_cols), sep=""),
    col=plate_cols, lwd=2, lty=c(1,2), 
    pch=19
    )
  dev.off()

  ## scatter plot showing correlation between rep1 and rep3

  qc_plot_file_3 <- paste(datapath, qc_plot_3, sep="")

  png(file=qc_plot_file_3, width=400, height=400)

    plot(
      NULL,
      NULL,
      xlim=c(min(rep1,na.rm=TRUE),max(rep1,na.rm=TRUE)),
      ylim=c(min(rep3,na.rm=TRUE),max(rep3,na.rm=TRUE)),
      main=paste("Rep1 vs Rep3 (r=",round(r13.cor, 2),")", sep=""),
      xlab="Normalised scores for Rep1",
      ylab="Normalised scores for Rep3",
      cex=1.5
    )
  
   for (plate in numberOfPlates) {
      points(
        rep1[which(sample_summary_data$plate == plate)], 
        rep3[which(sample_summary_data$plate == plate)],
        cex=1.5,
        pch=19,
        col=plate_cols[which(names(plate_cols) == plate)]
      )	
    }
    legend( 
      x="bottomright", 
      legend=paste("Plate",names(plate_cols), sep=""),
      col=plate_cols, lwd=2, lty=c(1,2), 
      pch=19
    )
    dev.off()
}
