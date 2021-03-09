data4TyersMake <- function(est_noctrls.annot, cond, fit.main, contrast){
  #est_noctrls.annot: Matrix est_nocontrols with anotations of genes annotatetableC
  #cond: vector of conditions 
  #fit.main: limma model used in DE analysis
  #contrast: vector of vectors length 2 with contrasts
  
  #Get contrast's values from limma
  ConList <- vector("list", length(contrast)) 
  for (i in 1:length(contrast)) {
    ConList[[i]] <- topTable(fit.main,n=Inf,coef=i, adjust="fdr")[,c("logFC","P.Value","adj.P.Val")]
    ConList[[i]] <- ConList[[i]][order(rownames(ConList[[i]])),]
  }
  
  #Scale values from est_noctrls for heatmap
  ia <- which(colnames(est_noctrls.annot) == "Chrom")
  est_noctrls <- est_noctrls.annot[, c(1:(ia-1))]
  est_noctrls_s<-est_noctrls[order(rownames(est_noctrls)),] 
  est_centered<-est_noctrls_s-apply(est_noctrls_s,1,mean)
  est_scaled<-est_centered/apply(est_noctrls_s,1,sd)
  #est_scaled.o <- est_scaled[, order(cond)]
  est_scaled.o <- est_scaled
  colnames(est_scaled.o) <- paste(colnames(est_scaled.o),"scaled",sep=".")
  est_scaled.o$AffyID <- rownames(est_noctrls_s)
  
  #Build matrix with mean expression per condition
  u.cond <- unique(cond)
  mean.matrix <- matrix(data= NA, nrow=nrow(est_noctrls), ncol=length(u.cond))
  for (ic in 1:length(u.cond)) {
    uc <- u.cond[ic]
    mean.matrix[,ic] <- apply(est_noctrls_s[,cond==uc],1,mean)
  }
  colnames(mean.matrix) <- paste("mean",u.cond, sep=".")
  
  #Compute FC using previously computed mean values
  FC.matrix <- matrix(data= NA, nrow=nrow(est_noctrls), ncol=length(contrast))
  col.FC.Names <- vector()
  for (ic in 1:length(contrast)) {
    FC.matrix[,ic] <- 2^abs(ConList[[ic]]$logFC) * sign(ConList[[ic]]$logFC)
    col.FC.Names <- c(col.FC.Names, paste("FC",contrast[[ic]][1], "vs", contrast[[ic]][2], sep="."))
  }
  colnames(FC.matrix) <- col.FC.Names
  
  #Build topDiff matrix with c("FC", "logFC","P.Value","adj.P.Val")
  topDiff.mat <- matrix(data= NA, nrow=nrow(est_noctrls), ncol=length(contrast)*4)
  col.topDiff.Names <- vector()
  for (ic in 1:length(contrast)) {
    icc <- 1+(4*(ic-1))
    cont.name <- paste(contrast[[ic]][1], "vs", contrast[[ic]][2], sep=".")
    topDiff.mat[,icc] <- FC.matrix[,ic]
    topDiff.mat[,icc+1] <- ConList[[ic]][,1]
    topDiff.mat[,icc+2] <- ConList[[ic]][,2]
    topDiff.mat[,icc+3] <- ConList[[ic]][,3]
    col.topDiff.Names <- c(col.topDiff.Names, colnames(FC.matrix)[ic],
                           paste(colnames(ConList[[ic]]), cont.name, sep="."))
  }
  colnames(topDiff.mat) <- col.topDiff.Names
  
  #Sort and build matrix with annotation
  topDiff.annot<-est_noctrls.annot[order(rownames(est_noctrls.annot)), 
                                   c("Symbol","mrna","UCSC_symbols","GO_biological_process",
                                     "GO_cellular_component", "GO_molecular_function", 
                                     "pathway", "Description","Chrom","Strand","Start","Stop")]
  colnames(topDiff.annot)[c(4:6)] <- c("GOBP", "GOCC", "GOMF")
  colnames(est_noctrls_s) <- paste(colnames(est_noctrls_s),"RMA",sep=".")
  
  #Build final matrix 
  data4Tyers<-data.frame(est_scaled.o, topDiff.annot,
                         topDiff.mat, mean.matrix,
                         est_noctrls_s)
  return(data4Tyers)
  #write.csv2(data4Tyers,file=paste(resultsDir,"Data4Tyers.csv",sep="/"),row.names=FALSE)
}
