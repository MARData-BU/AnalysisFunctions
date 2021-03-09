RNAseq.data4Tyers <- function(expr.mat, annot.mat, cond, fit.main, contrast, specie="human", GOndKEGG=TRUE) {
  #expr.mat: matrix of expression in CPM(E expression provided by voom), TMM or in counts, with rownames = Geneid
  #annot.mat: matrix with annotations obtained with featurecounts it has to have a column with gene symbols named "Geneid"
  #cond: vector with sample conditions used in the contrasts, same order as colnames(expr.mat)
  #fit.main: model to use in limma
  #contrast: List of vectors with each contrast to use
  #specie: specie used to make annotations
  #GOndKEGG: Annotation with GO and KEGG. Takes some time
  
  #Obtain contrasts with limma
  ConList <- vector("list", length(contrast)) 
  for (i in 1:length(contrast)) {
    ConList[[i]] <- topTable(fit.main,n=Inf,coef=i, adjust="fdr")[,c("logFC","P.Value","adj.P.Val")]
    ConList[[i]] <- ConList[[i]][order(rownames(ConList[[i]])),]
  }
  
  #Scale values from count matrix (normalized or not) for the heatmap
  est_noctrls_s<-expr.mat[order(rownames(expr.mat)),]
  est_centered<-est_noctrls_s-apply(est_noctrls_s,1,mean)
  est_scaled<-est_centered/apply(est_noctrls_s,1,sd)
  colnames(est_scaled) <- paste(colnames(est_scaled),"scaled",sep=".")
  
  #Build matrix with mean expression per condition
  u.cond <- unique(cond)
  mean.matrix <- matrix(data= NA, nrow=nrow(expr.mat), ncol=length(u.cond))
  for (ic in 1:length(u.cond)) {
    uc <- u.cond[ic]
    mean.matrix[,ic] <- apply(est_noctrls_s[,cond==uc],1,mean)
  }
  colnames(mean.matrix) <- paste("mean",u.cond, sep=".")
  rownames(mean.matrix) <- rownames(est_noctrls_s)
  
  #Compute FC using previously computed mean values
  FC.matrix <- matrix(data= NA, nrow=nrow(expr.mat), ncol=length(contrast))
  col.FC.Names <- vector()
  for (ic in 1:length(contrast)) {
    FC.matrix[,ic] <- 2^abs(ConList[[ic]]$logFC) * sign(ConList[[ic]]$logFC)
    col.FC.Names <- c(col.FC.Names, paste("FC",contrast[[ic]][1], "vs", contrast[[ic]][2], sep="."))
  }
  colnames(FC.matrix) <- col.FC.Names
  
  #Build topDiff matrix with c("FC", "logFC","P.Value","adj.P.Val")
  topDiff.mat <- matrix(data= NA, nrow=nrow(expr.mat), ncol=length(contrast)*4)
  col.topDiff.Names <- vector()
  for (ic in 1:length(contrast)) {
    icc <- 1+(4*(ic-1))
    cont.name <- paste(contrast[[ic]][1], "vs", contrast[[ic]][2], sep=".")
    topDiff.mat[,icc] <- FC.matrix[,ic]
    topDiff.mat[,icc+1] <- ConList[[ic]][,1]
    topDiff.mat[,icc+2] <- ConList[[ic]][,2]
    topDiff.mat[,icc+3] <- ConList[[ic]][,3]
    col.topDiff.Names <- c(col.topDiff.Names, col.FC.Names[ic],
                           paste(colnames(ConList[[ic]]), cont.name, sep="."))
  }
  colnames(topDiff.mat) <- col.topDiff.Names
  
  
  #If GOndKEGG annotations are requested, call function to annoatate
  if(GOndKEGG) {
    #Sort and build matrix with annotations
    if (specie =="human") {
      annot.mat.complt <- Complete.Human.GO.nd.KEGG(annot.mat)
      
    } else if (specie =="mouse") {
      annot.mat.complt <- Complete.Mouse.GO.nd.KEGG(annot.mat)
    }

  } else {
    #Sort annotation matrix
    annot.mat.complt <- annot.mat[match(rownames(est_scaled),annot.mat$Geneid),]
  }
  
  #Make sure everything is in same order before merging and returning
  a=all.equal(rownames(est_scaled), rownames(est_noctrls_s))
  b=all.equal(rownames(est_scaled), rownames(ConList[[1]]))
  c=all.equal(rownames(est_scaled), rownames(mean.matrix))
  d=all.equal(rownames(est_scaled), annot.mat.complt$Geneid)
  tryCatch(
    expr = {
      a&b&c&d #If all objects are in same order
      message("All objects are in same order. data4Tyers successfully created!")
      
      #Build final matrix only if above expression is true
      data4Tyers<-data.frame(est_scaled, annot.mat.complt,
                             topDiff.mat, mean.matrix,
                             est_noctrls_s)
      return(data4Tyers)
      
    },
    error = function(e){ #If  a&b&c&d generated error, pop up following message:
      message('Something is not in the correct order. Plese check that all inputs have same length and Geneid:')
      
      message('all.equal(rownames(est_scaled), rownames(est_noctrls_s): ',a)
      message('all.equal(rownames(est_scaled), rownames(ConList[[1]])): ',b)
      message('all.equal(rownames(est_scaled), rownames(mean.matrix)): ',c)
      message('all.equal(rownames(est_scaled), annot.mat.complt$Geneid): ',d)

    }
  ) 


}
