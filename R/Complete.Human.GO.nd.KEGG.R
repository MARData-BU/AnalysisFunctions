Complete.Human.GO.nd.KEGG <- function(annot.mat, GeneidCol = "Geneid", IDtype="geneSymb" ){
  #annot.mat: matrix with a column "Geneid" corresponding to gene symbols human
  require(gtools)
  annot.mat.s <- annot.mat[order(annot.mat[,GeneidCol]),]
  
  ###############################################################
  #Agafem la description de org.Mm.eg.db
  require(org.Hs.eg.db)
  columns(org.Hs.eg.db)
  # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"      
  # [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MGI"         
  # [15] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"      
  # [22] "SYMBOL"       "UNIGENE"      "UNIPROT"     
  keytypes(org.Hs.eg.db)
  metadata(org.Hs.eg.db)
  
  if(IDtype=="geneSymb") {
    GENENAME.hs <-  select(org.Hs.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GENENAME"), keytype="SYMBOL")
  } else if (IDtype=="ENSEMBLid") {
    GENENAME.hs <-  select(org.Hs.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GENENAME"), keytype="ENSEMBL")
    colnames(GENENAME.hs)[colnames(GENENAME.hs) == "ENSEMBL"] <- "SYMBOL"
  }
  
  
  dim(GENENAME.hs)#40053     2 (surt un warning, pero es ok)
  GENENAME.hs.agg <-aggregate(GENENAME.hs, by=list(GENENAME.hs$SYMBOL), FUN=function(x) paste(x, collapse="//"))
  GENENAME.hs.agg <- GENENAME.hs.agg[, c("Group.1", "GENENAME")]
  GENENAME.hs.agg.s <- GENENAME.hs.agg[order(GENENAME.hs.agg$Group.1),]
  
  ##################################################################
  #Agafem el GO de GO.db
  require(GO.db)
  columns(GO.db)
  #[1] "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"      
  keytypes(GO.db)
  metadata(GO.db)
  
  if(IDtype=="geneSymb") {
    GO.hs <-  select(org.Hs.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GO"), keytype="SYMBOL")#warning pero es ok
  } else if (IDtype=="ENSEMBLid") {
    GO.hs <-  select(org.Hs.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GO"), keytype="ENSEMBL")
    colnames(GO.hs)[colnames(GO.hs) == "ENSEMBL"] <- "SYMBOL"
  }
  
  dim(GO.hs)#268499      4
  GO.Term <-  select(GO.db, keys=GO.hs$GO, columns=c("TERM"), keytype="GOID")#warning pero es ok
  all.equal(GO.hs$GO, GO.Term$GOID) #TRUE
  GO.annot <- cbind(GO.hs, GO.Term)
  
  #Separem els GO en MF/CC i BP
  GO.BP <- vector(mode="character", length=nrow(GO.annot))
  GO.CC <- vector(mode="character", length=nrow(GO.annot))
  GO.MF <- vector(mode="character", length=nrow(GO.annot))
  for (r in c(1:nrow(GO.annot))) {
    Ont <- GO.annot$ONTOLOGY[r]
    if (!is.na(Ont)){
      if(Ont == "BP"){
        GO.BP[r] <- paste(GO.annot[r,c("ONTOLOGY","GOID", "TERM")], collapse=" ")
      } else if (Ont == "CC") {
        GO.CC[r] <- paste(GO.annot[r,c("ONTOLOGY","GOID", "TERM")], collapse=" ")
      } else if (Ont == "MF") {
        GO.MF[r] <- paste(GO.annot[r,c("ONTOLOGY","GOID", "TERM")], collapse=" ")
      }
    }
  }
  GO.annot.all <- data.frame(SYMBOL=GO.annot$SYMBOL, GO.BP, GO.CC, GO.MF)
  
  symb.vect=unique(GO.annot$SYMBOL)
  GO.BP.p <- vector(mode="character", length=length(symb.vect))
  GO.CC.p <- vector(mode="character", length=length(symb.vect))
  GO.MF.p <- vector(mode="character", length=length(symb.vect))
  for(si in  c(1:length(symb.vect))) {
    symb <- symb.vect[si]
    GO.mat <- GO.annot.all[GO.annot.all$SYMBOL == symb,]
    GO.BP.p[si] <- paste(unique(GO.mat$GO.BP), collapse="")
    GO.CC.p[si] <- paste(unique(GO.mat$GO.CC), collapse="")
    GO.MF.p[si] <- paste(unique(GO.mat$GO.MF), collapse="")
    
    GO.BP.p[si] <- gsub("BP", "//", GO.BP.p[si])
    GO.CC.p[si] <- gsub("CC", "//", GO.CC.p[si])
    GO.MF.p[si] <- gsub("MF", "//", GO.MF.p[si])
    
    GO.BP.p[si] <- sub("//", "", GO.BP.p[si])
    GO.CC.p[si] <- sub("//", "", GO.CC.p[si])
    GO.MF.p[si] <- sub("//", "", GO.MF.p[si])
  }
  
  GO.annot.desg <- data.frame(SYMBOL=symb.vect, GO.BP.p, GO.CC.p, GO.MF.p)
  GO.annot.agg.s <- GO.annot.desg[order(as.character(GO.annot.desg$SYMBOL)),]
  
  
  NEW.annot.mat <- cbind(annot.mat.s[,c(1:ncol(annot.mat.s))],
                         GENENAME.hs.agg.s$GENENAME,
                         GO.annot.agg.s[,c(2:ncol(GO.annot.agg.s))])
  
  colnames(NEW.annot.mat) <- c(colnames(annot.mat.s[,c(1:ncol(annot.mat.s))]),
                               "Description", "GO.BP", "GO.CC", "GO.MF")
  return(NEW.annot.mat)
  
}

