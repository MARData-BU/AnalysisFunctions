Complete.Mouse.GO.nd.KEGG <- function(annot.mat, GeneidCol = "Geneid", IDtype="geneSymb"){
  #annot.mat: matrix with a column "Geneid" corresponding to gene symbols mouse
  require(gtools)
  annot.mat.s <- annot.mat[order(annot.mat[,GeneidCol]),]
  
  ###############################################################
  #Agafem la description de org.Mm.eg.db
  library(org.Mm.eg.db)
  columns(org.Mm.eg.db)
  # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"      
  # [8] "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MGI"         
  # [15] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"      
  # [22] "SYMBOL"       "UNIGENE"      "UNIPROT"     
  keytypes(org.Mm.eg.db)
  metadata(org.Mm.eg.db)
  
  if(IDtype=="geneSymb") {
    GENENAME.Mm <- AnnotationDbi::select(org.Mm.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GENENAME"), keytype="SYMBOL")
  } else if (IDtype=="ENSEMBLid") {
    GENENAME.Mm <-AnnotationDbi::select(org.Mm.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GENENAME"), keytype="ENSEMBL")
    colnames(GENENAME.Mm)[colnames(GENENAME.Mm) == "ENSEMBL"] <- "SYMBOL"
  }
  
  dim(GENENAME.Mm)#16472     2
  GENENAME.Mm.agg <-aggregate(GENENAME.Mm, by=list(GENENAME.Mm$SYMBOL), FUN=function(x) paste(x, collapse="//"))
  GENENAME.Mm.agg <- GENENAME.Mm.agg[, c("Group.1", "GENENAME")]
  GENENAME.Mm.agg.s <- GENENAME.Mm.agg[order(GENENAME.Mm.agg$Group.1),]
  
  ##################################################################
  #Agafem el GO de GO.db
  library(GO.db)
  columns(GO.db)
  #[1] "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"      
  keytypes(GO.db)
  metadata(GO.db)
  
  if(IDtype=="geneSymb") {
    GO.Mm <- AnnotationDbi::select(org.Mm.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GO"), keytype="SYMBOL")
  } else if (IDtype=="ENSEMBLid") {
    GO.Mm <- AnnotationDbi::select(org.Mm.eg.db, keys=annot.mat.s[,GeneidCol], columns=c("GO"), keytype="ENSEMBL")
    colnames(GO.Mm)[colnames(GO.Mm) == "ENSEMBL"] <- "SYMBOL"
  }
  
  dim(GO.Mm)#244067      4 // 239172      4
  GO.Term <- AnnotationDbi::select(GO.db, keys=GO.Mm$GO, columns=c("TERM"), keytype="GOID")
  all.equal(GO.Mm$GO, GO.Term$GOID) #comprovem que l'ordre ?s el mateix abans de fer el cbind
  GO.annot <- cbind(GO.Mm, GO.Term)
  
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
                         GENENAME.Mm.agg.s$GENENAME,
                         GO.annot.agg.s[,c(2:ncol(GO.annot.agg.s))])
  
  colnames(NEW.annot.mat) <- c(colnames(annot.mat.s[,c(1:ncol(annot.mat.s))]),
                               "Description", "GO.BP", "GO.CC", "GO.MF")
  return(NEW.annot.mat)
  
}
