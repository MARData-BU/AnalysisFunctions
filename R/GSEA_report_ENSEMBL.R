#Magda Arnal
#02/07/2019
#Function to make the GSEA report and anotate genes with ENSEMBL mouse annotation

GSEA_report_ENSEMBL <- function(GSEA_result, GeneSetAnnot=F ,GeneSets_annot=NULL, TabCols=NULL, 
                             GeneStats_annot=NULL, StatsCols=NULL, GSEA_Dir) {
  #GSEA_result:Directory path with results of one GSEA analysis
  #GeneSetAnnot: TRUE When gene sets are described in a table and we want them to be included
  #GeneSets_annot: Table with gene sets annotation and rownames(GeneSets_annot) == gene set name
  #TabCols: Columns of the GeneSets_annot to be included
  #GeneStats_annot: Table with the statistics of genes and rownames == Genesymbol
  #StatsCols: columns of GeneStats_annot to be included
  require(org.Mm.eg.db)
  require(data.table)
  GSEA_res_files <- list.files(GSEA_result, full.names = T)
  #################################################################################################
  #Resultats Positius del GSEA
  #Definim la matriu per guardar el summary dels resultats del GSEA
  pos_gs_summ <- data.frame()
  #Llegim el report de Genesets positius
  File_pos_gs <- grep("gsea_report_for_na_pos.*.xls$", GSEA_res_files, perl=T ,value=T)
  pos_gs_tab <- as.data.frame(fread(File_pos_gs))
  #Ordenem els GS per p.Val ajustat amb FDR (es pot canviar el camp que utilitzem per ordenar els GS si es vol!!)
  pos_gs_tab.o <- pos_gs_tab[order(pos_gs_tab$`FDR q-val`),]
  #Obtenim el nom dels Gene sets
  pos_gs_names <- paste(pos_gs_tab.o$NAME, "xls", sep=".")
  #Mirem a veure de quins gene sets disposem de reports
  xls.pos.tab <- pos_gs_tab.o[pos_gs_names %in% list.files(GSEA_result),]
  xls.pos.names <- paste(xls.pos.tab$NAME, "xls", sep=".")
  #Per cada gene set extraiem la info dels gens
  for (p in 1:nrow(xls.pos.tab)){
    gs.tab <- as.data.frame(fread(file.path(GSEA_result, xls.pos.names[p])))
    gs.genes.tab <- gs.tab[gs.tab$`CORE ENRICHMENT` == "Yes", c("PROBE", "RANK IN GENE LIST")]
    rownames(gs.genes.tab) <- gs.genes.tab$PROBE
    gs.genes <- gs.genes.tab$PROBE
    #Del paquet R (.db) org.Mm.eg.db Agafem la descripcio del gen i la location
    GENENAME.hs <- select(org.Mm.eg.db, keys=gs.genes, columns=c("SYMBOL", "GENENAME"), keytype="ENSEMBL")
    
    #Comprobem que l'ordre es el mateix i afegim la ponderacio del gen
    #all.equal(gs.genes.tab$PROBE,GENENAME.hs$ENSEMBL)#TRUE
    GENENAME.hs$Gene_Rank <- gs.genes.tab[GENENAME.hs$ENSEMBL,"RANK IN GENE LIST"]
    
    #Afegim els estadistics dels gens TRUE
    stat.tab <- GeneStats_annot[GENENAME.hs$ENSEMBL, StatsCols]
    stat.tab$ENSEMBL <- rownames(stat.tab)
    GENETAB <- merge(GENENAME.hs, stat.tab, all.x = T, all.y=F)
    
    #Afefim els estadistics del gene set
    GENETAB$GENE_SET <- rep(gsub(".xls","",xls.pos.names[p]),nrow(GENETAB))
    GENETAB$NES_score <- rep(xls.pos.tab$NES[p],nrow(GENETAB))
    GENETAB$FDR_qVal <- rep(xls.pos.tab$`FDR q-val`[p],nrow(GENETAB))
    GENETAB$NOM_pVal <- rep(xls.pos.tab$`NOM p-val`[p],nrow(GENETAB))
    
    if (GeneSetAnnot) {
      real.gs.name <- gsub("_ENSEMBL", "", xls.pos.tab$NAME[p])
      GS.Annot <- GeneSets_annot[real.gs.name,TabCols]
      GS.Annot.tab <- data.frame()
      for (r in 1:nrow(GENETAB)){
        GS.Annot.tab <- rbind(GS.Annot.tab, GS.Annot)
      }
      GENENAME.Def <- cbind(GENETAB[,c("GENE_SET", "NES_score","FDR_qVal","NOM_pVal")],
                            GS.Annot.tab, 
                            GENETAB[,c("SYMBOL", "ENSEMBL", "GENENAME", "Gene_Rank", StatsCols)])
    } else {
      GENENAME.Def <- GENETAB[,c("GENE_SET", "NES_score","FDR_qVal","NOM_pVal",
                                 "SYMBOL", "ENSEMBL", "GENENAME",  "Gene_Rank", StatsCols)]
    }
    
    pos_gs_summ <- rbind(pos_gs_summ, GENENAME.Def)#Anem guardant a una matriu tots els resultats ordenats per FDR.qval
  }
  #Guardem els resultats
  #write.csv2(pos_gs_summ, file=file.path(GSEA_result, "Summary_PositiveEnriched_GS.csv"), row.names=F)
  write.csv2(pos_gs_summ, file=file.path(GSEA_Dir, paste(inDirFiles[i], "_PositiveEnriched",".csv", sep="")), row.names=F)
  
  #################################################################################################
  #Resultats Negatius del GSEA
  #Definim la matriu per guardar el summary dels resultats del GSEA
  neg_gs_summ <- data.frame()
  #Llegim el report de Genesets negatius
  File_neg_gs <- grep("gsea_report_for_na_neg.*.xls$", GSEA_res_files, perl=T ,value=T)
  neg_gs_tab <- as.data.frame(fread(File_neg_gs))
  #Ordenem els GS per p.Val ajustat amb FDR (es pot canviar el camp que utilitzem per ordenar els GS si es vol!!)
  neg_gs_tab.o <- neg_gs_tab[order(neg_gs_tab$`FDR q-val`),]
  #Obtenim el nom dels Gene sets
  neg_gs_names <- paste(neg_gs_tab.o$NAME, "xls", sep=".")
  #Mirem a veure de quins gene sets disposem de reports
  xls.neg.tab <- neg_gs_tab.o[neg_gs_names %in% list.files(GSEA_result),]
  xls.neg.names <- paste(xls.neg.tab$NAME, "xls", sep=".")
  #Per cada gene set extraiem la info dels gens
  for (p in 1:nrow(xls.neg.tab)){
    gs.tab <- as.data.frame(fread(file.path(GSEA_result, xls.neg.names[p])))
    gs.genes.tab <- gs.tab[gs.tab$`CORE ENRICHMENT` == "Yes", c("PROBE", "RANK IN GENE LIST")]
    rownames(gs.genes.tab) <- gs.genes.tab$PROBE
    gs.genes <- gs.genes.tab$PROBE
    #Del paquet R (.db) org.Mm.eg.db Agafem la descripcio del gen i la location
    GENENAME.hs <- select(org.Mm.eg.db, keys=gs.genes, columns=c("SYMBOL", "GENENAME"), keytype="ENSEMBL")
    
    #Comprobem que l'ordre es el mateix i afegim la ponderacio del gen
    #all.equal(gs.genes.tab$PROBE,GENENAME.hs$SYMBOL)#TRUE
    GENENAME.hs$Gene_Rank <- gs.genes.tab[GENENAME.hs$ENSEMBL,"RANK IN GENE LIST"]
    
    #Afegim els estadistics dels gens TRUE
    stat.tab <- GeneStats_annot[GENENAME.hs$ENSEMBL, StatsCols]
    stat.tab$ENSEMBL <- rownames(stat.tab)
    GENETAB <- merge(GENENAME.hs, stat.tab, all.x = T, all.y=F)
    
    #Afefim els estadistics del gene set
    GENETAB$GENE_SET <- rep(gsub(".xls","",xls.neg.names[p]),nrow(GENETAB))
    GENETAB$NES_score <- rep(xls.neg.tab$NES[p],nrow(GENETAB))
    GENETAB$FDR_qVal <- rep(xls.neg.tab$`FDR q-val`[p],nrow(GENETAB))
    GENETAB$NOM_pVal <- rep(xls.neg.tab$`NOM p-val`[p],nrow(GENETAB))
    
    if (GeneSetAnnot) {
      real.gs.name <- gsub("_ENSEMBL", "", xls.neg.tab$NAME[p])
      GS.Annot <- GeneSets_annot[real.gs.name,TabCols]
      GS.Annot.tab <- data.frame()
      for (r in 1:nrow(GENETAB)){
        GS.Annot.tab <- rbind(GS.Annot.tab, GS.Annot)
      }
      GENENAME.Def <- cbind(GENETAB[,c("GENE_SET", "NES_score","FDR_qVal","NOM_pVal")],
                            GS.Annot.tab, 
                            GENETAB[,c("SYMBOL", "ENSEMBL", "GENENAME", "Gene_Rank", StatsCols)])
    } else {
      GENENAME.Def <- GENETAB[,c("GENE_SET", "NES_score","FDR_qVal","NOM_pVal",
                                 "SYMBOL", "ENSEMBL", "GENENAME", "Gene_Rank", StatsCols)]
    }
    
    neg_gs_summ <- rbind(neg_gs_summ, GENENAME.Def)#Anem guardant a una matriu tots els resultats ordenats per FDR.qval
  }
  #Guardem els resultats
  #write.csv2(neg_gs_summ, file=file.path(GSEA_result, "Summary_NegativeEnriched_GS.csv"), row.names=F)
  write.csv2(neg_gs_summ, file=file.path(GSEA_Dir, paste(inDirFiles[i], "_NegativeEnriched",".csv", sep="")), row.names=F)
 
}
