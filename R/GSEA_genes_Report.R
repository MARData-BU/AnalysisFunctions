#Magda Arnal
#02/07/2019

#Function to make the report by gene, Getting the output of the other report make functions

############################################################################################################################
##################################REPORT SUMMARY FOR EACH GENE##############################################################

GSEA_genes_Report <- function(GSEA_summary, name.f, symbol.col="SYMBOL", GSEA_Dir=GSEADir) {
  #GSEA_summary: csv amb el summary dels resultats de GSEA generat amb la funcio GSEA_report_Make
  #name.f: name of the output file
  #symbol.col: Columna que conte els gene symbols
  symb <- unique(GSEA_summary[,symbol.col])
  nTimes <- vector()
  GeneSets <- vector()
  symbols <- vector()
  description <- vector()
  pval <- vector()
  padj <- vector()
  log2FC <- vector()
  for (s in symb) {
    tab_symb <- GSEA_summary[GSEA_summary[,symbol.col] == s,]
    nTimes <- c(nTimes, nrow(tab_symb))
    GeneSets <- c(GeneSets, paste(unique(tab_symb$GENE_SET), collapse="//"))
    description <- c(description, paste(unique(tab_symb$GENENAME), collapse="//"))
    pval <- c(pval, tab_symb[1,grep("P.Value", colnames(tab_symb))])
    padj <- c(padj, tab_symb[1,grep("adj.P.Val", colnames(tab_symb))])
    log2FC <- c(log2FC, tab_symb[1,grep("logFC", colnames(tab_symb))])
  }
  gene.df <- data.frame(symb, description, nTimes, GeneSets, pval, padj, log2FC)
  #gene.df <- gene.df[order(gene.df$pval),]
  gene.df <- gene.df[order(gene.df$nTimes, decreasing = T),]
  write.csv2(gene.df, file=file.path(GSEA_Dir, paste(name.f, "_GenesEnriched",".csv", sep="")), row.names=F)
  
}
