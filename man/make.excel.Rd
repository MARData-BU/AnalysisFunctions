\name{make.excel}
\alias{MAKE.EXCEL}
\title{Create excel of results}
\usage{
make.excel(pathinput = pathinput,fileinput = fileinput,contrast=contrast,pathoutput = pathoutput,filename = filename)
}
\arguments{
\item{pathinput}{The pathway where data4Tyers file is}
\item{fileinput}{The name of data4Tyers file (.rds)}
\item{contrast}{The contrast used in DE analysis. It has to be a list (ex: list(c("WT1","WT2") or list of lists (ex: list(c("WT_1", "KO_1"), c("WT_2", "KO_2"), c("WT1","WT2"), c("KO1","KO2"))})
\item{pathoutput}{The pathway for the excel file output}
\item{filename}{The name for the excel file (without extension)}
\item{pvalue}{The value of p.value which we want to filter data. pvalue = NULL by default}
\item{padj}{The value of adjusted p.value which we want to filter data. padj < 0.05 by default. If not all sheets can be filtered by padj < 0.05, it will be used pvalue < 0.05 to filter all data}
\item{logFC}{The value of logFoldChange which we want to filter data. logFC = 1 by default}

}
\description{
Create excel of results from data4tyers rds file
}
\examples{
?
  #' ###### FILTERS BY DEFAULT ######
#' contrast <- list(c("WT_1", "KO_1"), c("WT_2", "KO_2"), c("WT1","WT2"), c("KO1","KO2"))
#' fileinput <- "data4Tyers.rds"
#' filename <- "ProjectRNASeq"
#' pathinput <- "Project/Analysis"
#' pathoutput <- "Project/Analysis/Results
#'
#' make.excel(fileinput = fileinput, filename = filename, pathinput = pathinput, pathoutput = pathoutput, contrast = contrast)

  #' ###### FILTERS CHANGED  ######
#' contrast <- list(c("WT_1", "KO_1"), c("WT_2", "KO_2"), c("WT1","WT2"), c("KO1","KO2"))
#' fileinput <- "data4Tyers.rds"
#' filename <- "ProjectRNASeq"
#' pathinput <- "Project/Analysis"
#' pathoutput <- "Project/Analysis/Results
#'
#' make.excel(fileinput = fileinput, filename = filename, pathinput = pathinput, pathoutput = pathoutput, contrast = contrast, pvalue = 0.01, logFC = 0.585)
  
  ########################## LIMITATIONS ###########################
  #' - Contrast has to be <= 25 characters (EX: pos_treat.vs.neg_no.treat) 
  #' - Contrast has no to include "FC". (EX: pos_FC.vs.neg_FC) 
  #' - 20 possible contrasts (you can include more colors in colors4stats if you want to perform more than 20 contrasts) 
}
