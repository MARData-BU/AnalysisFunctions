\name{MAKE.EXCEL}
\alias{make.excel}
\title{Create excel}
\usage{
MAKE.EXCEL(pathinput,fileinput,contrast,filename, pvalue,padj, FC)
}
\arguments{
\item{pathinput}{The pathway where data4Tyers file is}
\item{fileinput}{The name of data4Tyers file (.rds)}
\item{contrast}{The contrast used in DE analysis}
\item{filename}{The name for the excel file}
\item{pvalue}{The value of p.value which we want to filter data. For filter by pvalue, padj = NULL}
\item{padj}{The value of adjusted p.value which we want to filter data.}
\item{FC}{The value of FoldChange which we want to filter data.}

}
\description{
Create excel
}
\examples{
?
#' contrast <- list(c("PyMT_WT", "PyMT_LysMCre"), c("PyMT_WT", "KO_WT"))
#' fileinput <- "data4Tyers.rds"
#' filename <- "ProjectRNASeq.xlsx"
#' pvalue <- 0.05
#' FC <- 1
#' pathinput <- "/bicoh/MARGenomics/Projecte/Analysis/Results"
#'
#'
#' MAKE.EXCEL(fileinput = "Data4Tyers.rds", contrast, pvalue , padj = NULL, FC )

  
  ########################## LIMITATIONS ###########################
  #' - Contrast has to be <= 25 characters (EX: pos_MV.tr.vs.neg_MV.no.tr) 
  #' - Contrast has no to include "FC". (EX: pos_FC.vs.neg_FC) 
  #' - 20 possible contrasts (you can include more colors in colors4stats if you want to perform more than 20) 
}
