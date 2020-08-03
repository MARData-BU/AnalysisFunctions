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
#' filename <- ""20200526_Excel_EGonzalez_MCiscar.xlsx"
#' pvalue <- 0.05
#' FC <- 1
#' pathinput <- "/bicoh/MARGenomics/20200526_EGonzalez_MCiscar/Analysis/Results"
#'
#'
#' MAKE.EXCEL(fileinput = "Data4Tyers.rds", contrast, pvalue , padj = NULL, FC )

}
