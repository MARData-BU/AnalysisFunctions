\name{GSEA_genes_Report}
\alias{GSEA_genes_Report}

\title{
Function to make a GSEA summary table with Mouse gene symbol annotation
}
\description{
In order fo the results of GSEA to be more comprehensive, we make a function to summarize the results
}
\usage{
GSEA_genes_Report(GSEA_summary, name.f, symbol.col="SYMBOL", GSEA_Dir)
}

\arguments{
  \item{GSEA_summary}{
csv with summary of GSEA results generated with function GSEA_report_Make
}
  \item{name.f}{
Name of the output file
}
  \item{symbol.col}{
Column with gene symbols
}
  \item{GSEA_Dir}{
Directory to save the results
}
}
\details{
This function makes a summary of the GSEA report, but this time by gene
}
\value{
data.frame object
}
\examples{
inDirFiles <- grep(".csv",list.files(GSEADir, full.names = F), value=T)

for (i in 1:length(inDirFiles)) {
  GSEA_summary <- as.data.frame(fread(file.path(GSEADir, inDirFiles[i]), dec=","))
  print(dim(GSEA_summary))#6354   11
  name_summary <- gsub(".csv", "", inDirFiles[i])
  GSEA_genes_Report(GSEA_summary=GSEA_summary, name.f = name_summary, symbol.col="SYMBOL", GSEA_Dir)
}
}

