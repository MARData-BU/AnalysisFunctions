# AnalysisFunctions
Package with basic analysis functions for genomic data.

- **createTargets**: Function to create phenodata of the samples in the gene expression matrix.
- **DeleteCtrlsGral**: Function to delete control probe IDs from the expression matrix.
- **gct.data4Tyers**: Function to create gct format files for GSEA analysis.
- **generate_cls**: Function to create .cls format files for GSEA analysis.
- **generate_gct**: Function to create .gct format files for GSEA analysis.
- **generate_rnk**: Function to create .rnk files for GSEA analysis.
- **MaxMin.Sym.Row**: Function to resolve GeneID duplicates chosing the minimum or maximum values.
- **Unique.Sym.Estimates**: Function to resolve GeneID duplicates computing the mean and the median.
- **data4TyersMake**: Function to make a data frame with a structure of data4Tyers from arrays (RMA data).
- **RNAseq.data4Tyers**: Function to make a data frame with a structure of data4Tyers from RNAseq data.
- **Complete.Human.GO.nd.KEGG**: Function to complete annotations from RNAseq (featurecounts function) with GO and KEGG databases, for each gene and human specie.
- **Complete.Mouse.GO.nd.KEGG**: Function to complete annotations from RNAseq (featurecounts function) with GO and KEGG databases, for each gene and mouse specie.
- **GSEA_report_HS**: Function to make the report of GSEA based on human symbols
- **GSEA_report_MS**: Function to make the report of GSEA based on mouse symbols
- **GSEA_report_ENSEMBL**: Function to make the report of GSEA based on ENSEMBL IDs
- **GSEA_genes_Report**: Function to make the report of GSEA sumarized by gene
- **make.excel**: Function to make the excel of results from data4tyers file.


