#' Example dataset of gene expression matrix
#'
#' This dataset contains normalized expression values
#' for demonstration purposes.
#' Normalized process was performed by DESeq2 as follows:
#' 

#' \preformatted{
#' dds <- DESeqDataSetFromMatrix(
#' 	countData=countData,
#' 	colData=sample_info,
#' 	design=~factor(level))
#' rld <- rlog(dds, blind=FALSE)
#' GSE143371_count <- assay(rld)
#' }
#'
#' @format A matrix with 57131 rows and 6 columns
#' @usage data(GSE143371_count)
#' @source https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143371
"GSE143371_count"
