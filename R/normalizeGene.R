#' normalizeGene:Normalize the gene expression data.
#' @description Gene expression matrix normalization is necessary to eliminate
#' technical biases and variabilities, ensuring accurate and comparable analysis
#' of gene expression data. Here we provide `normalizeGene()`to normalize the data.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param method Method for normalizing for gene expression dataframe, one of
#' "logNormalize", "CPM", "TPM", "DESeq" or "limma"
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom limma normalizeBetweenArrays
#'
#' @return A normalized gene expression matrix.
#'
#' @export
#' @examples
#' # Load test data for cisSceQTL
#' data(testGene)
#' # normalize the data
#' eqtl <- normalizeGene(eqtl, method = "logNormalize")

normalizeGene <- function(eQTLObject, method = "logNormalize") {
  expressionMatrix <- eQTLObject@rawData$rawExpMat
  method = method

  if (!method %in% c("logNormalize", "CPM", "TPM", "DESeq", "limma")) {
    stop("Invalid method.
         Please choose from 'logNormalize', 'CPM, 'TPM', 'DESeq' or 'limma' .")
  }

  rowsum = apply(expressionMatrix, 1, sum)
  expressionMatrix = expressionMatrix[rowsum!=0,]

  if (method == "logNormalize") {
    # logNormalize normalization for non-zero values
    normalized_data <- log1p(sweep(expressionMatrix,
                                   2,
                                   Matrix::colSums(expressionMatrix),
                                   FUN = "/") * 1e4)
  } else if (method == "CPM") {
    # CPM normalization
    total_counts <- colSums(expressionMatrix)
    normalized_data <- log1p(sweep(expressionMatrix, 2, total_counts, "/") * 1e6)

  } else if (method == "TPM") {
    # TPM normalization
    library_size <- colSums(expressionMatrix)
    normalized_data <- log1p(expressionMatrix / library_size * 1e6)

  } else if (method == "DESeq") {
    # DESeq normalization
    requireNamespace("DESeq2")
    sample.df <- colnames(expressionMatrix)
    sample.df <- as.data.frame(sample.df)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = expressionMatrix,
                                          colData = sample.df,
                                          design = ~1)
    dds <- DESeq2::DESeq(dds)
    normalized_data <- DESeq2::counts(dds, normalized = TRUE)
  } else if (method == "limma") {
    # limma normalization
    requireNamespace("limma")
    normalized_data <- limma::normalizeBetweenArrays(expressionMatrix,
                                                     method="quantile")
  }

  # Output standardized result information
  cat(paste("Normalization completed using method:", method, "\n"))
  cat("Dimensions of normalized data:", dim(normalized_data), "\n")

  eQTLObject@rawData$normExpMat = normalized_data
  return(eQTLObject)
}


