#' filterGeneSNP: Filter gene expression matrix and genotype matrix.
#'
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param snp.number.of.cells.percent Only SNPs where cells with each of the
#' different genotypes (REF and ALT, or AA, Aa, and aa) individually account for
#'  at least `snp.number.of.cells.percent`% of the total cells will be considered.
#'  The default valus is 10.
#' @param expression.min threshold for valid gene expression levels, utilized
#' alongside another parameter, expression.number.of.cells. The default value is 0.
#' @param expression.number.of.cells.percent Only genes with expression levels
#' exceeding `expression.min` in at least `expression.number.of.cells.percent`%
#' of cells are considered. The default value is 10.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' eqtl <- filterGeneSNP(eqtl,
#'                       snp.number.of.cells.percent = 10,
#'                       expression.min = 0,
#'                       expression.number.of.cells.percent = 10
#'                       )
filterGeneSNP <- function(
    eQTLObject,
    snp.number.of.cells.percent = 10,
    expression.min = 0,
    expression.number.of.cells.percent = 10,
    ...
){
  nor_expressionMatrix <- as.data.frame(eQTLObject@rawData$normExpMat)
  if(!is.null(nor_expressionMatrix)){
    expression.number.of.cells <- ceiling((expression.number.of.cells.percent / 100) * ncol(nor_expressionMatrix))
    valid_genes <- rownames(nor_expressionMatrix)[apply(nor_expressionMatrix >
                                                          expression.min, 1, sum) >=
                                                    expression.number.of.cells]
    filtered_expressionMatrix <- nor_expressionMatrix[valid_genes, ]
    filtered_expressionMatrix <- as.matrix(filtered_expressionMatrix)
  }else{
    stop("Please normalize the raw expression data first.")
  }


  snpMatrix <- as.data.frame(eQTLObject@rawData$snpMat)
  snp.list <- rownames(snpMatrix)
  snp.number.of.cells <- ceiling((snp.number.of.cells.percent / 100) * ncol(snpMatrix))

  biClassify = eQTLObject@biClassify

  if(biClassify  == TRUE){

    snpMatrix[snpMatrix == 3] <- 2

    snp_counts <- apply(snpMatrix, 1, function(row) {
      c(sum(row == 1), sum(row == 2))
    })

    snp_counts_df <- as.data.frame(t(snp_counts))
    names(snp_counts_df) <- c("count_ref", "count_alt")

    filtered_snp_ids <- rownames(snp_counts_df)[apply(snp_counts_df, 1, function(row) {
      all(row > snp.number.of.cells)
    })]

    filtered_snpMatrix <- snpMatrix[filtered_snp_ids, ]

  }else if(biClassify  == FALSE){

    snp_counts <- apply(snpMatrix, 1, function(row) {
      c(sum(row == 1), sum(row == 3), sum(row == 2))
    })

    snp_counts_df <- as.data.frame(t(snp_counts))
    names(snp_counts_df) <- c("count_AA", "count_Aa", "count_aa")

    filtered_snp_ids <- rownames(snp_counts_df)[apply(snp_counts_df, 1, function(row) {
      all(row > snp.number.of.cells)
    })]

    filtered_snpMatrix <- snpMatrix[filtered_snp_ids, ]
    filtered_snpMatrix <- as.matrix(filtered_snpMatrix)
  }

  eQTLObject@filterData <- list(expMat = filtered_expressionMatrix, snpMat = filtered_snpMatrix)
  return(eQTLObject)
}
