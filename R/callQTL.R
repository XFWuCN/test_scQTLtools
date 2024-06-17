#' callQTL: Uncover single-cell eQTLs exclusively using scRNA-seq data.A function
#' designed to identify eQTLs from scRNA-seq data.
#'
#' @param useModel Model for fitting dataframe, one of "possion", "zinb", or "linear".
#' @param p.adjust.Threshold Only SNP-Gene pairs with adjusted p-values meeting
#' the threshold will be displayed. The default value is 0.05.
#' @param p.adjust.method Methods for p-value adjusting, one of "bonferroni",
#' "holm", "hochberg", "hommel" or "BH". The default option is "BH".
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param gene_ids A gene ID or a list of gene IDS.
#' @param downstream Being used to match SNPs within a base range defined by the
#' start position of the specified gene(s) or all genes.
#' @param upstream Being used to match SNPs within a base range defined by the
#' end position of the specified gene(s) or all genes.
#'
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb fitdistr
#' @importFrom stringr str_split
#' @importFrom VGAM dzinegbin
#' @importFrom bbmle mle2
#' @importFrom gamlss gamlssML
#' @importFrom pscl zeroinfl
#' @importFrom stats p.adjust pchisq plogis lm
#' @importFrom glmmTMB glmmTMB
#' @importFrom dplyr mutate_all
#' @importFrom stats glm
#'
#' @return A dataframe, each row describes eQTL discovering result of a SNP-Gene
#' pair.
#' @export
#'
#' @examples
#' eqtl <- callQTL(eqtl,
#'                 gene_ids = NULL,
#'                 downstream = NULL,
#'                 upstream = NULL,
#'                 p.adjust.method = "bonferroni",
#'                 useModel = "zinb",
#'                 p.adjust.Threshold = 0.05,
#'                 logfc.threshold = 0.1
#'                 )


callQTL <- function(
    eQTLObject,
    gene_ids = NULL,
    downstream = NULL,
    upstream = NULL,
    p.adjust.method = "bonferroni",
    useModel = "zinb",
    p.adjust.Threshold = 0.05,
    logfc.threshold = 0.1){

  if(length(eqtl@filterData) == 0){
    cat("Please filter the data first.")
  }else{
    expressionMatrix = eQTLObject@filterData$expMat
    snpMatrix = eQTLObject@filterData$snpMat
  }
  biClassify = eQTLObject@biClassify
  species = eQTLObject@species

  if(is.null(gene_ids) && is.null(upstream) && is.null(downstream)){
    return(NULL)
  }else{
    if (!is.null(species) && species != ""){
      if (species == "human"){
        snp_dataset = "hsapiens_snp"
        gene_dataset = "hsapiens_gene_ensembl"
        OrgDb = org.Hs.eg.db
      } else if (species == "mouse"){
        snp_dataset = "mmusculus_snp"
        gene_dataset = "mmusculus_gene_ensembl"
        OrgDb = org.Mm.eg.db
      } else {
        stop("Please enter 'human' or 'mouse'.")
      }
    } else {
      stop("The 'species' variable is NULL or empty.")
    }
  }


  creat_snps_loc <- function(snp.list){

    # Get location for each SNP
    snp_mart <- useMart(biomart = "ENSEMBL_MART_SNP", dataset = snp_dataset)

    snps_loc <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                      filters = "snp_filter",
                      values = snp.list,
                      mart = snp_mart)
    colnames(snps_loc)[which(colnames(snps_loc) == "chrom_start")] <- "position"
    rownames(snps_loc) <- snps_loc[, 1]
    return(snps_loc)
  }

  creat_gene_loc <- function(gene.list){

    # Get location for each gene
    gene_mart = useEnsembl(biomart="ensembl",
                           dataset="hsapiens_gene_ensembl")

    ensembls <- mapIds(OrgDb,
                       keys = gene.list,
                       keytype = "SYMBOL",
                       column="ENSEMBL")
    ensembls = as.data.frame(ensembls)
    ensembls_id = unique(ensembls$ensembls)

    #ensembl_gene_id,external_gene_name
    gene_attributes=c('external_gene_name',
                      'chromosome_name',
                      'start_position',
                      'end_position')

    gene_loc <- getBM(attributes = gene_attributes,
                      filters = "external_gene_name",
                      values = gene.list,
                      mart = gene_mart)

    rownames(gene_loc) <- gene_loc[, 1]
    return(gene_loc)
  }

  check_snpList <- function(snp.list){
    if(grepl("^rs", snp.list[[1]][1])){  # chromosome start with "rs"
      creat_snps_loc(snp.list)
    } else if(grepl("\\d+:\\d+", snp.list[[1]][1])){  # chr:position
      # build an empty dataframe
      snps_df <- data.frame(refsnp_id = character(),
                            chr_name = character(),
                            position = numeric(),
                            stringsAsFactors = FALSE)
      snps_loc <- snps_df
      for (i in 1:length(snp.list)){
        snp_parts <- strsplit(snp.list[i], ":")[[1]]
        snps_loc <- rbind(snps_loc, data.frame(refsnp_id = snp.list[i],
                                               chr_name = snp_parts[1],
                                               position = as.numeric(snp_parts[2])))
      }
      return(snps_loc)
    } else{
      stop("Error: SNP does not match expected format.")
    }
  }

  snp.list <- rownames(eQTLObject@filterData$snpMat)
  gene.list <- rownames(eQTLObject@filterData$expMat)
  # When users do not perform cis-eQTL filtering
  if(is.null(gene_ids) && is.null(upstream) && is.null(downstream)){
    matched_gene <- gene.list
    matched_snps <- snp.list
  }else if(!is.null(gene_ids) && is.null(upstream) && is.null(downstream)){
    matched_snps <- snp.list
    if (all(gene_ids %in% gene.list)) {
      matched_gene = gene_ids
    } else {
      stop("The input gene_ids contain non-existent gene IDs. Please re-enter.")
    }
  }else if(is.null(gene_ids) && !is.null(upstream) && !is.null(downstream)){
    if (downstream >= upstream) {
      stop("downstream should be smaller than upstream.")
    }

    # 连接数据库，构建snps_loc和gene_loc

    # Get location for each SNP
    snps_loc <- check_snpList(snp.list)

    gene_loc <- creat_gene_loc(gene.list)

    # 遍历基因和SNP，匹配在基因范围内的SNP，输出gene-SNP pairs
    matched_gene <- c()
    matched_snps <- c()

    for (i in 1:nrow(gene_loc)) {
      # get range
      gene_start1 <- gene_loc$start_position[i] + downstream
      gene_end1 <- gene_loc$end_position[i] + upstream

      # 遍历每个SNP
      for (j in 1:nrow(snps_loc)) {
        # Retrieve the position information of the current SNP
        snp_chr <- snps_loc$chr_name[j]
        snp_pos <- snps_loc$position[j]

        # Determine if the current SNP is within the range of the gene
        if (snp_chr == gene_loc$chromosome_name[i] &&
            snp_pos >= gene_start1 &&
            snp_pos <= gene_end1) {
          # Add matching SNPs to the current gene's SNP list
          matched_snps <- c(matched_snps, snps_loc$refsnp_id[j])

          matched_gene <- c(matched_gene, gene_loc$ensembl_gene_id[i])

        }
      }
    }

    matched_snps <- unique(matched_snps)
    matched_gene <- unique(matched_gene)

  }else if(!is.null(gene_ids) && !is.null(upstream) && !is.null(downstream)){
    if (downstream >= upstream) {
      stop("downstream should be smaller than upstream.")
    }

    # 连接数据库，构建snps_loc和gene_loc
    snps_loc <- check_snpList(snp.list)

    if (all(gene_ids %in% gene.list)) {
      gene.list = gene_ids
    } else {
      stop("Invalid gene ID！Please modify the filtering condition in the
            previous step or replace it with another gene_ids.")
    }
    gene_loc <- creat_gene_loc(gene.list)

    matched_gene <- c()
    matched_snps <- c()

    for (i in 1:nrow(gene_loc)) {
      gene_start1 <- gene_loc$start_position[i] + downstream
      gene_end1 <- gene_loc$end_position[i] + upstream

      for (j in 1:nrow(snps_loc)) {

        snp_chr <- snps_loc$chr_name[j]
        snp_pos <- snps_loc$position[j]

        # 判断当前SNP是否在基因的范围内
        if (snp_chr == gene_loc$chromosome_name[i] &&
            snp_pos >= gene_start1 &&
            snp_pos <= gene_end1) {
          # 将匹配的SNP添加到当前基因的SNP列表中
          matched_snps <- c(matched_snps, snps_loc$refsnp_id[j])
          matched_gene <- c(matched_gene, gene_loc$ensembl_gene_id[i])

        }
      }

    }
    matched_snps <- unique(matched_snps)
    matched_gene <- unique(matched_gene)
  }else{
    stop("Please enter upstream and downstream simultaneously.")
  }

  possionModel <- function(
    eQTLObject,
    genelist,
    snplist,
    biClassify = FALSE,
    p.adjust.method = "bonferroni",
    p.adjust.Threshold = 0.05,
    logfc.threshold = 0.1){

    expressionMatrix <- eQTLObject@filterData$expMat
    expressionMatrix <- round(expressionMatrix * 1000)
    snpMatrix <- eQTLObject@filterData$snpMat
    unique_group <- unique(eQTLObject@groupBy$group)

    result_all <- data.frame()

    message("Start calling eQTL")
    for(k in unique_group){

      result <- data.frame(
        SNPid = character(),
        group = character(),
        Geneid = character(),
        pvalue = double(),
        adjusted_pvalue = double(),
        b = double(),
        abs_b = double(),
        Remark = character(),
        stringsAsFactors=FALSE)

      split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == k]
      expressionMatrix_split <- expressionMatrix[, split_cells]
      snpMatrix_split <- snpMatrix[, split_cells]

      if(biClassify == FALSE){
        message(k,':')
        message('0%   10   20   30   40   50   60   70   80   90   100%' )
        message('[----|----|----|----|----|----|----|----|----|----|')
        pb <- progress_bar$new(total = length(snplist),
                               format = "[:bar]",
                               clear = FALSE,
                               width = 51)
        for(i in 1:length(snplist)){
          snpid <- snplist[i]
          snp_mat <- snpMatrix_split[snpid, ]
          snp_mat <- as.data.frame(snp_mat)
          snp_mat$cells = rownames(snp_mat)

          replace_2_and_3 <- function(x) {
            ifelse(x == 2, 3, ifelse(x == 3, 2, x))
          }

          snp_mat_new <- snp_mat %>%
            mutate_all(funs(replace_2_and_3))

          genes <- genelist

          for(j in 1:length(genes)){
            gene_id <- genes[j]
            gene_mat <- expressionMatrix_split[gene_id, ]
            gene_mat <- as.data.frame(gene_mat)
            gene_mat$cells = rownames(gene_mat)

            combined_df <- merge(snp_mat_new, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 0)

            lmodel = stats::glm(combined_df$gene_mat ~ combined_df$snp_mat, family = poisson());

            lmout_pvalue = summary(lmodel)$coefficients[2, "Pr(>|z|)"]
            lmout_b = summary(lmodel)$coefficients[2, "Estimate"]

            new_row <- data.frame(SNPid = snpid,
                                  group = k,
                                  Geneid = genes[j],
                                  pvalue = lmout_pvalue,
                                  b = lmout_b)
            result <- rbind(result, new_row)
          }
          pb$tick()
        }
        message("finished!")

      }else if(biClassify == TRUE){

        snpMatrix_split[snpMatrix_split == 3] <- 2

        message(k,':')
        message('0%   10   20   30   40   50   60   70   80   90   100%' )
        message('[----|----|----|----|----|----|----|----|----|----|')
        pb <- progress_bar$new(total = length(snplist),
                               format = "[:bar]",
                               clear = FALSE,
                               width = 51)
        for(i in 1:length(snplist)){
          snpid <- snplist[i]
          snp_mat <- snpMatrix_split[snpid, ]
          snp_mat <- as.data.frame(snp_mat)
          snp_mat$cells = rownames(snp_mat)

          genes <- genelist

          for(j in 1:length(genes)){
            gene_id <- genes[j]
            gene_mat <- expressionMatrix_split[gene_id, ]
            gene_mat <- as.data.frame(gene_mat)
            gene_mat$cells = rownames(gene_mat)

            combined_df <- merge(snp_mat, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 0)

            lmodel = glm(combined_df$gene_mat ~ combined_df$snp_mat, family = poisson());

            lmout_pvalue = summary(lmodel)$coefficients[2, "Pr(>|z|)"]
            lmout_b = summary(lmodel)$coefficients[2, "Estimate"]

            new_row <- data.frame(SNPid = snpid,
                                  group = k,
                                  Geneid = genes[j],
                                  pvalue = lmout,
                                  b = lmout_b)

            result <- rbind(result, new_row)
          }
          pb$tick()
        }
        message("finished!")

      }else{
        stop("biClassify can only be selected as 'TRUE' or 'FALSE'")
      }


      if (!p.adjust.method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
        stop("Invalid p-adjusted method.
         Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel', or'fdr or BH'.")
      }

      # adjust p-value
      result[,"adjusted_pvalue"] <- p.adjust(result[,"pvalue"], method = "BH")
      result <- result[order(result[,"adjusted_pvalue"]),]
      rownames(result) <- NULL
      result <- result[result$adjusted_pvalue <= p.adjust.Threshold, ]

      # abs_b
      result <- result %>%
        mutate(abs_b = abs(b))

      result <- result[result$abs_b >= logfc.threshold, ]

      result_all <- rbind(result_all, result)

    }
    return(result_all)

  }

  zinbModel <- function(
    eQTLObject,
    genelist,
    snplist,
    biClassify = FALSE,
    p.adjust.method = "bonferroni",
    p.adjust.Threshold = 1e-5){

    # p-value correction methods
    if (!p.adjust.method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
      stop("Invalid p-adjusted method.
           Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel', or'fdr or BH'.")
    }

    expressionMatrix <- eQTLObject@filterData$expMat
    snpMatrix <- eQTLObject@filterData$snpMat
    expressionMatrix <- round(expressionMatrix * 1000)
    unique_group <- unique(eQTLObject@groupBy$group)

    result_all <- data.frame()

    message("Start calling eQTL")
    for(j in unique_group){
      if(biClassify == TRUE){

        snpMatrix_split[snpMatrix_split == 3] <- 2

        eQTLcalling <- function(i){
          if (i %% 100 == 0) {
            gc()
          }

          snpid <- snplist[i]
          ref_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid,] == 1]
          alt_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid,] == 2]

          if(length(ref_cells) > 0 && length(alt_cells) > 0){

            genes <- genelist
            gene.cnt <- 0

            # snp ----------------------------------------------------------------
            results_snp <- data.frame(
              SNPid = character(),
              group = character(),
              Geneid = character(),
              sample_size_1 = integer(),
              sample_size_2 = integer(),
              theta_1 = double(),
              theta_2 = double(),
              mu_1 = double(),
              mu_2 = double(),
              size_1 = double(),
              size_2 = double(),
              prob_1 = double(),
              prob_2 = double(),
              total_mean_1 = double(),
              total_mean_2 = double(),
              foldChange = double(),
              chi = double(),
              pvalue = double(),
              adjusted_pvalue = double(),
              Remark = character(),
              stringsAsFactors = FALSE)

            # gene ---------------------------------------------------------------
            for (gene in genes){
              gene.cnt <- gene.cnt + 1
              # gene expression for ref group
              counts_1 <- unlist(expressionMatrix_split[gene, ref_cells])
              # gene expression for alt group
              counts_2 <- unlist(expressionMatrix_split[gene, alt_cells])
              results_gene <- data.frame(group = j,
                                         SNPid = snpid,
                                         Geneid = gene,
                                         sample_size_1 = length(counts_1),
                                         sample_size_2 = length(counts_2),
                                         theta_1 = NA,
                                         theta_2 = NA,
                                         mu_1 = NA,
                                         mu_2 = NA,
                                         size_1 = NA,
                                         size_2 = NA,
                                         prob_1 = NA,
                                         prob_2 = NA,
                                         total_mean_1 = NA,
                                         total_mean_2 = NA,
                                         foldChange = NA,
                                         chi = NA,
                                         pvalue = NA,
                                         adjusted_pvalue = NA,
                                         Remark = NA,
                                         stringsAsFactors = FALSE)

              # calculate fold change --------------------------------------------
              totalMean_1 <- mean(counts_1)
              totalMean_2 <- mean(counts_2)
              foldChange  <- totalMean_1/totalMean_2

              # filling data -----------------------------------------------------
              results_gene[1,"total_mean_1"] <- totalMean_1
              results_gene[1,"total_mean_2"] <- totalMean_2
              results_gene[1,"foldChange"]   <- foldChange


              # build  model -----------------------------------------------------

              build_model <- function(counts) {

                options(show.error.messages = FALSE)
                zinb_gamlssML <- try(gamlssML(counts, family = "ZINBI"),
                                     silent = TRUE)
                options(show.error.messages = TRUE)

                if('try-error' %in% class(zinb_gamlssML)){
                  print("MLE of ZINB failed! Please choose another model")
                  return(list(theta = NA, mu = NA, size = NA, prob = NA))
                }else{

                  zinb <- zinb_gamlssML
                  theta <- zinb$nu
                  mu <- zinb$mu
                  size <- 1 / zinb$sigma
                  prob <- size / (size + mu)
                }
                return(list(theta = theta, mu = mu, size = size, prob = prob))
              }

              build_model(counts_1)
              build_model(counts_2)


              # estimate Ref params ----------------------------------------------
              params_1 <- build_model(counts_1)
              theta_1 <- params_1[['theta']]
              mu_1 <- params_1[['mu']]
              size_1 <- params_1[['size']]
              prob_1 <- params_1[['prob']]

              #estimate Alt params -----------------------------------------------
              params_2 <- build_model(counts_2)
              theta_2 <- params_2[['theta']]
              mu_2 <- params_2[['mu']]
              size_2 <- params_2[['size']]
              prob_2 <- params_2[['prob']]

              #combine Ref and Alt -----------------------------------------------
              params_combined <- build_model(c(counts_1, counts_2))
              theta_res <- params_combined[['theta']]
              mu_res <- params_combined[['mu']]
              size_res <- params_combined[['size']]
              prob_res <- params_combined[['prob']]


              # calculate p-value ---------------------------------------
              # 原假设是，两个模型拟合数据的效果没有显著差别，
              # 也就是ref组和alt组的基因表达量没有显著性差异
              # 计算两组数据的对数似然之和
              logL <- function(counts_1,
                               theta_1,
                               size_1,
                               prob_1,
                               counts_2,
                               theta_2,
                               size_2,
                               prob_2){
                # log-likelihood for count1 under parameter
                logL_1 <- sum(dzinegbin(counts_1,
                                        size = size_1,
                                        prob = prob_1,
                                        pstr0 = theta_1,
                                        log = TRUE))
                # log-likelihood for count2 under parameter
                logL_2 <- sum(dzinegbin(counts_2,
                                        size = size_2,
                                        prob = prob_2,
                                        pstr0 = theta_2,
                                        log = TRUE))
                logL <- logL_1 + logL_2
                logL
              }
              # 基于两个不同模型的似然值计算
              logL_A <- logL(counts_1,
                             theta_1,
                             size_1,
                             prob_1,
                             counts_2,
                             theta_2,
                             size_2,
                             prob_2)
              # 基于合并模型的似然值计算
              logL_B <- logL(counts_1,
                             theta_res,
                             size_res,
                             prob_res,
                             counts_2,
                             theta_res,
                             size_res,
                             prob_res)
              # 假设检验计算卡方
              chi <- logL_A - logL_B
              # 进而计算p值
              pvalue <- 1 - pchisq(2 * chi , df = 3)


              # filling data into data frame ------------------------------

              results_gene[1,"theta_1"] <- theta_1
              results_gene[1,"theta_2"] <- theta_2
              results_gene[1,"mu_1"] <- mu_1
              results_gene[1,"mu_2"] <- mu_2
              results_gene[1,"size_1"] <- size_1
              results_gene[1,"size_2"] <- size_2
              results_gene[1,"prob_1"] <- prob_1
              results_gene[1,"prob_2"] <- prob_2
              results_gene[1,"chi"] <- chi
              results_gene[1,"pvalue"] <- pvalue

              results_snp <- rbind(results_snp, results_gene)
            }

            # return
            return(results_snp)

          }else{
            results_snp <- data.frame()
          }
        }


        # final result
        result <- data.frame(
          SNPid = character(),
          group = character(),
          Geneid = character(),
          sample_size_1 = integer(),
          sample_size_2 = integer(),
          theta_1 = double(),
          theta_2 = double(),
          mu_1 = double(),
          mu_2 = double(),
          size_1 = double(),
          size_2 = double(),
          prob_1 = double(),
          prob_2 = double(),
          total_mean_1 = double(),
          total_mean_2 = double(),
          foldChange = double(),
          chi = double(),
          pvalue = double(),
          adjusted_pvalue = double(),
          Remark = character(),
          stringsAsFactors=FALSE)

        split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == j]
        expressionMatrix_split <- expressionMatrix[, split_cells]
        snpMatrix_split <- snpMatrix[, split_cells]

        message(j,':')
        message('0%   10   20   30   40   50   60   70   80   90   100%' )
        message('[----|----|----|----|----|----|----|----|----|----|')
        pb <- progress_bar$new(total = length(snplist),
                               format = "[:bar]",
                               clear = FALSE,
                               width = 51)
        for(i in 1:length(snplist)){
          result <- rbind(result, eQTLcalling(i))
          pb$tick()
        }

        #change column names
        colnames(result)[colnames(result) == 'sample_size_1'] <- 'sample_size_Ref'
        colnames(result)[colnames(result) == 'sample_size_2'] <- 'sample_size_Alt'
        colnames(result)[colnames(result) == 'theta_1'] <- 'theta_Ref'
        colnames(result)[colnames(result) == 'theta_2'] <- 'theta_Alt'
        colnames(result)[colnames(result) == 'mu_1'] <- 'mu_Ref'
        colnames(result)[colnames(result) == 'mu_2'] <- 'mu_Alt'
        colnames(result)[colnames(result) == 'size_1'] <- 'size_Ref'
        colnames(result)[colnames(result) == 'size_2'] <- 'size_Alt'
        colnames(result)[colnames(result) == 'prob_1'] <- 'prob_Ref'
        colnames(result)[colnames(result) == 'prob_2'] <- 'prob_Alt'
        colnames(result)[colnames(result) == 'total_mean_1'] <- 'total_mean_Ref'
        colnames(result)[colnames(result) == 'total_mean_2'] <- 'total_mean_Alt'
      }else if(biClassify == FALSE){
        eQTLcalling <- function(i){
          if (i %% 100 == 0) {
            gc()
          }

          snpid <- snplist[i]
          AA_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid,] == 1]
          Aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid,] == 3]
          aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[snpid,] == 2]

          if(length(AA_cells) > 0 && length(Aa_cells) > 0 && length(aa_cells) > 0){

            genes <- genelist
            gene.cnt <- 0

            #snp
            results_snp <- data.frame(SNPid = character(),
                                      group = character(),
                                      Geneid = character(),
                                      sample_size_1 = integer(),
                                      sample_size_2 = integer(),
                                      sample_size_3 = integer(),
                                      theta_1 = double(),
                                      theta_2 = double(),
                                      theta_3 = double(),
                                      mu_1 = double(),
                                      mu_2 = double(),
                                      mu_3 = double(),
                                      size_1 = double(),
                                      size_2 = double(),
                                      size_3 = double(),
                                      prob_1 = double(),
                                      prob_2 = double(),
                                      prob_3 = double(),
                                      total_mean_1 = double(),
                                      total_mean_2 = double(),
                                      total_mean_3 = double(),
                                      chi = double(),
                                      pvalue = double(),
                                      adjusted_pvalue = double(),
                                      Remark = character(),
                                      stringsAsFactors=FALSE)

            #gene
            for (gene in genes){
              gene.cnt <- gene.cnt + 1
              counts_1 <- unlist(expressionMatrix_split[gene, AA_cells])
              counts_2 <- unlist(expressionMatrix_split[gene, Aa_cells])
              counts_3 <- unlist(expressionMatrix_split[gene, aa_cells])
              results_gene <- data.frame(group = j,
                                         SNPid = snpid,
                                         Geneid = gene,
                                         sample_size_1 = length(counts_1),
                                         sample_size_2 = length(counts_2),
                                         sample_size_3 = length(counts_3),
                                         theta_1 = NA,
                                         theta_2 = NA,
                                         theta_3 = NA,
                                         mu_1 = NA,
                                         mu_2 = NA,
                                         mu_3 = NA,
                                         size_1 = NA,
                                         size_2 = NA,
                                         size_3 = NA,
                                         prob_1 = NA,
                                         prob_2 = NA,
                                         prob_3 = NA,
                                         total_mean_1 = NA,
                                         total_mean_2 = NA,
                                         total_mean_3 = NA,
                                         chi = NA,
                                         pvalue = NA,
                                         adjusted_pvalue = NA,
                                         Remark = NA,
                                         stringsAsFactors=FALSE)

              #calculate mean
              totalMean_1 <- mean(counts_1)
              totalMean_2 <- mean(counts_2)
              totalMean_3 <- mean(counts_3)

              #filling data
              results_gene[1,"total_mean_1"] <- totalMean_1
              results_gene[1,"total_mean_2"] <- totalMean_2
              results_gene[1,"total_mean_3"] <- totalMean_3

              #build ZINB model
              zinb <- function(counts) {
                if (sum(counts == 0) == length(counts)) {
                  theta <- 1
                  mu <- 0
                  size <- 1
                  prob <- size / (size + mu)
                } else {
                  options(show.error.messages = FALSE)
                  zinb_gamlssML <- try(gamlssML(counts,
                                                family = "ZINBI"),
                                       silent = TRUE)
                  options(show.error.messages = TRUE)
                  if ('try-error' %in% class(zinb_gamlssML)) {
                    zinb_zeroinfl <- try(zeroinfl(formula = counts ~ 1 | 1,
                                                  dist = "negbin"),
                                         silent = TRUE)
                    if ('try-error' %in% class(zinb_zeroinfl)) {
                      print("MLE of ZINB failed!")
                      return(list(theta = NA, mu = NA, size = NA, prob = NA))
                    } else {
                      zinb <- zinb_zeroinfl
                      theta <- plogis(zinb$coefficients$zero)
                      mu <- exp(zinb$coefficients$count)
                      size <- zinb$theta
                      prob <- size / (size + mu)
                    }
                  } else {
                    # 如果第一次建模没有失败，就从建模结果中得到参数
                    zinb <- zinb_gamlssML
                    theta <- zinb$nu
                    mu <- zinb$mu
                    size <- 1 / zinb$sigma
                    prob <- size / (size + mu)
                  }
                }
                return(list(theta = theta, mu = mu, size = size, prob = prob))
              }

              zinb(counts_1)
              zinb(counts_2)
              zinb(counts_3)

              #estimate AA params
              params_1 <- zinb(counts_1)
              theta_1 <- params_1[['theta']]
              mu_1 <- params_1[['mu']]
              size_1 <- params_1[['size']]
              prob_1 <- params_1[['prob']]

              #estimate Aa params
              params_2 <- zinb(counts_2)
              theta_2 <- params_2[['theta']]
              mu_2 <- params_2[['mu']]
              size_2 <- params_2[['size']]
              prob_2 <- params_2[['prob']]

              #estimate aa params
              params_3 <- zinb(counts_3)
              theta_3 <- params_3[['theta']]
              mu_3 <- params_3[['mu']]
              size_3 <- params_3[['size']]
              prob_3 <- params_3[['prob']]

              #combine AA,Aa and aa
              params_combined <- zinb(c(counts_1, counts_2, counts_3))
              theta_res <- params_combined[['theta']]
              mu_res <- params_combined[['mu']]
              size_res <- params_combined[['size']]
              prob_res <- params_combined[['prob']]

              #calculate p-value
              logL <- function(counts_1,
                               theta_1,
                               size_1,
                               prob_1,
                               counts_2,
                               theta_2,
                               size_2,
                               prob_2,
                               counts_3,
                               theta_3,
                               size_3,
                               prob_3){
                logL_1 <- sum(dzinegbin(counts_1, size = size_1, prob = prob_1,
                                        pstr0 = theta_1, log = TRUE))
                logL_2 <- sum(dzinegbin(counts_2, size = size_2, prob = prob_2,
                                        pstr0 = theta_2, log = TRUE))
                logL_3 <- sum(dzinegbin(counts_3, size = size_3, prob = prob_3,
                                        pstr0 = theta_3, log = TRUE))
                logL <- logL_1 + logL_2 + logL_3
                logL
              }

              logL_1 <- logL(counts_1, theta_1, size_1, prob_1, counts_2, theta_2,
                             size_2, prob_2, counts_3, theta_3, size_3, prob_3)
              logL_2 <- logL(counts_1, theta_res, size_res, prob_res, counts_2,
                             theta_res, size_res, prob_res, counts_3, theta_res,
                             size_res, prob_res)
              chi <- logL_1 - logL_2

              pvalue <- try(1 - pchisq(2 * chi , df = 3))
              if (class(pvalue) == "try-error") return(NA)

              #filling data into data frame
              results_gene[1,"theta_1"] <- theta_1
              results_gene[1,"theta_2"] <- theta_2
              results_gene[1,"theta_3"] <- theta_3
              results_gene[1,"mu_1"] <- mu_1
              results_gene[1,"mu_2"] <- mu_2
              results_gene[1,"mu_3"] <- mu_3
              results_gene[1,"size_1"] <- size_1
              results_gene[1,"size_2"] <- size_2
              results_gene[1,"size_3"] <- size_3
              results_gene[1,"prob_1"] <- prob_1
              results_gene[1,"prob_2"] <- prob_2
              results_gene[1,"prob_3"] <- prob_3
              results_gene[1,"chi"] <- chi
              results_gene[1,"pvalue"] <- pvalue

              results_snp <- rbind(results_snp, results_gene)
            }

            #return
            return(results_snp)

          }else{
            results_snp <- data.frame()
          }
        }

        #final result
        result <- data.frame(
          group = character(),
          SNPid = character(),
          Geneid = character(),
          sample_size_1 = integer(),
          sample_size_2 = integer(),
          sample_size_3 = integer(),
          theta_1 = double(),
          theta_2 = double(),
          theta_3 = double(),
          mu_1 = double(),
          mu_2 = double(),
          mu_3 = double(),
          size_1 = double(),
          size_2 = double(),
          size_3 = double(),
          prob_1 = double(),
          prob_2 = double(),
          prob_3 = double(),
          total_mean_1 = double(),
          total_mean_2 = double(),
          total_mean_3 = double(),
          chi = double(),
          pvalue = double(),
          adjusted_pvalue = double(),
          Remark = character(),
          stringsAsFactors=FALSE)

        split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == j]
        expressionMatrix_split <- expressionMatrix[, split_cells]
        snpMatrix_split <- snpMatrix[, split_cells]

        message(j,':')
        message('0%   10   20   30   40   50   60   70   80   90   100%' )
        message('[----|----|----|----|----|----|----|----|----|----|')
        pb <- progress_bar$new(total = length(snplist),
                               format = "[:bar]",
                               clear = FALSE,
                               width = 51)
        for(i in 1:length(snplist)){
          result <- rbind(result, eQTLcalling(i))
          pb$tick()
        }

        message("finished!")

        #change column names
        colnames(result)[colnames(result) == 'sample_size_1'] <- 'sample_size_AA'
        colnames(result)[colnames(result) == 'sample_size_2'] <- 'sample_size_Aa'
        colnames(result)[colnames(result) == 'sample_size_3'] <- 'sample_size_aa'
        colnames(result)[colnames(result) == 'theta_1'] <- 'theta_AA'
        colnames(result)[colnames(result) == 'theta_2'] <- 'theta_Aa'
        colnames(result)[colnames(result) == 'theta_3'] <- 'theta_aa'
        colnames(result)[colnames(result) == 'mu_1'] <- 'mu_AA'
        colnames(result)[colnames(result) == 'mu_2'] <- 'mu_Aa'
        colnames(result)[colnames(result) == 'mu_3'] <- 'mu_aa'
        colnames(result)[colnames(result) == 'size_1'] <- 'size_AA'
        colnames(result)[colnames(result) == 'size_2'] <- 'size_Aa'
        colnames(result)[colnames(result) == 'size_3'] <- 'size_aa'
        colnames(result)[colnames(result) == 'prob_1'] <- 'prob_AA'
        colnames(result)[colnames(result) == 'prob_2'] <- 'prob_Aa'
        colnames(result)[colnames(result) == 'prob_3'] <- 'prob_aa'
        colnames(result)[colnames(result) == 'total_mean_1'] <- 'total_mean_AA'
        colnames(result)[colnames(result) == 'total_mean_2'] <- 'total_mean_Aa'
        colnames(result)[colnames(result) == 'total_mean_3'] <- 'total_mean_aa'
      }else{
        stop("biClassify can only be selected as 'TRUE' or 'FALSE'")
      }


      # adjust p-value
      result[,"adjusted_pvalue"] <- p.adjust(result[,"pvalue"],
                                             method = p.adjust.method)
      result <- result[order(result[,"adjusted_pvalue"]),]
      rownames(result) <- NULL
      result <- result[result$adjusted_pvalue <= p.adjust.Threshold, ]
      result_all <- rbind(result_all, result)
    }
    return(result_all)

  }


  linearModel <- function(
    eQTLObject,
    genelist,
    snplist,
    biClassify = FALSE,
    p.adjust.method = "bonferroni",
    p.adjust.Threshold = 0.05,
    logfc.threshold = 0.1){

    expressionMatrix <- eQTLObject@filterData$expMat
    snpMatrix <- eQTLObject@filterData$snpMat
    unique_group <- unique(eQTLObject@groupBy$group)

    result_all <- data.frame()

    message("Start calling eQTL")
    for(k in unique_group){

      result <- data.frame(
        SNPid = character(),
        group = character(),
        Geneid = character(),
        pvalue = double(),
        adjusted_pvalue = double(),
        b = double(),
        abs_b = double(),
        Remark = character(),
        stringsAsFactors=FALSE)

      split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == k]
      expressionMatrix_split <- expressionMatrix[, split_cells]
      snpMatrix_split <- snpMatrix[, split_cells]

      if(biClassify == FALSE){

        message(k,':')
        message('0%   10   20   30   40   50   60   70   80   90   100%' )
        message('[----|----|----|----|----|----|----|----|----|----|')
        pb <- progress_bar$new(total = length(snplist),
                               format = "[:bar]",
                               clear = FALSE,
                               width = 51)
        for(i in 1:length(snplist)){
          snpid <- snplist[i]
          snp_mat <- snpMatrix_split[snpid, ]
          snp_mat <- as.data.frame(snp_mat)
          snp_mat$cells = rownames(snp_mat)

          replace_2_and_3 <- function(x) {
            ifelse(x == 2, 3, ifelse(x == 3, 2, x))
          }

          snp_mat_new <- snp_mat %>%
            mutate_all(funs(replace_2_and_3))

          genes <- genelist

          for(j in 1:length(genes)){
            gene_id <- genes[j]
            gene_mat <- expressionMatrix_split[gene_id, ]
            gene_mat <- as.data.frame(gene_mat)
            gene_mat$cells = rownames(gene_mat)

            combined_df <- merge(snp_mat, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 5)

            lmodel = lm(gene_mat ~ snp_mat, data = combined_df);

            lmout_pvalue = summary(lmodel)$coefficients[2, "Pr(>|t|)"]
            lmout_b = summary(lmodel)$coefficients[2, "Estimate"]

            new_row <- data.frame(SNPid = snpid,
                                  group = k,
                                  Geneid = genes[j],
                                  pvalue = lmout_pvalue,
                                  b = lmout_b)
            result <- rbind(result, new_row)
          }
          pb$tick()
        }
        message("finished!")

      }else if(biClassify == TRUE){

        snpMatrix_split[snpMatrix_split == 3] <- 2

        message(k,':')
        message('0%   10   20   30   40   50   60   70   80   90   100%' )
        message('[----|----|----|----|----|----|----|----|----|----|')
        pb <- progress_bar$new(total = length(snplist),
                               format = "[:bar]",
                               clear = FALSE,
                               width = 51)
        for(i in 1:length(snplist)){
          snpid <- snplist[i]
          snp_mat <- snpMatrix_split[snpid, ]
          snp_mat <- as.data.frame(snp_mat)
          snp_mat$cells = rownames(snp_mat)

          genes <- genelist

          for(j in 1:length(genes)){
            gene_id <- genes[j]
            gene_mat <- expressionMatrix_split[gene_id, ]
            gene_mat <- as.data.frame(gene_mat)
            gene_mat$cells = rownames(gene_mat)

            combined_df <- merge(snp_mat, gene_mat, by = "cells")
            combined_df <- subset(combined_df, snp_mat != 5)

            lmodel = lm(gene_mat ~ snp_mat, data = combined_df);

            lmout_pvalue = summary(lmodel)$coefficients[2, "Pr(>|t|)"]
            lmout_b = summary(lmodel)$coefficients[2, "Estimate"]

            new_row <- data.frame(SNPid = snpid,
                                  group = k,
                                  Geneid = genes[j],
                                  pvalue = lmout,
                                  b = lmout_b)

            result <- rbind(result, new_row)
          }
          pb$tick()
        }
        message("finished!")

      }else{
        stop("biClassify can only be selected as 'TRUE' or 'FALSE'")
      }


      if (!p.adjust.method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
        stop("Invalid p-adjusted method.
           Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel', or'fdr or BH'.")
      }

      # adjust p-value
      result[,"adjusted_pvalue"] <- p.adjust(result[,"pvalue"], method = "BH")
      result <- result[order(result[,"adjusted_pvalue"]),]
      rownames(result) <- NULL
      result <- result[result$adjusted_pvalue <= p.adjust.Threshold, ]

      # abs_b
      result <- result %>%
        mutate(abs_b = abs(b))

      result <- result[result$abs_b >= logfc.threshold, ]

      result_all <- rbind(result_all, result)

    }
    return(result_all)
  }


  if(useModel == "zinb"){
    result <- zinbModel(eQTLObject = eQTLObject,
                        genelist = matched_gene,
                        snplist = matched_snps,
                        biClassify = biClassify,
                        p.adjust.method = p.adjust.method,
                        p.adjust.Threshold = p.adjust.Threshold)
  }else if(useModel == "possion"){
    result <- possionModel(eQTLObject = eQTLObject,
                           genelist = matched_gene,
                           snplist = matched_snps,
                           biClassify = biClassify,
                           p.adjust.method = p.adjust.method,
                           p.adjust.Threshold = p.adjust.Threshold)
  }else if(useModel == "linear"){
    result <- linearModel(eQTLObject = eQTLObject,
                          genelist = matched_gene,
                          snplist = matched_snps,
                          biClassify = biClassify,
                          p.adjust.method = p.adjust.method,
                          p.adjust.Threshold = p.adjust.Threshold,
                          logfc.threshold = logfc.threshold)
  }else{
    stop("Invalid model Please choose from 'zinb', 'possion' , or 'linear'.")
  }

  eQTLObject@eQTLResult <- result
  return(eQTLObject)
}
