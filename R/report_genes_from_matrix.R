#' Report Genes from Matrix
#'
#' This function identifies underlying genes for a specific function (gene set) and sample of interest
#' from a preprocessed matrix, instead of from a Seurat object. It calculates residuals for each
#' gene, normalized by the reference leave-one-out (LOO) distribution, and returns a dataframe of
#' genes and their respective "normalized residuals" above the expression threshold.
#'
#' @importFrom stats lm pchisq sd p.adjust
#'
#' @param datH Reference matrix (rows = samples, columns = genes) for healthy/control data.
#' @param datD Query matrix (rows = samples, columns = genes), which should include the
#'   sample_of_interest for driver gene identification.
#' @param go_sets Matrix or data frame indicating membership of genes in gene sets. Rows should be
#'   gene names, and columns should be specific gene sets. A value of 1 indicates membership
#'   in that gene set.
#' @param gene_set_name Name of the gene set (function/pathway) of interest. This should match one
#'   of the column names of the gene set database.
#' @param sample_of_interest Specific sample (row in datD) to analyze.
#' @param expression_threshold Threshold to filter genes of low expression. Genes must exceed this
#'   average expression level in both datH and datD.
#'
#' @return A dataframe with genes and their normalized residuals.
#'
#' #' @examples
#' # gene_set_name = "Postsynaptic Membrane Organization (GO:0001941)"
#' # sample_of_interest = "sample_1"
#' #
#' # driver_genes_from_matrix <- report_genes_from_matrix(
#' #   datH = reference_matrix,
#' #   datD = query_matrix,
#' #   go_sets = go_sets,
#' #   gene_set_name = gene_set_name,
#' #   sample_of_interest = sample_of_interest,
#' #   expression_threshold = 0
#' # )
#'
#' @export
report_genes_from_matrix <- function(datH, datD, go_sets, gene_set_name,
                                     sample_of_interest,
                                     expression_threshold = 0) {
  CHISQ.MAX = 4
  # Subset datD to the sample of interest
  datD = datD[sample_of_interest, , drop = FALSE]

  # Select genes associated with the specified gene set
  func = gene_set_name
  selected_genes <- rownames(go_sets)[go_sets[, func] == 1]
  common_genes <- selected_genes[selected_genes %in% colnames(datH)]

  if (length(common_genes) < length(selected_genes)) {
    warning(paste("Not all genes found for", func, "- proceeding with",
                  length(common_genes), "genes"))
  }

  # Subset both matrices to the selected genes
  dat_subsetH <- datH[, common_genes, drop = FALSE]
  dat_subsetD <- datD[, common_genes, drop = FALSE]

  # Filter low-expressed genes by average expression across rows
  avg_expression_H <- colMeans(dat_subsetH)
  avg_expression_D <- colMeans(dat_subsetD)
  filtered_genes <- names(avg_expression_H)[
    avg_expression_H > expression_threshold & avg_expression_D > expression_threshold
  ]

  dat_subsetH <- dat_subsetH[, filtered_genes, drop = FALSE]
  dat_subsetD <- dat_subsetD[, filtered_genes, drop = FALSE]

  # Leave-One-Out (LOO) approach on reference data
  nH_subset <- nrow(dat_subsetH)
  svdepsilon_H_subset <- matrix(nrow = nH_subset, ncol = ncol(dat_subsetH), data = 0)

  for (i in seq_len(nH_subset)) {
    datH_loo <- dat_subsetH[-i, ]
    svdH_loo <- svd(datH_loo)
    V1_loo   <- svdH_loo$v[, 1]

    lmres <- lm(as.numeric(dat_subsetH[i, ]) ~ V1_loo)
    svdepsilon_H_subset[i, ] <- lmres$residuals
  }

  # Compute standard deviation of gene residuals from LOO
  svdsigmag <- apply(svdepsilon_H_subset, 2, sd)

  # Perform SVD on full reference set, then compute residuals for the sample_of_interest
  svdH <- svd(dat_subsetH)
  V1_H <- svdH$v[, 1]

  lmres_query <- lm(as.numeric(dat_subsetD[1, ]) ~ V1_H)
  normalized_residuals <- lmres_query$residuals / svdsigmag

  # Compile results
  driver_genes <- data.frame(
    Gene = filtered_genes,
    Normalized_Residual = normalized_residuals,
    stringsAsFactors = FALSE
  )

  return(driver_genes)
}
