#' Run Spectra Analysis from Matrix Data
#'
#' Performs spectra analysis on given matrices query (datD) samples compared to a reference (datH).
#' It calculates R² and p-values for each gene set and saves the results as CSV files.
#'
#' @param datH A normalized (e.g. tpm + log) matrix of gene expression data for healthy samples, with genes as columns.
#' @param datD A normalized (e.g. tpm + log) matrix of gene expression data for diseased samples, similar to datH.
#' @param go_sets A dataframe or matrix where columns represent gene sets. Can be generated using process_gene_sets_for_matrix function.
#' @param output_folder The directory path where the result files will be saved.
#' @param coordination_key An optional key for file naming (default is 'nokey'). Could be cell type or condition.
#' @param CHISQ.MAX An optional maximum chi-square value for robustness (default is 4).
#' @param expression_threshold Genes with an average expression in the reference
#'        OR query matrix below this threshold are removed from the analysis.
#'        Default 0. Set to -1 for no filtering applied.
#' @param gene_number_threshold Number of genes in a gene set after filtering needed to perform analysis. If less genes remaining, results will be NA.
#'        Default 10 genes.
#'
#' @return This function does not return a value but saves two CSV files containing R² and p-values.
#'
#' @examples
#'
#' # create dummy object and gene sets
#' set.seed(123)
#' genes <- paste0("Gene", 1:10)
#' datH <- matrix(rexp(50, rate = 1), nrow = 10, ncol = 5)
#' rownames(datH) <- genes
#' colnames(datH) <- paste0("Healthy", 1:5)
#' datH <- t(datH)  # samples as rows
#' datD <- matrix(rexp(10, rate = 1), nrow = 10, ncol = 2)
#' rownames(datD) <- genes
#' colnames(datD) <- paste0("Disease", 1:2)
#' datD <- t(datD)  # samples as rows
#' go_sets <- matrix(sample(0:1, 10 * 5, replace = TRUE), nrow = 10, ncol = 5)
#' rownames(go_sets) <- genes
#' colnames(go_sets) <- paste0("Pathway", 1:5)
#'
#' # Create temporary output folder
#' output_folder <- tempfile("example_output_")
#' dir.create(output_folder)
#'
#' # Run function
#' run_spectra_from_matrix(
#'   datH = datH,
#'   datD = datD,
#'   go_sets = go_sets,
#'   output_folder = output_folder,
#'   coordination_key = "ExampleCellType",
#'   expression_threshold = 0,
#'   gene_number_threshold = 0
#' )
#' @importFrom stats lm pchisq p.adjust
#' @importFrom utils write.csv
#' @export

run_spectra_from_matrix <- function(datH, datD, go_sets, output_folder, coordination_key = "nokey", CHISQ.MAX = 4, expression_threshold = 0, gene_number_threshold = 10) {

  # Check and create R2 and Pval subfolders if they don't exist
  r2_folder <- file.path(output_folder, "/R2")
  pval_folder <- file.path(output_folder, "/Pval")
  padj_folder <- file.path(output_folder, "/Padj")

  if (!dir.exists(r2_folder)) {
    dir.create(r2_folder)
  }
  if (!dir.exists(pval_folder)) {
    dir.create(pval_folder)
  }
  if (!dir.exists(padj_folder)) {
    dir.create(padj_folder)
  }

    # Initialize dataframes for R² and p-values
    r2_df <- data.frame(matrix(nrow = nrow(datD), ncol = length(colnames(go_sets))))
    pval_df <- data.frame(matrix(nrow = nrow(datD), ncol = length(colnames(go_sets))))
    svdpvals_df <- data.frame(matrix(nrow = nrow(datD), ncol = length(colnames(go_sets))))
    svdchisq_df <- data.frame(matrix(nrow = nrow(datD), ncol = length(colnames(go_sets))))
    deviant.pvalD_df <- data.frame(matrix(nrow = nrow(datD), ncol = length(colnames(go_sets))))

    rownames(r2_df) <- rownames(datD)
    rownames(pval_df) <- rownames(datD)
    colnames(r2_df) <- colnames(go_sets)
    colnames(pval_df) <- colnames(go_sets)

    rownames(svdpvals_df) <- rownames(datD)
    rownames(svdchisq_df) <- rownames(datD)
    rownames(deviant.pvalD_df) <- rownames(datD)
    colnames(svdpvals_df) <- colnames(go_sets)
    colnames(svdchisq_df) <- colnames(go_sets)
    colnames(deviant.pvalD_df) <- colnames(go_sets)

    # Loop through each gene set
    for (func in colnames(go_sets)) {
      # Select genes associated with the current function
      selected_genes <- rownames(go_sets)[go_sets[, func] == 1]
      common_genes <- selected_genes[selected_genes %in% colnames(datH)]

      if (length(common_genes) < length(selected_genes)) {
        warning(paste("Not all genes found for", func, "- proceeding with", length(common_genes), "genes"))
      }

      if (length(common_genes) < gene_number_threshold) {

        warning(paste("Less genes than threshold in ", func, " - skipping gene set"))
        svdpvals_df[, func] <- NA
        svdchisq_df[, func] <- NA
        r2_df[, func] <- NA
        deviant.pvalD_df[, func] <- NA
        next  # Skip to the next iteration of the loop
      }

      dat_subsetH <- datH[, common_genes, drop = FALSE]
      dat_subsetD <- datD[, common_genes, drop = FALSE]

      #filtering low expressed genes
      # Calculate the average expression for each gene in dat_subsetH
      avg_expression_H <- colMeans(dat_subsetH)

      # Calculate the average expression for each gene in dat_subsetD
      avg_expression_D <- colMeans(dat_subsetD)

      # Find genes that are above the threshold in both dat_subsetH and dat_subsetD
      filtered_genes <- names(avg_expression_H)[avg_expression_H > expression_threshold & avg_expression_D > expression_threshold]

      if (length(filtered_genes) < gene_number_threshold) {

        warning(paste("Less genes than threshold in ", func, " - skipping gene set"))
        svdpvals_df[, func] <- NA
        svdchisq_df[, func] <- NA
        r2_df[, func] <- NA
        deviant.pvalD_df[, func] <- NA
        next  # Skip to the next iteration of the loop
      }

      # Subset dat_subsetH and dat_subsetD to keep only the filtered genes
      dat_subsetH <- dat_subsetH[, filtered_genes, drop = FALSE]
      dat_subsetD <- dat_subsetD[, filtered_genes, drop = FALSE]

      # Calculate R^2 for the healthy subset with LOO
      nH_subset <- nrow(dat_subsetH)
      r2_H_subset <- rep(0, nH_subset)
      svdpvals_H_subset <- rep(0, nH_subset)
      svdepsilon_H_subset = matrix(nrow=nrow(dat_subsetH), ncol=ncol(dat_subsetH), data=0)
      for (i in 1:nH_subset) {
        datH_loo <- dat_subsetH[-i,]
        svdH_loo <- svd(datH_loo)
        V1_loo <- svdH_loo$v[,1]
        lmres <- lm(as.numeric(dat_subsetH[i,]) ~ V1_loo)
        svdepsilon_H_subset[i,] = lmres$residuals
        svdpvals_H_subset[i] = summary(lmres)$coefficients[2,4]
        r2_H_subset[i] <- summary(lmres)$r.squared
      }

      svdsigmag = apply(svdepsilon_H_subset, 2, sd)
      svdH <- svd(dat_subsetH)
      V1_H <- svdH$v[,1]
      # Calculate R^2 for the diseased subset
      nD_subset <- nrow(dat_subsetD)
      r2_D_subset <- rep(0, nD_subset)
      # Regress each disease sample on V1
      #svdepsilonD = matrix(nrow=nrow(dat_subsetD), ncol=ncol(dat_subsetD), data=0)
      svdpvalsD = rep(NA, nD_subset)
      svdchisqD = rep(NA, nD_subset)  # This records the chi-square statistic assessing deviation of sample from normal reference.
      deviant.pvalD = rep(NA, nD_subset)
      chisq.df = length(svdsigmag)-1 # degree of freedom of chi-sq (number of genes)
      for (i in 1:nD_subset) {
        lmres <- lm(as.numeric(dat_subsetD[i,]) ~ V1_H)
        #svdepsilonD[i,] = lmres$residuals
        svdpvalsD[i] = summary(lmres)$coefficients[2,4]
        svdchisqD[i] = sum(pmin((lmres$residuals/svdsigmag)^2, CHISQ.MAX))
        r2_D_subset[i] <- summary(lmres)$r.squared
        deviant.pvalD[i]=pchisq(svdchisqD[i], df=chisq.df, lower.tail=FALSE)
        # Adjust R² values based on p-values

      }
      # Store results in respective dataframes
      svdpvals_df[, func] <- svdpvalsD
      svdchisq_df[, func] <- svdchisqD
      r2_df[, func] <- r2_D_subset
      deviant.pvalD_df[, func] <- deviant.pvalD

    }
    # Apply FDR correction to p-values for each gene set
    pval_df <- apply(deviant.pvalD_df, 1, function(p) p.adjust(p, method = "fdr"))
    pval_df <- t(pval_df)
    # Save R² and p-values to CSV files for each cell type
    write.csv(r2_df, paste0(output_folder, "/R2/R2_", coordination_key, ".csv"), row.names = TRUE)
    write.csv(pval_df, paste0(output_folder, "/Padj/Pval_", coordination_key, ".csv"), row.names = TRUE)
    write.csv(deviant.pvalD_df, paste0(output_folder, "/Pval/Pval_", coordination_key, ".csv"), row.names = TRUE)

}


