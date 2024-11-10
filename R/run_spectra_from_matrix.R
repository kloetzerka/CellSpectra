#' Run Spectra Analysis from Matrix Data
#'
#' Performs spectra analysis on given matrices query (datD) samples compared to a reference (datH).
#' It calculates R² and p-values for each gene set and saves the results as CSV files.
#'
#' @param datH A normalized (e.g. tpm + log) matrix of gene expression data for healthy samples, with genes as rows.
#' @param datD A normalized (e.g. tpm + log) matrix of gene expression data for diseased samples, similar to datH.
#' @param go_sets A dataframe or matrix where columns represent gene sets. Can be generated using process_gene_sets_for_matrix function.
#' @param output_folder The directory path where the result files will be saved.
#' @param coordination_key An optional key for file naming (default is 'nokey'). Could be cell type or condition.
#' @param CHISQ.MAX An optional maximum chi-square value for robustness (default is 4).
#'
#' @return This function does not return a value but saves two CSV files containing R² and p-values.
#'
#' @examples
#' # datH and datD are matrices with the same genes but different samples
#' #(healthy/reference vs diseased/query)
#' #run_spectra_from_matrix(datH, datD, go_sets, "path/to/output", "key_example")
#'
#' @importFrom stats lm pchisq p.adjust
#' @importFrom utils write.csv
#' @export


run_spectra_from_matrix <- function(datH, datD, go_sets, output_folder, coordination_key = "nokey", CHISQ.MAX = 4) {

  # Check and create R2 and Pval subfolders if they don't exist
  r2_folder <- file.path(output_folder, "/R2")
  pval_folder <- file.path(output_folder, "/Pval")

  if (!dir.exists(r2_folder)) {
    dir.create(r2_folder)
  }
  if (!dir.exists(pval_folder)) {
    dir.create(pval_folder)
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

      dat_subsetH <- datH[, common_genes, drop = FALSE]
      dat_subsetD <- datD[, common_genes, drop = FALSE]

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
      svdepsilonD = matrix(nrow=nrow(dat_subsetD), ncol=ncol(dat_subsetD), data=0)
      svdpvalsD = rep(NA, nD_subset)
      svdchisqD = rep(NA, nD_subset)  # This records the chi-square statistic assessing deviation of sample from normal reference.
      deviant.pvalD = rep(NA, nD_subset)
      chisq.df = length(svdsigmag)-1 # degree of freedom of chi-sq (number of genes)
      for (i in 1:nD_subset) {
        lmres <- lm(as.numeric(dat_subsetD[i,]) ~ V1_H)
        svdepsilonD[i,] = lmres$residuals
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
    pval_df <- apply(deviant.pvalD_df, 2, function(p) p.adjust(p, method = "fdr"))

    # Save R² and p-values to CSV files for each cell type
    write.csv(r2_df, paste0(output_folder, "/R2/R2_", coordination_key, ".csv"), row.names = TRUE)
    write.csv(pval_df, paste0(output_folder, "/Pval/Pval_", coordination_key, ".csv"), row.names = TRUE)

}
