#' Run Spectra Analysis
#'
#' This function performs the actual scSpectra analysis on our created references from a seurat object, calculates R² and p-values for each gene set, and saves the results in
#' specified subfolders ('R2' and 'Pval') within the output directory. Only needs the output folder from the create_references function and the cell types to analyze.
#'
#' @importFrom stats lm pchisq sd p.adjust
#' @importFrom utils write.csv
#'
#' @param output_folder_base The base directory for saving output files. The default is set to the
#'        value of the 'output_folder_base' variable defined outside the function.
#' @param cell_types Vector of cell types to be processed, defaulting to the 'cell_types' variable
#'        defined outside the function.
#' @param CHISQ.MAX Maximum value for chi-square statistic used to cut-off the contribution of an individual gene of a pathway, with a default
#'        value of 4.
#' @param expression_threshold Genes with an average expression in the reference
#'        OR query matrix below this threshold are removed from the analysis.
#'        Default 0.
#' @param gene_number_threshold Number of genes in a gene set after filtering needed to perform analysis. If less genes remaining, results will be NA.
#'        Default 10 genes.
#' @param restrict_to_sample Set TRUE if you want to specify the samples to run
#'        scSpectra on. If FALSE, the tool will run analysis based on the
#'        valid_samples file created by the create_references function.
#'        Make sure to specify the vector_of_samples.
#'        Default: FALSE.
#' @param results_sample_key Needed if restrict_to_sample = TRUE. Will specify
#'        the samplekey to use when generating R2 and Pval results into the
#'        output folder. Default: "sampleofinterest".
#' @param vector_of_samples A vector of sample in the shape of
#'        c("sample1", "sample2"). Only needed if restrict_to_sample = TRUE.
#' @param QC_report If TRUE a QC report of the reference coordination for each
#'        sample / cell type and every pathway will be generated (including the
#'        mean, median, and sd of the reference LOO R2 distribution and the
#'        variance explained by V1 of the reference).
#'        Default: TRUE.
#'
#'
#' @return Does not return any value, but outputs files to the specified directory.
#'
#' @examples
#' # Assuming 'seurat_object' is defined and 'output_folder_base' and 'cell_types' are set:
#' #run_spectra(output_folder_base = output_folder_base,
#' #cell_types = cell_types, CHISQ.MAX = 4, expression_threshold = 0,
#' #gene_number_threshold = 10, restrict_to_sample = FALSE,
#' #vector_of_samples = c())
#'
#'
#' @export

run_spectra <- function(output_folder_base = output_folder_base, cell_types = cell_types, CHISQ.MAX = 4, expression_threshold = 0, gene_number_threshold = 10,
                        restrict_to_sample = FALSE, results_sample_key = "sampleofinterest", vector_of_samples = c(), QC_report = TRUE) {

  # Check and create R2 and Pval subfolders if they don't exist
  r2_folder <- file.path(output_folder_base, "R2")
  pval_folder <- file.path(output_folder_base, "Padj")
  pval_raw_folder <- file.path(output_folder_base, "Pval")

  if (!dir.exists(r2_folder)) {
    dir.create(r2_folder)
  }
  if (!dir.exists(pval_folder)) {
    dir.create(pval_folder)
  }
  if (!dir.exists(pval_raw_folder)) {
    dir.create(pval_raw_folder)
  }

  go_sets = readRDS(file = paste0(output_folder_base, "go_sets_modified.rds"))

  # Loop through each cell type
  for (cell_type in cell_types) {
    # Read in expression data for each cell type

    celltype_output_folder <- paste0(output_folder_base, cell_type, "/")

    if (restrict_to_sample == TRUE) {

      diseased_samples = vector_of_samples

    } else {

      valid_samples = readRDS(file = paste0(celltype_output_folder, "valid_samples.rds"))

      samples_of_interest = readRDS(file = paste0(output_folder_base, "samples_of_interest.rds"))

      diseased_samples <- intersect(valid_samples, samples_of_interest)

    }

    # Initialize dataframes for R² and p-values
    r2_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
    pval_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
    svdpvals_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
    svdchisq_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
    deviant.pvalD_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))

    rownames(r2_df) <- diseased_samples
    rownames(pval_df) <- diseased_samples
    colnames(r2_df) <- colnames(go_sets)
    colnames(pval_df) <- colnames(go_sets)

    rownames(svdpvals_df) <- diseased_samples
    rownames(svdchisq_df) <- diseased_samples
    rownames(deviant.pvalD_df) <- diseased_samples
    colnames(svdpvals_df) <- colnames(go_sets)
    colnames(svdchisq_df) <- colnames(go_sets)
    colnames(deviant.pvalD_df) <- colnames(go_sets)

    # Loop through each diseased sample
    for (sample in diseased_samples) {
      # Create pseudobulk for the query sample

      sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")

      datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
      datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))

      if (QC_report == TRUE) {

        #QC check df
        QCresults <- data.frame(
          Pathway = character(),
          R2mean = numeric(),
          R2median = numeric(),
          R2sd = numeric(),
          V1VE = numeric(),
          stringsAsFactors = FALSE
        )

      }

      # Loop through each gene set
      for (func in colnames(go_sets)) {
        # Select genes associated with the current function
        selected_genes <- rownames(go_sets)[go_sets[, func] == 1]
        common_genes <- selected_genes[selected_genes %in% colnames(datH)]

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
          svdpvals_df[sample, func] <- NA
          svdchisq_df[sample, func] <- NA
          r2_df[sample, func] <- NA
          deviant.pvalD_df[sample, func] <- NA
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
        chisq.df = length(svdsigmag)-1 # degree of freedom of chi-sq (number of genes) minus 1
        lmres <- lm(as.numeric(dat_subsetD) ~ V1_H)
        #svdepsilonD[sample,] = lmres$residuals
        svdpvalsD[sample] = summary(lmres)$coefficients[2,4]
        svdchisqD[sample] = sum(pmin((lmres$residuals/svdsigmag)^2, CHISQ.MAX))
        r2_D_subset[sample] <- summary(lmres)$r.squared
        deviant.pvalD[sample]=pchisq(svdchisqD[sample], df=chisq.df, lower.tail=FALSE)
        # Adjust R² values based on p-values

        svdpvals_df[sample, func] <- svdpvalsD[sample]
        svdchisq_df[sample, func] <- svdchisqD[sample]
        r2_df[sample, func] <- r2_D_subset[sample]
        deviant.pvalD_df[sample, func] <- deviant.pvalD[sample]

        if (QC_report == TRUE) {

        # Add the results as a new row in the data frame
        new_QCrow <- data.frame(
          Pathway = func,
          R2mean = mean(r2_H_subset),
          R2median = median(r2_H_subset),
          R2sd = sd(r2_H_subset),
          V1VE = svdH$d[1]^2 / sum(svdH$d^2),
          stringsAsFactors = FALSE
        )

        QCresults <- rbind(QCresults, new_QCrow)

        }

      }

      if (QC_report == TRUE) {
      write.csv(QCresults, paste0(sample_output_folder, "pathway_QC.csv"), row.names = FALSE)
      }

    }
    # Apply FDR correction to p-values for each gene set
    pval_df <- apply(deviant.pvalD_df, 1, function(p) p.adjust(p, method = "fdr"))
    pval_df <- t(pval_df)

    if (restrict_to_sample == TRUE) {

    # Save R² and p-values to CSV files in the respective subfolders
    r2_file_path <- file.path(r2_folder, paste0("R2_", results_sample_key, "_", cell_type, ".csv"))
    pval_file_path <- file.path(pval_folder, paste0("Padj_", results_sample_key, "_", cell_type, ".csv"))
    pval_file_raw_path <- file.path(pval_raw_folder, paste0("Pval_", results_sample_key, "_", cell_type, ".csv"))

    } else {

    # Save R² and p-values to CSV files in the respective subfolders
    r2_file_path <- file.path(r2_folder, paste0("R2_", cell_type, ".csv"))
    pval_file_path <- file.path(pval_folder, paste0("Padj_", cell_type, ".csv"))
    pval_file_raw_path <- file.path(pval_raw_folder, paste0("Pval_", cell_type, ".csv"))

    }

    write.csv(r2_df, r2_file_path, row.names = TRUE)
    write.csv(pval_df, pval_file_path, row.names = TRUE)
    write.csv(deviant.pvalD_df, pval_file_raw_path, row.names = TRUE)
  }

}
