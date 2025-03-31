#' Report Genes
#'
#' This function identifies driver genes for a specific cell type, function, and sample of interest.
#' It calculates residuals for each gene, normalized by the reference LOO distribution, and returns a
#' dataframe of genes and respective "normalized residuals" above the expression threshold
#'
#'
#' @importFrom stats lm pchisq sd p.adjust
#'
#' @param output_folder_base Directory for loading the reference and query data files.
#' @param gene_set_name Name of the gene set or pathway of interest.
#' @param cell_type Cell type to analyze.
#' @param sample_of_interest Specific sample to analyze.
#' @param expression_threshold Threshold to filter genes of low expression.
#'
#' @return A dataframe with genes and their normalized residuals.
#'
#' #' @examples
#' #gene_set_name = "Postsynaptic Membrane Organization (GO:0001941)"
#' #output_folder_base = "/.../output_human_extended/"
#' #drivinggenes <- report_genes(output_folder_base, gene_set_name, cell_type = "TAL",
#' #sample_of_interest = "sample")
#'
#' @export
report_genes <- function(output_folder_base, gene_set_name, cell_type, sample_of_interest, expression_threshold = 0) {

  CHISQ.MAX = 4

  cell_types = c(cell_type)

  soi = sample_of_interest

  gene_set_oi = gene_set_name

  celltype_output_folder <- paste0(output_folder_base, cell_type, "/")

  # Set your threshold for average expression
  threshold <- expression_threshold

  go_sets = readRDS(file = paste0(output_folder_base, "go_sets_modified.rds"))

  go_sets <- subset(go_sets, select = gene_set_oi)

  # Check if gene set is available
  if (is.null(go_sets) || ncol(go_sets) == 0) {
    stop("Specified gene set not found in the gene sets file.")
  }

  # Loop through each cell type
  for (cell_type in cell_types) {
    # Read in expression data for each cell type

    #####in the mod. version we specify only one sample
    diseased_samples = c(soi)

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

        # Assuming dat_subsetH and dat_subsetD are already defined
        # Calculate the average expression for each gene in dat_subsetH
        avg_expression_H <- colMeans(dat_subsetH)

        # Calculate the average expression for each gene in dat_subsetD
        avg_expression_D <- colMeans(dat_subsetD)

        # Find genes that are above the threshold in both dat_subsetH and dat_subsetD
        filtered_genes <- names(avg_expression_H)[avg_expression_H > threshold & avg_expression_D > threshold]


        if (length(filtered_genes) < length(common_genes)) {
          warning(paste("Genes below the threshold were removed for ", func, " - proceeding with", length(filtered_genes), "genes"))
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
        svdchisqD = rep(NA, nD_subset)
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

        normalized_residuals = lmres$residuals/svdsigmag

      }


    }

  }

  # Compile results into a dataframe
  driver_genes <- data.frame(
    Gene = filtered_genes,
    Normalized_Residual = normalized_residuals,
    stringsAsFactors = FALSE
  )

  return(driver_genes)
}
