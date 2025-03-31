#' Leave-One-Out Coordination from Matrix
#'
#' Calculates R² values using leave-one-out coordination for each gene set in the given matrix.
#' The matrix can be generated as a cell - type - pseudobulk - matrix from any single cell dataset.
#' Computes R² values for each sample and gene set, thereby estimating the coordination for the given matrix.
#' The matrix should contain the samples of interest of a given cell type or condition (pseudobulk).
#' Make sure to subset datD to only contain samples of interest (e.g. healthy vs diseased)
#' The results are then saved to a specified output file.
#'
#' @importFrom stats lm
#' @importFrom utils write.csv
#'
#' @param datD A normalized matrix containing the dataset from which R² values are calculated. Rows are the samples of interest. Columns are genes.
#' @param go_sets A dataframe or matrix where each column is a gene set and rows are genes.
#'        The entries should be binary, indicating gene inclusion in the gene set.
#'        Can be generated using the process_gene_sets_for_matrix function.
#' @param output_folder String specifying the path to the output folder where results will be saved.
#' @param coordination_key A string identifier used in naming the output file (e.g. the cell type or condition).
#'
#' @return Does not return anything explicitly, but saves the R² results as a CSV file in the output folder.
#'
#' @examples
#' # Create a dummy expression matrix
#' set.seed(123)
#' datD <- matrix(rexp(30, rate = 1), nrow = 3, ncol = 10)
#' colnames(datD) <- paste0("Gene", 1:10)
#' rownames(datD) <- paste0("Sample", 1:3)
#'
#' # Create random gene sets
#' go_sets <- matrix(sample(0:1, 10 * 5, replace = TRUE), nrow = 10, ncol = 5)
#' rownames(go_sets) <- colnames(datD)  # gene names as rownames
#' colnames(go_sets) <- paste0("Pathway", 1:5)
#'
#' # Create temporary output folder
#' output_folder <- tempfile("coordination_output_example")
#' dir.create(output_folder)
#'
#' # Run function
#' loo_coordination_from_matrix(
#'   datD = datD,
#'   go_sets = go_sets,
#'   output_folder = output_folder,
#'   coordination_key = "ExampleCellType"
#' )
#' @export

loo_coordination_from_matrix <- function(datD, go_sets, output_folder, coordination_key = "nokey") {

  # Initialize an empty results dataframe for R^2 values
  r2_results_df <- data.frame(
    GeneSet = character(),
    Sample = character(),
    R2Value = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop through each gene set
  for (func in colnames(go_sets)) {
    # Subset data for current biological function
    selected_genes <- rownames(go_sets)[go_sets[, func] == 1]
    common_genes <- selected_genes[selected_genes %in% colnames(datD)]
    if(length(common_genes) < length(selected_genes)) {
      warning(paste("Not all genes found for", func, "- proceeding with", length(common_genes), "genes"))
    }

    dat_subsetD <- datD[, common_genes, drop=FALSE]
    nD_subset <- nrow(dat_subsetD)
    r2_D_subset <- numeric(nD_subset)

    for (i in 1:nD_subset) {
      datD_loo <- dat_subsetD[-i,]
      svdD_loo <- svd(datD_loo)
      V1_loo <- svdD_loo$v[,1]
      lmres <- lm(as.numeric(dat_subsetD[i,]) ~ V1_loo)
      r2_D_subset[i] <- summary(lmres)$r.squared
    }

    r2_D_subset <- r2_D_subset[!is.nan(r2_D_subset)]


    for (i in 1:length(r2_D_subset)) {
      r2_results_df <- rbind(r2_results_df, data.frame(
        GeneSet = func,
        Sample = rownames(datD)[i],
        R2Value = r2_D_subset[i]
      ))
    }
  }

  # Save the R^2 results for the current cell type
  output_filename <- paste0(output_folder, "/", coordination_key, "_r2_results.csv")
  write.csv(r2_results_df, file = output_filename, row.names = FALSE)

}

