#' Create Pseudobulk References
#'
#' This function generates pseudobulk references by aggregating single-cell expression data into pseudobulk samples.
#' It does this by sampling cells with replacement to reach a target count, ensuring each pseudobulk sample has a similar total count of expression.
#' This can be crucial in cases with lower cell numbers in query samples compared to reference samples (e.g. biopsy samples).
#'
#'
#' @param seurat_subset A Seurat object subset for which pseudobulk replicates are to be created. Usually consists of a specific cell type and the reference samples.
#' @param target_count An integer representing the target total count for each pseudobulk sample (counts of the query sample we want to compare to the reference).
#' @param num_replicates The number of pseudobulk replicates to create. If query has more cells than most of the reference this can be set to 1. If the query count is smaller than the reference samples, increasing replicates can improve the estimation of V1.
#' @param seed An optional integer seed for reproducibility of random reference cell sampling.
#'
#' @return A list of pseudobulk replicates, each replicate being a vector of summed expression counts.
#'
#' @examples
#' # Assuming 'seurat_subset' is a subset of a Seurat object
#' #(reference samples + cell type of interest):
#' #creating a minimal seurat object first:
#' counts <- matrix(rpois(20, lambda = 5), nrow = 4)
#' rownames(counts) <- paste0("Gene", 1:4)
#' colnames(counts) <- paste0("Cell", 1:5)
#' seurat_subset <- SeuratObject::CreateSeuratObject(counts = counts)
#' #change to V4 structure
#' seurat_subset[["RNA"]] <- as(seurat_subset[["RNA"]], Class = "Assay")
#' # Run function
#' pseudobulks <- create_pseudobulk(seurat_subset,
#' target_count = 50,
#' num_replicates = 3,
#' seed = 10)
#'
#' @export

create_pseudobulk <- function(seurat_subset, target_count, num_replicates, seed = seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  cell_names <- colnames(seurat_subset@assays$RNA@counts)
  pseudobulk_replicates <- list()

  for (replicate in 1:num_replicates) {
    shuffled_cells <- sample(cell_names, replace = TRUE)  # Shuffle for each replicate
    pseudobulk_counts <- rep(0, nrow(seurat_subset@assays$RNA@counts))

    while (sum(pseudobulk_counts) < target_count) {
      for (cell in shuffled_cells) {
        cell_counts <- seurat_subset@assays$RNA@counts[, cell]
        pseudobulk_counts <- pseudobulk_counts + cell_counts

        if (sum(pseudobulk_counts) >= target_count) {
          break
        }
      }

      if (sum(pseudobulk_counts) < target_count) {
        shuffled_cells <- sample(cell_names, replace = TRUE)  # Reshuffle if target not met
      }
    }

    pseudobulk_replicates[[paste0("replicate_", replicate)]] <- pseudobulk_counts
  }

  return(pseudobulk_replicates)
}
