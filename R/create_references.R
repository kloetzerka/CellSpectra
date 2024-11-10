#' Create Pseudobulk References Files from Seurat Object matching the counts of query samples.
#'
#' This is the "on the fly" input for our spectra analysis.
#' It ensures that differences in cell numbers between query and reference samples drive an inflation of dyscoordination.
#'
#' This function processes a Seurat object to create pseudobulk references. It subsets the Seurat object
#' by cell type, identifies valid samples based on a cell number threshold, and generates pseudobulk
#' references for both healthy and diseased samples. Results are saved as RDS files in specified output folders.
#'
#' @importFrom dplyr group_by filter %>%
#' @importFrom Matrix colSums t
#' @importFrom Seurat subset
#' @importFrom methods as
#'
#' @param seurat_object A Seurat object to process.
#' @param output_folder_base A string specifying the base path where output files will be saved. Default: output_folder_base.
#' @param num_replicates The number of replicates to generate for each reference sample.
#' @param cell_types A vector of cell types to process. Default: cell_types.
#' @param cell_number_threshold Minimum number of cells required to consider a sample as valid. Default = 10.
#' @param seed Optional; an integer seed for reproducibility of random processes.
#'
#' @return The function does not return a value but saves output files (datH and datD matrices) in specified directories.
#'
#' @examples
#' # Assuming 'seurat_object' is your Seurat object and 'output_folder_base' is defined:
#' #create_references(seurat_object,
#' #output_folder_base,
#' #num_replicates = 1,
#' #cell_types = c("CellType1", "CellType2"),
#' #cell_number_threshold = 10,
#' #seed = NULL)
#'
#' @export


create_references <- function(seurat_object, output_folder_base = output_folder_base, num_replicates = 1, cell_types = cell_types, cell_number_threshold = 10, seed = NULL) {

  # Ensure the output base folder exists
  if (!dir.exists(output_folder_base)) {
    dir.create(output_folder_base, recursive = TRUE)
  }

  for (cell_type in cell_types) {

    #create output folder for cell type
    celltype_output_folder <- paste0(output_folder_base, cell_type, "/")
    if (!dir.exists(celltype_output_folder)) {
      dir.create(celltype_output_folder, recursive = TRUE)
    }


    # Subset the Seurat object by cell type
    celltype_subset <- subset(seurat_object, subset = celltype_annotation  == cell_type)

    # Get the list of samples with at least 10 cells
    valid_samples <- celltype_subset@meta.data %>%
      group_by(orig_ident) %>%
      filter(n() >= cell_number_threshold) %>%
      pull(orig_ident) %>%
      unique()

    saveRDS(valid_samples, file = paste0(celltype_output_folder, "valid_samples.rds"))

    # Separate healthy and diseased samples
    set_subsetH <- subset(celltype_subset, subset = query_control == "control")
    set_subsetD <- subset(celltype_subset, subset = query_control == "query")

    samples_of_interest <- unique(set_subsetD@meta.data$orig_ident)

    saveRDS(samples_of_interest, file = paste0(output_folder_base, "samples_of_interest.rds"))

    set_subsetD <- subset(set_subsetD, subset = orig_ident %in% valid_samples)

    diseased_samples <- intersect(valid_samples, unique(set_subsetD@meta.data$orig_ident))

    # Loop through each diseased sample
    for (sample in diseased_samples) {
      # Create pseudobulk for the query sample
      query_sample_subset <- subset(set_subsetD, subset = orig_ident == sample)
      query_counts <- Matrix::colSums(Matrix::t(query_sample_subset@assays$RNA@counts))
      query_total_count <- sum(query_counts)
      #print(sum(query_counts))

      # Initialize ref_pseudobulk dataframe
      num_genes <- nrow(seurat_object@assays$RNA@counts)  # Assuming 'seurat_subset' is any of the subset objects like 'set_subsetH'
      ref_sample_ids <- intersect(valid_samples, unique(set_subsetH@meta.data$orig_ident))
      num_ref_samples <- num_replicates * length(ref_sample_ids)
      ref_pseudobulk <- data.frame(matrix(nrow = num_ref_samples, ncol = num_genes))

      row_names <- unlist(sapply(ref_sample_ids, function(id) paste(id, 1:num_replicates, sep = "_")))

      rownames(ref_pseudobulk) <- row_names

      colnames(ref_pseudobulk) <- rownames(seurat_object@assays$RNA@counts)

      # When generating pseudobulk for reference samples:
      for (ref_sample in ref_sample_ids) {
        ref_sample_subset <- subset(set_subsetH, subset = orig_ident == ref_sample)
        pseudobulk_replicates <- create_pseudobulk(ref_sample_subset, query_total_count, num_replicates, seed = seed)

        for (replicate_num in 1:num_replicates) {
          replicate_name <- paste0(ref_sample, "_", replicate_num)
          ref_pseudobulk[replicate_name, ] <- pseudobulk_replicates[[paste0("replicate_", replicate_num)]]
        }
      }

      # Normalize total counts for each sample
      total_counts <- rowSums(ref_pseudobulk)
      #print(rowSums(ref_pseudobulk))
      ref_pseudobulk <- ref_pseudobulk / total_counts * 10000

      # Logarithmic transformation
      datH <- log1p(ref_pseudobulk)

      nH_subset <- nrow(datH)

      # Normalize total counts for each sample
      query_counts <- t(query_counts)

      # Normalize total counts for each sample
      total_counts <- rowSums(query_counts)
      query_counts <- query_counts / total_counts * 10000

      # Logarithmic transformation
      datD <- log1p(query_counts)


      # Save datH and datD matrices
      sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
      if (!dir.exists(sample_output_folder)) {
        dir.create(sample_output_folder, recursive = TRUE)
      }

      saveRDS(datH, file = paste0(sample_output_folder, "datH.rds"))
      saveRDS(datD, file = paste0(sample_output_folder, "datD.rds"))
    }

  }
}
