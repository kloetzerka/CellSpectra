#' Prepare Seurat Object for reference generation / pseudobulking / spectra analysis
#'
#' Modifies a Seurat object to include specific metadata: cell type annotation, sample/patient identifier,
#' and condition information (mapped to 'query' or 'control' based on provided lists).
#'
#'
#' @param seurat_obj A Seurat object to be modified.
#' @param celltype_col Name of the column in the Seurat object metadata that contains cell type annotations.
#' @param sample_id_col Name of the column in the Seurat object metadata that contains sample/patient identifiers.
#' @param condition_col Name of the column in the Seurat object metadata that contains condition information.
#' @param query_list List of values in the condition column to be mapped as 'query'.
#' @param control_list List of values in the condition column to be mapped as 'control'.
#'
#' @return A modified Seurat object with additional metadata columns: celltype_annotation, orig_ident, and query_control.
#'
#' @examples
#' # Create dummy count matrix
#' counts <- matrix(rpois(20, lambda = 5), nrow = 4)
#' rownames(counts) <- paste0("Gene", 1:4)
#' colnames(counts) <- paste0("Cell", 1:5)
#'
#' # Create a minimal Seurat object
#' seurat_object <- SeuratObject::CreateSeuratObject(counts = counts)
#'
#' # Downgrade assay to Seurat v4-compatible structure (if Seurat v5 is used)
#' seurat_object[["RNA"]] <- as(seurat_object[["RNA"]], Class = "Assay")
#'
#' # Add minimal metadata: cell type, sample, and condition
#' seurat_object$celltypes <- c("CellType1", "CellType2", "CellType1", "CellType2", "CellType1")
#' seurat_object$samples <- c("Sample1", "Sample1", "Sample2", "Sample2", "Sample1")
#' seurat_object$condition <- c("ConditionA", "ConditionA", "ConditionB", "ConditionB", "ConditionA")
#'
#' # Prepare Seurat object for analysis
#' seurat_object <- prepare_seurat_object(
#'   seurat_object,
#'   celltype_col = "celltypes",
#'   sample_id_col = "samples",
#'   condition_col = "condition",
#'   query_list = c("ConditionB"),
#'   control_list = c("ConditionA")
#' )
#' @export

prepare_seurat_object <- function(seurat_obj, celltype_col, sample_id_col, condition_col, query_list, control_list) {
  # Ensure that the provided columns and lists exist and are valid
  if (!(celltype_col %in% colnames(seurat_obj@meta.data))) {
    stop("celltype_col not found in Seurat object metadata.")
  }
  if (!(sample_id_col %in% colnames(seurat_obj@meta.data))) {
    stop("sample_id_col not found in Seurat object metadata.")
  }
  if (!(condition_col %in% colnames(seurat_obj@meta.data))) {
    stop("condition_col not found in Seurat object metadata.")
  }

  # Assign the celltype and orig_ident metadata
  seurat_obj@meta.data$celltype_annotation <- seurat_obj@meta.data[, celltype_col]
  seurat_obj@meta.data$orig_ident <- seurat_obj@meta.data[, sample_id_col]

  # Initialize the query_control column with NA
  seurat_obj@meta.data$query_control <- NA

  # Map specific conditions to 'query' or 'control'
  query_indices <- which(seurat_obj@meta.data[, condition_col] %in% query_list)
  control_indices <- which(seurat_obj@meta.data[, condition_col] %in% control_list)

  seurat_obj@meta.data$query_control[query_indices] <- "query"
  seurat_obj@meta.data$query_control[control_indices] <- "control"

  # Return the modified Seurat object
  return(seurat_obj)
}
