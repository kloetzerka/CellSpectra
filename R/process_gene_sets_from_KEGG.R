#' Process_gene_sets_from_KEGG
#'
#' This function retrieves all KEGG pathways for a given organism and creates a binary matrix
#' indicating the presence of genes (represented in a Seurat object) in these pathways. The
#' matrix has to be saved manually to the output base folder ("go_sets_modified.rds")
#'
#' @param organism Character. The KEGG organism code (e.g., "mmu" for mouse).
#' @param seurat_object Seurat object used for the spectra analysis (gene symbols to use)
#' @param use_genes_directly  If TRUE gene_vec is used as gene symbol input.
#' @param gene_vec Vector of genes considered in the gene set database.
#' @return A dataframe where rows are genes and columns are KEGG pathways.
#'         Each cell contains 1 if the gene is part of the pathway, otherwise 0.
#' @import KEGGREST
#' @import dplyr
#' @import tidyr
#' @export
#' @examples
#'
#' organism <- "mmu"
#'
#' #kegg_gene_matrix <- process_gene_sets_from_KEGG(organism, seurat_object = seurat_object)
#' #saveRDS(kegg_gene_matrix, "go_sets_modified.rds")
process_gene_sets_from_KEGG <- function(organism, seurat_object, use_genes_directly = FALSE, gene_vec = c()) {

  if (use_genes_directly == TRUE) {
    genes = gene_vec
  }

  if (use_genes_directly == FALSE) {
    genes = rownames(seurat_object@assays$RNA)
  }

  get_all_kegg_pathways <- function(organism) {
    pathways <- keggList("pathway", organism)
    pathway_ids <- names(pathways)
    pathway_descriptions <- gsub(paste0(" - ", organism, ".*"), "", pathways)
    return(setNames(pathway_descriptions, pathway_ids))
  }

  get_genes_from_kegg_pathway <- function(pathway_id) {
    pathway_info <- keggGet(pathway_id)[[1]]
    genes <- pathway_info$GENE
    if (is.null(genes)) {
      return(character(0))
    }
    gene_names <- genes[seq(2, length(genes), 2)]  # Extract gene names
    gene_names <- gsub("\\;.*", "", gene_names)  # Remove descriptions
    return(gene_names)
  }

  pathways <- get_all_kegg_pathways(organism)
  pathway_ids <- names(pathways)
  pathway_descriptions <- as.character(pathways)

  # Initialize binary matrix
  gene_kegg_matrix <- matrix(0, nrow = length(genes), ncol = length(pathway_ids), dimnames = list(genes, pathway_ids))

  for (pathway_id in pathway_ids) {
    pathway_genes <- get_genes_from_kegg_pathway(pathway_id)
    gene_kegg_matrix[genes %in% pathway_genes, pathway_id] <- 1
  }

  # Adjust column names to include pathway descriptions
  colnames(gene_kegg_matrix) <- sapply(colnames(gene_kegg_matrix), function(id) paste(id, pathway_descriptions[which(pathway_ids == id)], sep = " "))

  # Convert to dataframe
  gene_kegg_df <- as.data.frame(gene_kegg_matrix)

  # Remove columns where the sum is zero
  gene_kegg_df <- gene_kegg_df[, colSums(gene_kegg_df) > 0]

  return(gene_kegg_df)
}

