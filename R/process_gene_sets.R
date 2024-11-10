#' Process and Filter Gene Sets
#'
#' Reads an enrichR gene - gene set library file, processes each line to extract gene sets and their genes,
#' and filters these gene sets based on the gene count criteria. It integrates with a Seurat object
#' to match and filter the gene sets. Output can be used for a scSpectra analysis.
#'
#' https://maayanlab.cloud/Enrichr/#libraries
#'
#' @param filename A string specifying the path to the gene sets file.
#' @param seurat_object A Seurat object containing the RNA assays.
#' @param output_folder_base A string specifying the base path where the output file will be saved. Default: output_folder_base.
#' @param min_genes An integer giving the minimum number of genes for a gene set to be included. Default: 10.
#' @param max_genes An integer giving the maximum number of genes for a gene set to be included. Default: 100.
#'
#' @return A data frame of gene sets where each column represents a gene set, and each row
#' represents a gene. The data frame contains binary values indicating gene presence in each set.
#'
#' @examples
#' # Assuming you have a valid 'seurat_object' and a gene sets file:
#' #go_sets <- process_gene_sets(
#' #filename = "path/to/file.txt",
#' #seurat_object = object,
#' #output_folder_base = output_folder_base,
#' #min_genes = 10,
#' #max_genes = 100)
#'
#' @export

process_gene_sets <- function(filename, seurat_object, output_folder_base = output_folder_base, min_genes = 10, max_genes = 100) {
  lines <- readLines(filename)

  # Process each line to extract gene sets and their genes
  gene_sets <- list()
  for (line in lines) {
    elements <- unlist(strsplit(line, "\t"))
    gene_set <- elements[1]
    genes <- elements[-1]  # Remove the first element (gene set name)
    gene_sets[[gene_set]] <- genes
  }

  # Get genes from Seurat object
  genes <- rownames(seurat_object@assays$RNA)

  # Create an empty dataframe
  gene_set_names <- names(gene_sets)
  df <- as.data.frame(matrix(0, nrow = length(genes), ncol = length(gene_set_names)))
  rownames(df) <- genes
  colnames(df) <- gene_set_names

  # Fill the dataframe
  for (gene_set in gene_set_names) {
    associated_genes <- gene_sets[[gene_set]]
    for (gene in associated_genes) {
      if (gene %in% rownames(df)) {
        df[gene, gene_set] <- 1
      }
    }
  }

  go_sets <- df

  # Filter go_sets based on gene count criteria
  go_sets <- go_sets[, sapply(go_sets, function(col) {
    sum_genes <- sum(col == 1)
    return(sum_genes >= min_genes && sum_genes <= max_genes)
  })]

  # Save go_sets
  output_path <- file.path(output_folder_base, "go_sets_modified.rds")
  saveRDS(go_sets, file = output_path)

  return(go_sets)
}
