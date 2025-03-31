#' Process Gene Sets (for a Matrix)
#'
#' Reads a file containing gene sets, processes each line to extract gene sets and their genes,
#' then matches these with genes in the provided matrix. The function filters the gene sets
#' based on a minimum and maximum number of genes. The difference to the process_gene_sets function
#' is the matching to a (pseudo)bulk matrix of samples (rows) and genes (columns).
#'
#'https://maayanlab.cloud/Enrichr/#libraries
#'
#'or
#'
#'https://www.gsea-msigdb.org/gsea/msigdb
#'
#' @param filename A string, the path to the file containing gene sets.
#' @param datD A matrix, typically a pseudobulk matrix, where cols are genes.
#' @param min_genes An integer, the minimum number of genes for a gene set to be considered valid.
#' @param max_genes An integer, the maximum number of genes for a gene set to be considered valid.
#'
#' @return A data frame where each column represents a gene set, rows represent genes,
#'         and entries are binary (1 if the gene is in the gene set, 0 otherwise).
#'         The data frame is filtered based on the specified gene count criteria.
#'
#' @examples
#' #gene_sets_file <- "path/to/gene_sets.txt"
#' #pseudobulk_matrix <- matrix(data, nrow = number_of_samples, ncol = number_of_genes)
#' #go_sets <- process_gene_sets_for_matrix(gene_sets_file, datD = pseudobulk_matrix, 10, 100)
#'
#' @export

process_gene_sets_for_matrix <- function(filename, datD, min_genes = 10, max_genes = 100) {
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
  genes <- colnames(datD)

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


  return(go_sets)
}










