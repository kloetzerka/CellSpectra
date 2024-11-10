#' Process Gene Sets from GO Annotations
#'
#' Reads Gene Ontology (GO) annotations from a GAF file and GO term names from an OBO file. It then
#' creates a binary matrix indicating the presence (1) or absence (0) of genes in each GO set.
#' The matrix is filtered based on a specified minimum and maximum number of genes.
#'
#'latest files can be downloaded from: http://current.geneontology.org
#'
#' @param gaf_file_path String; the file path to the GAF file containing GO annotations. Make sure to download file for species of interest.
#' @param obo_file_path String; the file path to the OBO file containing GO term names.
#' @param genes Vector of strings; the list of gene symbols to include in the analysis. e.g. gene symbols in a single cell or gene expression matrix.
#' @param min_genes Integer; the minimum number of genes for a GO set to be included (default = 10).
#' @param max_genes Integer; the maximum number of genes for a GO set to be included (default = 100).
#'
#' @return A binary matrix where rows represent genes, columns represent GO sets (with GO IDs
#'         and names), and values are 1 (gene is in the GO set) or 0 (gene is not in the GO set).
#'         The matrix includes only GO sets with gene counts within the specified range.
#'
#' @examples
#' # Assuming you have the necessary files and a list of genes
#' # gaf_file <- "path/to/gaf_file.gaf"
#' # obo_file <- "path/to/obo_file.obo"
#' # genes <- c("gene1", "gene2", "gene3")
#' # go_sets <- process_gene_sets_from_GO(gaf_file, obo_file, genes, 10, 100)
#'
#' @export



process_gene_sets_from_GO <- function(gaf_file_path, obo_file_path, genes, min_genes = 10, max_genes = 100) {
  # Read GAF File
  gaf_data <- read.table(gaf_file_path, header = FALSE, sep = "\t", comment.char = "!", quote = "", fill = TRUE, stringsAsFactors = FALSE)
  unique_go_ids <- unique(gaf_data$V5)

  # Read OBO File and Create GO ID to Name Mapping
  obo_lines <- readLines(obo_file_path)
  go_terms <- list()
  current_id <- ""
  for (line in obo_lines) {
    if (grepl("^id: GO:", line)) {
      current_id <- sub("id: ", "", line)
    } else if (grepl("^name:", line)) {
      go_terms[[current_id]] <- sub("name: ", "", line)
      current_id <- ""  # Reset current_id
    }
  }

  # Create a Binary Matrix/Dataframe for GO Sets
  gene_go_matrix <- matrix(0, nrow = length(genes), ncol = length(unique_go_ids), dimnames = list(genes, unique_go_ids))
  for (row in 1:nrow(gaf_data)) {
    gene <- gaf_data$V3[row]
    go_id <- gaf_data$V5[row]
    if (gene %in% genes && go_id %in% colnames(gene_go_matrix)) {
      gene_go_matrix[gene, go_id] <- 1
    }
  }

  # Adjust Column Names to Include GO Term Names
  colnames(gene_go_matrix) <- sapply(colnames(gene_go_matrix), function(id) paste(id, go_terms[[id]]))

  # Filter Based on Gene Count Criteria
  valid_cols <- colSums(gene_go_matrix) >= min_genes & colSums(gene_go_matrix) <= max_genes
  filtered_gene_go_matrix <- gene_go_matrix[, valid_cols]

  return(filtered_gene_go_matrix)
}
