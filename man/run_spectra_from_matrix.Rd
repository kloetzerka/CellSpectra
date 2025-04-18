% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_spectra_from_matrix.R
\name{run_spectra_from_matrix}
\alias{run_spectra_from_matrix}
\title{Run Spectra Analysis from Matrix Data}
\usage{
run_spectra_from_matrix(
  datH,
  datD,
  go_sets,
  output_folder,
  coordination_key = "nokey",
  CHISQ.MAX = 4,
  expression_threshold = 0,
  gene_number_threshold = 10
)
}
\arguments{
\item{datH}{A normalized (e.g. tpm + log) matrix of gene expression data for healthy samples, with genes as columns.}

\item{datD}{A normalized (e.g. tpm + log) matrix of gene expression data for diseased samples, similar to datH.}

\item{go_sets}{A dataframe or matrix where columns represent gene sets. Can be generated using process_gene_sets_for_matrix function.}

\item{output_folder}{The directory path where the result files will be saved.}

\item{coordination_key}{An optional key for file naming (default is 'nokey'). Could be cell type or condition.}

\item{CHISQ.MAX}{An optional maximum chi-square value for robustness (default is 4).}

\item{expression_threshold}{Genes with an average expression in the reference
OR query matrix below this threshold are removed from the analysis.
Default 0. Set to -1 for no filtering applied.}

\item{gene_number_threshold}{Number of genes in a gene set after filtering needed to perform analysis. If less genes remaining, results will be NA.
Default 10 genes.}
}
\value{
This function does not return a value but saves two CSV files containing R² and p-values.
}
\description{
Performs spectra analysis on given matrices query (datD) samples compared to a reference (datH).
It calculates R² and p-values for each gene set and saves the results as CSV files.
}
\examples{

# create dummy object and gene sets
set.seed(123)
genes <- paste0("Gene", 1:10)
datH <- matrix(rexp(50, rate = 1), nrow = 10, ncol = 5)
rownames(datH) <- genes
colnames(datH) <- paste0("Healthy", 1:5)
datH <- t(datH)  # samples as rows
datD <- matrix(rexp(10, rate = 1), nrow = 10, ncol = 2)
rownames(datD) <- genes
colnames(datD) <- paste0("Disease", 1:2)
datD <- t(datD)  # samples as rows
go_sets <- matrix(sample(0:1, 10 * 5, replace = TRUE), nrow = 10, ncol = 5)
rownames(go_sets) <- genes
colnames(go_sets) <- paste0("Pathway", 1:5)

# Create temporary output folder
output_folder <- tempfile("example_output_")
dir.create(output_folder)

# Run function
run_spectra_from_matrix(
  datH = datH,
  datD = datD,
  go_sets = go_sets,
  output_folder = output_folder,
  coordination_key = "ExampleCellType",
  expression_threshold = 0,
  gene_number_threshold = 0
)
}
