% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_pseudobulk.R
\name{create_pseudobulk}
\alias{create_pseudobulk}
\title{Create Pseudobulk References}
\usage{
create_pseudobulk(seurat_subset, target_count, num_replicates, seed = seed)
}
\arguments{
\item{seurat_subset}{A Seurat object subset for which pseudobulk replicates are to be created. Usually consists of a specific cell type and the reference samples.}

\item{target_count}{An integer representing the target total count for each pseudobulk sample (counts of the query sample we want to compare to the reference).}

\item{num_replicates}{The number of pseudobulk replicates to create. If query has more cells than most of the reference this can be set to 1. If the query count is smaller than the reference samples, increasing replicates can improve the estimation of V1.}

\item{seed}{An optional integer seed for reproducibility of random reference cell sampling.}
}
\value{
A list of pseudobulk replicates, each replicate being a vector of summed expression counts.
}
\description{
This function generates pseudobulk references by aggregating single-cell expression data into pseudobulk samples.
It does this by sampling cells with replacement to reach a target count, ensuring each pseudobulk sample has a similar total count of expression.
This can be crucial in cases with lower cell numbers in query samples compared to reference samples (e.g. biopsy samples).
}
\examples{
# Assuming 'seurat_subset' is a subset of a Seurat object
#(reference samples + cell type of interest):
#creating a minimal seurat object first:
counts <- matrix(rpois(20, lambda = 5), nrow = 4)
rownames(counts) <- paste0("Gene", 1:4)
colnames(counts) <- paste0("Cell", 1:5)
seurat_subset <- SeuratObject::CreateSeuratObject(counts = counts)
#change to V4 structure
seurat_subset[["RNA"]] <- as(seurat_subset[["RNA"]], Class = "Assay")
# Run function
pseudobulks <- create_pseudobulk(seurat_subset,
target_count = 50,
num_replicates = 3,
seed = 10)

}
