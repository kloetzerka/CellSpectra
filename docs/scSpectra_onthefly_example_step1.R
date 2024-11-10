library(Seurat)
library(Matrix)
library(dplyr)


output_folder_base <- "/home/kloetzer/Atlas/scSpectra/species_healthy_disease/onthefly_output/"  # Specify your base output directory

# Assuming you have set the number of pseudoreplicates
num_replicates <- 3  # Adjust this number as needed

create_pseudobulk <- function(seurat_subset, target_count, num_replicates, seed = NULL) {
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



#we need a count matrix to pick the genes to do that 

seurat_object = readRDS("/home/kloetzer/Atlas/objects/Atlas6.6.rds")

#the following not needed anymore. We extract genes directly from seurat object

# Set working directory and read in data
#setwd("/home/kloetzer/Atlas/Atlas_Extension/")

# Read and subsets for gene names in data
#dat = read.table("PT_pseudobulk_normalized.csv", sep=",", header=TRUE)

# create geneset data sets for scSPECTRA from enrichR website


go_sets = read.table("/home/kloetzer/Atlas/MN/go_sets_subset.csv", sep=",", header=TRUE)
gene_names = go_sets[,1]
row.names(go_sets) = gene_names
go_sets = go_sets[,-1]


# Filter go_sets based on gene count criteria
go_sets <- go_sets[, sapply(go_sets, function(col) {
  sum_genes <- sum(col == 1)
  return(sum_genes >= 10 && sum_genes <= 100)
})]


# Save go_sets
saveRDS(go_sets, file = paste0(output_folder_base, "go_sets_modified.rds"))

#cell_types = c("Podo", "TAL", "DTL_ATL", "DCT_CNT_CD", "EC", "Stromal", "PT", "IC", "Immune", "PEC")

cell_types = c("Podo")

##################sample level and celltype results############

# Loop through each cell type
for (cell_type in cell_types) {
  
  
  #create output folder for cell type
  celltype_output_folder <- paste0(output_folder_base, cell_type, "/")
  if (!dir.exists(celltype_output_folder)) {
    dir.create(celltype_output_folder, recursive = TRUE)
  }
  
  
  # Subset the Seurat object by cell type
  celltype_subset <- subset(seurat_object, subset = annotation_final_level1 == cell_type)
  
  # Get the list of samples with at least 10 cells
  valid_samples <- celltype_subset@meta.data %>% 
    group_by(orig_ident) %>%
    filter(n() >= 10) %>%
    pull(orig_ident) %>%
    unique()
  
  saveRDS(valid_samples, file = paste0(celltype_output_folder, "valid_samples.rds"))
  
  # Separate healthy and diseased samples
  set_subsetH <- subset(celltype_subset, subset = disease == "healthy")
  set_subsetD <- subset(celltype_subset, subset = disease != "healthy")
  
  set_subsetD <- subset(celltype_subset, subset = orig_ident %in% valid_samples)
  
  diseased_samples <- intersect(valid_samples, unique(set_subsetD@meta.data$orig_ident))
  
  # Loop through each diseased sample
  for (sample in intersect(valid_samples, unique(set_subsetD@meta.data$orig_ident))) {
    # Create pseudobulk for the query sample
    query_sample_subset <- subset(set_subsetD, subset = orig_ident == sample)
    query_counts <- Matrix::colSums(t(query_sample_subset@assays$RNA@counts))
    query_total_count <- sum(query_counts)
    print(sum(query_counts))
    
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
      pseudobulk_replicates <- create_pseudobulk(ref_sample_subset, query_total_count, num_replicates, seed = 10)
      
      for (replicate_num in 1:num_replicates) {
        replicate_name <- paste0(ref_sample, "_", replicate_num)
        ref_pseudobulk[replicate_name, ] <- pseudobulk_replicates[[paste0("replicate_", replicate_num)]]
      }
    }
    
    # Normalize total counts for each sample
    total_counts <- rowSums(ref_pseudobulk)
    print(rowSums(ref_pseudobulk))
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
