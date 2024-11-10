setwd(".../Atlas/scSPECTRA/multispecies/")
go_sets = read.table("go_sets_subset.csv", sep=",", header=TRUE)
# Extract gene names and set row names
gene_names = go_sets[,1]
row.names(go_sets) = gene_names
go_sets = go_sets[,-1]

# Filter go_sets based on gene count criteria
go_sets <- go_sets[, sapply(go_sets, function(col) {
  sum_genes <- sum(col == 1)
  return(sum_genes >= 10 && sum_genes <= 100)
})]

#metadat
metadat = read.table(".../Atlas/scSPECTRA/multispecies/multispecies_metadata.csv", sep=",", header=TRUE)
rownames(metadat) <- metadat$orig_ident

# Define cell types
#cell_types <- c("Podo", "TAL", "DTL_ATL", "DCT", "DCT2", "CNT", "CD_PC", "MD", "PTS1", "PTS2", "PTS3", "injPT", "injTAL", "injDCT_CNT", "EC", "Stromal", "ICA", "ICB",  "Immune", "PEC")
#cell_types <- c("Podo", "TAL", "DTL_ATL", "DCT", "DCT2", "CNT", "CD_PC", "MD", "PTS1", "PTS2", "PTS3", "injPT", "injTAL", "injDCT_CNT", "EC", "Stromal", "ICA", "ICB",  "Immune")
cell_types <- c("PEC")



# Loop through each cell type
for (specific_cell_type in cell_types) {

  # Read the data for the specific cell type
  datD <- read.table(paste0(".../Atlas/scSPECTRA/multispecies/Preprocessing/", specific_cell_type, "_pseudobulk_normalized.csv"), sep = ",", header = TRUE)
  # Preserve sample names as rownames
  row.names(datD) <- datD[, 1]
  datD <- datD[, -1]  # Remove the first column as it's now set as rownames

  # Initialize an empty results dataframe for R^2 values
  r2_results_df <- data.frame(
    GeneSet = character(),
    Sample = character(),
    R2Value = numeric()
  )

  # Loop through each gene set
  for (func in colnames(go_sets)) {

    # Subset data to contain only the genes associated with the current biological function
    selected_genes <- rownames(go_sets)[go_sets[, func] == 1]

    # Ensure all selected genes are in datD
    common_genes <- selected_genes[selected_genes %in% colnames(datD)]
    if(length(common_genes) < length(selected_genes)) {
      warning(paste("Not all genes found for", func, "- proceeding with", length(common_genes), "genes"))
    }

    dat_subsetD <- datD[, common_genes, drop=FALSE]

    # Calculate R^2 for the subset with LOO
    nD_subset <- nrow(dat_subsetD)
    r2_D_subset <- rep(0, nD_subset)
    for (i in 1:nD_subset) {
      datD_loo <- dat_subsetD[-i,]
      svdD_loo <- svd(datD_loo)
      V1_loo <- svdD_loo$v[,1]
      lmres <- lm(as.numeric(dat_subsetD[i,]) ~ V1_loo)
      r2_D_subset[i] <- summary(lmres)$r.squared
    }
    r2_D_subset = r2_D_subset[!is.nan(r2_D_subset)]

    # Store R^2 results
    for (i in 1:length(r2_D_subset)) {
      r2_results_df <- rbind(r2_results_df, data.frame(
        GeneSet = func,
        Sample = rownames(datD)[i],
        R2Value = r2_D_subset[i]
      ))
    }
  }

  # Save the R^2 results for the current cell type
  output_filename <- paste0(".../Atlas/scSPECTRA/multispecies_marker_new/", specific_cell_type, "_r2_results.csv")
  write.csv(r2_results_df, file = output_filename, row.names = FALSE)
}

