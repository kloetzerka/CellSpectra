library(Seurat)
library(Matrix)
library(dplyr)


output_folder_base <- "/home/kloetzer/Atlas/scSpectra/species_healthy_disease/onthefly_output/"  # Specify your base output directory

go_sets = readRDS(file = paste0(output_folder_base, "go_sets_modified.rds"))

#cell_types <- c("Podo", "TAL", "DTL_ATL", "DCT", "DCT2", "CNT", "CD_PC", "MD", "PTS1", "PTS2", "PTS3", "injPT", "injTAL", "injDCT_CNT", "EC", "Stromal", "ICA", "ICB",  "Immune", "PEC")
cell_types = c("Podo")
##################sample level and celltype results############

# Loop through each cell type
for (cell_type in cell_types) {
  # Read in expression data for each cell type
  
  celltype_output_folder <- paste0(output_folder_base, cell_type, "/")
  
  valid_samples = readRDS(file = paste0(celltype_output_folder, "valid_samples.rds"))
  
  samples_of_interest = readRDS(file = paste0(output_folder_base, "samples_of_interest.rds"))
  
  diseased_samples <- intersect(valid_samples, samples_of_interest)
  
  # Initialize dataframes for R² and p-values
  r2_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
  pval_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
  svdpvals_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
  svdchisq_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
  deviant.pvalD_df <- data.frame(matrix(nrow = length(diseased_samples), ncol = length(colnames(go_sets))))
  
  rownames(r2_df) <- diseased_samples
  rownames(pval_df) <- diseased_samples
  colnames(r2_df) <- colnames(go_sets)
  colnames(pval_df) <- colnames(go_sets)
  
  rownames(svdpvals_df) <- diseased_samples
  rownames(svdchisq_df) <- diseased_samples
  rownames(deviant.pvalD_df) <- diseased_samples
  colnames(svdpvals_df) <- colnames(go_sets)
  colnames(svdchisq_df) <- colnames(go_sets)
  colnames(deviant.pvalD_df) <- colnames(go_sets)
  
  # Loop through each diseased sample
  for (sample in diseased_samples) {
    # Create pseudobulk for the query sample
    
    sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
    
    datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
    datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))
    
    # Loop through each gene set
    for (func in colnames(go_sets)) {
      # Select genes associated with the current function
      selected_genes <- rownames(go_sets)[go_sets[, func] == 1]
      common_genes <- selected_genes[selected_genes %in% colnames(datH)]
      
      if (length(common_genes) < length(selected_genes)) {
        warning(paste("Not all genes found for", func, "- proceeding with", length(common_genes), "genes"))
      }
      
      dat_subsetH <- datH[, common_genes, drop = FALSE]
      dat_subsetD <- datD[, common_genes, drop = FALSE]
      
      # Calculate R^2 for the healthy subset with LOO
      nH_subset <- nrow(dat_subsetH)
      r2_H_subset <- rep(0, nH_subset)
      svdpvals_H_subset <- rep(0, nH_subset)
      svdepsilon_H_subset = matrix(nrow=nrow(dat_subsetH), ncol=ncol(dat_subsetH), data=0)
      for (i in 1:nH_subset) {
        datH_loo <- dat_subsetH[-i,]
        svdH_loo <- svd(datH_loo)
        V1_loo <- svdH_loo$v[,1]
        lmres <- lm(as.numeric(dat_subsetH[i,]) ~ V1_loo)
        svdepsilon_H_subset[i,] = lmres$residuals
        svdpvals_H_subset[i] = summary(lmres)$coefficients[2,4]
        r2_H_subset[i] <- summary(lmres)$r.squared
      }
      
      svdsigmag = apply(svdepsilon_H_subset, 2, sd)
      svdH <- svd(dat_subsetH)
      V1_H <- svdH$v[,1]
      # Calculate R^2 for the diseased subset
      nD_subset <- nrow(dat_subsetD)
      r2_D_subset <- rep(0, nD_subset)
      # Regress each disease sample on V1
      #svdepsilonD = matrix(nrow=nrow(dat_subsetD), ncol=ncol(dat_subsetD), data=0)
      svdpvalsD = rep(NA, nD_subset)
      svdchisqD = rep(NA, nD_subset)  # This records the chi-square statistic assessing deviation of sample from normal reference.
      deviant.pvalD = rep(NA, nD_subset)
      CHISQ.MAX = 4 # Cap the contribution of each gene at this value for robustness
      chisq.df = length(svdsigmag)-1 # degree of freedom of chi-sq (number of genes) minus 1
      lmres <- lm(as.numeric(dat_subsetD) ~ V1_H)
      #svdepsilonD[sample,] = lmres$residuals
      svdpvalsD[sample] = summary(lmres)$coefficients[2,4]
      svdchisqD[sample] = sum(pmin((lmres$residuals/svdsigmag)^2, CHISQ.MAX))
      r2_D_subset[sample] <- summary(lmres)$r.squared
      deviant.pvalD[sample]=pchisq(svdchisqD[sample], df=chisq.df, lower.tail=FALSE)
      # Adjust R² values based on p-values
      
      svdpvals_df[sample, func] <- svdpvalsD[sample]
      svdchisq_df[sample, func] <- svdchisqD[sample]
      r2_df[sample, func] <- r2_D_subset[sample]
      deviant.pvalD_df[sample, func] <- deviant.pvalD[sample]
      
      
    }
    # Store results in respective dataframes
    
    
  }
  # Apply FDR correction to p-values for each gene set
  pval_df <- apply(deviant.pvalD_df, 2, function(p) p.adjust(p, method = "fdr"))
  
  # Save R² and p-values to CSV files for each cell type
  write.csv(r2_df, paste0("/home/kloetzer/Atlas/scSpectra/species_healthy_disease/R2/R2_", cell_type, ".csv"), row.names = TRUE)
  write.csv(pval_df, paste0("/home/kloetzer/Atlas/scSpectra/species_healthy_disease/Pval/Pval_", cell_type, ".csv"), row.names = TRUE)
}