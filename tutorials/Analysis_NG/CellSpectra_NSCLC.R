library(Matrix)
library(dplyr)
library(Seurat)
library(CellSpectra)

#select cell type - consider that we need to change the epithelial cells to cancer in the reference if we run it on cancer
ctypes = c(Sys.getenv("CELLTYPE"))

#define output folder
output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_nonTumorRef_Primary/"

#here we define if we run it on metastasis or primary tumor

query_samples <- c("tumor_primary")

ref_samples <- c("normal", "normal_adjacent")

ref_and_query_samples <- c("tumor_primary", "normal", "normal_adjacent")

#first we load the mouse atlas object
seurat_object <- readRDS("/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/cancer_core_CaEp_convert.rds")

#we need to rename our assay after the conversion to a seurat
seurat_object@assays$RNA <- seurat_object@assays$originalexp


cancer_samples <- c("lung adenocarcinoma", "squamous cell lung carcinoma", "non-small cell lung carcinoma")
non_cancer_samples <- c("normal", "chronic obstructive pulmonary disease")

seurat_object = subset(seurat_object, subset = origin %in% ref_and_query_samples)

# Create a new annotation vector based on the cell_type_major
seurat_object$disease_grouped_qr <- as.character(seurat_object$disease)

# Map cell types to their new groups
seurat_object$disease_grouped_qr[seurat_object$origin %in% query_samples] <- "query_samples"
seurat_object$disease_grouped_qr[seurat_object$origin %in% ref_samples] <- "reference_samples"

# Convert the grouped cell type back to a factor
seurat_object$disease_grouped_qr <- factor(seurat_object$disease_grouped_qr)


#define your sample
sample_to_change <- Sys.getenv("SAMPLE")
#sample_to_change <- "Chen_Zhang_2020_NSCLC-10"

# Define the sample identifier and the condition you want to change
sample_id_col <- "donor_id"
condition_col <- "disease_grouped_qr"

# Convert the condition column to character if it is a factor
seurat_object@meta.data[[condition_col]] <- as.character(seurat_object@meta.data[[condition_col]])

# Change the condition key for the selected sample directly in the Seurat object
seurat_object@meta.data[[condition_col]][seurat_object@meta.data[[sample_id_col]] == sample_to_change & seurat_object@meta.data[["disease_grouped_qr"]] == "query_samples"] <- "sample_of_interest"

# Convert the condition column back to factor if necessary
seurat_object@meta.data[[condition_col]] <- as.factor(seurat_object@meta.data[[condition_col]])

#we use the prepare function to get our seurat object in shape
seurat_object <- prepare_seurat_object(
  seurat_obj = seurat_object, #our object to prepare
  celltype_col = "cell_type_major", #column in which annotation is stored
  sample_id_col = sample_id_col, #sample identifier
  condition_col = condition_col, #column storing the condition information
  query_list = c("sample_of_interest"), #categories in the condition_col defining the query samples
  control_list = c("reference_samples") #categories in the condition_col defining the reference
)

#Since CellSpectra uses cell_type as a variabel in the subset function we have to remove the cell_type column
seurat_object$cell_type <- NULL
seurat_object$sample <- NULL


#Step 1
create_references(
  seurat_object = seurat_object,
  output_folder_base = output_folder_base,
  num_replicates = 3,
  cell_types = ctypes,
  cell_number_threshold = 10,
  seed = 123
)



output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_nonTumorRef_Primary/"

cell_types = c(Sys.getenv("CELLTYPE"))

soi = Sys.getenv("SAMPLE")

run_spectra(output_folder_base = output_folder_base, cell_types = cell_types, CHISQ.MAX = 4, expression_threshold = 0, gene_number_threshold = 10,
            restrict_to_sample = TRUE, results_sample_key = soi, vector_of_samples = c(soi), QC_report = TRUE)

