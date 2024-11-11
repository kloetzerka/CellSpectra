library(Matrix)
library(dplyr)
library(Seurat)
library(CellSpectra)

#select cell type
ctypes = c(Sys.getenv("CELLTYPE"))

#define output folder
output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_Ref_human_diseased_all/"

#here we define if we run it on metastasis or primary tumor

#first we load the mouse atlas object
seurat_object <- readRDS("/home/isilon/users/o_kloetzer/Atlas/Atlas6.6.rds")

# Create a new annotation vector based on the cell_type_major
seurat_object$disease_grouped_qr <- as.character(seurat_object$disease)

# Map cell types to their new groups
seurat_object$disease_grouped_qr[seurat_object$species == "human" & seurat_object$disease == "diseased"] <- "reference_samples"
seurat_object$disease_grouped_qr[seurat_object$disease_grouped_qr != "reference_samples"] <- "query_samples"

# Convert the grouped cell type back to a factor
seurat_object$disease_grouped_qr <- factor(seurat_object$disease_grouped_qr)

# Define the sample identifier and the condition you want to change
sample_id_col <- "orig_ident"
condition_col <- "disease_grouped_qr"


#we use the prepare function to get our seurat object in shape
seurat_object <- prepare_seurat_object(
  seurat_obj = seurat_object, #our object to prepare
  celltype_col = "annotation_final_level1", #column in which annotation is stored
  sample_id_col = sample_id_col, #sample identifier
  condition_col = condition_col, #column storing the condition information
  query_list = c("query_samples"), #categories in the condition_col defining the query samples
  control_list = c("reference_samples") #categories in the condition_col defining the reference
)

#not necessary in this case:
#Since scSpectra uses cell_type as a variabel in the subset function we have to remove the cell_type column
#seurat_object$cell_type <- NULL
#seurat_object$sample <- NULL

#Step 1
create_references(
  seurat_object = seurat_object,
  output_folder_base = output_folder_base,
  num_replicates = 3,
  cell_types = ctypes,
  cell_number_threshold = 10,
  seed = 123
)


#Step 2
output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_Ref_human_diseased_all/"

cell_types = c(Sys.getenv("CELLTYPE"))

run_spectra(output_folder_base = output_folder_base, cell_types = cell_types, CHISQ.MAX = 4, expression_threshold = -1, gene_number_threshold = 10,
            restrict_to_sample = FALSE, QC_report = TRUE)

