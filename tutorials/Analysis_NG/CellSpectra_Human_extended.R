library(Matrix)
library(dplyr)
library(Seurat)
library(CellSpectra)
library(KEGGREST)

#select cell type
ctypes = c(Sys.getenv("CELLTYPE"))
#ctypes = "Podo"

#define output folder
output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/human_extended/output_human_extended/"

query_samples <- c("AKI", "CKD")

ref_samples <- c("Control")

ref_and_query_samples <- c("AKI", "CKD", "Control")

#first we load the mouse atlas object
seurat_object <- readRDS("/home/isilon/users/o_kloetzer/Atlas/Atlas_human_extension_II_convert.rds")

kegg_gene_matrix <- process_gene_sets_from_KEGG("hsa", seurat_object = seurat_object)
saveRDS(kegg_gene_matrix, "/home/isilon/users/o_kloetzer/Atlas/Revision/human_extended/output_human_extended/go_sets_modified.rds")

#metadat
metadat = read.table("/home/isilon/users/o_kloetzer/Atlas/Atlas_Extended_II_Albuminuria.csv", sep=",", header=TRUE)
rownames(metadat) <- metadat$Sample

# Add Disease_level2 information to each cell in the Seurat object
# Assuming 'orig_ident' in Seurat object matches 'Sample' in metadat
sample_to_disease <- metadat %>% select(Sample, Disease_level1) %>% distinct()
rownames(sample_to_disease) <- sample_to_disease$Sample

# Use 'orig_ident' in Seurat object to map each cell to its 'Disease_level2'
seurat_object <- seurat_object %>%
  AddMetaData(metadata = sample_to_disease[seurat_object@meta.data$orig_ident, "Disease_level1"], col.name = "Disease_level1")



# Create a new annotation vector based on the cell_type_major
seurat_object$disease_grouped_qr <- as.character(seurat_object$Disease_level1)


seurat_object$disease_grouped_qr[seurat_object$Disease_level1 %in% query_samples] <- "query_samples"
seurat_object$disease_grouped_qr[seurat_object$Disease_level1 %in% ref_samples] <- "reference_samples"

# Convert the grouped cell type back to a factor
seurat_object$disease_grouped_qr <- factor(seurat_object$disease_grouped_qr)

#define your sample
sample_to_change <- Sys.getenv("SAMPLE")
#sample_to_change <- "32-2"

# Define the sample identifier and the condition you want to change
sample_id_col <- "orig_ident"
condition_col <- "disease_grouped_qr"

#many empty samples - remove
#seurat_object@meta.data$donor_id <- droplevels(seurat_object@meta.data$donor_id)

# Convert the condition column to character if it is a factor
seurat_object@meta.data[[condition_col]] <- as.character(seurat_object@meta.data[[condition_col]])

# Change the condition key for the selected sample directly in the Seurat object
seurat_object@meta.data[[condition_col]][seurat_object@meta.data[[sample_id_col]] == sample_to_change & seurat_object@meta.data[["disease_grouped_qr"]] == "query_samples"] <- "sample_of_interest"

# Convert the condition column back to factor if necessary
seurat_object@meta.data[[condition_col]] <- as.factor(seurat_object@meta.data[[condition_col]])

#we use the prepare function to get our seurat object in shape
seurat_object <- prepare_seurat_object(
  seurat_obj = seurat_object, #our object to prepare
  celltype_col = "C_scANVI_modified1", #column in which annotation is stored
  sample_id_col = sample_id_col, #sample identifier
  condition_col = condition_col, #column storing the condition information
  query_list = c("sample_of_interest"), #categories in the condition_col defining the query samples
  control_list = c("reference_samples") #categories in the condition_col defining the reference
)


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
output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/human_extended/output_human_extended/"

cell_types = c(Sys.getenv("CELLTYPE"))

soi = Sys.getenv("SAMPLE")

run_spectra(output_folder_base = output_folder_base, cell_types = cell_types, CHISQ.MAX = 4, expression_threshold = -1, gene_number_threshold = 10,
            restrict_to_sample = TRUE, results_sample_key = soi, vector_of_samples = c(soi), QC_report = TRUE)

