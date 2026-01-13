library(Seurat)
library(dplyr)
library(tibble)
library(stringr)

# Set parent path to sample folders
base_dir <- "/scratch/alpine/ceck@xsede.org/Kuhn_Data/Sarah_Data/New_Library_CSP"
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
print(sample_dirs)

# List to hold Seurat objects
seurat_list <- list()

# Function to process each sample
process_sample <- function(sample_path) {
  sample_name <- basename(sample_path)
  message("Processing sample: ", sample_name)

  # Define output and matrix paths
  output_folder <- file.path(sample_path, paste0(sample_name, "_output"))
  gex_path <- file.path(output_folder, "outs", "per_sample_outs", paste0(sample_name, "_output"), "count", "sample_filtered_feature_bc_matrix")
  vdj_path <- file.path(output_folder, "outs", "per_sample_outs", paste0(sample_name, "_output"), "vdj_t", "filtered_contig_annotations.csv")

  if (!dir.exists(gex_path)) stop("GEX path not found: ", gex_path)

  # Load GEX and ADT matrices
  gex <- Read10X(gex_path)

  if (is.list(gex) && "Gene Expression" %in% names(gex)) {
    rna <- gex[["Gene Expression"]]
  } else {
    stop("Gene Expression matrix not found")
  }

  adt <- if ("Antibody Capture" %in% names(gex)) gex[["Antibody Capture"]] else NULL

  # Create Seurat object
  seu <- CreateSeuratObject(counts = rna, project = sample_name)

  if (!is.null(adt)) {
    seu[["ADT"]] <- CreateAssayObject(counts = adt)
  }

  # Load and merge VDJ if available
  if (file.exists(vdj_path)) {
    vdj <- read.csv(vdj_path)
    vdj <- vdj %>%
      filter(productive == "true") %>%
      group_by(barcode) %>%
      slice(1) %>%
      ungroup()

    matching_barcodes <- vdj$barcode %in% colnames(seu)
    if (any(matching_barcodes)) {
      vdj <- vdj[matching_barcodes, ]
      rownames(vdj) <- vdj$barcode
      metadata <- seu@meta.data %>% rownames_to_column("barcode")
      metadata <- left_join(metadata, vdj, by = "barcode") %>% column_to_rownames("barcode")
      seu@meta.data <- metadata
    }
  }

  return(seu)
}

# Loop through all samples
for (sample_path in sample_dirs) {
  seurat_obj <- tryCatch({
    process_sample(sample_path)
  }, error = function(e) {
    message("Error processing ", sample_path, ": ", e$message)
    return(NULL)
  })

  if (!is.null(seurat_obj)) {
    seurat_list[[basename(sample_path)]] <- seurat_obj
  }
}

# Save each Seurat object
for (name in names(seurat_list)) {
  saveRDS(seurat_list[[name]], file = paste0("seurat_", name, ".rds"))
}

message("Script complete.")
