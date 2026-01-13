library(Seurat)
# Load Data
seurat_861_epi   <- readRDS("seurat_861_epithelium.rds")
seurat_861_spl   <- readRDS("seurat_861_spleen_resequenced.rds")
seurat_863_epi   <- readRDS("seurat_863_epithelium.rds")
seurat_863_spl   <- readRDS("seurat_863_spleen_resequenced.rds")
seurat_864_epi   <- readRDS("seurat_864_epithelium.rds")
seurat_864_spl   <- readRDS("seurat_864_spleen_resequenced.rds")
seurat_865_epi   <- readRDS("seurat_865_epithelium.rds")
seurat_865_spl   <- readRDS("seurat_865_spleen_resequenced.rds")
seurat_pooled_jln <- readRDS("seurat_Pooled_JLN.rds")
seurat_pooled_spl <- readRDS("seurat_Pooled_spleen.rds")


# Add sample ID
seurat_861_epi$dataset   <- "861_epithelium"
seurat_861_spl$dataset   <- "861_spleen"
seurat_863_epi$dataset   <- "863_epithelium"
seurat_863_spl$dataset   <- "863_spleen"
seurat_864_epi$dataset   <- "864_epithelium"
seurat_864_spl$dataset   <- "864_spleen"
seurat_865_epi$dataset   <- "865_epithelium"
seurat_865_spl$dataset   <- "865_spleen"
seurat_pooled_jln$dataset <- "Pooled_JLN"
seurat_pooled_spl$dataset <- "Pooled_spleen"

merged_obj <- merge(
  seurat_861_epi,
  y = list(
    seurat_861_spl,
    seurat_863_epi,
    seurat_863_spl,
    seurat_864_epi,
    seurat_864_spl,
    seurat_865_epi,
    seurat_865_spl,
    seurat_pooled_jln,
    seurat_pooled_spl
  ),
  add.cell.ids = c(
    "861Epi", "861Spl",
    "863Epi", "863Spl",
    "864Epi", "864Spl",
    "865Epi", "865Spl",
    "PooledJLN", "PooledSpl"
  ),
  project = "Sarah_Merged"
)

saveRDS(merged_obj, file = "seurat_Sarah_merged_new_library.rds")