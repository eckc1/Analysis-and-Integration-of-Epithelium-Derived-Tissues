library(dplyr)
library(tibble)
library(Matrix)
library(Seurat)
library(harmony)
#dat <- readRDS("UMAP_new_cellIdents.rds")
#dat <- readRDS("../merged_vdj_dat_batch_corrected.rds")
dat <- readRDS("Batch_Corrected_coords_removeFailedDemultiplex.rds")
head(dat)
dim(dat)
table(dat$tissue)

# Make sure cluster IDs are the active identity
Idents(dat) <- "seurat_clusters"

# Then plot
p1 <- DimPlot(
  dat,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE,
  label.size = 2.5
) + NoLegend() +
  ggtitle("UMAP of CSP and RNAseq data") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_reverse() +   # reverse x-axis
  scale_y_reverse()     # reverse y-axis

p1

# Make sure we're plotting ADT explicitly
Assays(dat)
DefaultAssay(dat) <- "ADT"

ADT_plot <- FeaturePlot(
  dat,
  features = c("Ms.TCR.Bchain","Ms.CD3","Ms.CD8a","Ms.CD45",
               "Ms.CD45.2","Ms.CD45RB","Ms.CD19","Ms.NK.1.1","Ms.F4-80"),
  reduction = "wnn.umap",
  min.cutoff = "q10",
  max.cutoff = "q99",
  cols = c("lightgrey","darkgreen"),
  order=TRUE,
  combine = FALSE
)
ADT_plot <- lapply(ADT_plot, function(p) p + scale_x_reverse() + scale_y_reverse())
ADT_plot <- patchwork::wrap_plots(ADT_plot)
ADT_plot

DefaultAssay(dat) <- "RNA"
RNA_test <- FeaturePlot(dat, features = c("Ptprc","Cd4", "Cd8a","Cd14", "Prf1", "Cd19", "Cd3e", "Foxp3"), reduction="wnn.umap",
                        min.cutoff="q10", max.cutoff="q90", order=TRUE, combine=TRUE)
RNA_test <- lapply(RNA_test, function(p) p + scale_x_reverse() + scale_y_reverse())
RNA_test <- patchwork::wrap_plots(RNA_test)


# Extract CD4 T cells only
clusters_to_keep <- c(18,16,8,0,12,13,27,11,34)

# Subset Seurat object to keep only those clusters
dat_filtered <- subset(dat, idents = clusters_to_keep)
colnames(dat_filtered@meta.data)
# Generate subsetted UMAP
p3 <- DimPlot(dat_filtered, reduction = 'wnn.umap', label = TRUE, group.by = "tissue", repel = TRUE, label.size = 2.5) + NoLegend()+
  ggtitle("CD4 T Cells Isolated")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_x_reverse() + 
  scale_y_reverse()  
p3



## Filter out low gene counts based on RNA assay

## Work on RNA assay
DefaultAssay(dat_filtered) <- "RNA"

## Join layers so we can access a single matrix
dat_filtered <- JoinLayers(dat_filtered, assay = "RNA")

## Use raw counts if available; otherwise use "data"
if ("counts" %in% Layers(dat_filtered[["RNA"]])) {
  mat <- GetAssayData(dat_filtered, assay = "RNA", layer = "counts")
} else {
  mat <- GetAssayData(dat_filtered, assay = "RNA", layer = "data")
}

## Compute fraction of cells expressing each gene (> 0)
n_cells   <- ncol(mat)
pct_cells <- Matrix::rowSums(mat > 0) / n_cells

## Keep genes expressed in at least 1% of cells
genes_keep <- names(pct_cells[pct_cells >= 0.01])

## (Optional sanity check)
message("Genes in genes_keep that are in RNA assay: ",
        sum(genes_keep %in% rownames(dat_filtered[["RNA"]])))

## Subset ONLY the RNA assay, leave ADT and others unchanged
rna_assay <- dat_filtered[["RNA"]]
rna_assay_sub <- subset(rna_assay, features = genes_keep)
dat_filtered[["RNA"]] <- rna_assay_sub

## (Optional: quick checks)
DefaultAssay(dat_filtered) <- "RNA"
message("Number of RNA features after filtering: ",
        length(Features(dat_filtered)))

if ("ADT" %in% names(dat_filtered@assays)) {
  DefaultAssay(dat_filtered) <- "ADT"
  message("Number of ADT features (should be unchanged): ",
          length(Features(dat_filtered)))
}


# Check assays
Assays(dat_filtered)

# check RNA genes
DefaultAssay(dat_filtered) <- "RNA"
length(Features(dat_filtered))

# Check ADT is still there and unchanged
DefaultAssay(dat_filtered) <- "ADT"
length(Features(dat_filtered))  # should match what you had before

#saveRDS(dat_filtered, "CD4_filtered_only.rds")


# Re-cluster CD4 T cells
dat_filtered <- readRDS("../CD4_filtered_only.rds")
head(dat_filtered)
dim(dat_filtered)
table(dat_filtered$tissue)

# RNA assay
DefaultAssay(dat_filtered) <- "RNA"
dat_filtered <- NormalizeData(dat_filtered)
dat_filtered <- FindVariableFeatures(dat_filtered)
dat_filtered <- ScaleData(dat_filtered)
dat_filtered <- RunPCA(dat_filtered)
# Then re-run Harmony / integration
# Harmony on RNA PCs (from 'pca')
dat_filtered <- RunHarmony(
  object = dat_filtered,
  group.by.vars = "sample_id",
  reduction = "pca",
  assay.use = "RNA",
  dims.use = 1:15,
  reduction.save = "harmony",
  verbose = FALSE
)

# ADT assay
DefaultAssay(dat_filtered) <- "ADT"
dat_filtered <- NormalizeData(dat_filtered, normalization.method = "CLR")
dat_filtered <- ScaleData(dat_filtered)
dat_filtered <- RunPCA(dat_filtered, reduction.name = "apca")
dat_filtered <- FindMultiModalNeighbors(
  dat_filtered,
  reduction.list = list("pca", "apca"),
  dims.list = list(1:30, 1:30)
)
# Harmony on ADT PCs (from 'apca')
dat_filtered <- RunHarmony(
  object = dat_filtered,
  group.by.vars = "sample_id",
  reduction = "apca",
  assay.use = "ADT",
  dims.use = 1:18,
  reduction.save = "aharmony",
  verbose = FALSE
)



dat_filtered <- RunUMAP(
  dat_filtered,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap"
)
dat_filtered <- FindClusters(
  dat_filtered,
  graph.name = "wsnn",
  resolution = 0.5   # adjustable
)
p3 <- DimPlot(
  dat_filtered,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE,
  label.size = 2.5
) + 
  NoLegend() +
  ggtitle("Reclustered UMAP – CD4 T Cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse() +
  scale_y_reverse()

p3

dat_filtered$tissue <- factor(
  dat_filtered$tissue,
  levels = c("GreenSpleen", "JLN", "RedSpleen", "Epithelium")
)

cb_palette <- c(
  "Epithelium" = "#E69F00",  # orange
  "JLN"        = "#009E73",  # bluish-green
  "RedSpleen"  = "#0072B2"  , # blue
  "GreenSpleen" = "#CC79A7"
)
p3 <- DimPlot(
  dat_filtered,
  reduction = "wnn.umap",
  group.by = "tissue",
  label = FALSE
) +
  ggtitle("CD4 T Cells (WNN UMAP, colored by tissue)") +
  #scale_color_manual(values = cb_palette) +
  scale_color_manual(values = cb_palette,
                     name = "Tissue",
                    breaks = c("GreenSpleen", "Epithelium", "JLN", "RedSpleen"),
                    labels = c(
                      "GreenSpleen" = "Epithelial Derived Spleen",
                      "Epithelium" = "Epithelium",
                      "JLN"        = "Non-Epithelial Derived JLN",
                      "RedSpleen"  = "Non-Epithelial Derived Spleen"
                    ) ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.key.height = unit(0.4, "cm")
  ) +
  scale_x_reverse() +
  scale_y_reverse()

p3

#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/VDJ_Analysis_Batch_Correcte/TCR_Figures/CD4_Filtered_only.png", plot = p3)

#saveRDS(dat_filtered, "CD4_filtered_only_clustered.rds")

p3 <- DimPlot(
  dat_filtered,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE,
  label.size = 2.5,
  group.by = "sample_id"
) + 
  NoLegend() +
  ggtitle("Reclustered UMAP – CD4 T Cells") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse() +
  scale_y_reverse()

p3
ggsave(filename = "CD4_Filtered_only_sample_ID.png", plot = p3)
