 library(Seurat)
library(dplyr)
library(ggplot2)
library(harmony)

# Read in QC-filtered object
dat <- readRDS("NoRibosomal_Genes_seurat_Sarah_merged_MitoFiltered_Demultiplexed_no_doublets.rds")
head(dat)
dim(dat)
dim(dat@meta.data$tissue)
colnames(dat@meta.data)
table(dat@meta.data$tissue)
table(dat@meta.data$sample_id)

#mem.maxVSize(vsize = Inf)
#write.csv(dat@meta.data, "dat_metadata.csv")
table(dat@meta.data$orig.ident)
dat@assays

# Remove all cells that were not demultiplexed properly
dat <- subset(dat, subset = !sample_id %in% c("Pooled_JLN", "Pooled_RedSpleen"))
dim(dat)

## ADT QC
DefaultAssay(dat) <- "ADT"

# Drop any ADT detected in < 1% of cells
adt_counts <- GetAssayData(dat, slot = "counts")
nz_frac <- Matrix::rowMeans(adt_counts > 0)

# Keep ADTs that appear in ≥1% of cells
keep_ADTs <- names(which(nz_frac >= 0.01))
VariableFeatures(dat[["ADT"]]) <- keep_ADTs

cat("Keeping", length(keep_ADTs), "ADT features:\n")
print(keep_ADTs)
adt_data <- GetAssayData(dat, slot = "data")
adt_var <- apply(as.matrix(adt_data), 1, var)
keep_ADTs <- names(which(nz_frac >= 0.01 & adt_var > 0))
VariableFeatures(dat[["ADT"]]) <- keep_ADTs


DefaultAssay(dat) <- 'RNA'
dat<- NormalizeData(dat) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(reduction.name = 'pca')

DefaultAssay(dat) <- 'ADT'
VariableFeatures(dat) <- rownames(dat[["ADT"]])
dat <- NormalizeData(dat, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


###############################################
######## Add harmony batch correction #########

# Harmony batch correction by sample_id + WNN/UMAP colored by sample_id
if (!"sample_id" %in% colnames(dat@meta.data)) {
  dat$sample_id <- dat$orig.ident
}

set.seed(4537)

# Harmony on RNA PCs (from 'pca')
dat <- RunHarmony(
  object = dat,
  group.by.vars = "sample_id",
  reduction = "pca",
  assay.use = "RNA",
  dims.use = 1:15,
  reduction.save = "harmony",
  verbose = FALSE
)

# Harmony on ADT PCs (from 'apca')
dat <- RunHarmony(
  object = dat,
  group.by.vars = "sample_id",
  reduction = "apca",
  assay.use = "ADT",
  dims.use = 1:18,
  reduction.save = "aharmony",
  verbose = FALSE
)

# Use harmonized reductions for WNN
dat <- FindMultiModalNeighbors(
  dat,
  reduction.list = list("harmony", "aharmony"),
  dims.list = list(1:15, 1:18),
  modality.weight.name = "RNA.weight",
  k.nn = 20
)

# UMAP + clustering
dat <- RunUMAP(dat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
dat <- FindClusters(dat, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)


# Keep clusters with >= 200 cells (as before)
clust_counts <- table(Idents(dat))
keep_clusters <- names(clust_counts)[clust_counts >= 200]
dat <- subset(dat, idents = keep_clusters)

saveRDS(dat, "Batch_Corrected_coords_removeFailedDemultiplex.rds")

#colnames(dat@meta.data)


table(dat$sample_id)
table(dat$tissue)

# Plots
p_by_sample <- DimPlot(
  dat, reduction = "wnn.umap", group.by = "sample_id", label = FALSE
) + ggtitle("WNN UMAP (Harmony-corrected), colored by sample_id") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_reverse() +   # reverse x-axis
  scale_y_reverse()     # reverse y-axis
ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/batch_by_sample.png", plot = p_by_sample)


p_by_tissue <- DimPlot(
  dat, reduction = "wnn.umap", group.by = "tissue", label = FALSE
) + ggtitle("WNN UMAP (Harmony-corrected), colored by sample_id") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_reverse() +   # reverse x-axis
  scale_y_reverse()     # reverse y-axis
ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/batch_by_tissue.png", plot = p_by_tissue)


p_by_cluster <- DimPlot(
  dat, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 2.5
) + NoLegend() + ggtitle("WNN UMAP (Harmony-corrected), clusters") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_reverse() +   # reverse x-axis
  scale_y_reverse()     # reverse y-axis
ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/batch_by_cluster.png", plot = p_by_cluster)

p_by_dataset <- DimPlot(
  dat, reduction = "wnn.umap", label = TRUE, group.by = "dataset", repel = TRUE, label.size = 2.5
) + NoLegend() + ggtitle("WNN UMAP (Harmony-corrected), clusters") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_reverse() +   # reverse x-axis
  scale_y_reverse()     # reverse y-axis
ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/batch_by_dataset.png", plot = p_by_dataset)


p_by_dataset <- DimPlot(
  dat, reduction = "wnn.umap", label = TRUE, group.by = "dataset",
  repel = TRUE, label.size = 2.5
) + 
  NoLegend() +
  ggtitle("WNN UMAP (Harmony-corrected), clusters") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_reverse() +   # reverse x-axis
  scale_y_reverse()     # reverse y-axis
ggsave(
  filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/batch_by_dataset.png",
  plot = p_by_dataset
)


p_by_sample
p_by_cluster
p_by_tissue

head(dat)
table(dat$seurat_clusters)
table(dat$orig.ident)
table(dat$tissue)
unique(dat$dataset)

###########################################################
# Make sure ADT assay is used for plotting
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
ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/CSP_featurePlotmultiplexed.png", plot = ADT_plot)


DefaultAssay(dat) <- "RNA"
RNA_test <- FeaturePlot(dat, features = c("Ptprc","Cd4", "Cd8a","Cd14", "Prf1", "Cd19", "Cd3e", "Foxp3"), reduction="wnn.umap",
                        min.cutoff="q10", max.cutoff="q90", order=TRUE, combine=TRUE)
RNA_test <- lapply(RNA_test, function(p) p + scale_x_reverse() + scale_y_reverse())
RNA_test <- patchwork::wrap_plots(RNA_test)

ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/RNA_featurePlotmultiplexed.png", plot = RNA_test)




print("End Weighted Nearest Neighbor Analysis Pipeline")








