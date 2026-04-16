library(Seurat)
library(SeuratObject)
library(Matrix)
library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
library(ggrepel)   

# Load in CD4 T cells only
#dat <- readRDS("Batch_Corrected_coords_removeFailedDemultiplex.rds")
dat <- readRDS("CD4_filtered_only.rds")
dat_filtered <- dat
head(dat)
table(dat$sample_id)
#write.csv(table(dat$sample_id), "Muliplex_sampleCounts.csv")

p3 <- DimPlot(dat, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()+
  ggtitle("UMAP of CSP and RNAseq data")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_x_reverse() + 
  scale_y_reverse()  
p3


# Get all sample IDs
tbl <- table(dat$sample_id)

# Convert to data.frame for manipulation
df <- as.data.frame(tbl)
colnames(df) <- c("sample_id", "count")

# Summarize by sample type keywords
summary_df <- df %>%
  mutate(type = case_when(
    grepl("_RedSpleen", sample_id) ~ "RedSpleen",
    grepl("_GreenSpleen", sample_id) ~ "GreenSpleen",
    grepl("_JLN", sample_id) ~ "JLN",
    grepl("_epithelium", sample_id) ~ "Epithelium",
    TRUE ~ "Other"
  )) %>%
  group_by(type) %>%
  summarise(total_cells = sum(count)) %>%
  arrange(desc(total_cells))

summary_df
#write.csv(summary_df,"Multiplexed_cellCounts.csv")

# Foxp3 visualization
DefaultAssay(dat) <- "RNA"
RNA_test <- FeaturePlot(dat, features = c("Foxp3"), reduction="wnn.umap",
                        min.cutoff="q10", max.cutoff="q90", order=TRUE, combine=TRUE)
RNA_test <- lapply(RNA_test, function(p) p + scale_x_reverse() + scale_y_reverse())
RNA_test <- patchwork::wrap_plots(RNA_test)
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/Foxp3_Deseq_Results/foxp3_expression.png")


# Define the clusters of interest
clusters_to_keep <- c(11)

# Subset Seurat object to keep only those clusters
Idents(dat) <- "seurat_clusters"
dat_filtered <- subset(dat, idents = clusters_to_keep)


p3 <- DimPlot(dat_filtered, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()+
  ggtitle("UMAP of CSP and RNAseq data CD4 T Cells Isolated")+
  theme(plot.title = element_text(hjust=0.5))+
  scale_x_reverse() + 
  scale_y_reverse()  
p3

DefaultAssay(dat_filtered) <- "RNA"
foxp3_test <- FeaturePlot(dat_filtered, features = c("Foxp3"), reduction="wnn.umap",
                        min.cutoff="q10", max.cutoff="q90", order=TRUE, combine=TRUE)
foxp3_test <- lapply(foxp3_test, function(p) p + scale_x_reverse() + scale_y_reverse())
foxp3_test <- patchwork::wrap_plots(foxp3_test)
foxp3_test
dim(dat_filtered)
#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Batch_Correction/Foxp3_Deseq_Results/foxp3_expression_CD4.png")


# Check dimensions
dim(dat_filtered)
table(dat_filtered$tissue)
table(dat_filtered$seurat_clusters)
table(dat_filtered$sample_id)
sum(table(dat_filtered$seurat_clusters))
#write.csv(table(dat_filtered$sample_id), "CD4_TcellsOnlyMultiplexed.csv")

###### DEG Analysis with DESEQ ######

# Get full counts across all per-sample layers
count_layers <- grep("^counts", Layers(dat_filtered, assay = "RNA"), value = TRUE)
stopifnot(length(count_layers) > 0)

cnt_list <- lapply(count_layers, function(ly) {
  LayerData(dat_filtered, assay = "RNA", layer = ly)  # returns dgCMatrix (numeric @x)
})

# ensure consistent gene order
genes <- rownames(cnt_list[[1]])
cnt_list <- lapply(cnt_list, function(m) m[genes, , drop = FALSE])

# bind all cells; align to the Seurat object columns just in case
counts_all <- do.call(cbind, cnt_list)
counts_all  <- counts_all[, colnames(dat_filtered), drop = FALSE]

# Define groups and keep cells of interest
Idents(dat_filtered) <- "sample_id"
dat_filtered$tissue_group <- case_when(
  grepl("epithelium", dat_filtered$sample_id, ignore.case = TRUE) ~ "Epithelium",
  grepl("RedSpleen",  dat_filtered$sample_id) ~ "RedSpleen",
  grepl("GreenSpleen",dat_filtered$sample_id) ~ "GreenSpleen",
  grepl("JLN",        dat_filtered$sample_id) ~ "JLN",
  .default = NA_character_
)

keep_cells <- which(dat_filtered$tissue_group %in% c("Epithelium","RedSpleen","GreenSpleen","JLN"))
counts_all  <- counts_all[, keep_cells, drop = FALSE]
cell_meta   <- dat_filtered@meta.data[keep_cells, , drop = FALSE]

# PSEUDOBULK: sum counts per sample_id 
sample_id_vec <- droplevels(factor(cell_meta$sample_id))
X <- Matrix::sparse.model.matrix(~ 0 + sample_id_vec)
colnames(X) <- sub("^sample_id_vec", "", colnames(X))

pb_counts <- counts_all %*% X   # genes x samples (dgCMatrix)
nz <- Matrix::colSums(pb_counts) > 0
pb_counts <- pb_counts[, nz, drop = FALSE]

# Sample-level metadata
samples <- colnames(pb_counts)
length(samples)
samp2group <- tapply(cell_meta$tissue_group, sample_id_vec, function(v) {
  u <- unique(na.omit(v)); if (length(u) == 1) u else u[1]
})
tissue_group <- as.character(samp2group[samples])

colData_df <- data.frame(
  tissue_group = factor(tissue_group, levels = c("Epithelium","RedSpleen","GreenSpleen","JLN")),
  row.names = samples
)

# Give small dense integer matrix to DESeq2
pb_mat <- as.matrix(round(pb_counts))    
storage.mode(pb_mat) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = pb_mat,
  colData   = colData_df,
  design    = ~ tissue_group
)

# keep genes seen in >= 2 samples with >=1 count
keep <- rowSums(counts(dds) >= 1) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)

# Pairwise contrasts and volcano plots 
dir.create("DESeq_Results_Batch_Corrected", showWarnings = FALSE)
groups <- c("Epithelium","RedSpleen","GreenSpleen","JLN")
comparisons <- combn(groups, 2, simplify = FALSE)

# helper to shrink LFCs with contrast, use prefer ashr first
shrink_with_contrast <- function(dds, fac, g1, g2) {
  if (requireNamespace("ashr", quietly = TRUE)) {
    lfcShrink(dds, contrast = c(fac, g1, g2), type = "ashr")
  } else {
    lfcShrink(dds, contrast = c(fac, g1, g2), type = "normal")
  }
}

for (comp in comparisons) {
  g1 <- comp[1]; g2 <- comp[2]
  cname <- paste0(g1, "_vs_", g2)
  message("Processing: ", cname)
  
  # Results for p-values / padj
  res <- results(dds, contrast = c("tissue_group", g1, g2), alpha = 0.05)
  
  # Shrunk LFCs for better volcano visualization
  res_shrunk <- shrink_with_contrast(dds, "tissue_group", g1, g2)
  
  # Merge into one data frame
  res_df <- data.frame(
    gene           = rownames(res),
    log2FoldChange = res_shrunk$log2FoldChange,
    pvalue         = res$pvalue,
    padj           = res$padj,
    stringsAsFactors = FALSE
  ) |> na.omit()
  
  # significance and direction (padj + LFC thresholds)
  res_df <- res_df |>
    mutate(
      neg_log10_p = -log10(pvalue),
      status = case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange >=  0.5 ~ "Up",
        !is.na(padj) & padj < 0.05 & log2FoldChange <= -0.5 ~ "Down",
        TRUE                                                 ~ "NS"
      )
    )
  
  # label top genes among significant Up/Down only
  up_lab <- res_df |>
    filter(status == "Up") |>
    arrange(padj, desc(log2FoldChange)) |>
    slice_head(n = 15)
  
  down_lab <- res_df |>
    filter(status == "Down") |>
    arrange(padj, log2FoldChange) |>
    slice_head(n = 15)
  
  lab_df <- bind_rows(up_lab, down_lab)
  
  # consistent colors: Down=blue, Up=red, NS=grey
  volc_cols <- c("Down" = "blue", "Up" = "red", "NS" = "grey70")
  
  # use padj on y when available; fall back to p
  y_thr <- -log10(0.05)
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_p)) +
    geom_point(aes(color = status), alpha = 0.7, size = 1.6) +
    scale_color_manual(values = volc_cols) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    geom_hline(yintercept = y_thr, linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.size = 0.2,
      min.segment.length = 0
    ) +
    labs(
      title = cname,
      x = "log2(Fold Change) (shrunk)",
      y = "-log10(p-value)",
      color = "Direction"
    ) +
    theme_minimal(base_size = 13)
  
  # save outputs (use the same folder you created)
  write.csv(res_df, file.path("DESeq_Results_Batch_Corrected", paste0("DESeq2_", cname, ".csv")), row.names = FALSE)
  ggsave(file.path("DESeq_Results_Batch_Corrected", paste0("Volcano_", cname, ".png")),
         plot = p, width = 7, height = 6, dpi = 300)
  

message("Pseudobulk DESeq2 complete: ", ncol(pb_mat), " samples; ", nrow(dds), " genes after filtering.")
}


print("Differential Gene Expression Complete")