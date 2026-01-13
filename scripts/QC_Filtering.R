library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# Ignore local RAM throttle 
mem.maxVSize(vsize = Inf)

# Load in merged Seurat object
dat <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Intermediate_RDS_Objects/seurat_Sarah_merged_new_library.rds")
dim(dat)
head(dat)
dat@assays

# Check versions
packageVersion("Seurat"); packageVersion("SeuratObject")
original_cell_counts <- ncol(dat)
# Summary Metrics
cell_count <-table(dat$orig.ident)
cell_count
#write.csv(cell_count, "original_cell_count.csv")
unique(dat$orig.ident)

## Rename Spleen Samples to Green_Spleen or red spleen
spleen_samples <- c("861_spleen", "863_spleen", "864_spleen", "865_spleen")

# Update only those specific samples
dat$orig.ident[
  dat$orig.ident %in% spleen_samples
] <- paste0(
  sub("_spleen", "", dat$orig.ident[
    dat$orig.ident %in% spleen_samples
  ]),
  "_GreenSpleen"
)

dat$orig.ident[dat$orig.ident == "Pooled_spleen"] <- "Pooled_RedSpleen"

# Check names and cell counts
unique(dat$orig.ident)

#######################################################################
### De multiplexing
# Demultiplex pooled samples separately (tissue-matched)
# Requirements: Seurat v4/5, object `dat` with orig.ident like "861_spleen", "861_epithelium",
# and pooled labels exactly "Pooled_spleen", "Pooled_JLN".

# confidence cutoff
thr <- 0.40

# Mark which samples are pooled vs not pooled
is_pooled  <- dat$orig.ident %in% c("Pooled_RedSpleen", "Pooled_JLN")
non_pooled <- !is_pooled

# specify donor = numeric prefix (861, 863, 864, 865)
dat$donor <- NA_character_
dat$donor[non_pooled] <- sub("^([0-9]+).*", "\\1", dat$orig.ident[non_pooled])

# tissue classification
dat$tissue <- NA_character_
dat$tissue[grepl("_GreenSpleen$", dat$orig.ident)] <- "GreenSpleen"
dat$tissue[grepl("_epithelium$", dat$orig.ident)]  <- "Epithelium"
dat$tissue[dat$orig.ident == "Pooled_RedSpleen"]   <- "RedSpleen"
dat$tissue[dat$orig.ident == "Pooled_JLN"]         <- "JLN"

# Check Labels
table(dat$donor, useNA="ifany")
table(dat$tissue, useNA="ifany")
table(dat$orig.ident)

# Build References and queries for redspleen
ref_greenspleen <- subset(dat, subset = grepl("_GreenSpleen$", orig.ident))
qry_redspleen   <- subset(dat, subset = orig.ident == "Pooled_RedSpleen")

# Build References and queries for JLN
ref_epithelium <- subset(dat, subset = grepl("_epithelium$", orig.ident))
qry_jln        <- subset(dat, subset = orig.ident == "Pooled_JLN")

# Helper function for demultiplexing
.demux_one <- function(ref_obj, qry_obj, thr = 0.40, dims = 1:50) {
  if (ncol(qry_obj) == 0 || ncol(ref_obj) == 0) return(qry_obj)
  # Normalize and stabilize with SCT, create SCT assay
  if (!"SCT" %in% Assays(ref_obj)) ref_obj <- SCTransform(ref_obj, verbose = FALSE)
  if (!"SCT" %in% Assays(qry_obj)) qry_obj <- SCTransform(qry_obj, verbose = FALSE)
  
  # Compute PCAs
  ref_obj <- RunPCA(ref_obj, verbose = FALSE)
  qry_obj <- RunPCA(qry_obj, verbose = FALSE)
  
  # Specify reference identities
  Idents(ref_obj) <- ref_obj$donor
  
  # Use PCA embeddings for matching and use SCT normalization
  anchors <- FindTransferAnchors(
    reference = ref_obj, query = qry_obj,
    normalization.method = "SCT", dims = dims,
    reference.reduction = "pca", verbose = FALSE
  )
  # Use anchors to predict a donor level label for each cell
  preds <- TransferData(anchorset = anchors, refdata = ref_obj$donor, dims = dims, verbose = FALSE)
  # Attatch predictions to metadata
  qry_obj <- AddMetaData(qry_obj, metadata = preds)
  
  # Call donors purely by threshold
  qry_obj$predicted_donor <- ifelse(qry_obj$prediction.score.max >= thr, qry_obj$predicted.id, NA_character_)
  qry_obj$demux_call      <- ifelse(is.na(qry_obj$predicted_donor), "ambiguous", qry_obj$predicted_donor)
  # Return query object with new metadata
  qry_obj
}

# Run DEMUX function
qry_redspleen <- .demux_one(ref_greenspleen, qry_redspleen, thr = thr)
qry_jln       <- .demux_one(ref_epithelium,   qry_jln,       thr = thr)

# Merge query metadata with full seurat object
merge_cols <- c("predicted_donor","prediction.score.max","demux_call")
if (ncol(qry_redspleen) > 0) {
  dat@meta.data[colnames(qry_redspleen), merge_cols] <-
    qry_redspleen@meta.data[colnames(qry_redspleen), merge_cols]
}
if (ncol(qry_jln) > 0) {
  dat@meta.data[colnames(qry_jln), merge_cols] <-
    qry_jln@meta.data[colnames(qry_jln), merge_cols]
}

# Keep all cells and filter pooled cells based on DEMUX threshold
dat$keep_for_deg <- TRUE
dat$keep_for_deg[is_pooled] <- (dat$demux_call[is_pooled] != "ambiguous")

# Keep GreenSpleen and RedSpleen separate
dat$sample_id <- dat$orig.ident
kept_pooled <- which(is_pooled & dat$keep_for_deg)
dat$sample_id[kept_pooled] <- paste0(dat$predicted_donor[kept_pooled], "_", dat$orig.ident[kept_pooled])

# DEMUX summary output
cat("\n=== DEMUX SUMMARY ===\n")
cat("\nPooled_RedSpleen calls:\n")
print(table(dat$demux_call[dat$orig.ident == "Pooled_RedSpleen"], useNA="ifany"))
cat("Prediction score summary (Pooled_RedSpleen):\n")
print(summary(dat$prediction.score.max[dat$orig.ident == "Pooled_RedSpleen"]))
cat("Predicted donors (Pooled_RedSpleen, confident only):\n")
print(table(dat$predicted_donor[dat$orig.ident == "Pooled_RedSpleen" & dat$keep_for_deg], useNA="ifany"))

cat("\nPooled_JLN calls:\n")
print(table(dat$demux_call[dat$orig.ident == "Pooled_JLN"], useNA="ifany"))
cat("Prediction score summary (Pooled_JLN):\n")
print(summary(dat$prediction.score.max[dat$orig.ident == "Pooled_JLN"]))
cat("Predicted donors (Pooled_JLN, confident only):\n")
print(table(dat$predicted_donor[dat$orig.ident == "Pooled_JLN" & dat$keep_for_deg], useNA="ifany"))

# Check Demux assignment
unique(dat$sample_id)
table(dat$sample_id)

#saveRDS(dat,"DEMUX.rds")


# Sesnsitivity analysis to determine correct theshold
obj <- dat
# Re-label pooled calls at multiple thresholds using existing predictions
relabel <- function(obj, thr) {
  ii <- obj$orig.ident %in% c("Pooled_RedSpleen","Pooled_JLN")
  calls <- ifelse(obj$prediction.score.max >= thr, obj$sample_id, NA)
  out <- table(sample = obj$orig.ident[ii],
               call   = ifelse(is.na(calls[ii]), "ambiguous", calls[ii]))
  list(thr = thr, counts = out)
}
res <- lapply(c(0.40, 0.50, 0.55, 0.60), \(t) relabel(dat, t))
# Print how assignments change vs threshold
res  


#######################################################################
### Doublet Removal wit Doublet Finder

dat <- readRDS("DEMUX.rds")
Assays(dat)
head(dat)
## DoubletFinder with Seurat v5 (per-sample)
suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(dplyr)
  library(tibble)
})

# 10x multiplet rates table
.tenx_multiplet_rates <- function() {
  data.frame(
    Multiplet_rate = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
    Loaded_cells   = c(  800,  1600,  3200,  4800,  6400,  8000,  9600, 11200, 12800, 14400, 16000),
    Recovered_cells= c(  500,  1000,  2000,  3000,  4000,  5000,  6000,  7000,  8000,  9000, 10000)
  )
}

# Safe paramSweep wrapper to call correct version of doublet finder
.param_sweep_safe <- function(seu, PCs, sct = FALSE) {
  if ("paramSweep" %in% ls(getNamespace("DoubletFinder"))) {
    return(paramSweep(seu, PCs = PCs, sct = sct))
  } else if ("paramSweep_v3" %in% ls(getNamespace("DoubletFinder"))) {
    return(paramSweep_v3(seu, PCs = PCs, sct = sct))
  } else {
    stop("No paramSweep{,_v3} found in DoubletFinder installation.")
  }
}

# Calculate PCs until cumulative variance > 90% and elbow-based addition is under 0.1%
.pick_min_pc <- function(seu) {
  stdv <- seu[["pca"]]@stdev
  if (is.null(stdv) || length(stdv) == 0L) stop("PCA stdev not found; did RunPCA succeed?")
  percent_var    <- (stdv^2 / sum(stdv^2)) * 100
  cumulative_var <- cumsum(percent_var)
  co1 <- which(cumulative_var > 90 & percent_var < 5)[1]
  co2 <- suppressWarnings( which(diff(percent_var) < 0.1)[1] + 1 )
  co1 <- ifelse(is.na(co1), length(stdv), co1)
  co2 <- ifelse(is.na(co2), length(stdv), co2)
  min(co1, co2)
}

#  Robust HVG step for small samples (compute variance directly if VST fails)
.safe_hvg <- function(seu, nfeatures = 2000) {
  DefaultAssay(seu) <- "RNA"
  out <- tryCatch(
    FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE),
    warning = function(w) FindVariableFeatures(seu, selection.method = "mean.var.plot",
                                               nfeatures = min(1000, max(300, nrow(seu))), verbose = FALSE),
    error   = function(e) FindVariableFeatures(seu, selection.method = "mean.var.plot",
                                               nfeatures = min(1000, max(300, nrow(seu))), verbose = FALSE)
  )
  out
}

# Per-sample runner for doublet finder
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL) {
  # Use SampleID if present, else orig.ident
  sample_label <- NULL
  if ("SampleID" %in% colnames(seu_sample_subset@meta.data)) {
    sample_label <- unique(seu_sample_subset$SampleID)
  } else if ("sample_id" %in% colnames(seu_sample_subset@meta.data)) {
    sample_label <- unique(seu_sample_subset$sample_id)
  } else if ("orig.ident" %in% colnames(seu_sample_subset@meta.data)) {
    sample_label <- unique(seu_sample_subset$orig.ident)
  } else sample_label <- "sample"
  
  message(sprintf(">> Running DoubletFinder on sample: %s (n=%d cells)", paste(sample_label, collapse=","), ncol(seu_sample_subset)))
  
  # Very small samples are unreliable; mark singlet and return
  if (ncol(seu_sample_subset) < 60) {
    warning(sprintf("Sample %s has <60 cells; skipping DF and marking all as Singlet.", paste(sample_label, collapse=",")))
    df <- data.frame(doublet_finder = "Singlet", row_names = colnames(seu_sample_subset))
    return(df)
  }
  
  DefaultAssay(seu_sample_subset) <- "RNA"
  if (!"data" %in% Layers(seu_sample_subset, assay = "RNA"))
    seu_sample_subset <- NormalizeData(seu_sample_subset, verbose = FALSE)
  
  # Standard Seurat preprocessing, but v5-safe and robust
  sample <- NormalizeData(seu_sample_subset, verbose = FALSE)
  sample <- .safe_hvg(sample, nfeatures = 2000)
  hv <- VariableFeatures(sample)
  if (length(hv) < 100) stop("Too few HVGs after filtering; check QC thresholds.")
  
  sample <- ScaleData(sample, features = hv, verbose = FALSE)
  
  # Choose a safe number of PCs to compute
  npcs_safe <- max(10, min(50, ncol(sample) - 1L, length(hv) - 1L))
  if (npcs_safe < 10) stop("Not enough cells/features to compute 10 PCs.")
  
  sample <- RunPCA(sample, features = hv, npcs = npcs_safe, verbose = FALSE)
  
  # Pick PCs (min_pc), but cap to computed PCs
  min_pc <- .pick_min_pc(sample)
  min_pc <- max(5, min(min_pc, npcs_safe))
  
  sample <- RunUMAP(sample, dims = 1:min_pc, verbose = FALSE)
  sample <- FindNeighbors(sample, dims = 1:min_pc, verbose = FALSE)
  sample <- FindClusters(sample, resolution = 0.1, verbose = FALSE)
  
  # pK sweep & selection using paramSweep -> summarizeSweep -> find.pK)
  sweep_list  <- .param_sweep_safe(sample, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  
  # Robust pK extraction
  optimal.pk <- tryCatch({
    pkv <- bcmvn$pK[which.max(bcmvn$BCmetric)]
    as.numeric(as.character(pkv))
  }, error = function(e) NA_real_)
  
  if (!is.finite(optimal.pk)) {
    warning("find.pK returned NA; falling back to pK = 0.005")
    optimal.pk <- 0.005
  }
  
  # Homotypic proportion from clusters
  annotations <- sample$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # Multiplet rate: estimate from 10x table by recovered cells
  if (is.null(multiplet_rate)) {
    rates <- .tenx_multiplet_rates()
    # pick the largest Recovered_cells < current nCells
    multiplet_rate <- rates %>%
      dplyr::filter(Recovered_cells < ncol(sample)) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::pull(Multiplet_rate)
    if (length(multiplet_rate) == 0) multiplet_rate <- 0.008  # mild default
    message(sprintf("  Using estimated multiplet_rate = %.3f", multiplet_rate))
  } else {
    message(sprintf("  Using user-specified multiplet_rate = %.3f", multiplet_rate))
  }
  
  nExp.poi <- round(multiplet_rate * ncol(sample))
  nExp.poi.adj <- max(1L, round(nExp.poi * (1 - homotypic.prop)))
  
  # Run DoubletFinder (single pass)
  sample <- doubletFinder(
    seu = sample, PCs = 1:min_pc, pK = optimal.pk, nExp = nExp.poi.adj
  )
  
  # Standardize the output column name to "doublet_finder"
  df_cols <- grep("^DF.classifications", colnames(sample@meta.data), value = TRUE)
  if (length(df_cols) == 0) stop("DoubletFinder did not generate DF.classifications column.")
  colnames(sample@meta.data)[colnames(sample@meta.data) == tail(df_cols, 1)] <- "doublet_finder"
  
  out <- sample@meta.data["doublet_finder"] |> rownames_to_column("row_names")
  out
}


# Apply to dat object

# Choose column sample_id or else orig.ident
split_col <- if ("sample_id" %in% colnames(dat@meta.data)) "sample_id" else "orig.ident"

samp_list <- SplitObject(dat, split.by = split_col)

# Run per-sample
set.seed(123)
df_calls_list <- lapply(samp_list, run_doubletfinder_custom)

# Merge per-sample results to one data.frame, add to metadata
df_calls <- bind_rows(df_calls_list)
rownames(df_calls) <- df_calls$row_names
df_calls$row_names <- NULL

dat <- AddMetaData(dat, df_calls, col.name = "doublet_finder")

# Summary table
doublets_summary <- dat@meta.data %>%
  mutate(.sample = .data[[split_col]]) %>%
  group_by(.sample, doublet_finder) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  group_by(.sample) %>%
  mutate(countT = sum(total_count),
         percent = paste0(round(100 * total_count / countT, 2), "%")) %>%
  select(.sample, doublet_finder, total_count, percent) %>%
  arrange(.sample, doublet_finder)

print(doublets_summary)
write.csv(doublets_summary, "doublets_summary.csv")

# Remove doublets (returns singlets only)
dat_singlets <- subset(dat, subset = doublet_finder == "Singlet")

# Save both 
 #saveRDS(dat, "seurat_with_doubletfinder.rds")
 #saveRDS(dat_singlets, "seurat_after_doubletfinder_singlets_only.rds")
length(colnames(dat_singlets))
length(colnames(dat))
table(dat_singlets$orig.ident)
table(dat_singlets$tissue)
table(dat$tissue)

table(dat_singlets$sample_id)
table(dat$sample_id)

#write.csv(dat_singlets@meta.data, "metadata.csv")

#######################################################################
### QC CONTROL

# Load in data file after doublet removal and demultiplexing

dat <- readRDS("/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Intermediate_RDS_Objects/seurat_after_doubletfinder_singlets_only.rds")
#dat <- dat_singlets
dim(dat)
head(dat)
dat@assays
# Check versions
packageVersion("Seurat"); packageVersion("SeuratObject")
original_cell_counts <- ncol(dat)

# Summary Metrics
cell_count <-table(dat$tissue)
cell_count

#write.csv(cell_count, "NoDoublets_cell_count.csv")
unique(dat$orig.ident)
length(colnames(dat))
table(dat$sample_id)

# RNA filtering QC

# Calculate percent mitochondrial per cell
#dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
# Mouse
dat[["percent.mt"]] <- PercentageFeatureSet(
  dat,
  #pattern = "(?i)^mt-",
  pattern = "mt-",
  assay="RNA"
)


# Normalize and set default layer

# Visualize distribution
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0.1) +  theme_minimal()

#ggsave("QC_violin.pdf", plot = p, width = 6, height = 4)


# Visualize distribution for only mt
p <- VlnPlot(dat, features = c("percent.mt"), 
             ncol = 3, pt.size = 0.1)+  theme_minimal()
#ggsave("percent_mt_violin.pdf", plot = p, width = 6, height = 4)

# Visualize distribution for only rp
#p <- VlnPlot(dat, features = c("percent.rp"), 
             #ncol = 3, pt.size = 0.1)+  theme_minimal()
#ggsave("percent_mt_violin.pdf", plot = p, width = 6, height = 4)



# Apply QC filter: keep cells with 200–2500 genes and <5% mito
dat <- subset(dat, 
                       subset = nFeature_RNA > 200 & 
                         #nFeature_RNA < 2500 & 
                         percent.mt < 5)
dat@assays
# Visualize filtered distribution
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
             ncol = 3, pt.size = 0.1) +
  theme_minimal()
#ggsave("QC_violin.pdf_filtered.pdf", plot = p, width = 6, height = 4)

# Visualize mito only filtered distribution
p <- VlnPlot(dat, features = c("percent.mt"), 
             ncol = 3, pt.size = 0.1) +
  theme_minimal()
#ggsave("percent_mt_violin_filtered.pdf", plot = p, width = 6, height = 4)




# Remove MT genes from the RNA assay
DefaultAssay(dat) <- "RNA"
rna_feats <- rownames(dat[["RNA"]])
mt_genes <- grep("(?i)^mt-", rna_feats, value = TRUE)   # case-insensitive
length(mt_genes)
if (length(mt_genes) > 0) {
  keep_feats <- setdiff(rna_feats, mt_genes)
  # Subset just the RNA assay's features
  dat[["RNA"]] <- subset(dat[["RNA"]], features = keep_feats)
}
dim(dat)

if ("RNA" %in% names(Assays(dat))) {
  cat("RNA features:", nrow(dat[["RNA"]]), 
      "| Cells:", ncol(dat[["RNA"]]), "\n")
}
if ("ADT" %in% names(Assays(dat))) {
  cat("ADT features:", nrow(dat[["ADT"]]), 
      "| Cells:", ncol(dat[["ADT"]]), "\n")
}
dim(dat)

###### REMOVE RIBOSOMAL GENES ##############

# Identify ribosomal genes (case-insensitive)
ribo_genes <- grep("(?i)^RPS|^RPL", rna_feats, value = TRUE)
length(ribo_genes)

# Subset to remove ribosomal genes
if (length(ribo_genes) > 0) {
  keep_feats <- setdiff(rna_feats, ribo_genes)
  dat[["RNA"]] <- subset(dat[["RNA"]], features = keep_feats)
}
# Print RNA and ADT dimensions
if ("RNA" %in% names(Assays(dat))) {
  cat("RNA features:", nrow(dat[["RNA"]]), 
      "| Cells:", ncol(dat[["RNA"]]), "\n")
}
if ("ADT" %in% names(Assays(dat))) {
  cat("ADT features:", nrow(dat[["ADT"]]), 
      "| Cells:", ncol(dat[["ADT"]]), "\n")
}

dim(dat)


# Cell Count and Gene Lists
cell_count <-table(dat$sample_id)
length(colnames(dat))
#write.csv(cell_count, "cell_count_filtered_demulitplexed.csv")
cell_count <-table(dat$tissue)
#write.csv(cell_count, "tissues_cell_count_filtered_demulitplexed.csv")

has_rna <- dat@assays$RNA
length(rownames(has_rna))
write.csv(rownames(has_rna), "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/Final_Figures/RNA_Genes.csv")

has_adt <- dat@assays$ADT
#write.csv(rownames(has_adt), "CSP_Genes.csv")


# Check cell counts
cat("Original cells:", original_cell_counts, "\n")
cat("Filtered cells:", ncol(dat), "\n")

# ADT filtering QC
saveRDS(dat, "NoRibosomal_Genes_seurat_Sarah_merged_MitoFiltered_Demultiplexed_no_doublets.rds")



### END QC PIPELINE ###






