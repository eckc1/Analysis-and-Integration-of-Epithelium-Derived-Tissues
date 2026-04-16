library(tibble)
library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(Seurat)

dat <- readRDS("../CD4_filtered_only_clustered.rds")
original_metadata <- read.csv("/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/TCR_Data/Original_TCR_meta.csv")

## Re-merge VDJ data ##

# Start from Seurat metadata
meta <- dat@meta.data %>%
  rownames_to_column("cell_id")

# Join on barcodes
merged_meta <- meta %>%
  left_join(original_metadata, by = c("cell_id" = "X"))

#Put cell_id back as rownames and drop the column
merged_meta <- merged_meta %>%
  column_to_rownames("cell_id")

# Replace Seurat metadata
dat@meta.data <- merged_meta

dim(dat@meta.data)
colnames(dat@meta.data)
table(dat@meta.data$cdr3, useNA = "ifany")

sum(rownames(dat@meta.data) %in% original_metadata$X)
length(intersect(rownames(dat@meta.data), original_metadata$X))

# Create a table of the top 20 clonotypes, meaning the most abundant values in the cdr3 column

# Extract cdr3 column
cdr3_vals <- dat$cdr3.y

# Remove missing or empty values
cdr3_vals <- cdr3_vals[!is.na(cdr3_vals) & cdr3_vals != ""]

# Count frequencies
cdr3_table <- sort(table(cdr3_vals), decreasing = TRUE)

# Convert to data frame
cdr3_df <- as.data.frame(cdr3_table)
colnames(cdr3_df) <- c("cdr3", "count")
#write.csv(cdr3_df, "cdr3_df.csv")

#write.csv(cdr3_df, "cdr3_counts.csv")
# Take top 20
top20_df <- head(cdr3_df, 20)

top20_df


### Clonotype Bar plot

# Start from metadata
meta <- dat@meta.data %>% 
  dplyr::select(clonotype = cdr3.y,
         group     = tissue) %>%    
  filter(!is.na(clonotype))

# Keep only the tissue groups you want
keep_groups <- c("Epithelium", "JLN", "RedSpleen")

# Count cells per clonotype × group, after dropping GreenSpleen
df <- meta %>%
  filter(group %in% keep_groups) %>%  
  group_by(clonotype, group) %>%
  summarise(n = n(), .groups = "drop")

# Convert to percentages within each clonotype
df <- df %>%
  group_by(clonotype) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# Identify expanded clonotypes 
#expanded_clones <- df %>% 
  #group_by(clonotype) %>% 
  #summarise(total = sum(n), .groups = "drop") %>%
  #filter(total >= 2) %>%      # threshold
  #pull(clonotype)

# Remove clones in only one group
expanded_clones <- df %>% 
  group_by(clonotype) %>%
  filter(n_distinct(group) > 1) %>%
  pull(clonotype) %>%
  unique()


# Subset to expanded clonotypes
df_sub <- df %>% 
  filter(clonotype %in% expanded_clones)


head(df_sub)
dim(df_sub)
table(df_sub$group)

# Create clonotype key with CT labels
clonotype_key <- df_sub %>%
  distinct(clonotype) %>%
  arrange(clonotype) %>%    
  mutate(clonotype_id = paste0("CT", row_number())) %>%
  dplyr::select(clonotype_id, clonotype)

# Join CT labels into plotting dataframe
df_plot <- df_sub %>%
  left_join(clonotype_key, by = "clonotype")

#write.csv(clonotype_key, "clonotype_key_CT1_CT15.csv", row.names = FALSE)

cb_palette <- c(
  "Epithelium" = "#E69F00",  # orange
  "JLN"        = "#009E73",  # bluish-green
  "RedSpleen"  = "#0072B2"   # blue
)
plot <- ggplot(
  df_plot,
  aes(
    x = factor(clonotype_id, levels = clonotype_key$clonotype_id),
    y = pct,
    fill = group
  )
) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = cb_palette,
                    name= "Tissue",
                    breaks = c("Epithelium", "JLN", "RedSpleen"),
                    labels = c(
                      "Epithelium" = "Epithelium",
                      "JLN"        = "Non-Epithelial Derived JLN",
                      "RedSpleen"  = "Non-Epithelial Derived Spleen"
                   ) ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Clonotype",
    y = "Percent Cells Within Clonotype"
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title   = element_text(size = 12),
    legend.position = "right"
  )
plot

#ggsave(filename = "/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/VDJ_Analysis_Batch_Correcte/TCR_Figures/Top_clonotypes_multiCell_nofilter_formatted.png", plot = plot)

# Create plot with ggplot
#plot <- ggplot(df_sub, aes(x = factor(clonotype), y = pct, fill = group)) +
  #geom_bar(stat = "identity", width = 0.8) +
  #scale_y_continuous(expand = c(0,0)) +
  #labs(x = "Clonotype", y = "% cells among clonotype") +
  #theme_classic() +
 # theme(
   # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
   # legend.title = element_blank(),
   # axis.title = element_text(size = 12),
   # legend.position = "right"
  #)

#ggsave(filename = "Top_clonotypes_multiCell_nofilter.png", plot = plot)


## Dot plot of gene expression 

## Clonotypes to plot 
clones_use <- unique(df_sub$clonotype)

DefaultAssay(dat) <- "RNA"

## Gene list

gene_list <- c("Foxp3", "Cd44", "Nt5e", "Izumo1r")

# Get RNA expression matrix
rna_data <- GetAssayData(dat, assay = "RNA", slot = "data")
#write.csv(rownames(rna_data), "genes_test.csv")
# Keep only genes that actually exist in the assay
valid_genes   <- intersect(gene_list, rownames(rna_data))
missing_genes <- setdiff(gene_list, valid_genes)

if (length(valid_genes) == 0) {
  stop("None of the genes in gene_list are present in the RNA assay.")
}
if (length(missing_genes) > 0) {
  message("These genes were not found and will be skipped: ",
          paste(missing_genes, collapse = ", "))
}

expr_mat <- rna_data[valid_genes, , drop = FALSE]

## Tissues to plot
keep_tissues <- c("RedSpleen", "JLN", "Epithelium")

## Build dataframe with metadata + expression (long format)
expr_df <- expr_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "cell",
    values_to = "expr"
  ) %>%
  left_join(
    dat@meta.data %>% 
      tibble::rownames_to_column("cell") %>%
      dplyr::select(
        cell,
        clonotype = cdr3.y,
        tissue    = tissue   
      ),
    by = "cell"
  ) %>%
  filter(
    !is.na(clonotype),
    clonotype %in% clones_use,
    tissue %in% keep_tissues
  )

## Summarize per tissue × clonotype × gene
dot_df <- expr_df %>%
  group_by(tissue, clonotype, gene) %>%
  summarise(
    pct_exp = mean(expr > 0) * 100,
    avg_exp = mean(expr),
    .groups = "drop"
  ) %>%
  # Include all clones, regardless of presence in tissue
  #complete(
    #tissue = keep_tissues,
    #clonotype = clones_use,
    #gene = valid_genes,
    #fill = list(
     # pct_exp = 0,
     # avg_exp = 0
    #)
  #) %>%
  mutate(
    clonotype = factor(clonotype, levels = rev(clones_use)),
    tissue    = factor(tissue, levels = keep_tissues)
  )

## Function to make one plot per tissue
make_dotplot <- function(tissue_name) {
  ggplot(
    dot_df %>% filter(tissue == tissue_name),
    aes(x = gene, y = clonotype)
  ) +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    scale_size(range = c(1, 8), breaks = c(0, 25, 50, 75)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      title = tissue_name,
      x = "",
      y = "Clonotype",
      size  = "Percent expressing",
      color = "Average Expression"
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 7),
      legend.position = "right"
    )
}

## Create separate plots
plot_red   <- make_dotplot("RedSpleen")
plot_jln   <- make_dotplot("JLN")
plot_epi   <- make_dotplot("Epithelium")

plot_red
plot_jln
plot_epi

## Save to files
#ggsave("Dotplot_RedSpleen.png", plot_red, width = 5, height = 7, dpi = 300)
#ggsave("Dotplot_JLN.png",       plot_jln, width = 5, height = 7, dpi = 300)
#ggsave("Dotplot_Epithelium.png",plot_epi, width = 5, height = 7, dpi = 300)
## Dotplot: tissue = outline color, expression = fill gradient, fixed dot size

# Ensure consistent ordering
dot_df <- dot_df %>%
  mutate(
    clonotype = factor(clonotype, levels = rev(clones_use)),
    gene      = factor(gene, levels = valid_genes),
    tissue    = factor(tissue, levels = keep_tissues)
  )

# Position dodge so each tissue becomes a separate dot
pd <- position_dodge(width = 0.7)

p_all_tissues <- ggplot(dot_df, aes(x = gene, y = clonotype)) +
  geom_point(
    aes(
      fill  = avg_exp,   # expression gradient
      color = tissue     # tissue identity
    ),
    position = pd,
    size = 3,
    shape = 21,          # allows fill + outline
    stroke = 0.8
  ) +
  scale_fill_gradient2(
    low  = "blue",
    mid  = "white",
    high = "red",
    name = "Average expression"
  ) +
  scale_color_manual(
    values = c(
      RedSpleen  = "#D55E00",
      JLN        = "#0072B2",
      Epithelium = "#009E73"
    ),
    name = "Tissue"
  ) +
  labs(
    title = "Gene expression by clonotype (dots = tissues)",
    x = "",
    y = "Clonotype"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

p_all_tissues

#ggsave(p_all_tissues, file="DOT_p_all_tissues.png")


##### Bar chart
# Shared colorblind-friendly palette (Okabe–Ito)
cb_palette <- c(
  "Epithelium" = "#E69F00",  # orange
  "JLN"        = "#009E73",  # bluish-green
  "RedSpleen"  = "#0072B2"   # blue
)

# Attach clonotype key + use CT labels (clonotype_id)
dot_df <- dot_df %>%
  left_join(clonotype_key, by = "clonotype") %>%   
  filter(!is.na(clonotype_id)) %>%                 
  mutate(
    clonotype_id = factor(clonotype_id, levels = rev(clonotype_key$clonotype_id)),
    gene         = factor(gene, levels = valid_genes),
    tissue       = factor(tissue, levels = names(cb_palette))  # enforce same order
  )

pd <- position_dodge(width = 0.7)

# Numeric positions of clonotype rows (now CT labels)
clono_y <- seq_along(levels(dot_df$clonotype_id))
sep_lines <- clono_y[-length(clono_y)] + 0.5

p_all_tissues <- ggplot(
  dot_df,
  aes(y = clonotype_id, x = avg_exp, fill = tissue)  
) +
  geom_col(
    position = pd,
    width = 0.6,
    color = "black",
    linewidth = 0.2
  ) +
  facet_wrap(~ gene, nrow = 1) +
  #scale_fill_manual(values = cb_palette, breaks = names(cb_palette), name = "Tissue") +
  labs(
    title = "Gene Expression By Clonotype and Tissue",
    x = "Average Expression",
    y = "Clonotype"
  ) +
  theme_classic() +
  scale_fill_manual(values = cb_palette,
                    name = "Tissue",
                    breaks = c("Epithelium", "JLN", "RedSpleen"),
                    labels = c(
                      "Epithelium" = "Epithelium",
                      "JLN"        = "Non-Epithelial Derived JLN",
                      "RedSpleen"  = "Non-Epithelial Derived Spleen"
                    ) ) +
  theme(
    axis.text.y  = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_all_tissues
ggsave(p_all_tissues, file="clonotype_bar.png")

## Subset Bar Plot

# Helper to keep CT order
make_subset_bar <- function(df, ct_keep, plot_title) {
  df_sub <- df %>%
    filter(clonotype_id %in% ct_keep) %>%
    mutate(
      clonotype_id = factor(clonotype_id, levels = rev(ct_keep))  # controls order on y-axis
    )
  
  ggplot(df_sub, aes(y = clonotype_id, x = avg_exp, fill = tissue)) +
    geom_col(
      position = pd,
      width = 0.6,
      color = "black",
      linewidth = 0.2
    ) +
    facet_wrap(~ gene, nrow = 1) +
    scale_fill_manual(values = cb_palette, breaks = names(cb_palette), name = "Tissue") +
    labs(
      title = plot_title,
      x = "Average Expression",
      y = "Clonotype"
    ) +
    scale_fill_manual(values = cb_palette,
                      name = "Tissue",
                      breaks = c("Epithelium", "JLN", "RedSpleen"),
                      labels = c(
                        "Epithelium" = "Epithelium",
                        "JLN"        = "Non-Epithelial Derived JLN",
                        "RedSpleen"  = "Non-Epithelial Derived Spleen"
                      ) ) +
    theme_classic() +
    theme(
      axis.text.y  = element_text(size = 7),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

# Plot 1: CT3, CT7, CT8
p_ct3_ct7_ct8 <- make_subset_bar(
  dot_df,
  ct_keep = c("CT3", "CT7", "CT8"),
  plot_title = "Gene Expression By Clonotype and Tissue (CT3, CT7, CT8)"
)

# Plot 2: CT1, CT12
p_ct1_ct12 <- make_subset_bar(
  dot_df,
  ct_keep = c("CT1", "CT12"),
  plot_title = "Gene Expression By Clonotype and Tissue (CT1, CT12)"
)

p_ct3_ct7_ct8
p_ct1_ct12

# Save plots
ggsave("/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/VDJ_Analysis_Batch_Correcte/TCR_Figures/bar_CT3_CT7_CT8.png", p_ct3_ct7_ct8)
ggsave("/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/VDJ_Analysis_Batch_Correcte/TCR_Figures/bar_CT1_CT12.png",    p_ct1_ct12)


## UMAP Plotting code
Reductions(dat)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

clone_of_interest <- unique(df_sub$clonotype)
keep_tissues <- c("RedSpleen", "JLN", "Epithelium")

join_key <- "cell_id"
meta_cols <- colnames(dat@meta.data)

clono_col <- dplyr::case_when(
  "cdr3"   %in% meta_cols ~ "cdr3",
  "cdr3.y" %in% meta_cols ~ "cdr3.y",
  "cdr3.x" %in% meta_cols ~ "cdr3.x",
  TRUE ~ NA_character_
)

if (is.na(clono_col)) {
  stop("No clonotype column found. I looked for: cdr3, cdr3.y, cdr3.x. Available columns include: ",
       paste(head(meta_cols, 50), collapse = ", "))
}

# Build UMAP df
umap_df <- Embeddings(dat, reduction = "wnn.umap") %>%
  as.data.frame() %>%
  rownames_to_column(join_key)

colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")

# Join metadata to get tissue + raw clonotype
meta_df <- dat@meta.data %>%
  rownames_to_column(join_key) %>%
  dplyr::select(
    !!join_key,
    tissue,
    clonotype_raw = all_of(clono_col)
  ) %>%
  mutate(
    clonotype_raw = as.character(clonotype_raw),
    tissue        = as.character(tissue)
  )

umap_df <- umap_df %>%
  left_join(meta_df, by = join_key) %>%
  dplyr::rename(clonotype = clonotype_raw) %>%
  filter(is.na(tissue) | tissue %in% keep_tissues)

# attach CT labels from clonotype_key
umap_df <- umap_df %>%
  left_join(clonotype_key, by = "clonotype") %>%  # adds clonotype_id (CT1..)
  mutate(
    clonotype_id = as.character(clonotype_id),
    clonotype_plot = ifelse(!is.na(clonotype_id), clonotype_id, "Other"),
    tissue_plot    = ifelse(is.na(tissue), "Unknown", tissue)
  )

# Update clone_of_interest to be CT labels
clone_of_interest <- clonotype_key %>%
  filter(clonotype %in% unique(df_sub$clonotype)) %>%
  pull(clonotype_id) %>%
  unique()

# Factors for stable legends
umap_df$clonotype_plot <- factor(umap_df$clonotype_plot, levels = c("Other", clone_of_interest))
umap_df$tissue_plot    <- factor(umap_df$tissue_plot, levels = c(keep_tissues, "Unknown"))

bg_df <- umap_df %>% filter(clonotype_plot == "Other")
hi_df <- umap_df %>% filter(clonotype_plot != "Other")

# palette keyed to CT labels (not raw cdr3)
clone_pal <- setNames(scales::hue_pal()(length(clone_of_interest)), clone_of_interest)

shape_map <- c(
  RedSpleen  = 16,
  JLN        = 17,
  Epithelium = 15,
  Unknown    = 4
)

p_umap_clones <- ggplot() +
  geom_point(
    data = bg_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "grey84",
    size = 0.25
  ) +
  geom_point(
    data = hi_df,
    aes(x = UMAP_1, y = UMAP_2, color = clonotype_plot, shape = tissue_plot),
    size = 1.8,          # <-- slightly bigger highlighted clonotype dots (was 0.7)
    alpha = 0.9
  ) +
  scale_color_manual(values = clone_pal, name = "Clonotype") +
  scale_shape_manual(values = shape_map, name = "Tissue",
                     breaks = c("Epithelium", "JLN", "RedSpleen"),
                     labels = c(
                       "Epithelium" = "Epithelium",
                       "JLN"        = "Non-Epithelial Derived JLN",
                       "RedSpleen"  = "Non-Epithelial Derived Spleen")) +
  labs(
    title = paste0("Tissue-distributed Clonotypes in CD4 T-cells"),
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right") +
  scale_x_reverse() +
  scale_y_reverse() +
  coord_fixed()
p_umap_clones <- p_umap_clones +
  guides(
    shape = guide_legend(order = 1),   
    color = guide_legend(order = 2)    
  )
# Trim X-axis
p_umap_clones <- p_umap_clones +
  coord_cartesian(xlim = c(1.5, 6.5))

p_umap_clones
ggsave(p_umap_clones, file="/Users/eckco/Desktop/Kuhn_Lab/Sarah_Single_Cell_Analysis/Use_CSP_Data/VDJ_Analysis_Batch_Correcte/TCR_Figures/p_umap_clones_supplementary.png")
