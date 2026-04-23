library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(scales) #to look at color palettes #show_col(pal)
library(AnnotationDbi) #for converting between gene symbols and Entrez IDs 
library(org.Mm.eg.db) #mouse database for conversions between gene symbols and Entrez IDs
library(ChIPseeker) #for getting gene names and other gene annotations from the marker peaks
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
library(ggplot2)
library(ggrepel)
library(ggtext)
library(colorblindr)

CD4subsetproj <- loadArchRProject(path = "./save_CD4subset")

#marker peaks red spleen vs IEL----------------- 
markerpeaks <- getMarkerFeatures(
  ArchRProj = CD4subsetproj,
  groupBy = "TissueType", 
  useMatrix = "PeakMatrix", 
  bias = c("TSSEnrichment", "log10(nFrags)"), 
  testMethod = "wilcoxon", 
  useGroups = "IEL", 
  bgdGroups = "epithelial_spleen"
)

markerpeakslist <- getMarkers(
  seMarker = markerpeaks, 
  cutOff = "FDR <= 0.05 & abs(Log2FC) >= 2", 
  returnGR = TRUE
)

#motif enrichment----------------
#homer
CD4subsetproj <- addMotifAnnotations(
 ArchRProj = CD4subsetproj, 
motifSet = "homer", 
annoName = "motif")

motifEnrichmentUp_homer <- peakAnnoEnrichment(
  seMarker = markerpeaks, 
  ArchRProj = CD4subsetproj, 
  peakAnnotation = "motif", 
  cutOff = "FDR <= 0.05 & Log2FC >= 2"
)

motifEnrichmentDown_homer <- peakAnnoEnrichment(
  seMarker = markerpeaks, 
  ArchRProj = CD4subsetproj, 
  peakAnnotation = "motif", 
  cutOff = "FDR <= 0.05 & Log2FC <= -2"
)

#plots initial 
up_df <- data.frame(TF = rownames(motifEnrichmentUp_homer), 
                    mlog10Padj = assay(motifEnrichmentUp_homer)[,1])
up_df <- up_df[order(up_df$mlog10Padj, decreasing = TRUE),]
up_df$rank <- seq_len(nrow(up_df))

down_df <- data.frame(TF = rownames(motifEnrichmentDown_homer), 
                      mlog10Padj = assay(motifEnrichmentDown_homer)[,1])
down_df <- down_df[order(down_df$mlog10Padj, decreasing = TRUE),]
down_df$rank <- seq_len(nrow(down_df))

upplot <- ggplot(up_df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = up_df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  ggtitle("Up in cIELs vs red splenocytes")

downplot <- ggplot(down_df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = down_df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  ggtitle("Down in cIELs vs red splenocytes")

plotPDF(downplot, upplot, 
        name = "enriched_homer_motifs_redvsIEL",
        ArchRProj = CD4subsetproj, 
        width = 9, height = 9, 
        addDOC = FALSE)

#adjusted plots for figure 
TF_label <- c("MYB.HTH_168", "E2F.E2F_56", "Nkx6.1.Homeobox_186")
rownumbers <- which(up_df$TF%in% TF_label)
up_df$labels <- ""
up_df$labels[rownumbers] <- up_df$TF[rownumbers]

upplot <- ggplot(up_df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = up_df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, 
                                          label = labels), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  ggtitle("Up in CD4<sup>+</sup> cIELs vs. epithelial-derived splenocytes") + 
  theme(plot.title = element_markdown(face = "bold", hjust = 0.5), 
        legend.position = c(0.7, 0.2), 
        #legend.position = "bottom",
        legend.direction = "horizontal", 
        axis.title = element_text(face = "bold"))

ggsave(filename = "motif_enrichment_UP.png", plot = upplot, 
       width = 5.5, height = 4, units = "in")

#choose labels to plot 
TF_label <- c("c.Jun.CRE.bZIP_142", "E2A.bHLH..near_PU.1_57", "Atf2.bZIP_11", "GATA3.Zf..DR4_103", 
              "Atf7.bZIP_14", "JunD.bZIP_143")
rownumbers <- which(down_df$TF%in% TF_label)
down_df$labels <- ""
down_df$labels[rownumbers] <- down_df$TF[rownumbers]

downplot <- ggplot(down_df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2) +
  ggrepel::geom_label_repel(
    data = down_df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = labels), 
    size = 3,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  ggtitle("Down in CD4<sup>+</sup> cIELs vs. epithelial-derived splenocytes") +
  theme(plot.title = element_markdown(face = "bold", hjust = 0.5), 
        #legend.position = "bottom", 
        legend.position = c(0.7, 0.2), 
        legend.direction = "horizontal", 
        axis.title = element_text(face = "bold"))

ggsave(filename = "motif_enrichment_DOWN.png", plot = downplot, 
       width = 5.5, height = 4, units = "in")



