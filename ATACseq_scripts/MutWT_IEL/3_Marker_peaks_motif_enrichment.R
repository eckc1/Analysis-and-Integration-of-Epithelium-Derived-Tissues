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

#Identify marker peaks----------------
#get and list marker peaks 
markerpeaks <- getMarkerFeatures(
  ArchRProj = CD4subsetproj,
  groupBy = "Genotype", 
  useMatrix = "PeakMatrix", 
  bias = c("TSSEnrichment", "log10(nFrags)"), 
  testMethod = "wilcoxon"
)

markerpeakslist <- getMarkers(
  seMarker = markerpeaks, 
  cutOff = "FDR <= 0.05 & abs(Log2FC) >= 2", 
  returnGR = TRUE
)

#annotate marker peaks 
MHCIIannotatedpeaks <- annotatePeak(
  markerpeakslist$MHCII, 
  tssRegion = c(-3000, 3000), 
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
  level = "transcript", 
  assignGenomicAnnotation = TRUE, 
  annoDb = "org.Mm.eg.db",
  addFlankGeneInfo = TRUE, 
  flankDistance = 5000, 
  sameStrand = FALSE,
  columns = c("SYMBOL", "GENENAME")
)

WTannotatedpeaks <- annotatePeak(
  markerpeakslist$WT, 
  tssRegion = c(-3000, 3000), 
  TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
  level = "transcript", 
  assignGenomicAnnotation = TRUE, 
  annoDb = "org.Mm.eg.db",
  addFlankGeneInfo = TRUE, 
  flankDistance = 5000, 
  sameStrand = FALSE,
  columns = c("SYMBOL", "GENENAME")
)
write.table(MHCIIannotatedpeaks, 
            file = "./MHCIImarkerpeaks_anno.csv",
            sep = ",", dec = ".")


write.table(WTannotatedpeaks, 
            file = "./WTmarkerpeaks_anno.csv",
            sep = ",", dec = ".")



#Motif enrichment in marker peaks---------
#motif enrichment 
#repeat, making new motif annotation, with other databases 
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

#plot enrichment- initial 
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
  ggtitle("Up in Mut CD4 IELs vs WT IELs full sample")

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
  ggtitle("Down in Mut CD4 IELs vs WT IELs full sample")

plotPDF(downplot, upplot, 
        name = "enriched_motifs_fullsample",
        ArchRProj = CD4subsetproj, 
        width = 9, height = 9,
        addDOC = FALSE)

#adjusted plot for figure 
downplot <- ggplot(down_df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 3) +
  ggrepel::geom_label_repel(
    data = down_df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 4,
    nudge_x = 2,
    color = "black", 
    max.overlaps = 1
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
  ggtitle("Less accesible TF motifs in CD4<sup>+</sup> MHCII<sup>ΔIEC</sup> cIELs") + 
  theme(plot.title = element_markdown(face = "bold", hjust = 0, size = 14), 
        #legend.position = c(0.7, 0.15), 
        legend.position = "right", 
        #legend.direction = "horizontal", 
        legend.direction = "vertical", 
        panel.border = element_blank(), 
        axis.line.x = element_line(color = "black"), 
        axis.line.y = element_line(color = "black"), 
        axis.title = element_text(face = "bold", size = 12))
ggsave("TF_enrichment.tiff", plot = downplot, height = 4, width = 6)

#ChromVar deviations enrichment------------
CD4subsetproj <- addBgdPeaks(CD4subsetproj)
CD4subsetproj <- addDeviationsMatrix(
  ArchRProj = CD4subsetproj, 
  peakAnnotation = "motif"
)

#Violin plot of NFAT enrichment score--------
#make default plot 
nfat_violin_default <- plotGroups(ArchRProj = CD4subsetproj, groupBy = "Genotype", colorBy = "motifMatrix", 
                                  name = "z:NFAT.RHD_178", plotAs = "violin", pal = pal)
#extract data
nfat_violin_df <- nfat_violin_default@data

#rename columns 
colnames(nfat_violin_df) <- c("Genotype", "Zscore")

#make plot 
pal <- c("WT" = "#107F80",
         "MHCII" = "#40007F")

nfat_violin_scratch <- ggplot(data = nfat_violin_df, 
                              aes(x = Genotype, y = Zscore, fill = Genotype, colour = Genotype)) + 
  geom_violin(alpha = 0.5) + geom_boxplot(alpha = 0) +
  scale_x_discrete(limits = rev(levels(nfat_violin_df$Genotype)), 
                   labels = c("WT" = "MHCII<sup>WT</sup>", "MHCII" = "MHCII<sup>ΔIEC</sup>")) +
  theme_classic() +
  scale_fill_manual(values = pal) + 
  scale_color_manual(values = pal) + 
  ylab("Z-score") + 
  labs(title = "NFAT binding motif accessibility in CD4<sup>+</sup> cIELs") + 
  theme(legend.position = "none",
        plot.title = element_markdown(hjust = 0.5, face = "bold", size = 12), 
        axis.text.x = element_markdown(size = 10, face = "bold"), 
        axis.title = element_text(face = "bold"), 
        axis.title.x = element_blank())
ggsave(filename = "NFAT_violin.tiff", plot = nfat_violin_scratch, 
       height = 3, width = 4.2)
