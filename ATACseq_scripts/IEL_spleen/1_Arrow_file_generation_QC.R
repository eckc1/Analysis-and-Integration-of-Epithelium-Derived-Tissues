library(ArchR)
library(pheatmap)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRGenome("mm10")

inputFiles_spleen <- c("~/Downloads/ATAC_analysis/resequence/frags/red1_fragments.tsv.gz", 
                "~/Downloads/ATAC_analysis/resequence/frags/red2_fragments.tsv.gz", 
                "~/Downloads/ATAC_analysis/resequence/frags/green1_fragments.tsv.gz", 
                "~/Downloads/ATAC_analysis/resequence/frags/green2_fragments.tsv.gz")


inputFiles_IEL_1 <- c("~/Downloads/ATAC_analysis/resequence/frags/IEL1_fragments.tsv.gz", 
                      "~/Downloads/ATAC_analysis/resequence/frags/IEL2_fragments.tsv.gz")

inputFiles_IEL_2 <- c("~/Downloads/ATAC_analysis/all_colon_ATAC/05_24_frags/IEL3_fragments.tsv.gz")

names(inputFiles_spleen) <- c("red1", "red2", "green1", "green2")
names(inputFiles_IEL_1) <- c("IEL1", "IEL2")
names(inputFiles_IEL_2) <- c("IEL3")

#make arrow files---------------
spleen_arrow_files <- createArrowFiles(
  inputFiles = inputFiles_spleen, 
  sampleNames = names(inputFiles_spleen), 
  minTSS = 20, 
  minFrags = 10000, 
  addTileMat = TRUE, 
  addGeneScoreMat = TRUE
)
IEL_arrowFiles_1 <- createArrowFiles(
  inputFiles = inputFiles_IEL_1, 
  sampleNames = names(inputFiles_IEL_1), 
  minTSS = 20, 
  minFrags = 5600, 
  addTileMat = TRUE, 
  addGeneScoreMat = TRUE
)
IEL_arrowFiles_2 <- createArrowFiles(
  inputFiles = inputFiles_IEL_2, 
  sampleNames = names(inputFiles_IEL_2), 
  minTSS = 22, 
  minFrags = 10000, 
  maxFrags = 39800, 
  addTileMat = TRUE, 
  addGeneScoreMat = TRUE
)

#add doublet scores 
doublet_scores <- addDoubletScores(
  input = spleen_arrow_files, 
  k=10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)
doublet_scores_IEL_1 <- addDoubletScores(
  input = IEL_arrowFiles_1, 
  k=10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)
doublet_scores_IEL_2 <- addDoubletScores(
  input = IEL_arrowFiles_2, 
  k=10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

#Make ArchR project-----------
arrowfiles <- c("~/Downloads/ATAC_analysis/resequence/colon_initial/IEL1.arrow", 
                "~/Downloads/ATAC_analysis/resequence/colon_initial/IEL2.arrow", 
                "~/Downloads/ATAC_analysis/all_colon_ATAC/IEL3.arrow",
                "~/Downloads/ATAC_analysis/resequence/spleen/green1.arrow", 
                "~/Downloads/ATAC_analysis/resequence/spleen/green2.arrow", 
                "~/Downloads/ATAC_analysis/resequence/spleen/red1.arrow", 
                "~/Downloads/ATAC_analysis/resequence/spleen/red2.arrow"
)
proj1 <- ArchRProject(
  ArrowFiles = arrowfiles, 
  outputDirectory = "upstreamoutput", 
  copyArrows = TRUE
)

#add metadata-----------------
#date of sac 
sacdates <- gsub("IEL1", "11_17_23",
                 gsub("IEL2", "11_21_23", 
                      gsub("IEL3", "04_16_24", 
                           gsub("green1", "11_17_23", 
                                gsub("green2", "11_21_23", 
                                     gsub("red1", "11_17_23", 
                                          gsub("red2", "11_21_23", proj1$Sample)))))))
proj1$SacDates <- sacdates
#tissue type 
tissuetype <- gsub("IEL1", "IEL",
                   gsub("IEL2", "IEL",
                        gsub("IEL3", "IEL",
                             gsub("green1", "spleen",
                                  gsub("green2", "spleen", 
                                       gsub("red1", "epithelial_spleen", 
                                            gsub("red2", "epithelial_spleen", proj1$Sample)))))))
proj1$TissueType <- tissuetype

#QC plots before filtering 
violin1 <- plotGroups(ArchRProj = proj1, 
                      groupBy = "Sample", 
                      colorBy = "cellColData", 
                      plotAs = "violin", 
                      name = "TSSEnrichment",
                      alpha = 0.4, 
                      addBoxPlot = TRUE)

violin2 <- plotGroups(ArchRProj = proj1, 
                      groupBy = "SacDates", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)

violin3 <- plotGroups(ArchRProj = proj1, 
                      groupBy = "TissueType", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)


plotPDF(violin1, violin2, violin3, name = "TSS_violin_plots_prefilter.pdf", ArchRProj = proj1, addDOC = FALSE, height = 7, width = 7)
rm(violin1, violin2, violin3)


fragsize <- plotFragmentSizes(ArchRProj = proj1)
TSS <- plotTSSEnrichment(ArchRProj = proj1)
plotPDF(fragsize, TSS, name = "fragsize_TSSenrichment_prefilter.pdf", height = 7, width = 7, addDOC = FALSE, ArchRProj = proj1)
rm(fragsize, TSS)

#filter doublets--------------
proj2 <- filterDoublets(proj1)

#post-filtering QC plots---------------
violin1 <- plotGroups(ArchRProj = proj2, 
                      groupBy = "Sample", 
                      colorBy = "cellColData", 
                      plotAs = "violin", 
                      name = "TSSEnrichment",
                      alpha = 0.4, 
                      addBoxPlot = TRUE)

violin2 <- plotGroups(ArchRProj = proj2, 
                      groupBy = "SacDates", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)

violin3 <- plotGroups(ArchRProj = proj2, 
                      groupBy = "TissueType", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)


plotPDF(violin1, violin2, violin3, name = "TSS_violin_plots_prefilter.pdf", ArchRProj = proj2, addDOC = FALSE, height = 7, width = 7)
rm(violin1, violin2, violin3)


fragsize <- plotFragmentSizes(ArchRProj = proj2)
TSS <- plotTSSEnrichment(ArchRProj = proj2)
plotPDF(fragsize, TSS, name = "fragsize_TSSenrichment_prefilter.pdf", height = 7, width = 7, addDOC = FALSE, ArchRProj = proj2)
rm(fragsize, TSS)

#save ArchR project------------
proj2 <- saveArchRProject(ArchRProj = proj2, 
                               outputDirectory = "save_proj2", 
                               load = TRUE)

























