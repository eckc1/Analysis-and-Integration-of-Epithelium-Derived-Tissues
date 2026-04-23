library(ArchR)
library(pheatmap)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
addArchRGenome("mm10")

inputFiles_IEL_1 <- c("~/Downloads/ATAC_analysis/resequence/frags/IEL1_fragments.tsv.gz", 
                    "~/Downloads/ATAC_analysis/resequence/frags/IEL2_fragments.tsv.gz", 
                    "~/Downloads/ATAC_analysis/resequence/frags/IELmut1_fragments.tsv.gz")


inputFiles_IEL_2 <- c("~/Downloads/ATAC_analysis/all_colon_ATAC/05_24_frags/IELmut2_fragments.tsv.gz", 
                      "~/Downloads/ATAC_analysis/all_colon_ATAC/05_24_frags/IEL3_fragments.tsv.gz")

names(inputFiles_IEL_1) <- c("IEL1", "IEL2", "IELmut1")
names(inputFiles_IEL_2) <- c("IELmut2", "IEL3")

#create Arrow files-----------------------
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

#Add doublet scores-------------------------------
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
                "~/Downloads/ATAC_analysis/resequence/colon_initial/IELmut1.arrow", 
                "~/Downloads/ATAC_analysis/all_colon_ATAC/IELmut2.arrow", 
                "~/Downloads/ATAC_analysis/all_colon_ATAC/IEL3.arrow")

colonproj <- ArchRProject(
  ArrowFiles = arrowfiles, 
  outputDirectory = "colonprojOutput", 
  copyArrows = TRUE
)

#Add metadata to ArchR project----------
#date of takedown
sacdates <- gsub("IEL1", "11_17_23",
                 gsub("IEL2", "11_21_23", 
                      gsub("IEL3", "04_16_24", 
                           gsub("IELmut1", "12_04_23",
                                gsub("IELmut2", "04_16_24", mutwtIEL_proj1$Sample)))))
mutwtIEL_proj1$SacDates <- sacdates

#date of sequencing
seqdates <- gsub("IEL1", "12_15_23",
                 gsub("IEL2", "12_15_23",
                      gsub("IEL3", "05_10_24",
                           gsub("IELmut1", "12_15_23",
                                gsub("IELmut2", "05_10_24", mutwtIEL_proj1$Sample)))))
mutwtIEL_proj1$SeqDates <- seqdates

#genotype
genotype <- gsub("IEL1", "WT",
                 gsub("IEL2", "WT",
                      gsub("IEL3", "WT",
                           gsub("IELmut1", "MHCII",
                                gsub("IELmut2", "MHCII", mutwtIEL_proj1$Sample)))))
mutwtIEL_proj1$Genotype <- genotype

#sex
sex <- gsub("IEL1", "F",
            gsub("IEL2", "M",
                 gsub("IEL3", "F",
                      gsub("IELmut1", "F",
                           gsub("IELmut2", "M", mutwtIEL_proj3$Sample)))))
mutwtIEL_proj3$Sex <- sex

#age
age <- gsub("IEL1", "11",
            gsub("IEL2", "29",
                 gsub("IEL3", "12",
                      gsub("IELmut1", "13",
                           gsub("IELmut2", "11", mutwtIEL_proj3$Sample)))))
mutwtIEL_proj3$Age <- age

#QC plots before doublet filtering--------------
violin1 <- plotGroups(ArchRProj = colonproj, 
                      groupBy = "Sample", 
                      colorBy = "cellColData", 
                      plotAs = "violin", 
                      name = "TSSEnrichment",
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
violin2 <- plotGroups(ArchRProj = colonproj, 
                      groupBy = "SacDates", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
violin3 <- plotGroups(ArchRProj = colonproj, 
                      groupBy = "SeqDates", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
violin4 <- plotGroups(ArchRProj = colonproj, 
                      groupBy = "Genotype", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
plotPDF(violin1, violin2, violin3, violin4, name = "TSS_violin_plots_prefilter.pdf", ArchRProj = colonproj, addDOC = FALSE, height = 7, width = 7)
rm(violin1, violin2, violin3, violin4)


fragsize <- plotFragmentSizes(ArchRProj = colonproj)
TSS <- plotTSSEnrichment(ArchRProj = colonproj)
plotPDF(fragsize, TSS, name = "fragsize_TSSenrichment_prefilter.pdf", height = 7, width = 7, addDOC = FALSE, ArchRProj = colonproj)
rm(fragsize, TSS)

#Filter doublets----------------------
colonproj2 <- filterDoublets(colonproj)

#QC plots after doublet filtering
violin1 <- plotGroups(ArchRProj = colonproj2, 
                      groupBy = "Sample", 
                      colorBy = "cellColData", 
                      plotAs = "violin", 
                      name = "TSSEnrichment",
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
violin2 <- plotGroups(ArchRProj = colonproj2, 
                      groupBy = "SacDates", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
violin3 <- plotGroups(ArchRProj = colonproj2, 
                      groupBy = "SeqDates", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
violin4 <- plotGroups(ArchRProj = colonproj2, 
                      groupBy = "Genotype", 
                      colorBy = "cellColData",
                      name = "TSSEnrichment",
                      plotAs = "violin", 
                      alpha = 0.4, 
                      addBoxPlot = TRUE)
plotPDF(violin1, violin2, violin3, violin4, name = "TSS_violin_plots_postfilter.pdf", ArchRProj = colonproj2, addDOC = FALSE, height = 7, width = 7)
rm(violin1, violin2, violin3, violin4)

fragsize <- plotFragmentSizes(ArchRProj = colonproj2)
TSS <- plotTSSEnrichment(ArchRProj = colonproj2)
plotPDF(fragsize, TSS, name = "fragsize_TSSenrichment_postfilter.pdf", height = 7, width = 7, addDOC = FALSE, ArchRProj = colonproj2)
rm(fragsize, TSS)


#save ArchR project------------
colonproj2 <- saveArchRProject(ArchRProj = colonproj2, 
                               outputDirectory = "save_colonProj", 
                               load = TRUE)


