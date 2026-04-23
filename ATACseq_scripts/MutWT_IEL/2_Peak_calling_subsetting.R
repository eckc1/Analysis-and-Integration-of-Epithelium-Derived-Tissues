library(ArchR)
library(pheatmap)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)

#load ArchR project 
mutwtIEL_proj2 <- loadArchRProject(path = "./save_colonProj")

#Add iterative LSI--------
mutwtIEL_proj2 <- addIterativeLSI(
  ArchRProj = mutwtIEL_proj2, 
  useMatrix = "TileMatrix",
  name = "LSI", 
  iterations = 3,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.8, 1), 
    samplecells = 1000,  
    n.start = 10
  ), 
  firstSelection = "Top", 
  depthCol = "nFrags", 
  varFeatures = 25000,
  dimsToUse = 1:30, 
  seed = 1, 
  force = TRUE
)

#Add UMAP clusters-----------
mutwtIEL_proj2 <- addClusters(
  input = mutwtIEL_proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster02", 
  resolution = 0.2
)

mutwtIEL_proj2 <- addClusters(
  input = mutwtIEL_proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster04", 
  resolution = 0.4
)

mutwtIEL_proj2 <- addClusters(
  input = mutwtIEL_proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster08", 
  resolution = 0.8
)

mutwtIEL_proj2 <- addClusters(
  input = mutwtIEL_proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster1", 
  resolution = 1
)

#Add UMAP--------
mutwtIEL_proj2 <- addUMAP(
  ArchRProj = mutwtIEL_proj2, 
  reducedDims = "LSI", 
  name = "UMAP", 
  nNeighbors = 250, 
  minDist = 0.4, 
  metric = "cosine", 
  force = TRUE)

#plot UMAP-------------------
umap1 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "Sample", size = 0.1)
umap2 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "SacDates", size = 0.1)

pal = c("05_10_24" = "#DC494C", "12_15_23" = "#371377")
umap3 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "SeqDates", 
                       pal = pal, 
                       size = 0.1)
pal = c("WT" = "#4c34b9", "MHCII" = "#c0692c")
umap4 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "Genotype", 
                       pal = pal, 
                       size = 0.1)
umap5 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "Sex", 
                       size = 0.1)
umap6 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "Age", 
                       size = 0.1)
umap7 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "cluster02", 
                       size = 0.1)
umap8 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "cluster04", 
                       size = 0.1)
umap9 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                       embedding = "UMAP", 
                       colorBy = "cellColData", 
                       name = "cluster08", 
                       size = 0.1)
umap10 <- plotEmbedding(ArchRProj = mutwtIEL_proj2, 
                        embedding = "UMAP", 
                        colorBy = "cellColData", 
                        name = "cluster1", 
                        size = 0.1)
plotPDF(umap1, umap2, umap3, umap4, umap5, umap6, umap7, umap8, umap9, umap10,
        name = "initialUMAPs", addDOC = FALSE, ArchRProj = mutwtIEL_proj2,
        width = 7, height = 7)
rm(umap1, umap2, umap3, umap4, umap5, umap6, umap7, umap8, umap9, umap10, pal)

#Add gene scores and plot T cell subset marker genes-----
mutwtIEL_proj2 <- addImputeWeights(ArchRProj = mutwtIEL_proj2, reducedDims = "LSI")


markergenes <- c("Cd4", "Cd8a")
umapMarkersWeighted <- plotEmbedding(
  ArchRProj = mutwtIEL_proj2, 
  colorBy = "GeneScoreMatrix", 
  name = markergenes, 
  embedding = "UMAP", 
  quantCut = NULL, 
  imputeWeights = getImputeWeights(mutwtIEL_proj2)
)

plotPDF(umapMarkersWeighted$Cd4, 
        umapMarkersWeighted$Cd8a, 
        name = "mutGenescoreMarkersUMAP.pdf", 
        height = 7, width = 7, 
        ArchRProj = mutwtIEL_proj2,
        addDOC = FALSE
)
rm(umapMarkersWeighted, markergenes)

#Call peaks for CD4+ T cells-----------
#subset CD4 T cells into new ArchR project
cells <- BiocGenerics::which(mutwtIEL_proj2$cluster04 %in% c("C1", "C2"))
cellnames <- mutwtIEL_proj2$cellNames[cells]
CD4subsetproj <- subsetArchRProject(ArchRProj = mutwtIEL_proj2, cells = cellnames, 
                                    dropCells = TRUE, 
                                    outputDirectory = "./CD4subset", 
                                    force = TRUE)
rm(cellnames, cells)

#Prepare for pseudobulk replicates
IELmut1 <- BiocGenerics::which(CD4subsetproj$Sample %in% "IELmut1")
IELmut1_1 <- CD4subsetproj$Genotype[IELmut1]

IELmut2 <- BiocGenerics::which(CD4subsetproj$Sample %in% "IELmut2")
IELmut2_1 <- CD4subsetproj$Genotype[IELmut2]

IEL1 <- BiocGenerics::which(CD4subsetproj$Sample %in% "IEL1")
IEL1_1 <- CD4subsetproj$Genotype[IEL1]

IEL2 <- BiocGenerics::which(CD4subsetproj$Sample %in% "IEL2")
IEL2_1 <- CD4subsetproj$Genotype[IEL2]

IEL3 <- BiocGenerics::which(CD4subsetproj$Sample %in% "IEL3")
IEL3_1 <- CD4subsetproj$Genotype[IEL3]

table(IELmut1_1)
table(IELmut2_1)
table(IEL1_1)
table(IEL2_1)
table(IEL3_1)

rm(IELmut1, IELmut1_1, IELmut2, IELmut2_1, IEL1, IEL1_1, IEL2, IEL2_1, IEL3, IEL3_1)

#call peaks
path2MACS2 <- findMacs2()
CD4subsetproj <- addGroupCoverages(ArchRProj = CD4subsetproj, groupBy = "Genotype", 
                                   minCells = 200,
                                   maxCells = 500, 
                                   minReplicates = 3, 
                                   maxReplicates = 4, 
                                   sampleRatio = 0.8,
                                   logFile = createLogFile("addGroupCoverages"))
CD4subsetproj <- addReproduciblePeakSet(ArchRProj = CD4subsetproj, 
                                        groupBy = "Genotype",
                                        pathToMacs2 = path2MACS2)
getPeakSet(CD4subsetproj)
CD4subsetproj <- addPeakMatrix(CD4subsetproj)

#save ArchR projects------------
mutwtIEL_proj2 <- saveArchRProject(ArchRProj = mutwtIEL_proj2, 
                                  outputDirectory = "save_mutwtIEL", 
                                  load = TRUE)

CD4subsetproj <- saveArchRProject(ArchRProj = CD4subsetproj, 
                               outputDirectory = "save_CD4subset", 
                               load = TRUE)













