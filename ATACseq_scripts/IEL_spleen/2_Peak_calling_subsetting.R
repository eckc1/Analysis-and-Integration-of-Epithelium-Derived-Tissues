library(ArchR)
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)

#Load ArchR project
proj2 <- loadArchRProject(path = "./save_proj2")

#Add iterative LSI--------------
proj2 <- addIterativeLSI(
  ArchRProj = proj2, 
  useMatrix = "TileMatrix",
  name = "LSI", 
  iterations = 4,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.6, 0.8, 1), 
    samplecells = 30000,
    n.start = 10
  ), 
  firstSelection = "Top", 
  depthCol = "nFrags", 
  varFeatures = 25000,
  dimsToUse = 1:30, 
  seed = 1, 
  force = TRUE
)

#add UMAP clusters-------------
proj2 <- addClusters(
  input = proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster02", 
  resolution = 0.2)

proj2 <- addClusters(
  input = proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster04", 
  resolution = 0.4)

proj2 <- addClusters(
  input = proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster06", 
  resolution = 0.6)

proj2 <- addClusters(
  input = proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster08", 
  resolution = 0.8)

proj2 <- addClusters(
  input = proj2, 
  reducedDims = "LSI", 
  method = "Seurat", 
  name = "cluster1", 
  resolution = 1)

#add UMAP-----------------
proj2 <- addUMAP(
  ArchRProj = proj2, 
  reducedDims = "LSI", 
  name = "UMAP", 
  nNeighbors = 15, 
  minDist = 0.4,
  metric = "cosine", 
  force = TRUE)

#UMAP plots--------------------
umaps <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "Sample")
umapc <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "TissueType")
umap1 <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "cluster02")
umap2 <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "cluster04")
umap3 <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "cluster06")
umap4 <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "cluster08")
umap5 <- plotEmbedding(ArchRProj = proj2, embedding = "UMAP", colorBy = "cellColData", name = "cluster1")

plotPDF(umaps, umapc, umap1, umap2, umap3, umap4, umap5,
        name = "UMAPs.pdf", addDOC = FALSE, 
        ArchRProj = proj2, 
        width = 7, height = 7)
rm(umaps, umapc, umap1, umap2, umap3, umap4, umap5)

#Impute weights--------------
proj2 <- addImputeWeights(ArchRProj = proj2, reducedDims = "LSI")

#Immune cell subsets marker genes----------------
markergenes <- c("Cd3e", "Cd4", "Cd8a", "Cd8b1")

umapMarkersWeighted <- plotEmbedding(
  ArchRProj = proj2, 
  colorBy = "GeneScoreMatrix", 
  name = markergenes, 
  embedding = "UMAP", 
  quantCut = NULL, 
  imputeWeights = getImputeWeights(proj2)
)

plotPDF(
  umapMarkersWeighted$Cd3e,
  umapMarkersWeighted$Cd4, 
  umapMarkersWeighted$Cd8a, 
  umapMarkersWeighted$Cd8b1,
  name = "genescoreMarkersWeighted.pdf", height = 7, width = 7, 
  ArchRProj = proj2, addDOC = FALSE
)

rm(umapMarkersWeighted, markergenes)

#Call peaks for CD4+ T cells-----------------
#subset CD4 T cells
cells <- BiocGenerics::which(proj2$cluster02 %in% "C8")
cellnames <- proj2$cellNames[cells]
CD4subsetproj <- subsetArchRProject(ArchRProj = proj2, cells = cellnames, dropCells = TRUE, force = TRUE,
                                    outputDirectory = "./CD4_subset")

#prepare for pseudobulk replicates- identify min and max # of cells per group 
tissuecheck <- CD4subsetproj@cellColData$TissueType
tissuecheck
rm(tissuecheck)

#call peaks
path2MACS2 <- findMacs2()
CD4subsetproj <- addGroupCoverages(ArchRProj = CD4subsetproj, groupBy = "TissueType", 
                                   minCells = 200, 
                                   maxCells = 600,  
                                   minReplicates = 2, 
                                   maxReplicates = 3, 
                                   sampleRatio = 0.8,
                                   logFile = createLogFile("addGroupCoverages"))

CD4subsetproj <- addReproduciblePeakSet(ArchRProj = CD4subsetproj, 
                                        groupBy = "TissueType",
                                        pathToMacs2 = path2MACS2)
getPeakSet(CD4subsetproj)
CD4subsetproj2 <- addPeakMatrix(CD4subsetproj)

#save projects 
proj2 <- saveArchRProject(ArchRProj = proj2, 
                                   outputDirectory = "save_proj2", 
                                   load = TRUE)

CD4subsetproj <- saveArchRProject(ArchRProj = CD4subsetproj, 
                                  outputDirectory = "save_CD4subset", 
                                  load = TRUE)











