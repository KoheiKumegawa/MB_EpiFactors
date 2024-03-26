#----------------------------------------------------------------------------
# scMultiomics-Tumor.R
#----------------------------------------------------------------------------
# ArchR-based analysis
library(ArchR)
addArchRThreads(threads = 8) 
addArchRGenome("mm10")

#-------- 01 make ArchR project --------#
#sample assignment
sampleName <- list.files("data/", pattern = "_fragments.tsv.gz", full.names = T)
names(sampleName) <- c("RS03060_T")
outFile <- as.character(sampleName)

#make arrow files
ArrowFiles = character(length(sampleName))
ArrowFiles <- createArrowFiles(inputFiles = outFile,
                               sampleNames = names(sampleName), 
                               minTSS = 4, minFrags = 1000, 
                               addTileMat = TRUE, addGeneScoreMat = TRUE)

#doublet detection
doubScores <- addDoubletScores(ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

#make ArchR project and filter doublets!
pre_arc <- ArchRProject(ArrowFiles, outputDirectory = "output", copyArrows = F)
pre_arc <- filterDoublets(pre_arc)

#quality filter
arc <- pre_arc[which(pre_arc$TSSEnrichment > 6 & pre_arc$nFrags > 2500)]
saveRDS(arc, "rds/01_arc.rds")

#-------- 02 combine scRNA and clustering --------#
#read RNA count matrix
inputRNA <- "data/RS03060_filtered_feature_bc_matrix.h5"
seRNA <- import10xFeatureMatrix(input = inputRNA, names = "RS03060_T")

#add gene expression to ArchR project
arc <- addGeneExpressionMatrix(input = arc, seRNA = seRNA)
arc <- arc[!is.na(arc$Gex_nUMI)]
arc <- arc[which(arc$Gex_nGenes > 400 & arc$Gex_nGenes < 9000)]

#LSI-ATAC
arc <- addIterativeLSI(
  ArchRProj = arc, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

#LSI-RNA
arc <- addIterativeLSI(
  ArchRProj = arc, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

#Combined Dims
arc <- addCombinedDims(arc, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

#UMAPs
arc <- addUMAP(arc, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
arc <- addUMAP(arc, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
arc <- addUMAP(arc, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

#Add Clusters
arc <- addClusters(arc, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)

#Plot Embedding
p1 <- plotEmbedding(arc, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(arc, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(arc, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p1, p2, p3, name = "02_UMAP-scATAC-scRNA-Combined.pdf", ArchRProj = arc, addDOC = FALSE)

saveRDS(arc, "rds/02_arc.rds")

#-------- 03 marker vis --------#
markerGenes <- c("Col3a1", #Vascular fibroblast
                 "C1qb", "Aif1", #Microglia
                 "Sox10", "Olig2", #Oligodendrocytes
                 "Aqp4", #Astrocytes
                 "Pax2", "Ascl1", #interneuron progenitors/stem
                 "Gli1", "Gli2", "Ccnd2", "Barhl1", "Atoh1", #GNP or GNP-like tumor
                 "Cntn2", "Grin2b", #differentiated granule cells/tumor cells
                 "Mki67",  "Mcm2", "Ccnd1", "Ccnb1")

arc <- addImputeWeights(arc, reducedDims = "LSI_Combined")
p4 <- plotEmbedding(arc, 
                    colorBy = "GeneScoreMatrix", 
                    name = markerGenes, 
                    embedding = "UMAP_Combined", 
                    imputeWeights = getImputeWeights(arc), 
                    plotAs = "points", size = 1.5, rastr = T)
p5 <- plotEmbedding(arc, 
                    colorBy = "GeneExpressionMatrix", 
                    name = markerGenes, 
                    embedding = "UMAP_Combined", log2Norm = T, 
                    imputeWeights = getImputeWeights(arc),
                    plotAs = "points", size = 1.5, rastr = T)
pdf("output/Plots/03_UMAP_linegeMarkers.pdf", width = 4, height = 4)
p4
p5
dev.off()

#-------- 04 peak call and p2g linkage--------#
#peak call and add peak matrix to arrow files
pathToMacs2 <- findMacs2()
arc <- addGroupCoverages(arc, groupBy = "Clusters")
arc <- addReproduciblePeakSet(arc, groupBy = "Clusters", pathToMacs2 = pathToMacs2, method = "q", cutOff = 0.1)
arc <- addPeakMatrix(arc)
arc <- addPeak2GeneLinks(arc,
                         reducedDims = "LSI_Combined", 
                         useMatrix = "GeneExpressionMatrix", 
                         maxDist = 500000)
p2g_all <- getPeak2GeneLinks(arc, corCutOff = 0, FDRCutOff = 1, resolution = 500, returnLoops = F) #492420 peak-gene associations
p2g_500k <- getPeak2GeneLinks(arc, corCutOff = 0.2, FDRCutOff = 0.01, resolution = 500, returnLoops = F) #28423 peak-gene associations

#df <- data.frame(Correlation = p2g_all$Correlation, FDR = -log10(p2g_all$FDR))
df <- data.frame(Rank = c(1:492420), Correlation = p2g_all$Correlation)
p5.2 <- ggplot(df, aes(x=Correlation)) + 
  geom_density(color = "lightgray", fill="red",alpha=.2) + theme_ArchR() +
  geom_vline(xintercept = 0.2, lty = "dotted") + xlim(c(0,1))
pdf("output/Plots/04_density_p2g_all.pdf", width = 5, height = 5)
p5.2
dev.off()

#Neurod1
idx1 <- which(metadata(p2g_500k)$geneSet$name == "Neurod1") 
idx2 <- p2g_500k$idxATAC[which(p2g_500k$idxRNA == idx1)]
neurod1_rd <- metadata(p2g_500k)$peakSet[idx2]

d <- data.frame(seqnames = seqnames(neurod1_rd), start = start(neurod1_rd)-1, end = end(neurod1_rd))
write.table(d, "output/output_bed/04_Neurod1_RD.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#Vis
neurod1_gr <- metadata(p2g_500k)$geneSet[idx1]
neurod1_gr <- resize(neurod1_gr, width = 1, fix = "start")

p2g_500k_gr <- getPeak2GeneLinks(arc, corCutOff = 0.2, FDRCutOff = 0.01, resolution = 1, returnLoops = T)
p2g_500k_gr <- p2g_500k_gr$Peak2GeneLinks

start_gr <- resize(p2g_500k_gr, width = 1, fix = "start")
end_gr <- resize(p2g_500k_gr, width = 1, fix = "end")

n1 <- findOverlaps(start_gr, neurod1_gr) %>% queryHits()
n2 <- findOverlaps(end_gr, neurod1_gr) %>% queryHits()
gr <- p2g_500k_gr[unique(sort(c(n1,n2)))]
#gr$FDR2 <- -log10(gr$FDR)
#gr$value <- gr$FDR2

width(gr)

p6 <- plotBrowserTrack(
  ArchRProj = arc,
  groupBy = "Clusters",
  region = GRanges(seqnames = "chr2", ranges = IRanges(start = 79390001, end = 79540000)),
  features = arc@peakSet,
  loops = gr)
grid::grid.newpage()
grid::grid.draw(p6)

pdf("output/Plots/04_GenomeTrack_Neurod1_v3.pdf", width = 10, height = 8)
grid::grid.draw(p6)
dev.off()

p6.2 <- plotBrowserTrack(
  ArchRProj = arc, 
  groupBy = "Clusters",
  region = GRanges(seqnames = "chr2", ranges = IRanges(start = 79390001, end = 79540000)),
  features = neurod1_rd,
  loops = gr)
grid::grid.newpage()
grid::grid.draw(p6.2)

pdf("output/Plots/04_GenomeTrack_Neurod1_v4.pdf", width = 10, height = 8)
grid::grid.draw(p6.2)
dev.off()

saveRDS(arc, "rds/04_arc.rds")

#-------- 05 cell type markers  --------#
# arc <- readRDS("rds/04_arc.rds")
GSMarkerClusters <- getMarkerFeatures(arc, groupBy = "Clusters", useMatrix = "GeneScoreMatrix")
GEMarkerClusters <- getMarkerFeatures(arc, groupBy = "Clusters", useMatrix = "GeneExpressionMatrix")

GSMarkers <- getMarkers(GSMarkerClusters, cutOff = "FDR < 0.05 & Log2FC > 0.5")
GEMarkers <- getMarkers(GEMarkerClusters, cutOff = "FDR < 0.05 & Log2FC > 0.5")

GSMarkersTable <- lapply(paste0("C", c(1:7)), function(i){
  out <- GSMarkers[[i]]
  if(nrow(out) != 0){
    out$clusters <- i
    return(out)
  } else {
    return(NULL)
    }
})
GSMarkersTable <- do.call(rbind, GSMarkersTable)

GEMarkersTable <- lapply(paste0("C", c(1:7)), function(i){
  out <- GEMarkers[[i]]
  if(nrow(out) != 0){
    out$clusters <- i
    return(out)
  } else {
    return(NULL)
  }
})
GEMarkersTable <- do.call(rbind, GEMarkersTable)

write.csv(GSMarkersTable, "output/Tables/05_ClusterMarkers_GeneScore.csv")
write.csv(GEMarkersTable, "output/Tables/05_ClusterMarkers_GeneExpression.csv")

GSMarkersTable <- read.csv("output/Tables/05_ClusterMarkers_GeneScore.csv")
GEMarkersTable <- read.csv("output/Tables/05_ClusterMarkers_GeneExpression.csv")

#gene_list
library(clusterProfiler)
GESigGenes <- lapply(paste0("C", c(1:7)), function(i) GEMarkersTable$name[GEMarkersTable$clusters == i & GEMarkersTable$FDR < 0.01 & GEMarkersTable$Log2FC > 1])
names(GESigGenes) <- paste0("C", c(1:7))

eGO <- lapply(paste0("C", c(1:7)), function(i) enrichGO(GESigGenes[[i]], OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.2))
p6.2 <- lapply(eGO, function(x) barplot(x, showCategory=20))
p6.2 <- p6.2[-3]
pdf("output/Plots/05_clusterProfiler_GO.pdf", width = 10, height = 5)
p6.2
dev.off()

#-------- 06 Quality plot  --------#
p7 <- plotGroups(arc, groupBy = "Sample", colorBy = "colData", 
                 name = c("TSSEnrichment", "nFrags", "Gex_nUMI", "Gex_nGenes", "Gex_MitoRatio"), 
                 plotAs = "violin", alpha = .2)
p8 <- plotTSSEnrichment(arc)
p9 <- plotFragmentSizes(arc)

pdf("output/Plots/06_QualityPlots.pdf", width = 5, height = 5)
p7
p8
p9
dev.off()
