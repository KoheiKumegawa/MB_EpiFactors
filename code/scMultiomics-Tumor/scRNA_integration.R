#scRNA integration
library(Seurat)
library(dplyr)
library(ArchR)
library(scales)
library(viridis)

d1 <- Read10X_h5("data/RS03060_filtered_feature_bc_matrix.h5")
d2 <- Read10X("data/data_fileterd/")

seu1 <- CreateSeuratObject(counts = d1$`Gene Expression`, project = "Ptch1KO")
seu2 <- CreateSeuratObject(counts = d2, project = "DKO")

seu <- merge(seu1, seu2)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-") # mm

p1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
p2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("output/Plots/GEx_QCplots.pdf")
p1
p2
p3
dev.off()

#quality filter
seu <- subset(seu, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 25)
table(seu$orig.ident)
# DKO Ptch1KO 
# 1645    6264 

saveRDS(seu, "rds/20241101_seu_preInteg.rds")

# split the dataset into a list of two seurat objects
seu.list <- SplitObject(seu, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu.list)
anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)
seu_integrated <- IntegrateData(anchorset = anchors)

# Run the standard workflow for visualization and clustering
DefaultAssay(seu_integrated) <- "integrated"
seu_integrated <- ScaleData(seu_integrated, verbose = FALSE)
seu_integrated <- RunPCA(seu_integrated, npcs = 50, verbose = FALSE)
seu_integrated <- RunUMAP(seu_integrated, reduction = "pca", dims = 1:50)
seu_integrated <- FindNeighbors(seu_integrated, reduction = "pca", dims = 1:50)
seu_integrated <- FindClusters(seu_integrated, resolution = 0.2)

cluster_colors <- ArchRPalettes$kelly[c(1:11)] %>% `names<-`(.,paste0("IC", c(1:11)))
sample_colors <- ArchRPalettes$calm[c(1:2)] %>% `names<-`(., c("Ptch1KO","DKO"))

seu_integrated$Clusters <- paste0("IC", as.numeric(seu_integrated$seurat_clusters))
p4 <- DimPlot(seu_integrated, reduction = "umap", label = T, group.by = "seurat_clusters", raster = T) + theme_ArchR()
p5 <- DimPlot(seu_integrated, reduction = "umap", cols = sample_colors, label = F, group.by = "orig.ident", raster = T) + theme_ArchR()
p6 <- DimPlot(seu_integrated, reduction = "umap", cols = cluster_colors, label = T, group.by = "Clusters", raster = T) + theme_ArchR()
pdf("output/Plots/I1_UMAP_Integrated.pdf", width = 4, height = 5)
p4
p5
p6
dev.off()

cM <- table(seu_integrated$Clusters, seu_integrated$orig.ident)
write.csv(cM, "output/Tables/I1_ClusterSampleMatrix.csv")

cM2 <- cM[c(1:3),]
t(cM2) / colSums(cM2) * 100
# IC1      IC2      IC3
# DKO     28.64516 41.29032 30.06452
# Ptch1KO 59.94511 29.30255 10.75234

cM2 <- reshape2::melt(cM2)
cM2$Var2 <- factor(cM2$Var2, levels = c("Ptch1KO", "DKO"))
p7 <- ggplot(cM2, aes(x = Var2, y = value, fill = Var1)) +
  geom_bar(stat = "identity", position = "fill") + scale_y_continuous(labels = percent) +
  scale_fill_manual(values = cluster_colors) + theme_ArchR() +
  labs(x = "Sample", y = "%Cluster", fill = "Cluster") +
  theme(legend.key.size = unit(0.1, 'cm'), axis.text.x = element_text(angle = 90, hjust = 1))
pdf("output/Plots/I1_Barplot_TumorClusterSample.pdf", width = 4, height = 5)
p7
dev.off()

#UMAP overlay
markerGenes <- c("Gli1", "Gli2", "Ccnd2", "Barhl1", "Atoh1", "Cntn2", "Grin2b", "Mki67",  "Mcm2", "Ccnd1", "Ccnb1", "Neurod1")
p8 <- lapply(markerGenes, function(x) FeaturePlot(seu_integrated, features = x, pt.size = 2, order = T, raster = T) + 
               scale_colour_gradientn(colours = viridis(256, option = "A")) + theme_ArchR())
pdf("output/Plots/I1_UMAPOL_markers1.pdf", height = 5, width = 4)
p8
dev.off()

#Violin
p9 <- lapply(markerGenes, function(x) VlnPlot(seu_integrated, group.by = "Clusters", cols = cluster_colors, assay = "RNA", features = x, pt.size = 0) + theme_ArchR())
VlnPlot(seu_integrated, features = markerGenes, pt.size = 0)
pdf("output/Plots/I1_VlnPlot_Integrated_MarkerGene1.pdf", width = 5, height = 4)
p9
dev.off()

#DEG
Idents(seu_integrated) <- "Clusters"
DEG <- FindAllMarkers(seu_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DEG, "output/Tables/I1_DEGClusters.csv")

#SKO vs DKO in IC3
seu_tmp <- seu_integrated[,seu_integrated$Clusters == "IC3"]
seu_tmp$Sample <- factor(seu_tmp$orig.ident, levels = c("Ptch1KO","DKO"))

p10 <- VlnPlot(seu_tmp, assay = "RNA", group.by = "Sample", cols = sample_colors, features = "Neurod1", pt.size = 0) + theme_ArchR()
p11 <- VlnPlot(seu_tmp, assay = "RNA", group.by = "Sample", cols = sample_colors, features = "Mki67", pt.size = 0) + theme_ArchR()

pdf("output/Plots/I1_VlnPlot_IC3_Mki67_Neurod1.pdf", width = 3, height = 4)
p10
p11
dev.off()

saveRDS(seu_integrated, "rds/I1_seu_integrated.rds")
