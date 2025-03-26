#------------------------------------------------------------------------------
# 02_ATAC-analysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
source("code/edgeR_PairwiseFunction.R")
se <- readRDS("rds/ATAC_se.rds")

#--------- Identify differential accessible regions ---------#
DiffTest <- edgeR_pairwise(se, compareCol = "sampleType", topGroup = "Chd7KO", bottomGroup = "Tumor")
Chd7KD_DN  <- rownames(DiffTest)[which(assay(DiffTest)[,"log2FC"] < -1 & assay(DiffTest)[,"FDR"] < 0.25)]
Chd7KD_UP  <- rownames(DiffTest)[which(assay(DiffTest)[,"log2FC"] > 1 & assay(DiffTest)[,"FDR"] < 0.25)]

#export
peaks <- rowRanges(se)
gr <- GRangesList(Chd7KD_DN = peaks[Chd7KD_DN], Chd7KD_UP = peaks[Chd7KD_UP])
lapply(names(gr), function(x){
  g <- gr[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed/02r2_", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#visualize
df <- data.frame(Log2FC = assay(DiffTest)[,"log2FC"], FDR = -log10(assay(DiffTest)[,"FDR"]))
df$Sig <- "X"
df$Sig[which(df$Log2FC <  -1 & df$FDR > -log10(0.25))] <- "DN"
df$Sig[which(df$Log2FC > 1 & df$FDR > -log10(0.25))] <- "UP"
df <- df[order(df$Sig, decreasing = T),]

p1 <- ggplot(df, aes(x = Log2FC, y = FDR, color = Sig)) + ggrastr::geom_point_rast(size = 0.5) + ArchR::theme_ArchR() + 
  geom_vline(xintercept = 0, lty = "dotted") + ggtitle("ATAC") +
  scale_color_manual(values = c("DN" = "blue", "UP" = "red", "X" = "darkgray")) + labs(x = "Log2FC", y = "-log10(FDR)") 

pdf("output/Plots/02_VolcanoPlot_ATAC_v2.pdf", width = 4.5, height = 5)
p1
dev.off()
