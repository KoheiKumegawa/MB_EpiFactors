#------------------------------------------------------------------------------
# 02_DiffTest.R
#------------------------------------------------------------------------------
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(Rcpp)
library(Rsamtools)
library(data.table)
library(ggplot2)
library(rtracklayer)
library(SummarizedExperiment)
library(ggrepel)
source("code/edgeR_PairwiseFunction.R")
'%ni%' <- Negate("%in%")

#--------- DiffPeaks ---------#
se_K4me1 <- readRDS("rds/01r_se_K4me1.rds")
se_K4me3 <- readRDS("rds/01r_se_K4me3.rds")
se_K27ac <- readRDS("rds/01r_se_K27ac.rds")
se_K27me3<- readRDS("rds/01r_se_K27me3.rds")

diff_ls <- list(Kmt2cKO = list("Kmt2cKO", "NC"), Chd7KO = list("Chd7KO", "NC"))
DT_K4me1 <- lapply(diff_ls, function(x) edgeR_pairwise(se_K4me1, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DT_K4me3 <- lapply(diff_ls, function(x) edgeR_pairwise(se_K4me3, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DT_K27ac <- lapply(diff_ls, function(x) edgeR_pairwise(se_K27ac, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DT_K27me3 <- lapply(diff_ls, function(x) edgeR_pairwise(se_K27me3, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))

#--------- |Log2FC| > 1 & P-value < 0.01 ---------#
#H3K4me1
K4me1_Kmt2cKO <- rowRanges(se_K4me1)[which(assay(DT_K4me1[[1]])[,"log2FoldChange"] < -1 & assay(DT_K4me1[[1]])[,"pval"] < 0.01)]
mcols(K4me1_Kmt2cKO) <- NULL
d <- data.frame(seqnames = seqnames(K4me1_Kmt2cKO), start = start(K4me1_Kmt2cKO)-1, end = end(K4me1_Kmt2cKO))
write.table(d, "output/output_bed/02r_K4me1_DN_Kmt2cKO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

K4me1_Chd7KO <- rowRanges(se_K4me1)[which(assay(DT_K4me1[[2]])[,"log2FoldChange"] < -1 & assay(DT_K4me1[[2]])[,"pval"] < 0.01)]
mcols(K4me1_Chd7KO) <- NULL
d <- data.frame(seqnames = seqnames(K4me1_Chd7KO), start = start(K4me1_Chd7KO)-1, end = end(K4me1_Chd7KO))
write.table(d, "output/output_bed/02r_K4me1_DN_Chd7KO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#H3K4me3
K4me3_Kmt2cKO <- rowRanges(se_K4me3)[which(assay(DT_K4me3[[1]])[,"log2FoldChange"] < -1 & assay(DT_K4me3[[1]])[,"pval"] < 0.01)]
mcols(K4me3_Kmt2cKO) <- NULL
d <- data.frame(seqnames = seqnames(K4me3_Kmt2cKO), start = start(K4me3_Kmt2cKO)-1, end = end(K4me3_Kmt2cKO))
write.table(d, "output/output_bed/02r_K4me3_DN_Kmt2cKO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

K4me3_Chd7KO <- rowRanges(se_K4me3)[which(assay(DT_K4me3[[2]])[,"log2FoldChange"] < -1 & assay(DT_K4me3[[2]])[,"pval"] < 0.01)]
mcols(K4me3_Chd7KO) <- NULL
d <- data.frame(seqnames = seqnames(K4me3_Chd7KO), start = start(K4me3_Chd7KO)-1, end = end(K4me3_Chd7KO))
write.table(d, "output/output_bed/02r_K4me3_DN_Chd7KO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#H3K27ac
K27ac_Kmt2cKO <- rowRanges(se_K27ac)[which(assay(DT_K27ac[[1]])[,"log2FoldChange"] < -1 & assay(DT_K27ac[[1]])[,"pval"] < 0.01)]
mcols(K27ac_Kmt2cKO) <- NULL
d <- data.frame(seqnames = seqnames(K27ac_Kmt2cKO), start = start(K27ac_Kmt2cKO)-1, end = end(K27ac_Kmt2cKO))
write.table(d, "output/output_bed/02r_K27ac_DN_Kmt2cKO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

K27ac_Chd7KO <- rowRanges(se_K27ac)[which(assay(DT_K27ac[[2]])[,"log2FoldChange"] < -1 & assay(DT_K27ac[[2]])[,"pval"] < 0.01)]
mcols(K27ac_Chd7KO) <- NULL
d <- data.frame(seqnames = seqnames(K27ac_Chd7KO), start = start(K27ac_Chd7KO)-1, end = end(K27ac_Chd7KO))
write.table(d, "output/output_bed/02r_K27ac_DN_Chd7KO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

#H3K27me3
K27me3_Kmt2cKO <- rowRanges(se_K27me3)[which(assay(DT_K27me3[[1]])[,"log2FoldChange"] < -1 & assay(DT_K27me3[[1]])[,"pval"] < 0.01)]
mcols(K27me3_Kmt2cKO) <- NULL
d <- data.frame(seqnames = seqnames(K27me3_Kmt2cKO), start = start(K27me3_Kmt2cKO)-1, end = end(K27me3_Kmt2cKO))
write.table(d, "output/output_bed/02r_K27me3_DN_Kmt2cKO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

K27me3_Chd7KO <- rowRanges(se_K27me3)[which(assay(DT_K27me3[[2]])[,"log2FoldChange"] < -1 & assay(DT_K27me3[[2]])[,"pval"] < 0.01)]
mcols(K27me3_Chd7KO) <- NULL
d <- data.frame(seqnames = seqnames(K27me3_Chd7KO), start = start(K27me3_Chd7KO)-1, end = end(K27me3_Chd7KO))
write.table(d, "output/output_bed/02r_K27me3_DN_Chd7KO_p001.bed", row.names = F, col.names = F, quote = F, sep = "\t")

DT_list <- c(DT_K4me1, DT_K4me3, DT_K27ac, DT_K27me3)
names(DT_list) <- paste0(c("K4me1","K4me1","K4me3","K4me3","K27ac","K27ac","K27me3","K27me3"), c("_Kmt2cKO", "_Chd7KO"))

p1 <- lapply(names(DT_list), function(i){
  df <- data.frame(Log2FC = assay(DT_list[[i]])[,"log2FoldChange"], P = -log10(assay(DT_list[[i]])[,"pval"]))
  df$Sig <- "X"
  df$Sig[which(df$Log2FC <  -1 & df$P > 2)] <- "DN"
  df$Sig[which(df$Log2FC > 1 & df$P > 2)] <- "UP"
  df <- df[order(df$Sig, decreasing = T),]
  
  out <- ggplot(df, aes(x = Log2FC, y = P, color = Sig)) + ggrastr::geom_point_rast(size = 0.5) + ArchR::theme_ArchR() + 
    geom_vline(xintercept = 0, lty = "dotted") + ggtitle(i) +
    scale_color_manual(values = c("DN" = "blue", "UP" = "red", "X" = "darkgray")) + labs(x = "Log2FC", y = "-log10(P-value)") 
  return(out)
})
pdf("output/Plots/02r_VolcanoPlot_DT_Histones.pdf", width = 4.5, height = 5)
p1
dev.off()
