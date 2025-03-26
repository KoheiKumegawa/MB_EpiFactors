#------------------------------------------------------------------------------
# 02r2_HistoneMarksDiff
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

#pseudo-count
assays(se_K4me1)$counts <- assays(se_K4me1)$counts+1
assays(se_K4me3)$counts <- assays(se_K4me3)$counts+1
assays(se_K27ac)$counts <- assays(se_K27ac)$counts+1
assays(se_K27me3)$counts <- assays(se_K27me3)$counts+1


diff_ls <- list(Kmt2cKO = list("Kmt2cKO", "NC"), Chd7KO = list("Chd7KO", "NC"))
DT_K4me1 <- lapply(diff_ls, function(x) edgeR_pairwise(se_K4me1, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DT_K4me3 <- lapply(diff_ls, function(x) edgeR_pairwise(se_K4me3, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DT_K27ac <- lapply(diff_ls, function(x) edgeR_pairwise(se_K27ac, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DT_K27me3 <- lapply(diff_ls, function(x) edgeR_pairwise(se_K27me3, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))

DT_list <- c(DT_K4me1, DT_K4me3, DT_K27ac, DT_K27me3)
names(DT_list) <- paste0(c("K4me1","K4me1","K4me3","K4me3","K27ac","K27ac","K27me3","K27me3"), c("_Kmt2cKO", "_Chd7KO"))
saveRDS(DT_list, "rds/02r2_DT_list.rds")

UP_peaks <- lapply(names(DT_list), function(i){
  se_tmp <- DT_list[[i]]
  out <- rowRanges(se_tmp)[which(assay(se_tmp)[,"log2FoldChange"] > 1 & assay(se_tmp)[,"FDR"] < 0.25)]
  names(out) <- NULL
  return(out)
}) %>% unlist
names(UP_peaks) <- paste0(names(DT_list), "_", "UP")
DN_peaks <- lapply(names(DT_list), function(i){
  se_tmp <- DT_list[[i]]
  out <- rowRanges(se_tmp)[which(assay(se_tmp)[,"log2FoldChange"] < -1 & assay(se_tmp)[,"FDR"] < 0.25)]
  names(out) <- NULL
  return(out)
}) %>% unlist
names(DN_peaks) <- paste0(names(DT_list), "_", "DN")
DiffPeaks <- c(UP_peaks,DN_peaks)

lapply(names(DiffPeaks), function(x){
  gr <- DiffPeaks[[x]]
  mcols(gr) <- NULL
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0("output/output_bed/02r2_", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

p1 <- lapply(names(DT_list), function(i){
  df <- data.frame(Log2FC = assay(DT_list[[i]])[,"log2FoldChange"], FDR = -log10(assay(DT_list[[i]])[,"FDR"]))
  df$Sig <- "X"
  df$Sig[which(df$Log2FC < -1 & df$FDR > -log10(0.25))] <- "DN"
  df$Sig[which(df$Log2FC > 1 & df$FDR > -log10(0.25))] <- "UP"
  df <- df[order(df$Sig, decreasing = T),]
  
  out <- ggplot(df, aes(x = Log2FC, y = FDR, color = Sig)) + ggrastr::geom_point_rast(size = 0.5) + ArchR::theme_ArchR() + 
    geom_vline(xintercept = 0, lty = "dotted") + ggtitle(i) +
    scale_color_manual(values = c("DN" = "blue", "UP" = "red", "X" = "darkgray")) + labs(x = "Log2FC", y = "-log10(FDR)") 
  return(out)
})
pdf("output/Plots/02r2_VolcanoPlot_DT_Histones.pdf", width = 4.5, height = 5)
p1
dev.off()
