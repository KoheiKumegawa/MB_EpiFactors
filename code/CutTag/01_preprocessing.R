#------------------------------------------------------------------------------
# 01_preprocessing
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
'%ni%' <- Negate("%in%")

#------------------
# functions
#------------------
bamToFragmentGR <- function(
  bamPATH = NULL,
  bamNAME = NULL,
  offsetPlus = 4,
  offsetMinus = -5,
  bamFlag = NULL
){
  if(is.null(bamPATH)){
    stop("Please set PATH to bam files")
  }
  if(is.null(bamNAME)){
    stop("No input bamNAME; please recheck your input")
  }
  if(is.null(bamFlag)){
    stop("Please set bamFlag using Rsamtools's scanBamFlag!")
  }
  
  #1. Read In Bam File
  sF <- scanBam(bamPATH, param = ScanBamParam(flag = bamFlag, what = c("rname","pos", "isize")))[[1]]
  
  #2. Make Fragment Table
  dt <- data.table(seqnames = sF$rname, start = sF$pos + offsetPlus, end = sF$pos + abs(sF$isize) - 1 + offsetMinus)
  
  #3. Make fragment Granges and remove unwanted chromosomes
  gr <- GRanges(seqnames = dt$seqnames, IRanges(start = dt$start, end = dt$end))
  idy = which(seqnames(gr) %in% seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  gr <- gr[-idy]
  gr <- dropSeqlevels(gr, seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  mcols(gr) <- DataFrame(sample = bamNAME)
  
  #4. output Granges List
  return(gr)
}

FragmentGRToBED <- function(
  gr = NULL,
  name = NULL,
  outputDir = NULL
){
  d <- data.frame(seqnames = seqnames(gr), start = start(gr)-1, end = end(gr))
  write.table(d, paste0(outputDir, "/", name, "-fragments.bed"), row.names = F, col.names = F, quote = F, sep = "\t")
  return(NULL)
}

RunMacs2 <- function(
  inputBamPATH = NULL,
  inputBamName = NULL,
  genome = NULL,
  outputDir = NULL,
  method = c("p", "q"),
  cutoff = 0.05
){
  if(genome %ni% c("hg19", "hg38", "mm10")){
    stop("Please set genome as hg19, hg38 and mm10!")
  }
  if(genome %in% c("hg19", "hg38")){
    gen <- "hs"
  }
  if(genome == "mm10"){
    gen <- "mm"
  }
  
  commandPeaks <- sprintf("macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BAM --call-summits --keep-dup auto", 
                          gen, inputBamName, inputBamPATH, outputDir)
  
  if (tolower(method) == "p") {
    commandPeaks <- sprintf("%s -p %s", commandPeaks, cutoff)
  } else {
    commandPeaks <- sprintf("%s -q %s", commandPeaks, cutoff)
  }
  
  message("Running Macs2...")
  message(commandPeaks)
  system(commandPeaks, intern = TRUE)
  
  return(NULL)
}

MakeSamplePeakSet <- function(gr, by = "score"){
  #nonOverlappingGRanges
  stopifnot(by %in% colnames(mcols(gr)))
  
  #function for picking up most significant peaks
  clusterGRanges <- function(gr, by = "score"){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = TRUE),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  #iteration of filtering overlapping peaks
  i <-  0
  gr_converge <- gr
  while(length(gr_converge) > 0){
    i <-  i + 1
    gr_selected <- clusterGRanges(gr = gr_converge, by = by)
    gr_converge <- subsetByOverlaps(gr_converge, gr_selected, invert=TRUE) #blacklist selected gr
    if(i == 1){ #if i=1 then set gr_all to clustered
      gr_all <- gr_selected
    }else{
      gr_all <- c(gr_all, gr_selected)
    }
  }
  gr_all <- sort(sortSeqlevels(gr_all))
  return(gr_all)
}

MakeCTSummarizedExperiment <- function(
  fragmentGRangesList = NULL,
  unionPeaks = NULL,
  blacklist = NULL,
  sampleName = NULL,
  prior.count = 1,
  by = "sample"
){
  fragments <- unlist(as(fragmentGRangesList, "GRangesList"))
  overlapDF <- DataFrame(findOverlaps(unionPeaks, fragments, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(fragments)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(unionPeaks), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  sparseM <- sparseM[, sampleName]
  rownames(sparseM) <- unionPeaks$name
  sparseM.cpm <- edgeR::cpm(sparseM, log = TRUE, prior.count = prior.count)
  rownames(sparseM.cpm) <- unionPeaks$name
  #SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = sparseM, log2cpm = sparseM.cpm),
                                                   rowRanges = unionPeaks, 
                                                   colData = DataFrame(sample = sampleName))
  return(se)
}

#------------------
# analysis
#------------------

#--------- Importing fragments ---------#
bamFiles <- read.csv("data/ct_samples_v3.csv", header = T)
bamFiles$AnalysisName <- paste0(gsub("_", "", bamFiles$Condition), "_", bamFiles$Antibody, "_", bamFiles$Rep)
fragments <- parallel::mclapply(seq_along(bamFiles$bamFile),
                                function(x){
                                  out <- bamToFragmentGR(bamPATH = paste0("data/", bamFiles$bamFile[x]), 
                                                         bamNAME = bamFiles$AnalysisName[x],
                                                         bamFlag = scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE))
                                  return(out)
                                }, mc.cores = 8)
fragments <- GRangesList(fragments)
names(fragments) <- bamFiles$AnalysisName

#--------- Calculating and Visualizing Quality ---------#
#fragment size
fragmentsSizeCT <- lapply(fragments, function(x) width(x))
df <- lapply(names(fragmentsSizeCT), function(i) data.frame(l = fragmentsSizeCT[[i]], sample = i)) %>% do.call(rbind, .)
p1 <- ggplot(df, aes(x = l, color = sample)) + geom_line(stat = "density", size = 0.5) + ArchR::theme_ArchR() + 
  xlim(0, 600) + ylab("Density") + xlab("Size of fragments (bp)") + 
  scale_color_manual(values = c(ArchR::ArchRPalettes$stallion, ArchR::ArchRPalettes$calm) %>% `names<-`(.,unique(df$sample)))
pdf("output/Plots/01r_FragmentWidth.pdf", width = 4, height = 6)
p1
dev.off()

#output
parallel::mclapply(names(fragments), 
                   function(x) FragmentGRToBED(gr = fragments[[x]], name = x, outputDir = "output/fragments_bed_v2/"), 
                   mc.cores = 16)

#--------- Peak call ---------#
parallel::mclapply(seq_along(bamFiles$bamFile), function(x){
  RunMacs2(inputBamPATH = paste0("data/", bamFiles$bamFile[x]),
           inputBamName = bamFiles$AnalysisName[x],
           genome = "mm10",
           outputDir = "output/sample_peaks_v2/",
           method = "q",
           cutoff = 0.1)
}, mc.cores = 16)

#--------- Making peak set per sample ---------#
BSgenome   <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
blacklist  <- rtracklayer::import.bed("ref/mm10-blacklist.v2.bed")

gr_ls <- lapply(bamFiles$AnalysisName, function(x) import.bed(paste0("output/sample_peaks_v2/",x, "_summits.bed")))
names(gr_ls) <- bamFiles$AnalysisName

gr_ls_proc <- parallel::mclapply(gr_ls, function(x){
  gr <- resize(x, width = 501, fix = "center") %>%
    subsetByOverlaps(., chromSizes, type = "within") %>%
    subsetByOverlaps(., blacklist, invert=TRUE) %>%
    MakeSamplePeakSet(., by = "score")
  mcols(gr)$scorePerMillion <- mcols(gr)$score / (sum(mcols(gr)$score) / 1000000)
  return(gr)
}, mc.cores = 16) %>% GenomicRanges::GRangesList(.)
saveRDS(gr_ls_proc, "rds/01r_gr_ls_proc.rds")

#--------- Making Consensus Peak Set ---------#
#H3K4me1
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K4me1")]
tmp <- MakeSamplePeakSet(unlist(gr_ls_proc[idx]), by = "scorePerMillion")
mcols(tmp)$sampleOverlap <- countOverlaps(tmp, gr_ls_proc[idx])
K4me1_peaks <- tmp[which(mcols(tmp)$sampleOverlap >= 2 & seqnames(tmp) %ni% c("chrY", "chrM"))]
mcols(K4me1_peaks)$name2 <- paste0("H3K4me1_", seq_along(K4me1_peaks))

#H3K4me3
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K4me3")]
tmp <- MakeSamplePeakSet(unlist(gr_ls_proc[idx]), by = "scorePerMillion")
mcols(tmp)$sampleOverlap <- countOverlaps(tmp, gr_ls_proc[idx])
K4me3_peaks <- tmp[which(mcols(tmp)$sampleOverlap >= 2 & seqnames(tmp) %ni% c("chrY", "chrM"))]
mcols(K4me3_peaks)$name2 <- paste0("H3K4me3_", seq_along(K4me3_peaks))

#H3K27ac
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K27ac")]
tmp <- MakeSamplePeakSet(unlist(gr_ls_proc[idx]), by = "scorePerMillion")
mcols(tmp)$sampleOverlap <- countOverlaps(tmp, gr_ls_proc[idx])
K27ac_peaks <- tmp[which(mcols(tmp)$sampleOverlap >= 2 & seqnames(tmp) %ni% c("chrY", "chrM"))]
mcols(K27ac_peaks)$name2 <- paste0("H3K27ac_", seq_along(K27ac_peaks))

#H3K27me3
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K27me3")]
tmp <- MakeSamplePeakSet(unlist(gr_ls_proc[idx]), by = "scorePerMillion")
mcols(tmp)$sampleOverlap <- countOverlaps(tmp, gr_ls_proc[idx])
K27me3_peaks <- tmp[which(mcols(tmp)$sampleOverlap >= 2 & seqnames(tmp) %ni% c("chrY", "chrM"))]
mcols(K27me3_peaks)$name2 <- paste0("H3K27me3_", seq_along(K27me3_peaks))

gr <- GRangesList(K4me1_peaks = K4me1_peaks,
                  K4me3_peaks = K4me3_peaks,
                  K27ac_peaks = K27ac_peaks,
                  K27me3_peaks = K27me3_peaks)
saveRDS(gr, "rds/01r_AllPeaks_GRL.rds")

#--------- Constructing NFI SE ---------#
#H3K4me1
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K4me1")]
fragments_tmp <- fragments[idx]
fragments_named <- lapply(names(fragments_tmp), function(x){
  fr <- fragments_tmp[[x]]
  mcols(fr)$sample <- x
  return(fr)
})
names(fragments_named) <- names(fragments_tmp)
fragments_named <- GRangesList(fragments_named)
se_K4me1 <- MakeCTSummarizedExperiment(fragmentGRangesList = fragments_named,
                                       unionPeaks = K4me1_peaks,
                                       blacklist = blacklist,
                                       sampleName = names(fragments_named),
                                       by = "sample",
                                       prior.count = 5)
se_K4me1$sampleType <- stringr::str_split(se_K4me1$sample, "_", simplify = T)[,1]
se_K4me1$rep <- stringr::str_split(se_K4me1$sample, "_", simplify = T)[,3]

#H3K4me3
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K4me3")]
fragments_tmp <- fragments[idx]
fragments_named <- lapply(names(fragments_tmp), function(x){
  fr <- fragments_tmp[[x]]
  mcols(fr)$sample <- x
  return(fr)
})
names(fragments_named) <- names(fragments_tmp)
fragments_named <- GRangesList(fragments_named)
se_K4me3 <- MakeCTSummarizedExperiment(fragmentGRangesList = fragments_named,
                                       unionPeaks = K4me3_peaks,
                                       blacklist = blacklist,
                                       sampleName = names(fragments_named),
                                       by = "sample",
                                       prior.count = 5)
se_K4me3$sampleType <- stringr::str_split(se_K4me3$sample, "_", simplify = T)[,1]
se_K4me3$rep <- stringr::str_split(se_K4me3$sample, "_", simplify = T)[,3]

#H3K27ac
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K27ac")]
fragments_tmp <- fragments[idx]
fragments_named <- lapply(names(fragments_tmp), function(x){
  fr <- fragments_tmp[[x]]
  mcols(fr)$sample <- x
  return(fr)
})
names(fragments_named) <- names(fragments_tmp)
fragments_named <- GRangesList(fragments_named)
se_K27ac <- MakeCTSummarizedExperiment(fragmentGRangesList = fragments_named,
                                       unionPeaks = K27ac_peaks,
                                       blacklist = blacklist,
                                       sampleName = names(fragments_named),
                                       by = "sample",
                                       prior.count = 5)
se_K27ac$sampleType <- stringr::str_split(se_K27ac$sample, "_", simplify = T)[,1]
se_K27ac$rep <- stringr::str_split(se_K27ac$sample, "_", simplify = T)[,3]

#H3K27me3
idx <- bamFiles$AnalysisName[which(bamFiles$Antibody == "H3K27me3")]
fragments_tmp <- fragments[idx]
fragments_named <- lapply(names(fragments_tmp), function(x){
  fr <- fragments_tmp[[x]]
  mcols(fr)$sample <- x
  return(fr)
})
names(fragments_named) <- names(fragments_tmp)
fragments_named <- GRangesList(fragments_named)
se_K27me3 <- MakeCTSummarizedExperiment(fragmentGRangesList = fragments_named,
                                        unionPeaks = K27me3_peaks,
                                        blacklist = blacklist,
                                        sampleName = names(fragments_named),
                                        by = "sample",
                                        prior.count = 5)
se_K27me3$sampleType <- stringr::str_split(se_K27me3$sample, "_", simplify = T)[,1]
se_K27me3$rep <- stringr::str_split(se_K27me3$sample, "_", simplify = T)[,3]

saveRDS(se_K4me1, "rds/01r_se_K4me1.rds")
saveRDS(se_K4me3, "rds/01r_se_K4me3.rds")
saveRDS(se_K27ac, "rds/01r_se_K27ac.rds")
saveRDS(se_K27me3, "rds/01r_se_K27me3.rds")
