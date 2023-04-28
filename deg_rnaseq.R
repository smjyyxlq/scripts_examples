setwd("~/Desktop/cam_rnaseq/")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2)) 
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(pheatmap))
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)

# ---- ReadsPerGene to raw counts table ----

filepath <- "~/Desktop/cam_rnaseq/Raw_/"
countfiles <- grep(".ReadsPerGene.out.tab",list.files(filepath), value = T)
l <- lapply(countfiles, function(x){ read.table(paste0(filepath, "/", file.path(x)), fill = T)[-(1:4),]})
# check rows
if( ! all( sapply( l, function(x) all( x$V1 == l[[1]]$V1 ) ) ) ){
  stop( "Gene IDs differ between files" )}

col_use <- 4  # -s reverse or stranded == "+-,-+"
counts_raw <- as.data.frame(sapply( l, function(x) x[, col_use] ))
rownames(counts_raw) <- l[[1]]$V1
colnames <- sub(".ReadsPerGene.out.tab", "", countfiles)
colnames(counts_raw) <- colnames

orders <- c("ESC_B", "ESC_C", "ESC_D", 
            "EpiLC_A", "EpiLC_B", "EpiLC_C", 
            "PGC_A",	"PGC_B",	"PGC_C")

counts_raw <- counts_raw[, orders, drop=F]

# ---- sample table ----

celltype <- c(rep("ESC", 3), rep("EpiLC", 3), rep("PGC", 3))
reps <- rep(c("rep1", "rep2", "rep3"), 3)
samplename <- paste(celltype, reps, sep = "_")
sampleTable <- data.frame(cellType = factor(celltype, levels = c("ESC", "EpiLC", "PGC")))
rownames(sampleTable) <- samplename
colnames(counts_raw) <- samplename
saveRDS(counts_raw, "cam_rnaseq.counts_raw.rds")

# all(rownames(sampleTable) %in% colnames(counts_raw))
set.seed(314)

# ---- DEG ----
ds_cam <- DESeqDataSetFromMatrix(countData = counts_raw, 
                                 colData = sampleTable,
                                 design = ~cellType)
dds_cam <- DESeq(ds_cam, test = "Wald")
saveRDS(dds_cam, "cam_rnaseq.dds.rds")

# annotate gene ID to symbol
id2symbol <- read.delim("~/Desktop/mm10_v93/mm10_v93.id2symbol",
                        header = T, stringsAsFactors = F)

counts_raw <- mutate(counts_raw, geneID = rownames(counts_raw))
counts_raw$geneSymbol <- id2symbol[
  match(counts_raw$geneID, counts_raw$geneID),
]$geneSymbol
write.table(counts_raw, "cam_rnaseq_xxxxxx.counts_raw.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# output counts normalized by size factor
counts_norm <- as.data.frame(counts(de_cam, normalized=T))
counts_norm <- mutate(counts_norm, geneID = rownames(counts_norm))
counts_norm$geneSymbol <- id2symbol[
  match(counts_norm$geneID, counts_norm$geneID),
]$geneSymbol
write.table(counts_norm, "cam_rnaseq_xxxxxx.counts_normalized.txt", sep = "\t", quote = F, row.names = T, col.names = T)

# ---- 1. EpiLC vs ESC -----
comparison <- "EpiLC_vs_ESC"
res1 <- results(dds_cam, contrast = c("cellType", "EpiLC", "ESC"), 
                alpha = 0.05, lfcThreshold = 0.58)  # to test whether |LFC| > 0.58 
sres1 <- lfcShrink(dds_cam, coef = "cellType_EpiLC_vs_ESC",
                   res = res1, type = "apeglm")  # update posterior estimates of LFC only
# range(abs(sres1$log2FoldChange))
sig_res1 <- as.data.frame(subset(sres1, padj < 0.05))
sig_res1$geneSymbol <- id2symbol[match(rownames(sig_res1), id2symbol$geneID),]$geneSymbol # add gene symbol
sig_res1_up <- filter(sig_res1, log2FoldChange > 0)
sig_res1_dn <- filter(sig_res1, log2FoldChange < 0)
write.table(sig_res1_up, paste0("DE_cam.", comparison, ".lfc058padj005_apeglm_up.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)
write.table(sig_res1_dn, paste0("DE_cam.", comparison, ".lfc058padj005_apeglm_dn.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# ---- 3. PGCLC vs ESC -----
comparison <- "PGC_vs_ESC"
res3 <- results(dds_cam, contrast = c("cellType", "PGC", "ESC"), 
                alpha = 0.05, lfcThreshold = 0.58)  
sres3 <- lfcShrink(dds_cam, coef = "cellType_PGC_vs_ESC",
                   res = res3, type = "apeglm")  
# range(abs(sres3$log2FoldChange))
sig_res3 <- as.data.frame(subset(sres3, padj < 0.05))
sig_res3$geneSymbol <- id2symbol[match(rownames(sig_res3), id2symbol$geneID),]$geneSymbol
sig_res3_up <- filter(sig_res3, log2FoldChange > 0)
sig_res3_dn <- filter(sig_res3, log2FoldChange < 0)
write.table(sig_res3_up, paste0("DE_cam.", comparison, ".lfc058padj005_apeglm_up.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)
write.table(sig_res3_dn, paste0("DE_cam.", comparison, ".lfc058padj005_apeglm_dn.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)

# ---- 2. PGC vs EpiLC -----
comparison <- "PGC_vs_EpiLC"

# To get coef for "cellType_PGC_vs_EpiLC" for lfcShrink using apeglm
dds_cam$cellType <- relevel(dds_cam$cellType, ref = "EpiLC")
dds_cam <- nbinomWaldTest(dds_cam)
resultsNames(dds_cam)
## [1] "Intercept"  "cellType_ESC_vs_EpiLC"     "cellType_PGC_vs_EpiLC"
res2 <- results(dds_cam, contrast = c("cellType", "PGC", "EpiLC"), 
                alpha = 0.05, lfcThreshold = 0.58)
sres2 <- lfcShrink(dds_cam, coef = "cellType_PGC_vs_EpiLC",
                   res = res2, type = "apeglm")
# range(abs(sres2$log2FoldChange))
sig_res2 <- as.data.frame(subset(sres2, padj < 0.05))
sig_res2$geneSymbol <- id2symbol[match(rownames(sig_res2), id2symbol$geneID),]$geneSymbol
sig_res2_up <- filter(sig_res2, log2FoldChange > 0)
sig_res2_dn <- filter(sig_res2, log2FoldChange < 0)
write.table(sig_res2_up, paste0("DE_cam.", comparison, ".lfc058padj005_apeglm_up.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)
write.table(sig_res2_dn, paste0("DE_cam.", comparison, ".lfc058padj005_apeglm_dn.txt"),
            row.names = T, col.names = T, sep = "\t", quote = F)


# ---- some visualization ----
source("~/Desktop/x_visul/x_theme.R")
# PCA
rld_cam <- rlog(dds_cam, blind = T)
saveRDS(rld_cam, "cam_rnaseq.rld.rds")
topVarGenes <- head(order(rowVars(assay(rld_cam)), decreasing = T), 500) # use top 500 varient genes 
rld_top <- assay(rld_cam)[ topVarGenes, ]
pca <- prcomp(t(rld_top))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
#keepdf <- as.data.frame(colData(rld_cam)[, "cellType", drop=F])
df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sampleTable, name=colnames(rld_cam))

pdf("PCA.cam_xxxx.3cells.pdf")
ggplot(data=df, aes_string(x="PC1", y="PC2", label = "name")) +
  geom_point(
    size=3,
    alpha = 0.8,
    position = position_jitter(width = -0.5, height = 0.5), colour = "black", shape = 21) +
  xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
  coord_fixed() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm")) +  
  x_theme() +
  geom_text_repel(
    #vjust = 1.5,
    #hjust = -0.5,
    size = 2.5,
    #colour = c(rep(c("#2b5466", "#DD2D4A"), 4), rep("#FFCC00", 4),
    #           rep(c("#56b4e9", "#FFD9DA"), 2), rep("#FFCC00", 2)
    )
dev.off()



# heatmap for genelist from clustering 
cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
k = 10
set.seed(314)
pattern <- pheatmap(rld_cam, kmeans_k = k, 
                    scale="row", color = cols,
                    border_color = NA, 
                    cellwidth = 20, cellheight = 8, cluster_cols=F, display_numbers = F)
pdf(paste0("heatmap.cam_rnaseq.xxxx.k", k, ".pdf"))
pattern
dev.off()

for(i in 1:k){
  july_i <- pattern$kmeans$cluster[which(pattern$kmeans$cluster == i)]
  july_i <- names(july_i)
  i_count <- rld_cam[july_i]
  i_count <- as.data.frame(i_count)
  i_count$geneID <- rownames(i_count)
  i_count$geneSymbol <- id2symbol[match(
    i_count$geneID, id2symbol$geneID),]$geneSymbol
  write.table(i_count, paste0("cam_rnaseq.xxxx.k", k, "_", i, "_list.txt"), 
              row.names = F, col.names = T, sep = "\t", quote = F)
  rownames(i_count) <- i_count$geneSymbol
  i_count <- i_count[,1:9]
  pdf(paste0("heatmap.cam_de_rld.k", k, "_", i, ".pdf"), 7, 20)
  pheatmap(i_count,
           scale="row",
           color = cols,
           cellwidth = 8, cellheight = 4,
           border_color=NA,
           fontsize_row = 4,
           fontsize_col = 6,
           fontsize = 6, legend = T,
           annotation_legend = T,
           cluster_rows= T,
           cluster_cols = F
  )
  dev.off()
}