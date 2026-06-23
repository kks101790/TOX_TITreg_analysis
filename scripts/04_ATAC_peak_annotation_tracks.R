# --- explicit package loading (added for standalone reproducibility; analysis logic UNCHANGED;
#     original scripts assumed an interactive session with these already attached) ---
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))

# =============================================================================
# ATAC-seq peak annotation (ChIPseeker) + locus tracks (trackViewer)
# Manuscript: "TOX enforces the immunosuppressive program of tumor-infiltrating
#              regulatory T cells" (Park et al., Nature Immunology)
# Figures: Fig 7b, 7c
# Environment: ATAC env (R 4.4.1)   (exact versions in ../envs/ and README.md)
# Provenance: original project script "Analyzing_ATACseq_2_Visualize.R" — analysis logic UNCHANGED; only this
#             header added. Configure input/output paths per README directory layout.
# =============================================================================

#!/usr/bin/Rscript
#usage : ./plotting_tsne.R [matrix file] [output prefix]

cKO <- read.delim("./../Data/ATACseq/result_ATACseq/2.Peak_Annotation/Tox_cKO/Tox_cKO_annotatedPeakList_blfiltered.csv", sep='\t')
FF <- read.delim("./../Data/ATACseq/result_ATACseq/2.Peak_Annotation/Tox_FF/Tox_FF_annotatedPeakList_blfiltered.csv", sep='\t')


feat <- c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Satb1", "Tnf", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

library(ChIPseeker)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

peaks.cKO <- readPeakFile("Tox_cKO_peaks.filt.broadPeak")
peaks.FF <- readPeakFile("Tox_FF_peaks.filt.broadPeak")

txdb <- makeTxDbFromGFF("mm10.refGene.gtf")
txdb.id <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnno.cKO <- annotatePeak("Tox_cKO_peaks.filt.broadPeak", tssRegion=c(-3000, 3000), TxDb=txdb)
peakAnno.FF <- annotatePeak("Tox_FF_peaks.filt.broadPeak", tssRegion=c(-3000, 3000), TxDb=txdb)
peakAnno.list <- list(peakAnno.cKO, peakAnno.FF)

files <- list(WT="Tox_FF_peaks.filt.broadPeak", KO="Tox_cKO_peaks.filt.broadPeak")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

pdf("./../Results/ATACseq/Whl_PCF_TSS_Merged.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

pdf("./../Results/ATACseq/Whl_PCF_TSS_Split.pdf")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
dev.off()

pdf("./../Results/ATACseq/Whl_PCF_GeneBody_Merged.pdf")
plotPeakProf2(files, upstream = rel(0.2), downstream = rel(0.2), conf = 0.95, by = "gene", type = "body", TxDb = txdb, nbin = 800)
dev.off()

pdf("./../Results/ATACseq/Whl_PCF_GeneBody_Split.pdf")
plotPeakProf2(files, upstream = rel(0.2), downstream = rel(0.2), conf = 0.95, by = "gene", type = "body", TxDb = txdb, facet = "row", nbin = 800)
dev.off()

pdf("./../Results/ATACseq/Whl_Peak_TSS_Heatmap.pdf")
tagHeatmap(tagMatrixList)
dev.off()


peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList.tmp <- peakAnnoList
for(a in c("WT", "KO")){
  peakAnnoList.tmp[[a]]@annoStat$Feature <- as.character(peakAnnoList.tmp[[a]]@annoStat$Feature)
  for(term in c("Promoter", "Exon", "Intron")){
    termsum.tmp <- sum(as.numeric(peakAnnoList.tmp[[a]]@annoStat$Frequency[grepl(term, peakAnnoList.tmp[[a]]@annoStat$Feature)]))
    peakAnnoList.tmp[[a]]@annoStat <- peakAnnoList.tmp[[a]]@annoStat[!grepl(term, peakAnnoList.tmp[[a]]@annoStat$Feature),]
    peakAnnoList.tmp[[a]]@annoStat <- rbind(peakAnnoList.tmp[[a]]@annoStat, c(term, termsum.tmp))
  }
}
peakAnnoList.tmp[["WT"]]@annoStat$Sample <- "WT"
peakAnnoList.tmp[["KO"]]@annoStat$Sample <- "KO"
peakAnnot <- rbind(peakAnnoList.tmp[["WT"]]@annoStat, peakAnnoList.tmp[["KO"]]@annoStat)
peakAnnot$Feature <- as.factor(peakAnnot$Feature)
peakAnnot$Frequency <- as.numeric(peakAnnot$Frequency)
peakAnnot$Sample <- as.factor(peakAnnot$Sample)

pdf("./../Results/ATACseq/Whl_PeakStat_Stats_Merged.pdf")
ggplot(peakAnnot, aes(x=Frequency, y=Sample, fill=Feature)) + geom_bar(stat="identity") + scale_y_discrete(limits = c("WT", "KO")) + theme_bw()
dev.off()

pdf("./../Results/ATACseq/Whl_PeakStat_Stats.pdf")
plotAnnoBar(peakAnnoList)
dev.off()

p <- list()
p[[1]] <- ggplot(peakAnnoList$WT@annoStat, aes(x='', y=Frequency, fill=Feature))+
  geom_bar(stat='identity')+
  theme_void()+
  coord_polar('y', start=0) +
	ggtitle("WT")
p[[2]] <- ggplot(peakAnnoList$KO@annoStat, aes(x='', y=Frequency, fill=Feature))+
  geom_bar(stat='identity')+
  theme_void()+
  coord_polar('y', start=0) +
	ggtitle("KO")
  #geom_text(aes(label=paste0(round(Frequency,1), '%')), position=position_stack(vjust=0.5))
pdf("./../Results/ATACseq/Whl_PeakStat_Stats_PieChart.pdf", width=14)
print(grid.arrange(grobs = p, ncol=2))
dev.off()

p <- list()
p[[1]] <- ggplot(peakAnnoList$WT@annoStat, aes(x='', y=Frequency, fill=Feature))+
  geom_bar(stat='identity')+
  theme_void()+
  coord_polar('y', start=0) +
	ggtitle("WT") +
	theme(text = element_text(size=0))
p[[2]] <- ggplot(peakAnnoList$KO@annoStat, aes(x='', y=Frequency, fill=Feature))+
  geom_bar(stat='identity')+
  theme_void()+
  coord_polar('y', start=0) +
	ggtitle("KO") +
	theme(text = element_text(size=0))
  #geom_text(aes(label=paste0(round(Frequency,1), '%')), position=position_stack(vjust=0.5))
pdf("./../Results/ATACseq/Whl_PeakStat_Stats_PieChart_nofont.pdf", width=14)
print(grid.arrange(grobs = p, ncol=2))
dev.off()


tmp <- cbind(peakAnnoList$WT@annoStat[,2], peakAnnoList$KO@annoStat[,2])
colnames(tmp) <- c("WT", "KO")
row.names(tmp) <- peakAnnoList$WT@annoStat[,1]


Mosaic <- vcd::structable(t(tmp))
#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
pdf("./../Results/ATACseq/Whl_PeakStat_Stats_MosaicPlot.pdf", width=8)
vcd::mosaic(Mosaic, 
  shade=TRUE, 
  legend=TRUE, 
  dir = c("v", "h"), 
  labeling_args = list(rot_labels = c(top = 90, left = 0), 
    offset_varnames = c(top = 10, left=10), 
    offset_labels = c(top = 0),
    just_labels = c("left", "center", "center", "right")),
  margins = c(top = 8),
  )
dev.off()



  

pdf("./../Results/ATACseq/Whl_PeakStat_Distribution.pdf")
plotDistToTSS(peakAnnoList)
dev.off()

peakAnnoList.id <- lapply(files, annotatePeak, TxDb=txdb.id, tssRegion=c(-3000, 3000), verbose=FALSE)
saveRDS(peakAnnoList.id, "./../Results/ATACseq/Whole_peakAnnoList.rds")
genes = lapply(peakAnnoList.id, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
saveRDS(genes, "./../Results/ATACseq/Whole_genes.rds")

# compGO <- compareCluster(geneCluster = genes, fun = "enrichGO", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, pAdjustMethod = "BH", OrgDb = "org.Mm.eg.db", ont = "BP",)
# tmp <- compGO@compareClusterResult[,c(1,3,10)]
# tmp.sig <- tmp[tmp$p.adjust < 0.05,]
# GO.diff <- names(which(table(tmp.sig$Description) == 1))
# compGO@compareClusterResult <- compGO@compareClusterResult[compGO@compareClusterResult$Description %in% GO.diff,]
# compGO@compareClusterResult$p.adjust <- -log(compGO@compareClusterResult$p.adjust)
# pdf("./../Results/ATACseq/Whl_CompKEGG_Dotplot.pdf", width=6, height=10)
# dotplot(compGO, showCategory = nrow(compGO)/2, title = "GO Pathway Enrichment Analysis") + scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = -log(0.05)) + guides(colour=guide_legend(title="-log(p.adjust)"))
# dev.off()

dotplot(compGO, showCategory = 15, title = "GO Pathway Enrichment Analysis")

compKEGG <- compareCluster(geneCluster = genes, fun = "enrichKEGG", pvalueCutoff  = 1, qvalueCutoff  = 1, pAdjustMethod = "BH", organism="mmu")
tmp <- compKEGG@compareClusterResult[,c(1,5,9)]
tmp.sig <- tmp[tmp$p.adjust < 0.05,]
KEGG.diff <- names(which(table(tmp.sig$Description) == 1))
compKEGG@compareClusterResult <- compKEGG@compareClusterResult[compKEGG@compareClusterResult$Description %in% KEGG.diff,]
compKEGG@compareClusterResult$Description <- gsub(" - Mus musculus ", "", gsub("\\(.*?\\)", "", compKEGG@compareClusterResult$Description))
compKEGG@compareClusterResult$p.adjust <- -log(compKEGG@compareClusterResult$p.adjust)
saveRDS(compKEGG, "./../Results/ATACseq/Whole_CompKEGG.rds")
pdf("./../Results/ATACseq/Whl_CompKEGG_Dotplot.pdf", width=6, height=10)
dotplot(compKEGG, showCategory = nrow(compKEGG)/2, title = "KEGG Pathway Enrichment Analysis") + scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = -log(0.05)) + guides(colour=guide_legend(title="-log(p.adjust)"))
dev.off()

pdf("./../Results/ATACseq/Whl_CompKEGG_Dotplot_nofont.pdf", width=6, height=10)
dotplot(compKEGG, showCategory = nrow(compKEGG)/2, title = "KEGG Pathway Enrichment Analysis", font.size=0) + scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = -log(0.05)) + guides(colour=guide_legend(title="-log(p.adjust)")) + theme(text = element_text(size=0), axis.ticks.x=element_blank())
dev.off()

compPathway <- compareCluster(geneCluster = genes, fun = "enrichPathway", pvalueCutoff  = 1, qvalueCutoff  = 1, pAdjustMethod = "BH", organism="mouse")
tmp <- compPathway@compareClusterResult[,c(1,3,10)]
tmp.sig <- tmp[tmp$p.adjust < 0.05,]
Pathway.diff <- names(which(table(tmp.sig$Description) == 1))
compPathway@compareClusterResult <- compPathway@compareClusterResult[compPathway@compareClusterResult$Description %in% Pathway.diff,]
compPathway@compareClusterResult$p.adjust <- -log(compPathway@compareClusterResult$p.adjust)
pdf("./../Results/ATACseq/Whl_CompPathway_Dotplot.pdf", width=6, height=15)
dotplot(compPathway, showCategory = nrow(compPathway)/2, title = "Reactome Pathway Enrichment Analysis") + scale_fill_gradient2(low = "blue", mid = "gray", high = "red", midpoint = -log(0.05)) + guides(colour=guide_legend(title="-log(p.adjust)"))
dev.off()


genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
pdf("./../Results/ATACseq/Whl_PeakStat_Venndiagram.pdf")
vennplot(genes)
dev.off()


feat <- c("Cd4","Foxp3","Gata3","Nfil3","Batf","Lef1","Tcf7","Id2","Il2ra","Tox","Pdcd1","Havcr2","Tigit","Ctla4","Lag3","Tnfrsf4","Tnfrsf9","Tnfrsf18","Entpd1","Il1rl1","Ccr8","Cd44","Il10","Ccr7","Il7r","Sell")

Genes.cKO <- table(peakAnno.cKO@anno$geneId)
Genes.FF <- table(peakAnno.FF@anno$geneId)

Genes.df <- as.data.frame(cbind(Gene= feat, cKO=Genes.cKO[feat], FF=Genes.FF[feat]))
Genes.df$cKO <- as.numeric(Genes.df$cKO)
Genes.df$FF <- as.numeric(Genes.df$FF)

pdf("./../Results/ATACseq/Whl_PeakStat_ScatterPlot.pdf")
ggplot(Genes.df, aes(x=FF, y=cKO, label=Gene)) +
 geom_point() +
 geom_text_repel(hjust=0, vjust=0) +
 geom_abline(slope=1, intercept = 0) +
 ggtitle("Number of ATACseq peaks") + 
 theme_bw()
dev.off()

WT.High <- c('Foxp3','Gata3','Nfil3','Batf','Il2ra','Tox','Pdcd1','Havcr2','Tigit','Ctla4','Lag3','Tnfrsf4','Tnfrsf9','Tnfrsf18','Entpd1','Il1rl1','Ccr8','Cd44','Il10','Il7r','Icos','Cd27','Cd28','Layn','Mageh1','Gzmb','Cd83')
KO.High <- c('Lef1','Tcf7','Ccr4','Ccr7','Sell','Anxa1','Anxa2','Cxcr4','Il17a','Cd52')

# tmp <- read.delim("mm10.refGene.gtf")
# genes.gtf <- unique(gsub("gene_id ", "", str_split_fixed(tmp[,9], "; ", 2)[,1]))

peak.WT.WT.High <- peakAnnoList$WT@anno$V4[peakAnnoList$WT@anno$geneId %in% WT.High]
peak.KO.WT.High <- peakAnnoList$KO@anno$V4[peakAnnoList$KO@anno$geneId %in% WT.High]
peak.WT.KO.High <- peakAnnoList$WT@anno$V4[peakAnnoList$WT@anno$geneId %in% KO.High]
peak.KO.KO.High <- peakAnnoList$KO@anno$V4[peakAnnoList$KO@anno$geneId %in% KO.High]

WT.BP <- read.delim("Tox_FF_peaks.filt.broadPeak", header=F)
write.table(WT.BP[WT.BP$V4 %in% peak.WT.WT.High,], "Tox_FF_peaks.WTHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)
write.table(WT.BP[WT.BP$V4 %in% peak.WT.KO.High,], "Tox_FF_peaks.KOHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)
KO.BP <- read.delim("Tox_cKO_peaks.filt.broadPeak", header=F)
write.table(KO.BP[KO.BP$V4 %in% peak.KO.WT.High,], "Tox_cKO_peaks.WTHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)
write.table(KO.BP[KO.BP$V4 %in% peak.KO.KO.High,], "Tox_cKO_peaks.KOHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)

tagMatrixList.Int <- lapply(list(WT_WTHigh="Tox_FF_peaks.WTHigh.broadPeak", WT_KOHigh="Tox_FF_peaks.KOHigh.broadPeak", KO_WTHigh="Tox_cKO_peaks.WTHigh.broadPeak", KO_KOHigh="Tox_cKO_peaks.KOHigh.broadPeak"), getTagMatrix, windows=promoter)

pdf("./../Results/ATACseq/Int_PCF_TSS_Merged.pdf")
plotAvgProf(tagMatrixList.Int, xlim=c(-3000, 3000))
dev.off()

pdf("./../Results/ATACseq/Int_PCF_GeneBody_Merged.pdf")
plotPeakProf2(list(WT_WTHigh="Tox_FF_peaks.WTHigh.broadPeak", WT_KOHigh="Tox_FF_peaks.KOHigh.broadPeak", KO_WTHigh="Tox_cKO_peaks.WTHigh.broadPeak", KO_KOHigh="Tox_cKO_peaks.KOHigh.broadPeak"), upstream = rel(0.2), downstream = rel(0.2), conf = 0.95, by = "gene", type = "body", TxDb = txdb, nbin = 800)
dev.off()


peakAnnoList.Int <- lapply(list(WT_WTHigh="Tox_FF_peaks.WTHigh.broadPeak", WT_KOHigh="Tox_FF_peaks.KOHigh.broadPeak", KO_WTHigh="Tox_cKO_peaks.WTHigh.broadPeak", KO_KOHigh="Tox_cKO_peaks.KOHigh.broadPeak"), annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("./../Results/ATACseq/Int_PeakStat_Stats.pdf")
plotAnnoBar(peakAnnoList.Int)
dev.off()

pdf("./../Results/ATACseq/Int_PeakStat_Distribution.pdf")
plotDistToTSS(peakAnnoList.Int)
dev.off()



DEGs <- readRDS("./../Results/Seurat/Treg_DEG_1vs0.rds")
WT.High <- row.names(DEGs)[DEGs$avg_log2FC < 0 & DEGs$p_val_adj < 0.001]
KO.High <- row.names(DEGs)[DEGs$avg_log2FC > 0 & DEGs$p_val_adj < 0.001]


peak.WT.WT.High <- peakAnnoList$WT@anno$V4[peakAnnoList$WT@anno$geneId %in% WT.High]
peak.KO.WT.High <- peakAnnoList$KO@anno$V4[peakAnnoList$KO@anno$geneId %in% WT.High]
peak.WT.KO.High <- peakAnnoList$WT@anno$V4[peakAnnoList$WT@anno$geneId %in% KO.High]
peak.KO.KO.High <- peakAnnoList$KO@anno$V4[peakAnnoList$KO@anno$geneId %in% KO.High]

WT.BP <- read.delim("Tox_FF_peaks.filt.broadPeak", header=F)
write.table(WT.BP[WT.BP$V4 %in% peak.WT.WT.High,], "Tox_FF_peaks.WTHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)
write.table(WT.BP[WT.BP$V4 %in% peak.WT.KO.High,], "Tox_FF_peaks.KOHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)
KO.BP <- read.delim("Tox_cKO_peaks.filt.broadPeak", header=F)
write.table(KO.BP[KO.BP$V4 %in% peak.KO.WT.High,], "Tox_cKO_peaks.WTHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)
write.table(KO.BP[KO.BP$V4 %in% peak.KO.KO.High,], "Tox_cKO_peaks.KOHigh.broadPeak", sep='\t', quote=F, col.names=F, row.names=F)

tagMatrixList.Int <- lapply(list(WT_WTHigh="Tox_FF_peaks.WTHigh.broadPeak", WT_KOHigh="Tox_FF_peaks.KOHigh.broadPeak", KO_WTHigh="Tox_cKO_peaks.WTHigh.broadPeak", KO_KOHigh="Tox_cKO_peaks.KOHigh.broadPeak"), getTagMatrix, windows=promoter)

pdf("./../Results/ATACseq/Int.scRNAseq_PCF_TSS_Merged.pdf")
plotAvgProf(tagMatrixList.Int, xlim=c(-3000, 3000))
dev.off()

pdf("./../Results/ATACseq/Int.scRNAseq_PCF_GeneBody_Merged.pdf")
plotPeakProf2(list(WT_WTHigh="Tox_FF_peaks.WTHigh.broadPeak", WT_KOHigh="Tox_FF_peaks.KOHigh.broadPeak", KO_WTHigh="Tox_cKO_peaks.WTHigh.broadPeak", KO_KOHigh="Tox_cKO_peaks.KOHigh.broadPeak"), upstream = rel(0.2), downstream = rel(0.2), conf = 0.95, by = "gene", type = "body", TxDb = txdb, nbin = 800)
dev.off()


peakAnnoList.Int <- lapply(list(WT_WTHigh="Tox_FF_peaks.WTHigh.broadPeak", WT_KOHigh="Tox_FF_peaks.KOHigh.broadPeak", KO_WTHigh="Tox_cKO_peaks.WTHigh.broadPeak", KO_KOHigh="Tox_cKO_peaks.KOHigh.broadPeak"), annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

pdf("./../Results/ATACseq/Int.scRNAseq_PeakStat_Stats.pdf")
plotAnnoBar(peakAnnoList.Int)
dev.off()

pdf("./../Results/ATACseq/Int.scRNAseq_PeakStat_Distribution.pdf")
plotDistToTSS(peakAnnoList.Int)
dev.off()





# The function importScore is used to import BED, WIG, bedGraph or BigWig files. The function importBam is employed to import the bam files. Here is an example.

library(trackViewer)
WT <- importBam("Tox_FF.fixed.bam")
KO <- importBam("Tox_cKO_sub.fixed.bam")

# The function coverageGR could be used to calculate the coverage after the data is imported.

# dat.WT <- coverageGR(WT.raw$dat)
# dat.KO <- coverageGR(KO.raw$dat)

## We can split the data by strand into two different track channels
## Here, we set the dat2 slot to save the negative strand info. 
 
# WT$dat <- dat.WT[strand(dat.WT)=="+"]
# WT$dat2 <- dat.WT[strand(dat.WT)=="-"]

# KO$dat <- dat.KO[strand(dat.KO)=="+"]
# KO$dat2 <- dat.KO[strand(dat.KO)=="-"]

#Step 2. Build the gene model
#The gene model can be built for a given genomic range using geneModelFromTxdb function which uses the TranscriptDb object as the input.

#gr <- GRanges("chr11", IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))
# trs <- geneModelFromTxdb(TxDb.Mmusculus.UCSC.mm10.knownGene,org.Hs.eg.db, gr=gr)

# Users can generate a track object with the geneTrack function by inputting a TxDb and a list of gene Entrez IDs. Entrez IDs can be obtained from other types of gene IDs such as gene symbol by using the ID mapping function. For example, to generate a track object given gene FMR1 and human TxDb, refer to the code below.
WT.High <- c('Foxp3','Gata3','Nfil3','Batf','Il2ra','Tox','Pdcd1','Havcr2','Tigit','Ctla4','Lag3','Tnfrsf4','Tnfrsf9','Tnfrsf18','Entpd1','Il1rl1','Ccr8','Cd44','Il10','Il7r','Icos','Cd27','Cd28','Layn','Mageh1','Gzmb','Cd83')
KO.High <- c('Lef1','Tcf7','Ccr4','Ccr7','Sell','Anxa1','Anxa2','Cxcr4','Il17a','Cd52')

for(gene in WT.High){
	pdf(paste0("./../Results/ATACseq/Int.WTHigh_trackViewer_", gene, ".pdf"), height=3)
	theTrack <- geneTrack(gene,txdb)[[1]]

	gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

	viewerStyle <- trackViewerStyle()
	setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
	trackList <- trackList(WT, KO, theTrack)
	names(trackList) <- c("WT", "KO", gene)
	print(viewTracks(trackList, gr=gr, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE))
	dev.off()
}

for(gene in KO.High){
	pdf(paste0("./../Results/ATACseq/Int.KOHigh_trackViewer_",gene, ".pdf"), height=3)
	theTrack <- geneTrack(gene,txdb)[[1]]

	gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

	viewerStyle <- trackViewerStyle()
	setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
	trackList <- trackList(WT, KO, theTrack)
	names(trackList) <- c("WT", "KO", gene)
	print(viewTracks(trackList, gr=gr, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE))
	dev.off()
}




gene <- "Foxp3"
theTrack <- geneTrack(gene,txdb)[[1]]

gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
trackList <- trackList(KO, WT, theTrack)
names(trackList) <- c("KO", "WT", gene)

setTrackStyleParam(trackList[[1]], "color", c("red", "black"))
setTrackStyleParam(trackList[[2]], "color", c("black", "black"))

setTrackStyleParam(trackList[[1]], "ylim", c(0, 82))
setTrackStyleParam(trackList[[2]], "ylim", c(0, 82))

pdf(paste0("./../Results/ATACseq/Int.WTHigh_trackViewer_", gene, ".pdf"), height=3)
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE)
addGuideLine(c(7579204, 7579623), vp=vp)
addGuideLine(c(7581695, 7582035), vp=vp)
addGuideLine(c(7583960, 7584385), vp=vp)
addGuideLine(c(7586562, 7586795), vp=vp)
dev.off()



gene <- "Foxp3"
theTrack <- geneTrack(gene,txdb)[[1]]

gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
trackList <- trackList(KO, WT, theTrack)
names(trackList) <- c("KO", "WT", gene)

setTrackStyleParam(trackList[[1]], "color", c("red", "black"))
setTrackStyleParam(trackList[[2]], "color", c("black", "black"))

setTrackStyleParam(trackList[[1]], "ylim", c(0, 82))
setTrackStyleParam(trackList[[2]], "ylim", c(0, 82))

setTrackYaxisParam(trackList[[1]], "gp", list(cex=0))
setTrackYaxisParam(trackList[[2]], "gp", list(cex=0))
setTrackYaxisParam(trackList[[3]], "gp", list(cex=0))

trackList[[1]]@style@ylabgp <- list(cex=0)
trackList[[2]]@style@ylabgp <- list(cex=0)
trackList[[3]]@style@ylabgp <- list(cex=0)


pdf(paste0("./../Results/ATACseq/Int.WTHigh_trackViewer_", gene, "_noFont.pdf"), height=3)
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
addGuideLine(c(7579204, 7579623), vp=vp)
addGuideLine(c(7581695, 7582035), vp=vp)
addGuideLine(c(7583960, 7584385), vp=vp)
addGuideLine(c(7586562, 7586795), vp=vp)
dev.off()




KO.High <- c("Lef1","Tcf7","Ccr7","Sell","Il17a","Cd52")
max.KO.High <- list(122,72,66,52,31,65)
names(max.KO.High) <- KO.High

for(gene in KO.High){
	theTrack <- geneTrack(gene,txdb)[[1]]
	
	gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

	viewerStyle <- trackViewerStyle()
	setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
	trackList <- trackList(KO, WT, theTrack)
	names(trackList) <- c("KO", "WT", gene)

	setTrackStyleParam(trackList[[1]], "color", c("red", "black"))
	setTrackStyleParam(trackList[[2]], "color", c("black", "black"))

	setTrackStyleParam(trackList[[1]], "ylim", c(0, max.KO.High[[gene]]))
	setTrackStyleParam(trackList[[2]], "ylim", c(0, max.KO.High[[gene]]))

	pdf(paste0("./../Results/ATACseq/Int.KOHigh_trackViewer_", gene, ".pdf"), height=3)
	vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE)
	dev.off()
}

WT.High <- c("Gata3","Nfil3","Batf","Pdcd1","Havcr2","Lag3","Tnfrsf4","Il1rl1","Ccr8","Il7ra","Icos","Cd27")
max.WT.High <- list(96,63,82,203,39,77,120,35,86,74,71,59)
names(max.WT.High) <- WT.High


WT.High <- c("Il7r","Icos","Cd27")
max.WT.High <- list(74,71,59)
names(max.WT.High) <- WT.High

for(gene in WT.High){
	theTrack <- geneTrack(gene,txdb)[[1]]
	
	gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

	viewerStyle <- trackViewerStyle()
	setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
	trackList <- trackList(KO, WT, theTrack)
	names(trackList) <- c("KO", "WT", gene)

	setTrackStyleParam(trackList[[1]], "color", c("red", "black"))
	setTrackStyleParam(trackList[[2]], "color", c("black", "black"))

	setTrackStyleParam(trackList[[1]], "ylim", c(0, max.WT.High[[gene]]))
	setTrackStyleParam(trackList[[2]], "ylim", c(0, max.WT.High[[gene]]))

	pdf(paste0("./../Results/ATACseq/Int.WTHigh_trackViewer_", gene, ".pdf"), height=3)
	vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle, autoOptimizeStyle=TRUE)
	dev.off()
}

final <-c("Tcf7","Ccr7","Ccr8")
max.final <- list(72,66,86)

for(gene in final){
	theTrack <- geneTrack(gene,txdb)[[1]]

	gr <- GRanges(unique(seqnames(theTrack@dat)), IRanges(min(start(ranges(theTrack@dat))) - 1000, max(end(ranges(theTrack@dat)))) + 10000, strand=unique(strand(theTrack@dat)))

	viewerStyle <- trackViewerStyle()
	setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
	setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
	trackList <- trackList(KO, WT, theTrack)
	names(trackList) <- c("KO", "WT", gene)

	setTrackStyleParam(trackList[[1]], "color", c("red", "black"))
	setTrackStyleParam(trackList[[2]], "color", c("black", "black"))

	setTrackStyleParam(trackList[[1]], "ylim", c(0, max.final[[gene]]))
	setTrackStyleParam(trackList[[2]], "ylim", c(0, max.final[[gene]]))

	setTrackYaxisParam(trackList[[1]], "gp", list(cex=0))
	setTrackYaxisParam(trackList[[2]], "gp", list(cex=0))
	setTrackYaxisParam(trackList[[3]], "gp", list(cex=0))

	trackList[[1]]@style@ylabgp <- list(cex=0)
	trackList[[2]]@style@ylabgp <- list(cex=0)
	trackList[[3]]@style@ylabgp <- list(cex=0)

	pdf(paste0("./../Results/ATACseq/Int.final_trackViewer_", gene, ".pdf"), height=3)
	vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
	dev.off()
}




getYlim <- function(tl, op){
    yscales <- lapply(tl, function(.ele){
        ylim <- .ele@style@ylim
        if(length(ylim)!=2){
            if(.ele@type %in% c("data", "lollipopData")){
                if(length(.ele@dat)>0){
                    ylim <- unique(round(range(.ele@dat$score)))
                }else{
                    ylim <- c(0, 0)
                }
                if(length(.ele@dat2)>0 && is.null(op)){
                    ylim2 <- unique(round(range(.ele@dat2$score)))
                    ylim <- c(ylim, -1*ylim2)
                }
                ylim <- range(c(0, ylim))
            }else{
              if(.ele@type == "interactionData"){
                ## max interaction height
                ylim <- c(0, 1)
              }else{
                ylim <- c(0, 0)
              }
            }            
        }
        ylim
    })
    
    yscaleR <- range(unlist(yscales))
    if(diff(yscaleR)==0) yscaleR <- c(0, 1)
    yscales <- lapply(yscales, function(.ele){
        if(diff(.ele)==0){
            if(all(.ele>=0)){
                .ele <- c(0, yscaleR[2])
            }else{
                .ele <- c(yscaleR[1], 0)
            }
        }
        if(.ele[1]>.ele[2]) .ele[1] <- 0
        .ele
    })
    names(yscales) <- names(tl)
    yscales
}

getYheight <- function(tl){
    yHeights <- sapply(tl, function(.ele){
        yh <- .ele@style@height
        if(length(yh)==0) yh <- -1
        yh[1]
    })
    noY <- yHeights == -1
    yHeightsT <- sum(yHeights[!noY])
    if(yHeightsT>1)
        stop("total heights of data tracks is greater than 1.")
    if(length(yHeights[noY]) > 0){
        yHeights[noY] <- 
            (1 - yHeightsT) / length(yHeights[noY])
    }
    names(yHeights) <- names(tl)
    yHeights
}

getYheight(trackList)
getYlim(trackList, "+")

