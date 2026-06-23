# --- explicit package loading (added for standalone reproducibility; analysis logic UNCHANGED;
#     original scripts assumed an interactive session with these already attached) ---
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggvenn))

# =============================================================================
# scRNA-seq + scVDJ-seq analysis of TI Tregs (TOXWT vs TOXKO)
# Manuscript: "TOX enforces the immunosuppressive program of tumor-infiltrating
#              regulatory T cells" (Park et al., Nature Immunology)
# Figures: Fig 6; Extended Data 7-8
# Environment: scRNA env (R 4.0.5)   (exact versions in ../envs/ and README.md)
# Provenance: original project script "Analyzing_scRNAseq.R" — analysis logic UNCHANGED; only this
#             header added. Configure input/output paths per README directory layout.
# =============================================================================

#!/usr/bin/Rscript
#usage : ./plotting_tsne.R [matrix file] [output prefix]

library(Seurat)
library(dplyr)

#Multiprocess
library(future)
options(future.globals.maxSize = 10000 * 1024^2)
plan("multicore", workers = 20)
plan("sequential")

################################################################################
# Analysis starts!
################################################################################

# Loading object.
Data.raw <- Read10X("./../Data/scRNAseq/HN00196315/HN00199157_result_10X/YFPTIL-Treg_GEX/filtered_feature_bc_matrix/")
Data.GEX <- Data.raw[[1]]
Data.AB <- Data.raw[[2]]
# Filtering wanted cells.



# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(Data.GEX), colnames(Data.AB))

# Subset RNA and HTO counts by joint cell barcodes
Data.GEX <- Data.GEX[, joint.bcs]
Data.AB <- as.matrix(Data.AB[, joint.bcs])

# Setup Seurat object
SeuratObj <- CreateSeuratObject(counts = Data.GEX)

# Normalize RNA data with log normalization
SeuratObj <- NormalizeData(SeuratObj)
# Find and scale variable features
SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "mean.var.plot")
SeuratObj <- ScaleData(SeuratObj, features = VariableFeatures(SeuratObj))

# Add HTO data as a new assay independent from RNA
SeuratObj[["HTO"]] <- CreateAssayObject(counts = Data.AB)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
SeuratObj <- NormalizeData(SeuratObj, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using
# the default settings
SeuratObj <- HTODemux(SeuratObj, assay = "HTO", positive.quantile = 0.95)

# Group cells based on the max HTO signal
Idents(SeuratObj) <- "HTO_maxID"
pdf("./../Results/Seurat/HTO_RidgePlot.pdf")
RidgePlot(SeuratObj, assay = "HTO", features = rownames(SeuratObj[["HTO"]])[1:2], ncol = 2)
dev.off()

pdf("./../Results/Seurat/HTO_FeatureScatter.pdf")
FeatureScatter(SeuratObj, feature1 = "C0301", feature2 = "C0302")
dev.off()

# Manual clustering...
SeuratObj.HTO.df <- as.data.frame(t(SeuratObj@assays$HTO@data))

tmp <- rep("NA", ncol(SeuratObj))
tmp[SeuratObj.HTO.df[,1] < 2 & SeuratObj.HTO.df[,2] < 0.2] <- "Negative"
tmp[SeuratObj.HTO.df[,1] > 2.2 & SeuratObj.HTO.df[,2] > 0.3] <- "Doublet"
tmp[SeuratObj.HTO.df[,1] > 2.2 & SeuratObj.HTO.df[,2] < 0.2] <- "C0301"
tmp[SeuratObj.HTO.df[,1] < 1.4 & SeuratObj.HTO.df[,2] > 0.3] <- "C0302"
#tmp[SeuratObj.HTO.df[,1] > 2 & SeuratObj.HTO.df[,1] < 0.7 & SeuratObj.HTO.df[,2] > 0.4 & SeuratObj.HTO.df[,2] < 0.6] <- "Doublet"
tmp <- as.factor(tmp)

SeuratObj <- AddMetaData(SeuratObj, metadata = tmp, col.name = "hash.ID")

Idents(SeuratObj) <- "hash.ID"
pdf("./../Results/Seurat/HTO_RidgePlot.pdf")
RidgePlot(SeuratObj, assay = "HTO", features = rownames(SeuratObj[["HTO"]])[1:2], ncol = 2)
dev.off()
pdf("./../Results/Seurat/HTO_FeatureScatter.pdf")
FeatureScatter(SeuratObj, feature1 = "C0301", feature2 = "C0302")
dev.off()
pdf("./../Results/Seurat/HTO_nCount.RNA.pdf")
VlnPlot(SeuratObj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# Extract the singlets
SeuratObj.singlet <- subset(SeuratObj, idents = c("C0301", "C0302"))

SeuratObj.singlet[["percent.mt"]] <- PercentageFeatureSet(SeuratObj.singlet, pattern = "^mt-")
pdf("./../Results/Seurat/Singlet_QC.pdf", width = 14)
par(mfrow = c(1,2))
plot(SeuratObj.singlet$nCount_RNA, SeuratObj.singlet$percent.mt)
abline(h = 5)
plot(SeuratObj.singlet$nCount_RNA, SeuratObj.singlet$nFeature_RNA)
abline(h = 200)
abline(h = 8000)
dev.off()

SeuratObj.singlet <- subset(x = SeuratObj.singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)


# Select the top 1000 most variable features
SeuratObj.singlet <- FindVariableFeatures(SeuratObj.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
SeuratObj.singlet <- ScaleData(SeuratObj.singlet, features = row.names(SeuratObj.singlet))

# Run PCA
SeuratObj.singlet <- RunPCA(SeuratObj.singlet, npcs=100, features = row.names(SeuratObj.singlet))

feat <- c("Cd4", "Cd8", "Tox", "Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Ifit3", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

for(dim in c(10,15,20,25,30,35,40,45,50)){
	print(dim)
	# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
	SeuratObj.singlet <- FindNeighbors(SeuratObj.singlet, reduction = "pca", dims = 1:dim)
	SeuratObj.singlet <- RunUMAP(SeuratObj.singlet, reduction = "pca", dims = 1:dim)

	png(paste0("./../Results/Seurat/Singlet_UMAP", ".", dim, "_Feature.png"), width = 6400, height = 4000, res = 144)
	print(FeaturePlot(SeuratObj.singlet, reduction = "umap", features = feat, ncol = 8))
	dev.off()

	pdf(paste0("./../Results/Seurat/Singlet_UMAP", ".", dim, "_Hash.ID.pdf"))
	print(DimPlot(SeuratObj.singlet, group.by = "hash.ID"))
	dev.off()
}


dim <- 40
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
SeuratObj.singlet <- FindNeighbors(SeuratObj.singlet, reduction = "pca", dims = 1:dim)
SeuratObj.singlet <- FindClusters(SeuratObj.singlet, resolution = seq(0.1, 1, 0.1), verbose = FALSE)
SeuratObj.singlet <- RunUMAP(SeuratObj.singlet, reduction = "pca", dims = 1:dim)

png("./../Results/Seurat/Singlet_UMAP_Feature.png", width = 6400, height = 4000, res = 144)
print(FeaturePlot(SeuratObj.singlet, reduction = "umap", features = feat, ncol = 8))
dev.off()

p <- list()
p[[1]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Cd4"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.1,0.2))
p[[2]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Foxp3"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.1,0.2))
p[[3]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Il2ra"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.1,0.2))
p[[4]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Tox"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.1,0.2))
png("./../Results/Seurat/Singlet_UMAP_Feature.png", width = 1600, height = 1600, res = 144)
print(grid.arrange(grobs = p, ncol = 2))
dev.off()

p <- list()
p[[1]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Cd4"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.8,0.85))
p[[2]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Foxp3"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.8,0.85))
p[[3]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Il2ra"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.8,0.85))
p[[4]] <- FeaturePlot(SeuratObj.singlet, reduction = "umap", features = c("Tox"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.8,0.85))
png("./../Results/Seurat/Singlet_UMAP_Feature.png", width = 1600, height = 1600, res = 144)
print(grid.arrange(grobs = p, ncol = 2))
dev.off()


pdf("./../Results/Seurat/Singlet_UMAP_Hash.ID.pdf")
print(DimPlot(SeuratObj.singlet, group.by = "hash.ID"))
dev.off()

pdf("./../Results/Seurat/Singlet_UMAP_Clusters.pdf")
for(j in seq(0.1, 1, 0.1)){
	print(DimPlot(SeuratObj.singlet, reduction = "umap", label = T, pt.size = 0.5, group.by = paste0("RNA_snn_res.", j)))
}
dev.off()

SeuratObj.singlet <- FindClusters(SeuratObj.singlet, resolution = 0.9)

pdf("./../Results/Seurat/Singlet_UMAP_Cluster.pdf")
DimPlot(SeuratObj.singlet, reduction = "umap", label = T, pt.size = 0.9)
dev.off()

pdf("./../Results/Seurat/Singlet_UMAP_Cluster_NoLegend.pdf")
DimPlot(SeuratObj.singlet, reduction = "umap", label = F, pt.size = 0.9) + NoLegend()
dev.off()

# Reading VDJ data.
add_clonotype <- function(seurat_obj, tcr_location){
	tcr <- read.csv(paste0(tcr_location, "/filtered_contig_annotations.csv"))
  
	# Remove the -1 at the end of each barcode.
	# Subsets so only the first line of each barcode is kept,
	# as each entry for given barcode will have same clonotype.
	#tcr$barcode <- gsub("-1", "", tcr$barcode)
	tcr <- tcr[!duplicated(tcr$barcode), ]
  
	# Only keep the barcode and clonotype columns. 
	# We'll get additional clonotype info from the clonotype table.
	tcr <- tcr[,c("barcode", "raw_clonotype_id")]
	names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
	# Clonotype-centric info.
	clono <- read.csv(paste0(tcr_location, "/clonotypes.csv"))
  
	# Slap the AA sequences onto our original table by clonotype_id.
	tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa", "cdr3s_nt")])
  
	# Reorder so barcodes are first column and set them as rownames.
	rownames(tcr) <- tcr$barcode
	tcr <- tcr[, -2]
  
	# Add to the Seurat object's metadata.
	clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
	return(clono_seurat)
}
SeuratObj.singlet <- add_clonotype(SeuratObj.singlet, "./../Data/scRNAseq/HN00196315/HN00199157_result_10X/YFPTIL-Treg_TCR/")
clonotype_id.tmp <- paste0(SeuratObj.singlet$hash.ID, "-", SeuratObj.singlet$clonotype_id)
clonotype_id.tmp[is.na(SeuratObj.singlet$clonotype_id)] <- NA
SeuratObj.singlet$clonotype_id <- clonotype_id.tmp

SeuratObj.singlet <- div_clonotype(SeuratObj.singlet)

pdf(paste0("./../Results/Seurat/Singlet_UMAP_clonotype.pdf"))
print(DimPlot(SeuratObj.singlet, group.by = "clonotype"))
dev.off()

png(paste0("./../Results/Seurat/Singlet_VlnPlot_Feature.png"), width = 3200, height = 1600, res = 144)
VlnPlot(SeuratObj.singlet, features = c("Cd4", "Cd8a", "Cd14", "Lyz", "Tox", "Pdcd1", "Foxp3", "Il2ra"), ncol=4)
dev.off()


saveRDS(SeuratObj.singlet, "./../Results/Seurat/Singlet_SeuratObj.rds")

########################
# Filtering only Tregs #
########################
Treg <- subset(SeuratObj.singlet, idents = c(1,2,4,6,7,8,9,10,12,14,16))

Treg <- Treg[,Treg$clonotype != "NA"]

Treg <- ScaleData(Treg, features = row.names(Treg))

# Run PCA
Treg <- RunPCA(Treg, npcs=100, features = row.names(Treg))

feat <- c("Cd4", "Cd8", "Tox", "Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Ifit3", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

dim <- 20

Treg <- FindNeighbors(Treg, reduction = "pca", dims = 1:dim)
Treg <- FindClusters(Treg, resolution = seq(0.1, 0.2, 0.01), verbose = FALSE)
Treg <- RunUMAP(Treg, reduction = "pca", dims = 1:dim)

div_clonotype <- function(seurat_obj){
	# Extracting clonotype information
	clonotable <- table(seurat_obj$clonotype_id)
	clonotable <- clonotable[clonotable >= 1]
	multiclone <- names(clonotable[clonotable > 1])
	singleclone <- names(clonotable[clonotable == 1])

	# Making object
	clonotype <- rep("NA", ncol(seurat_obj))
	clonotype[which(seurat_obj$clonotype_id %in% multiclone)] <- "multi-clonotype"
	clonotype[which(seurat_obj$clonotype_id %in% singleclone)] <- "single-clonotype"

	# Calculating the number of each clones
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotable[as.character(seurat_obj$clonotype_id)], col.name = "cloneNum")
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotype, col.name = "clonotype")
	return(seurat_obj)
}
Treg <- div_clonotype(Treg)

png("./../Results/Seurat/Treg_UMAP_Feature.png", width = 6400, height = 4000, res = 144)
print(FeaturePlot(Treg, reduction = "umap", features = feat, ncol = 8))
dev.off()

pdf("./../Results/Seurat/Treg_UMAP_Feature.Tox.pdf")
FeaturePlot(Treg, reduction = "umap", features = c("Tox"), ncol = 1, order=T) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20), legend.position = c(0.8,0.2))
dev.off()

umap.max <- apply(Treg@reductions$umap@cell.embeddings,2,max)
umap.min <- apply(Treg@reductions$umap@cell.embeddings,2,min)

Treg$hash.name <- Treg$hash.ID
levels(Treg$hash.name) <- c("YFP+ Treg", "YFP-CD25+ Treg", "Doublet", "NA", "Negative")
pdf("./../Results/Seurat/Treg_UMAP_Hash.Name.pdf")
Treg.301 <- Treg[,Treg$hash.name == "YFP+ Treg"]
print(DimPlot(Treg.301, group.by = "hash.name") + NoLegend() + ggtitle("YFP+ Treg") + xlim(umap.min[1], umap.max[1]) + ylim(umap.min[2], umap.max[2]) + theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20)))
Treg.302 <- Treg[,Treg$hash.name == "YFP-CD25+ Treg"]
print(DimPlot(Treg.302, group.by = "hash.name") + NoLegend() + ggtitle("YFP-CD25+ Treg") + xlim(umap.min[1], umap.max[1]) + ylim(umap.min[2], umap.max[2])+ theme(axis.text = element_text(size = 20), plot.title = element_text(size = 20)))
dev.off()

Treg.tmp <- Treg
Treg.tmp$clonotype <- as.factor(Treg.tmp$clonotype)
levels(Treg.tmp$clonotype) <- c("expanded-clonotype", "unique-clonotype")
pdf("./../Results/Seurat/Treg_UMAP_clonotype.pdf")
print(DimPlot(Treg.tmp, group.by = "clonotype", order = c("expanded-clonotype")) + scale_colour_manual(values=c("Gray", "Red")) + theme(legend.position = c(0.7,0.1)))
print(DimPlot(Treg.tmp, group.by = "clonotype", order = c("expanded-clonotype")) + scale_colour_manual(values=c("Gray", "Red")) + theme(legend.position = ""))
dev.off()

pdf("./../Results/Seurat/Treg_UMAP_Clusters.pdf")
for(j in seq(0.1, 0.2, 0.01)){
	print(DimPlot(Treg, reduction = "umap", label = T, pt.size = 0.5, group.by = paste0("RNA_snn_res.", j)))
}
dev.off()

Treg@active.ident <- Treg$RNA_snn_res.0.2

pdf("./../Results/Seurat/Treg_UMAP_Percent.MT.pdf")
print(FeaturePlot(Treg, reduction = "umap", features = "percent.mt"))
dev.off()

Treg$Cluster <- Treg@active.ident
levels(Treg$Cluster) <- c("0", "1", "2", "3", "4.5", "4.5", "6")

saveRDS(Treg, "./../Results/Seurat/Treg_SeuratObj.rds")


################################################################################
# DEG analysis
# Finding DEGs of clusters.
DEGs <- FindAllMarkers(Treg, logfc.threshold = 0, min.pct = 0)
DEGs$threshold <- DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > 0.25
write.table(DEGs, "./../Results/Seurat/Treg_DEG_Cluster.tsv", quote = F, sep = '\t', col.names=NA)
DEGs <- read.delim("./../Results/Seurat/Treg_DEG_Cluster.tsv", sep='\t', header=T, row.names=1)

PlottingVolcano <- function(WholeDEG, Output, feat){
	pdf(Output)
	WholeDEG$cluster <- as.factor(as.character(WholeDEG$cluster))
	for(cluster in levels(WholeDEG$cluster)){
		DEG.list <- WholeDEG[WholeDEG$cluster == cluster,]
		DEG.list$threshold <- DEG.list$p_val_adj < 0.05 & abs(DEG.list$avg_log2FC) > 0.25
		DEG.list$pct.AbsLogRatio <- abs(log2(DEG.list$pct.1 / DEG.list$pct.2))
		DEG.list <- DEG.list[order(DEG.list$avg_log2FC, decreasing = T),]

		DEG.list$genelabels <- "FALSE"
		SigMark <- DEG.list[DEG.list$threshold,]
		DEG.list$genelabels[row.names(DEG.list) %in% head(row.names(SigMark), 20)] <- "TRUE"
		DEG.list$genelabels[row.names(DEG.list) %in% tail(row.names(SigMark), 20)] <- "TRUE"
		DEG.list$genelabels[DEG.list$gene %in% feat & row.names(DEG.list) %in% row.names(SigMark)] <- "TRUE"


		gg <- ggplot(DEG.list) +
			geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = threshold)) +
			ggrepel::geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = ifelse(genelabels == T, as.character(gene),"")), size=5, max.overlaps=Inf, force=40, segment.color="grey") +
      ggtitle(paste0("VolcanoPlot of cluster ", cluster)) +
			xlab("Average log2 fold change") +
			ylab("-log10 adjusted p-value") +
			theme(legend.position = "none",
			plot.title = element_text(size = rel(1.5), hjust = 0.5),
			axis.title = element_text(size = rel(1.25)))
		print(gg)
	}
	dev.off()
}

PlottingVolcano.single <- function(WholeDEG, Output, feat){
	pdf(Output)
	DEG.list <- WholeDEG
	DEG.list$gene <- row.names(DEG.list)
	DEG.list$threshold <- DEG.list$p_val_adj < 0.05 & abs(DEG.list$avg_log2FC) > 0.25
	DEG.list$pct.AbsLogRatio <- abs(log2(DEG.list$pct.1 / DEG.list$pct.2))
	DEG.list <- DEG.list[order(DEG.list$avg_log2FC, decreasing = T),]

	DEG.list$genelabels <- "FALSE"
	SigMark <- DEG.list[DEG.list$threshold,]
	DEG.list$genelabels[row.names(DEG.list) %in% head(row.names(SigMark), 20)] <- "TRUE"
	DEG.list$genelabels[row.names(DEG.list) %in% tail(row.names(SigMark), 20)] <- "TRUE"
	DEG.list$genelabels[DEG.list$gene %in% feat & row.names(DEG.list) %in% row.names(SigMark)] <- "TRUE"


	gg <- ggplot(DEG.list) +
		geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = threshold)) +
		ggrepel::geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = ifelse(genelabels == T, as.character(gene),"")), max.overlaps=Inf, force=20, segment.color="grey") +
		xlab("Average log2 fold change") +
		ylab("-log10 adjusted p-value") +
		theme(legend.position = "none",
		plot.title = element_text(size = rel(1.5), hjust = 0.5),
		axis.title = element_text(size = rel(1.25)))
	print(gg)
	dev.off()
}

PlottingVolcano.special <- function(WholeDEG, Output, feat){
	pdf(Output)
	DEG.list <- WholeDEG
	DEG.list$gene <- row.names(DEG.list)
	DEG.list$threshold <- DEG.list$p_val_adj < 0.05 & abs(DEG.list$avg_log2FC) > 0.25
	DEG.list$pct.AbsLogRatio <- abs(log2(DEG.list$pct.1 / DEG.list$pct.2))
	DEG.list <- DEG.list[order(DEG.list$avg_log2FC, decreasing = T),]

	DEG.list$genelabels <- "FALSE"
	DEG.list$genelabels[DEG.list$gene %in% c("Ikzf2", "Pdcd1", "Ccr8", "Tnfrsf9", "Dgat2", "Tnfrsf4", "Tigit", "Tox", "Nr4a1", "Klrg1", "Id2", "Ccl5", "Batf", "Cd83", "Foxp3", "Areg", "Gata3", "Gzmb", "Havcr2", "Il1rl1", "Il2ra", "Tbx21", "Cd5", "Cd44", "Tmem176b", "Foxp1", "Saraf", "Tmem176a", "S1pr1", "Tbc1d4", "Tcf7", "B4galt1", "Ccr4", "Klf2", "Il10")] <- "TRUE"


	gg <- ggplot(DEG.list) +
		geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = threshold)) +
		ggrepel::geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = ifelse(genelabels == T, as.character(gene),"")), max.overlaps=Inf, force=20, segment.color="grey") +
		xlab("Average log2 fold change") +
		ylab("-log10 adjusted p-value") +
		theme(legend.position = "none",
		plot.title = element_text(size = rel(1.5), hjust = 0.5),
		axis.title = element_text(size = rel(1.25)))
	print(gg)
	dev.off()
}


feats <- c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1", "Tox")
PlottingVolcano(DEGs, "./../Results/Seurat/Treg_DEG_Volcano_Cluster.pdf", feats)
WholeDEG <- DEGs
WholeDEG$cluster <- as.factor(as.character(WholeDEG$cluster))
cluster <- "1"
DEG.list <- WholeDEG[WholeDEG$cluster == cluster,]
DEG.list$threshold <- DEG.list$p_val_adj < 0.05 & abs(DEG.list$avg_log2FC) > 0.25
DEG.list$pct.AbsLogRatio <- abs(log2(DEG.list$pct.1 / DEG.list$pct.2))
DEG.list <- DEG.list[order(DEG.list$avg_log2FC, decreasing = T),]

DEG.list$genelabels <- "FALSE"
SigMark <- DEG.list[DEG.list$threshold,]
DEG.list$genelabels[row.names(DEG.list) %in% head(row.names(SigMark), 20)] <- "TRUE"
DEG.list$genelabels[row.names(DEG.list) %in% tail(row.names(SigMark), 20)] <- "TRUE"
DEG.list$genelabels[DEG.list$gene %in% feat & row.names(DEG.list) %in% row.names(SigMark)] <- "TRUE"

pdf("./../Results/Seurat/Treg_DEG_Volcano_Cluster1_reverse.pdf")
gg <- ggplot(DEG.list) +
	geom_point(aes(x = -avg_log2FC, y = -log10(p_val_adj), colour = threshold)) +
	ggrepel::geom_text_repel(aes(x = -avg_log2FC, y = -log10(p_val_adj), label = ifelse(genelabels == T, as.character(gene),"")), size=5, max.overlaps=Inf, force=40, segment.color="grey") +
	ggtitle(paste0("VolcanoPlot of cluster ", cluster)) +
	xlab("Average log2 fold change") +
	ylab("-log10 adjusted p-value") +
	theme_bw() +
	theme(legend.position = "none",
	plot.title = element_text(size = rel(1.5), hjust = 0.5),
	axis.title = element_text(size = rel(1.25)))
print(gg)
dev.off()


DEGs.5vs4 <- FindMarkers(Treg, ident.1="5", ident.2="4", logfc.threshold = 0, min.pct = 0)
DEGs.1vs0 <- FindMarkers(Treg, ident.1="1", ident.2="0", logfc.threshold = 0, min.pct = 0)
PlottingVolcano.single(DEGs.5vs4, "./../Results/Seurat/Treg_DEG.5vs4_Volcano_Cluster.pdf", feats)
PlottingVolcano.special(DEGs.1vs0, "./../Results/Seurat/Treg_DEG.1vs0_Volcano_Cluster.pdf", feats)

# DEGs between hashtags
for(cl in levels(Treg@active.ident)){
	print(cl)
	Seurat.tmp <- Treg[,Treg@active.ident == cl]
	Seurat.tmp@active.ident <- Seurat.tmp$hash.ID

	DEGs <- FindMarkers(Seurat.tmp, ident.1="C0302", ident.2="C0301", logfc.threshold = 0, min.pct = 0)
	DEGs$threshold <- DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > 0.25
	write.table(DEGs, paste0("./../Results/Seurat/Treg_DEG.hashID_Cluster.", cl, ".tsv"), quote = F, sep = '\t', col.names=NA)
	PlottingVolcano.single(DEGs, paste0("./../Results/Seurat/Treg_DEG.hashID_Volcano_Cluster.", cl, ".pdf"), feats)
}

for(cl in levels(Treg@active.ident)){
	DEGs <- read.delim(paste0("./../Results/Seurat/Treg_DEG.hashID_Cluster.", cl, ".tsv"), row.names=1)
	PlottingVolcano.single(DEGs, paste0("./../Results/Seurat/Treg_DEG.hashID_Volcano_Cluster.", cl, ".pdf"), feats)
}


tmp <- Treg
tmp@active.ident <- tmp$hash.ID
DEGs <- FindMarkers(tmp, ident.1="C0302", ident.2="C0301", logfc.threshold = 0, min.pct = 0)
DEGs$threshold <- DEGs$p_val_adj < 0.05 & abs(DEGs$avg_log2FC) > 0.25
write.table(DEGs, "./../Results/Seurat/Treg_DEG.hashID_Whole.tsv", quote = F, sep = '\t', col.names=NA)
DEGs <- read.delim("./../Results/Seurat/Treg_DEG.hashID_Whole.tsv", row.names=1)
PlottingVolcano.single(DEGs, paste0("./../Results/Seurat/Treg_DEG.hashID_Volcano_Whole.pdf"), feats)



GS.MSigDB <- getGmt("./../../Treg_MJ/Data/Signature/MSigDB/Hallmark.symbols_mouse.gmt")

# GSVA of whole clusters
gsva.tmp <- gsva(as.matrix(Treg@assays$RNA@data), GS.MSigDB, parallel.sz = 28)
saveRDS(gsva.tmp, "./../Results/Seurat/Treg_Whole_GSVA_MSigDB.norm.rds")

p <- list()
i <- 1
volcano.df <- data.frame()
for(term in row.names(gsva.tmp)){
	term.df <- as.data.frame(cbind(GSVA = gsva.tmp[term, ], hash.ID = as.character(Treg$hash.ID)))
	term.df$GSVA <- as.numeric(as.character(term.df$GSVA))

	p[[i]] <- ggplot(term.df, aes(x = hash.ID, y = GSVA)) +
		geom_violin(adjust = 3) +
		geom_boxplot() +
		ggtitle(term) +
		xlab("hash.ID") +
		ylab("GSVA score") +
		stat_compare_means(method="wilcox.test")

	i <- i + 1
	volcano.df <- rbind(volcano.df, c(term, mean(term.df$GSVA[term.df$hash.ID == "C0302"]) - mean(term.df$GSVA[term.df$hash.ID == "C0301"]), wilcox.test(term.df$GSVA[term.df$hash.ID == "C0301"], term.df$GSVA[term.df$hash.ID == "C0302"])$p.value))
}
pdf("./../Results/Seurat/Treg_Whole_GSVA_MSigDB.norm_Violin.pdf", width = 36, height = 24)
print(grid.arrange(grobs = p, ncol = 9))
dev.off()

colnames(volcano.df) <- c("Term", "Diff", "pvalue")
volcano.df$Diff <- as.numeric(volcano.df$Diff)
volcano.df$pvalue <- as.numeric(volcano.df$pvalue)
volcano.df$Term <- str_split_fixed(volcano.df$Term, "[_]", 2)[,-1]

pdf("./../Results/Seurat/Treg_Whole_GSVA_MSigDB.norm_Volcano.pdf", width=11, height=11)
print(ggplot(volcano.df, aes(x=Diff, y=-log10(pvalue))) +
	geom_point() +
	xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
	geom_text_repel(aes(x = Diff, y = -log10(pvalue), label = ifelse(pvalue < 0.05, as.character(Term),"")), max.overlaps=Inf, force=20, segment.color="grey") +
	geom_hline(yintercept=-log10(0.05)))
dev.off()


for(cl in levels(Treg@active.ident)){
	print(cl)
	Seurat.tmp <- Treg[,Treg@active.ident == cl]
	gsva.tmp <- gsva(as.matrix(Seurat.tmp@assays$RNA@data), GS.MSigDB, parallel.sz = 28)
	saveRDS(gsva.tmp, paste0("./../Results/Seurat/Treg_Cluster.", cl, "_GSVA_MSigDB.norm.rds"))

	p <- list()
	i <- 1
	volcano.df <- data.frame()
	for(term in row.names(gsva.tmp)){
		term.df <- as.data.frame(cbind(GSVA = gsva.tmp[term, ], hash.ID = as.character(Seurat.tmp$hash.ID)))
		term.df$GSVA <- as.numeric(as.character(term.df$GSVA))

		p[[i]] <- ggplot(term.df, aes(x = hash.ID, y = GSVA)) +
			geom_violin(adjust = 3) +
			geom_boxplot() +
			ggtitle(term) +
			xlab("hash.ID") +
			ylab("GSVA score") +
			stat_compare_means(method="wilcox.test")

		i <- i + 1
		volcano.df <- rbind(volcano.df, c(term, mean(term.df$GSVA[term.df$hash.ID == "C0302"]) - mean(term.df$GSVA[term.df$hash.ID == "C0301"]), wilcox.test(term.df$GSVA[term.df$hash.ID == "C0301"], term.df$GSVA[term.df$hash.ID == "C0302"])$p.value))
	}
	pdf(paste0("./../Results/Seurat/Treg_Cluster.", cl, "_GSVA_MSigDB.norm_Violin.pdf"), width = 36, height = 24)
	print(grid.arrange(grobs = p, ncol = 9))
	dev.off()

	colnames(volcano.df) <- c("Term", "Diff", "pvalue")
	volcano.df$Diff <- as.numeric(volcano.df$Diff)
	volcano.df$pvalue <- as.numeric(volcano.df$pvalue)
	volcano.df$Term <- str_split_fixed(volcano.df$Term, "[_]", 2)[,-1]

	pdf(paste0("./../Results/Seurat/Treg_Cluster.", cl, "_GSVA_MSigDB.norm_Volcano.pdf"), width=11, height=11)
	print(ggplot(volcano.df, aes(x=Diff, y=-log10(pvalue))) +
		geom_point() +
		xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
		geom_text_repel(aes(x = Diff, y = -log10(pvalue), label = ifelse(pvalue < 0.05, as.character(Term),"")), max.overlaps=Inf, force=20, segment.color="grey") +
		geom_hline(yintercept=-log10(0.05)))
	dev.off()
}

# Mosaic plot between clonotype and hashID from each clusters
MosaicData <- Treg
MosaicData$clonotype <- as.factor(as.character(MosaicData$clonotype))
MetaInfo <- as.factor(as.character(MosaicData$hash.ID))
ClusterInfo <- MosaicData$clonotype
MosaicTable <- table(MetaInfo, ClusterInfo)

Mosaic <- vcd::structable(ClusterInfo ~ MetaInfo)
#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
pdf("./../Results/Seurat/Treg_Whole_MosaicPlot_clonotype.hashID.pdf", width=8)
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

for(cl in levels(Treg@active.ident)){
	Seurat.tmp <- Treg[,Treg@active.ident == cl]
	MosaicData <- Seurat.tmp
	MosaicData$clonotype <- as.factor(as.character(MosaicData$clonotype))
	MetaInfo <- as.factor(as.character(MosaicData$hash.ID))
	ClusterInfo <- MosaicData$clonotype
	MosaicTable <- table(MetaInfo, ClusterInfo)

	Mosaic <- vcd::structable(ClusterInfo ~ MetaInfo)
	#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
	pdf(paste0("./../Results/Seurat/Treg_Cluster.", cl, "_MosaicPlot_clonotype.hashID.pdf"), width=8)
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
}


# Mosaic plot between clonotype and clusters
MosaicData <- Treg
MosaicData$clonotype <- as.factor(as.character(MosaicData$clonotype))
MetaInfo <- as.factor(as.character(MosaicData@active.ident))
ClusterInfo <- MosaicData$clonotype
MosaicTable <- table(MetaInfo, ClusterInfo)

Mosaic <- vcd::structable(ClusterInfo ~ MetaInfo)
#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
pdf("./../Results/Seurat/Treg_Whole_MosaicPlot_clonotype.cluster.pdf", width=8)
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

# Clonotype barplot
Treg.tmp <- Treg[,Treg@active.ident %in% c(0,1)]
Treg.df <- as.data.frame(cbind(Group = as.character(Treg@active.ident), clonotype = Treg$clonotype))
Treg.tmp.df <- as.data.frame(cbind(Group = as.character(Treg.tmp@active.ident), clonotype = Treg.tmp$clonotype))

pdf("./../Results/Seurat/Treg_Whole_clonotype_Barplot.pdf")
ggplot(Treg.df, aes(Group, fill = clonotype)) +
  geom_bar(position = 'fill') + 
  #scale_x_discrete(labels = c('Iso PD-1+ Treg' = 'Iso\nPD-1+\nTreg', 'aPD-1 PD-1+ Treg' = 'aPD-1\nPD-1+\nTreg', 'Iso PD-1- Treg' = 'Iso\nPD-1-\nTreg', 'aPD-1 PD-1- Treg' = 'aPD-1\nPD-1-\nTreg')) +
  ylab("Proportion") +
  #guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme_cowplot() +
  #theme(text = element_blank()) +
  #theme(legend.position = "none") +
  scale_fill_manual(values=c('red','grey'))
dev.off()

pdf("./../Results/Seurat/Treg_Cluster.0.1_clonotype_Barplot.pdf")
ggplot(Treg.tmp.df, aes(Group, fill = clonotype)) +
  geom_bar(position = 'fill') + 
  #scale_x_discrete(labels = c('Iso PD-1+ Treg' = 'Iso\nPD-1+\nTreg', 'aPD-1 PD-1+ Treg' = 'aPD-1\nPD-1+\nTreg', 'Iso PD-1- Treg' = 'Iso\nPD-1-\nTreg', 'aPD-1 PD-1- Treg' = 'aPD-1\nPD-1-\nTreg')) +
  ylab("Proportion") +
  #guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme_cowplot() +
  #theme(text = element_blank()) +
  #theme(legend.position = "none") +
  scale_fill_manual(values=c('red','grey'))
dev.off()


# Gini Index

Gini.df <- data.frame()
for(ct in levels(Treg@active.ident)){
	Gini.df <- rbind(Gini.df, c(ct, ineq::ineq(table(Treg$clonotype_id[Treg@active.ident == ct]), type='Gini')))
}
colnames(Gini.df) <- c("Cluster", "Gini")
Gini.df$Cluster <- as.factor(Gini.df$Cluster)
Gini.df$Gini <- as.numeric(Gini.df$Gini)

pdf("./../Results/Seurat/Treg_GiniIndex_Whole_Cluster.pdf")
ggplot(Gini.df, aes(x = Cluster, y = Gini)) +
	geom_boxplot() +
	#stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
	xlab("Cell types") +
	ylab("Gini Index")
dev.off()


Gini.df <- data.frame()
for(hi in c("C0301", "C0302")){
	Gini.df <- rbind(Gini.df, c(hi, ineq::ineq(table(Treg$clonotype_id[Treg$hash.ID == hi]), type='Gini')))
}
colnames(Gini.df) <- c("hash.ID", "Gini")
Gini.df$hash.ID <- as.factor(Gini.df$hash.ID)
Gini.df$Gini <- as.numeric(Gini.df$Gini)

pdf("./../Results/Seurat/Treg_GiniIndex_Whole_hashID.pdf")
ggplot(Gini.df, aes(x = hash.ID, y = Gini)) +
	geom_boxplot() +
	#stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
	xlab("Cell types") +
	ylab("Gini Index")
dev.off()


for(cl in levels(Treg@active.ident)){
	Seurat.tmp <- Treg[,Treg@active.ident == cl]
	Gini.df <- data.frame()
	for(hi in c("C0301", "C0302")){
		Gini.df <- rbind(Gini.df, c(hi, ineq::ineq(table(Seurat.tmp$clonotype_id[Seurat.tmp$hash.ID == hi]), type='Gini')))
	}
	colnames(Gini.df) <- c("hash.ID", "Gini")
	Gini.df$hash.ID <- as.factor(Gini.df$hash.ID)
	Gini.df$Gini <- as.numeric(Gini.df$Gini)

	pdf(paste0("./../Results/Seurat/Treg_GiniIndex_Cluster.", cl, "_hashID.pdf"))
	print(ggplot(Gini.df, aes(x = hash.ID, y = Gini)) +
		geom_boxplot() +
		#stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
		xlab("Cell types") +
		ylab("Gini Index"))
	dev.off()
}

# GSVA between cluster 0 and 1
Seurat.tmp <- Treg[,Treg@active.ident %in% c("0", "1")]
gsva.tmp <- gsva(as.matrix(Seurat.tmp@assays$RNA@data), GS.MSigDB, parallel.sz = 28)
saveRDS(gsva.tmp, "./../Results/Seurat/Treg_Cluster.0.1_GSVA_MSigDB.norm.rds")

p <- list()
i <- 1
volcano.df <- data.frame()
for(term in row.names(gsva.tmp)){
	term.df <- as.data.frame(cbind(GSVA = gsva.tmp[term, ], cluster = as.character(Seurat.tmp@active.ident)))
	term.df$GSVA <- as.numeric(as.character(term.df$GSVA))

	p[[i]] <- ggplot(term.df, aes(x = cluster, y = GSVA)) +
		geom_violin(adjust = 3) +
		geom_boxplot() +
		ggtitle(term) +
		xlab("cluster") +
		ylab("GSVA score") +
		stat_compare_means(method="wilcox.test")

	i <- i + 1
	volcano.df <- rbind(volcano.df, c(term, mean(term.df$GSVA[term.df$cluster == "0"]) - mean(term.df$GSVA[term.df$cluster == "1"]), wilcox.test(term.df$GSVA[term.df$cluster == "0"], term.df$GSVA[term.df$cluster == "1"])$p.value))
}
pdf("./../Results/Seurat/Treg_Cluster.0.1_GSVA_MSigDB.norm_Violin.pdf", width = 36, height = 24)
print(grid.arrange(grobs = p, ncol = 9))
dev.off()

colnames(volcano.df) <- c("Term", "Diff", "pvalue")
volcano.df$Diff <- as.numeric(volcano.df$Diff)
volcano.df$pvalue <- as.numeric(volcano.df$pvalue)
volcano.df$Term <- str_split_fixed(volcano.df$Term, "[_]", 2)[,-1]

pdf("./../Results/Seurat/Treg_Cluster.0.1_GSVA_MSigDB.norm_Volcano.pdf", width=11, height=11)
print(ggplot(volcano.df, aes(x=Diff, y=-log10(pvalue))) +
	geom_point() +
	xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
	geom_text_repel(aes(x = Diff, y = -log10(pvalue), label = ifelse(pvalue < 0.05, as.character(Term),"")), max.overlaps=Inf, force=20, segment.color="grey") +
	geom_hline(yintercept=-log10(0.05)))
dev.off()

volcano.df$p.adj <- p.adjust(volcano.df$pvalue)
volcano.df <- volcano.df[volcano.df$p.adj < 0.05,]
volcano.df <- volcano.df[order(volcano.df$Diff),]
volcano.df.cut <- volcano.df[-log(volcano.df$p.adj) > 250,]
pdf("./../Results/Seurat/Treg_Cluster.0.1_GSVA_MSigDB.norm_Barplot.pdf", width=11, height=11)
ggplot(volcano.df.cut, aes(x=Term, y=Diff, fill=-log(p.adj))) + geom_bar(stat="identity") + scale_x_discrete(limits=volcano.df.cut$Term) + coord_flip() + theme(axis.text.y = element_text(size = 15))
dev.off()

# Similar GSVA as in Treg_MJ
Treg.tmp <- Treg[,Treg@active.ident %in% c(0,1)]
GS.Interest.gmt <- getGmt("./../../Treg_MJ/Data/Signature/Custom/Hallmark_Custom.gmt")
GSVA.Interest.norm <- gsva(as.matrix(Treg@assays$RNA@data), GS.Interest.gmt, parallel.sz = 28)
GSVA.Interest.tmp <- gsva(as.matrix(Treg.tmp@assays$RNA@data), GS.Interest.gmt, parallel.sz = 28)

# With font
p <- list()
i <- 1
for(term in row.names(GSVA.Interest.norm)){
    term.df <- as.data.frame(cbind(GSVA = GSVA.Interest.norm[term, ], celltype = as.character(Treg@active.ident)))
    term.df$GSVA <- as.numeric(as.character(term.df$GSVA))

		groups <- levels(Treg@active.ident)
		combs <- expand.grid(groups, groups)
		combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

		my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
		colnames(my.comps) <- NULL
		rownames(my.comps) <- NULL
		my.comps <- as.list(my.comps)


    p[[i]] <- ggplot(term.df, aes(x = celltype, y = GSVA)) +
      geom_violin(adjust=3) +
			geom_boxplot() +
			stat_compare_means(comparisons = my.comps, label="p.signif") +
      ggtitle(term) +
      xlab("Cell types") +
      ylab("GSVA score") +
      theme_cowplot()
    i <- i + 1
}
pdf("./../Results/Seurat/Treg_Whole_GSVA_Custom_Boxplot.pdf", width=21, height=7)
grid.arrange(grobs = p, ncol = 3)
dev.off()

p <- list()
i <- 1
for(term in row.names(GSVA.Interest.tmp)){
    term.df <- as.data.frame(cbind(GSVA = GSVA.Interest.tmp[term, ], celltype = as.character(Treg.tmp@active.ident)))
    term.df$GSVA <- as.numeric(as.character(term.df$GSVA))

		groups <- levels(Treg.tmp@active.ident)
		combs <- expand.grid(groups, groups)
		combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

		my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
		colnames(my.comps) <- NULL
		rownames(my.comps) <- NULL
		my.comps <- as.list(my.comps)


    p[[i]] <- ggplot(term.df, aes(x = celltype, y = GSVA)) +
      geom_violin(adjust=3) +
			geom_boxplot() +
			stat_compare_means(comparisons = my.comps, label="p.signif") +
      ggtitle(term) +
      xlab("Cell types") +
      ylab("GSVA score") +
      theme_cowplot()
    i <- i + 1
}
pdf("./../Results/Seurat/Treg_Cluster.0.1_GSVA_Custom_Boxplot.pdf", width=21, height=7)
grid.arrange(grobs = p, ncol = 3)
dev.off()


GS.Interest.gmt <- read.gmt("./../../Treg_MJ/Data/Signature/Custom/Hallmark_Custom.gmt")

ranks <- DEGs.1vs0$avg_log2FC
names(ranks) <- row.names(DEGs.1vs0)
ranks <- ranks[order(ranks, decreasing=T)]

GSEA <- GSEA(ranks, TERM2GENE = GS.Interest.gmt, verbose=FALSE, pvalueCutoff = 1)
saveRDS(GSEA, "./../Results/Seurat/Treg_DEG.1vs0_GSEA_Custom.rds")

library(enrichplot)
p <- list()
p[[1]] <- gseaplot2(GSEA, geneSetID = 1, title = GSEA$Description[1], color="chartreuse3", pvalue_table=F)
p[[2]] <- gseaplot2(GSEA, geneSetID = 2, title = GSEA$Description[2], color="chartreuse3", pvalue_table=F)
p[[3]] <- gseaplot2(GSEA, geneSetID = 3, title = GSEA$Description[3], color="chartreuse3", pvalue_table=F)
pdf("./../Results/Seurat/Treg_DEG.1vs0_GSEA_Custom_Enrichplot.pdf", width=21, height=7)
grid.arrange(grobs = p, ncol = 3)
dev.off()

# Visualizing chemokines
Chemokines <- read.delim("./../../Treg_MJ/Data/Signature/FunctionalGeneset/Chemokines", sep='\t', header=F)[,1]
ChemokineReceptors <- read.delim("./../../Treg_MJ/Data/Signature/FunctionalGeneset/ChemokineReceptors", header=F, sep='\t')[,1]

Chemokine.list <- list(Chemokines, ChemokineReceptors)
names(Chemokine.list) <- c("Chemokines", "ChemokineReceptors")

groups <- levels(Treg$Cluster)
combs <- expand.grid(groups, groups)
combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
colnames(my.comps) <- NULL
rownames(my.comps) <- NULL
my.comps <- as.list(my.comps)

p <- list()
for(geneset in names(Chemokine.list)){
  genes <- Chemokine.list[[geneset]]
  
  i <- 1
  q <- list()
  for(gene in genes){
    if(!(gene %in% row.names(Treg@assays$RNA@data))){next}
    VlnObj <- Treg[row.names(Treg) == gene,]
    VlnMat <- VlnObj@assays$RNA@data

		Vln.df <- as.data.frame(t(as.matrix(VlnMat)))
		Vln.df$Cluster <- Treg$Cluster
		
		q[[i]] <- ggplot(Vln.df, aes_string(x="Cluster", y=gene)) + geom_violin() + geom_boxplot() +
		stat_compare_means(comparisons=my.comps, label = "p.signif", method='wilcox.test') +
		ggtitle(gene) +
		theme_bw() +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank())

		i <- i + 1
  }
	p[[geneset]] <- q
}
pdf(paste0("./../Results/Seurat/Chemokines_Chemokines.pdf"), width=28, height=16)
print(grid.arrange(grobs = p[["Chemokines"]], ncol=7, left=grid::textGrob(gene)))
dev.off()
pdf(paste0("./../Results/Seurat/Chemokines_ChemokineReceptors.pdf"), width=24, height=16)
print(grid.arrange(grobs = p[["ChemokineReceptors"]], ncol=6, left=grid::textGrob(gene)))
dev.off()

# Visualizing cytokines.
Cytokines <- c("Il12a", "Ebi3", "Il10", "Il6", "Il4", "Il17a", "Il2", "Il7", "Tgfb1", "Ifng", "Gzmb", "Il23")

groups <- levels(Treg$Cluster)
combs <- expand.grid(groups, groups)
combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
colnames(my.comps) <- NULL
rownames(my.comps) <- NULL
my.comps <- as.list(my.comps)

genes <- Cytokines

i <- 1
q <- list()
for(gene in genes){
	if(!(gene %in% row.names(Treg@assays$RNA@data))){next}
	VlnObj <- Treg[row.names(Treg) == gene,]
	VlnMat <- VlnObj@assays$RNA@data

	Vln.df <- as.data.frame(t(as.matrix(VlnMat)))
	Vln.df$Cluster <- Treg$Cluster
	
	q[[i]] <- ggplot(Vln.df, aes_string(x="Cluster", y=gene)) + geom_violin() + geom_boxplot() +
	stat_compare_means(comparisons=my.comps, label = "p.signif", method='wilcox.test') +
	ggtitle(gene) +
	theme_bw() +
	theme(axis.title.x = element_blank(), axis.title.y = element_blank())

	i <- i + 1
}
pdf(paste0("./../Results/Seurat/Cytokines.pdf"), width=16, height=12)
print(grid.arrange(grobs = q, ncol=4))
dev.off()



Cytokines <- c("Il12a", "Ebi3", "Il10", "Il6", "Il4", "Il17a", "Il2", "Il7", "Tgfb1", "Ifng", "Gzmb", "Il23")
Chemokines <- Chemokine.list[[1]]
ChemokineReceptors <- Chemokine.list[[2]]

genes <- c(Cytokines, Chemokines, ChemokineReceptors)
geneset <- c(rep("Cytokines", length(Cytokines)),rep("Chemokines", length(Chemokines)),rep("ChemokineReceptors", length(ChemokineReceptors)))

geneset <- geneset[genes %in% row.names(DEGs.1vs0)]
genes <- genes[genes %in% row.names(DEGs.1vs0)]

data <- DEGs.1vs0[genes,]
data$geneset <- as.factor(geneset)
data <- data[data$p_val < 0.05,]
data$gene <- row.names(data)

empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$geneset), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$geneset <- rep(levels(data$geneset), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(geneset, avg_log2FC)
data$id <- seq(1, nrow(data))
data$color <- ifelse(data$avg_log2FC >= 0, "Red", "Blue")

# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data <- data

# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
# ----- ------------------------------------------- ---- #

base_data <- data %>% 
  group_by(geneset) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
base_data$pos <- c(-3.2, -2.6,-2.55)

pdf(paste0("./../Results/Seurat/CircularBarplot_noFont.pdf"), width=16, height=12)
ggplot(data, aes(x=as.factor(id), y=avg_log2FC, fill=color)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
	geom_bar(stat="identity") +
	ylim(2*min(na.exclude(data$avg_log2FC)), 2*max(na.exclude(data$avg_log2FC))) +
	theme_minimal() +
	theme(
		axis.text = element_blank(),
		axis.title = element_blank(),
		legend.position = 'none',
		plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
	) +
	coord_polar(start = 0) +
#	geom_text(data=label_data, aes(x=id, y=1.5*max(na.exclude(data$avg_log2FC)), label=gene, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
  geom_segment(data=base_data, aes(x = start, y = min(na.exclude(data$avg_log2FC))-0.1, xend = end, yend = min(na.exclude(data$avg_log2FC))- 0.1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )
 # geom_text(data=base_data, aes(x = title, y = pos, label=geneset), hjust=c(0.5,0.5,0.5), angle=c(mean(angle[data$geneset == "ChemokineReceptors"& !is.na(data$color)])-90,mean(angle[data$geneset == "Chemokines"& !is.na(data$color)])+90,mean(angle[data$geneset == "Cytokines"& !is.na(data$color)])-90), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
dev.off()

pdf(paste0("./../Results/Seurat/CircularBarplot_noFontGrid.pdf"), width=16, height=12)
ggplot(data, aes(x=as.factor(id), y=avg_log2FC, fill=color)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
	geom_bar(stat="identity") +
	ylim(2*min(na.exclude(data$avg_log2FC)), 2*max(na.exclude(data$avg_log2FC))) +
	theme_minimal() +
	theme(
		axis.text = element_blank(),
		axis.title = element_blank(),
		legend.position = 'none',
		plot.margin = unit(rep(-1,4), "cm"),      # Adjust the margin to make in sort labels are not truncated!
		panel.grid.major = element_blank(), panel.grid.minor = element_blank()
	) +
	coord_polar(start = 0) +
  geom_segment(data=base_data, aes(x = start, y = min(na.exclude(data$avg_log2FC))-0.1, xend = end, yend = min(na.exclude(data$avg_log2FC))- 0.1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )
dev.off()




##########################################################################################
# Integration analysis
#Int <- readRDS("./../../Treg_MJ/Results/NewTreg/Seurat/Int_SeuratObj.rds")

# Loading object.
Int.list <- readRDS("./../../Treg_MJ/Results/Seurat/Int.list_Raw.rds")

# Filtering wanted cells.
Treg_MJ.list <- list()
for(l in 1:3){
	Treg_MJ.list[[l]] <- Int.list[[l]][,(Int.list[[l]]$Celltype == "Treg_Norm" | Int.list[[l]]$Celltype == "Treg_KO") & Int.list[[l]]$orig.ident == "Iso"]
}
Treg_MJ.list[[1]]$orig.ident <- "Treg.PD1.1"
Treg_MJ.list[[2]]$orig.ident <- "Treg.PD1.2"
Treg_MJ.list[[3]]$orig.ident <- "Treg.PD1.3"

Treg_MJ.list[[1]]$clonotype_id <- gsub("Iso", "Treg.PD1.1", Treg_MJ.list[[1]]$clonotype_id)
Treg_MJ.list[[2]]$clonotype_id <- paste0("Treg.PD1.2_", Treg_MJ.list[[2]]$clonotype_id)
Treg_MJ.list[[3]]$clonotype_id <- paste0("Treg.PD1.3_", Treg_MJ.list[[3]]$clonotype_id)

saveRDS(Treg_MJ.list, "./../Results/Treg_Int/Treg.MJ.list_SeuratObj.rds")

Treg$Treg.TOX_cluster <- Treg@active.ident
tmp <- Treg[,Treg$Treg.TOX_cluster != 3]
tmp$orig.ident <- "Treg.TOX"

Treg.list <- list.append(tmp, Treg_MJ.list)
names(Treg.list) <- c("Treg.TOX", "Treg.PD1.1", "Treg.PD1.2", "Treg.PD1.3")

saveRDS(Treg.list, "./../Results/Treg_Int/Int.list_SeuratObj.rds")


# Integration starts!!
Treg.list <- lapply(X = Treg.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Treg.list)

Treg.anchors <- FindIntegrationAnchors(object.list = Treg.list, anchor.features = features)

Treg.intd <- IntegrateData(anchorset = Treg.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Treg.intd) <- "integrated"

# Run the standard workflow for visualization and clustering
Treg.intd <- ScaleData(Treg.intd, verbose = FALSE)
Treg.intd <- RunPCA(Treg.intd, npcs = 30, verbose = FALSE)
Treg.intd <- RunUMAP(Treg.intd, reduction = "pca", dims = 1:30)
Treg.intd <- FindNeighbors(Treg.intd, reduction = "pca", dims = 1:30)
Treg.intd <- FindClusters(Treg.intd, resolution = seq(0.1, 1, 0.1))

feat <- c("Cd4", "Cd8", "Tox", "Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Ifit3", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

png(paste0("./../Results/Treg_Int/Intd.2000_UMAP_Feature.png"), width = 6400, height = 4000, res = 144)
print(FeaturePlot(Treg.intd, reduction = "umap", features = feat, ncol = 8))
dev.off()

pdf(paste0("./../Results/Treg_Int/Intd.2000_UMAP_Hash.ID.pdf"))
print(DimPlot(Treg.intd, group.by = "hash.ID"))
dev.off()

pdf(paste0("./../Results/Treg_Int/Intd.2000_UMAP_Celltype.pdf"))
print(DimPlot(Treg.intd, group.by = "Celltype"))
dev.off()

pdf(paste0("./../Results/Treg_Int/Intd.2000_UMAP_Cluster.TOX.pdf"))
print(DimPlot(Treg.intd, group.by = "Treg.TOX_cluster"))
dev.off()

saveRDS(Treg.intd, "./../Results/Treg_Int/Intd.2000_SeuratObj.rds")


# Integrating using whole geens
for(x in 1:length(Treg.list)){
    Treg.list[[x]] <- NormalizeData(Treg.list[[x]])
    VariableFeatures(Treg.list[[x]]) <- row.names(Treg.list[[x]])
}

features <- SelectIntegrationFeatures(object.list = Treg.list, nfeatures = 200000) #13686

Treg.anchors <- FindIntegrationAnchors(object.list = Treg.list, anchor.features = features)

Treg.intd <- IntegrateData(anchorset = Treg.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Treg.intd) <- "integrated"

# Run the standard workflow for visualization and clustering
Treg.intd <- ScaleData(Treg.intd, verbose = FALSE)
Treg.intd <- RunPCA(Treg.intd, npcs = 30, verbose = FALSE)
Treg.intd <- RunUMAP(Treg.intd, reduction = "pca", dims = 1:30)
Treg.intd <- FindNeighbors(Treg.intd, reduction = "pca", dims = 1:30)
Treg.intd <- FindClusters(Treg.intd, resolution = seq(0.1, 1, 0.1))

feat <- c("Cd4", "Cd8", "Tox", "Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Ifit3", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

png(paste0("./../Results/Treg_Int/Intd.13686_UMAP_Feature.png"), width = 6400, height = 4000, res = 144)
print(FeaturePlot(Treg.intd, reduction = "umap", features = feat, ncol = 8))
dev.off()

pdf(paste0("./../Results/Treg_Int/Intd.13686_UMAP_Hash.ID.pdf"))
print(DimPlot(Treg.intd, group.by = "hash.ID"))
dev.off()

pdf(paste0("./../Results/Treg_Int/Intd.13686_UMAP_Celltype.pdf"))
print(DimPlot(Treg.intd, group.by = "Celltype"))
dev.off()

pdf(paste0("./../Results/Treg_Int/Intd.13686_UMAP_Cluster.TOX.pdf"))
print(DimPlot(Treg.intd, group.by = "Treg.TOX_cluster"))
dev.off()

saveRDS(Treg.intd, "./../Results/Treg_Int/Intd.13686_SeuratObj.rds")

Treg.intd <- readRDS("./../Results/Treg_Int/Intd.2000_SeuratObj.rds")

pdf("./../Results/Treg_Int/Intd_UMAP_Clusters.pdf")
for(j in seq(0.1, 1, 0.1)){
	print(DimPlot(Treg.intd, reduction = "umap", label = T, pt.size = 0.5, group.by = paste0("integrated_snn_res.", j)))
}
dev.off()

Treg.intd <- FindClusters(Treg.intd, resolution = 0.7)

pdf("./../Results/Treg_Int/Intd_UMAP_Cluster.pdf")
DimPlot(Treg.intd, reduction = "umap", label = T, pt.size = 0.9)
dev.off()

Treg.intd$Cluster.int <- Treg.intd@active.ident
levels(Treg.intd$Cluster.int) <- c("WT_1", "WT_1", "TOX_KO", "WT_1", "ISG_High", "PD1_KO", "WT_1", "WT_2", "TOX_High", "WT_1")

pdf("./../Results/Treg_Int/Intd_UMAP_Cluster.pdf")
print(DimPlot(Treg.intd, reduction = "umap", label = T, pt.size = 0.9))
print(DimPlot(Treg.intd, reduction = "umap", label = T, pt.size = 0.9, group.by="Cluster.int"))
dev.off()

# DEG analysis
# Finding DEGs of clusters.
DEGs.intd <- FindAllMarkers(Treg.intd, logfc.threshold = 0, min.pct = 0)
DEGs.intd$threshold <- DEGs.intd$p_val_adj < 0.05 & abs(DEGs.intd$avg_log2FC) > 0.25
write.table(DEGs.intd, "./../Results/Treg_Int/Intd_DEG_Cluster.tsv", quote = F, sep = '\t', col.names=NA)

PlottingVolcano <- function(WholeDEG, Output, feat){
	pdf(Output)
	WholeDEG$cluster <- as.factor(as.character(WholeDEG$cluster))
	for(cluster in levels(WholeDEG$cluster)){
		DEG.list <- WholeDEG[WholeDEG$cluster == cluster,]
		DEG.list$threshold <- DEG.list$p_val_adj < 0.05 & abs(DEG.list$avg_log2FC) > 0.25
		DEG.list$pct.AbsLogRatio <- abs(log2(DEG.list$pct.1 / DEG.list$pct.2))
		DEG.list <- DEG.list[order(DEG.list$avg_log2FC, decreasing = T),]

		DEG.list$genelabels <- "FALSE"
		SigMark <- DEG.list[DEG.list$threshold,]
		DEG.list$genelabels[row.names(DEG.list) %in% head(row.names(SigMark), 20)] <- "TRUE"
		DEG.list$genelabels[row.names(DEG.list) %in% tail(row.names(SigMark), 20)] <- "TRUE"
		DEG.list$genelabels[DEG.list$gene %in% feat & row.names(DEG.list) %in% row.names(SigMark)] <- "TRUE"


		gg <- ggplot(DEG.list) +
			geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj), colour = threshold)) +
			ggrepel::geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = ifelse(genelabels == T, as.character(gene),"")), max.overlaps=Inf, force=20, segment.color="grey") +
      ggtitle(paste0("VolcanoPlot of cluster ", cluster)) +
			xlab("Average log2 fold change") +
			ylab("-log10 adjusted p-value") +
			theme(legend.position = "none",
			plot.title = element_text(size = rel(1.5), hjust = 0.5),
			axis.title = element_text(size = rel(1.25)))
		print(gg)
	}
	dev.off()
}
feats <- c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1", "Tox")
PlottingVolcano(DEGs.intd, "./../Results/Treg_Int/Intd_DEG_Volcano_Cluster.pdf", feats)

Treg.intd$Sample <- Treg.intd$orig.ident
Treg.intd$Sample[Treg.intd$hash.ID %in% "C0301"] <- "Treg.TOX.1"
Treg.intd$Sample[Treg.intd$hash.ID %in% "C0302"] <- "Treg.TOX.2"

Treg.intd$clonotype_id[Treg.intd$Sample == "Treg.TOX.1"] <- gsub("C0301-", "Treg.TOX.1_", Treg.intd$clonotype_id[Treg.intd$Sample == "Treg.TOX.1"])
Treg.intd$clonotype_id[Treg.intd$Sample == "Treg.TOX.2"] <- gsub("C0302-", "Treg.TOX.2_", Treg.intd$clonotype_id[Treg.intd$Sample == "Treg.TOX.2"])

div_clonotype <- function(seurat_obj){
	# Extracting clonotype information
	clonotable <- table(seurat_obj$clonotype_id)
	clonotable <- clonotable[clonotable >= 1]
	multiclone <- names(clonotable[clonotable > 1])
	singleclone <- names(clonotable[clonotable == 1])

	# Making object
	clonotype <- rep("NA", ncol(seurat_obj))
	clonotype[which(seurat_obj$clonotype_id %in% multiclone)] <- "multi-clonotype"
	clonotype[which(seurat_obj$clonotype_id %in% singleclone)] <- "single-clonotype"

	# Calculating the number of each clones
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotable[as.character(seurat_obj$clonotype_id)], col.name = "cloneNum")
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotype, col.name = "clonotype")
	return(seurat_obj)
}
Treg.intd <- div_clonotype(Treg.intd)

pdf("./../Results/Treg_Int/Intd_UMAP_clonotype.pdf")
print(DimPlot(Treg.intd, group.by = "clonotype") + scale_colour_manual(values=c("Red", "Gray")))
dev.off()


# Gini Index
Gini.df <- data.frame()
for(ct in levels(Treg.intd$Cluster.int)){
	Gini.df <- rbind(Gini.df, c(ct, ineq::ineq(table(Treg.intd$clonotype_id[Treg.intd$Cluster.int == ct]), type='Gini')))
}
colnames(Gini.df) <- c("Cluster", "Gini")
Gini.df$Cluster <- as.factor(Gini.df$Cluster)
Gini.df$Gini <- as.numeric(Gini.df$Gini)

# pdf("./../Results/Treg_Int/Intd_GiniIndex_Whole_Cluster.pdf")
# ggplot(Gini.df, aes(x = Cluster, y = Gini)) +
# 	geom_boxplot() +
# 	#stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
# 	xlab("Cell types") +
# 	ylab("Gini Index")
# dev.off()

pdf("./../Results/Treg_Int/Intd_GiniIndex_Whole_Cluster_Barplot.pdf")
ggplot(Gini.df, aes(x = Cluster, y = Gini)) +
	geom_bar(stat="identity") +
	xlab("Cell types") +
	ylab("Gini Index") +
	theme_bw()
dev.off()

clonotype.df <- as.data.frame(cbind(Treg.intd$clonotype, as.character(Treg.intd$Cluster.int)))
colnames(clonotype.df) <- c("Clonotype", "Cluster")
clonotype.df$Cluster <- as.factor(clonotype.df$Cluster)

pdf("./../Results/Treg_Int/Intd_Whole_BarPlot_clonotype.cluster.pdf")
ggplot(clonotype.df, aes(x = Cluster, fill=Clonotype)) +
	geom_bar(position="fill") +
	xlab("Cell types") +
	ylab("Proportion") +
	theme_classic()
dev.off()


# Mosaic plot between clonotype and clusters
MosaicData <- Treg.intd
MosaicData$clonotype <- as.factor(as.character(MosaicData$clonotype))
MetaInfo <- as.factor(as.character(MosaicData$Cluster.int))
ClusterInfo <- MosaicData$clonotype
MosaicTable <- table(MetaInfo, ClusterInfo)

Mosaic <- vcd::structable(ClusterInfo ~ MetaInfo)
#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
pdf("./../Results/Treg_Int/Intd_Whole_MosaicPlot_clonotype.cluster.pdf", width=8)
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

MosaicData$Sample <- as.factor(as.character(MosaicData$Sample))
MetaInfo <- as.factor(as.character(MosaicData$Cluster.int))
ClusterInfo <- MosaicData$Sample
MosaicTable <- table(MetaInfo, ClusterInfo)

Mosaic <- vcd::structable(ClusterInfo ~ MetaInfo)
#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
pdf("./../Results/Treg_Int/Intd_Whole_MosaicPlot_Sample.cluster.pdf", width=8)
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

Treg.intd@active.ident <- Treg.intd$Cluster.int
saveRDS(Treg.intd, "./../Results/Treg_Int/Intd.13686_SeuratObj.rds")

# DEGs between clusters
DEGs.WTvsPD1KO <- FindMarkers(Treg.intd, logfc.threshold = 0, min.pct = 0, ident.1="WT_1", ident.2="PD1_KO")
DEGs.WTvsTOXKO <- FindMarkers(Treg.intd, logfc.threshold = 0, min.pct = 0, ident.1="WT_1", ident.2="TOX_KO")

DEGs.WTvsPD1KO$threshold <- DEGs.WTvsPD1KO$p_val_adj < 0.05 & abs(DEGs.WTvsPD1KO$avg_log2FC) > 0.25
DEGs.WTvsTOXKO$threshold <- DEGs.WTvsTOXKO$p_val_adj < 0.05 & abs(DEGs.WTvsTOXKO$avg_log2FC) > 0.25

write.table(DEGs.WTvsPD1KO, "./../Results/Treg_Int/Intd_DEG.Intd2000_WTvsPD1KO.tsv", quote = F, sep = '\t', col.names=NA)
write.table(DEGs.WTvsTOXKO, "./../Results/Treg_Int/Intd_DEG.Intd2000_WTvsTOXKO.tsv", quote = F, sep = '\t', col.names=NA)
DEGs.WTvsPD1KO <- read.delim("./../Results/Treg_Int/Intd_DEG.Intd2000_WTvsPD1KO.tsv", row.names=1)
DEGs.WTvsTOXKO <- read.delim("./../Results/Treg_Int/Intd_DEG.Intd2000_WTvsTOXKO.tsv", row.names=1)

PlottingVolcano.single(DEGs.WTvsPD1KO, "./../Results/Treg_Int/Treg_DEG.Intd2000_WTvsPD1KO_Volcano.pdf", feats)
PlottingVolcano.single(DEGs.WTvsTOXKO, "./../Results/Treg_Int/Treg_DEG.Intd2000_WTvsTOXKO_Volcano.pdf", feats)

DEGs.WTvsPD1KO.Pos <- row.names(DEGs.WTvsPD1KO)[DEGs.WTvsPD1KO$avg_log2FC > 0 & DEGs.WTvsPD1KO$p_val_adj < 0.05]
DEGs.WTvsPD1KO.Neg <- row.names(DEGs.WTvsPD1KO)[DEGs.WTvsPD1KO$avg_log2FC < 0 & DEGs.WTvsPD1KO$p_val_adj < 0.05]
DEGs.WTvsTOXKO.Pos <- row.names(DEGs.WTvsTOXKO)[DEGs.WTvsTOXKO$avg_log2FC > 0 & DEGs.WTvsTOXKO$p_val_adj < 0.05]
DEGs.WTvsTOXKO.Neg <- row.names(DEGs.WTvsTOXKO)[DEGs.WTvsTOXKO$avg_log2FC < 0 & DEGs.WTvsTOXKO$p_val_adj < 0.05]

WTHigh <- list(WTHigh_PD1KO = DEGs.WTvsPD1KO.Pos, WTHigh_TOXKO = DEGs.WTvsTOXKO.Pos)
WTLow <- list(WTLow_PD1KO = DEGs.WTvsPD1KO.Neg, WTLow_TOXKO = DEGs.WTvsTOXKO.Neg)
# columns = c("Stage 1", "Stage 2")
pdf("./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_Pos_VennDiagram.pdf",height=5)
ggvenn(WTHigh)
dev.off()
pdf("./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_Neg_VennDiagram.pdf",height=5)
ggvenn(WTLow)
dev.off()

library(grid)
library(pheatmap)
library(circlize)

Expression <- t(Treg.intd@assays$integrated@scale.data)
genes.intersect <- intersect(WTHigh[[1]], WTHigh[[2]])
genes <- c(setdiff(WTHigh[[1]], genes.intersect), genes.intersect, setdiff(WTHigh[[2]], genes.intersect))
Expression.cut <- Expression[,genes]
Expression.cut <- Expression.cut[order(Treg.intd$Cluster.int),]

annot_col <- as.data.frame(c(rep("WTHigh_PD1KO",length(setdiff(WTHigh[[1]], genes.intersect))), rep("WTHigh_Both", length(genes.intersect)), rep("WTHigh_TOXKO", length(setdiff(WTHigh[[2]], genes.intersect)))))
row.names(annot_col) <- genes
colnames(annot_col) <- "DEGs"
annot_col$DEGs <- as.factor(annot_col$DEGs)
annot_col$DEGs <- factor(annot_col$DEGs, levels=c("WTHigh_PD1KO", "WTHigh_Both", "WTHigh_TOXKO"))

annot_row <- as.data.frame(sort(Treg.intd$Cluster.int))
colnames(annot_row) <- "Cluster"

Expression.cut <- Expression.cut[annot_row$Cluster %in% c("WT_1", "TOX_KO", "PD1_KO"),]
annot_row <- as.data.frame(annot_row[annot_row$Cluster %in% c("WT_1", "TOX_KO", "PD1_KO"),])
colnames(annot_row) <- "Cluster"
row.names(annot_row) <- row.names(Expression.cut)
annot_row$Cluster <- as.factor(as.character(annot_row$Cluster))
annot_row$Cluster <- factor(annot_row$Cluster, levels=c("WT_1","TOX_KO", "PD1_KO"))

Expression.Cellmean <- data.frame()
for(group in levels(annot_row$Cluster)){
	Expression.Cellmean <- rbind(Expression.Cellmean, colMeans(Expression.cut[annot_row$Cluster == group,]))
}
colnames(Expression.Cellmean) <- colnames(Expression.cut)
row.names(Expression.Cellmean) <- levels(annot_row$Cluster)

pdf("./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_Pos_Heatmap.pdf", width=10, height=10)
# png("tmp.png", width=1000, height=1000)
print(pheatmap(Expression.cut,
	color=colorRampPalette(c("green", "black", "red"))(102),
	legend_breaks =  c(-2, -1, 0, 1, 2),
	breaks =  c(-Inf, seq(-2, 2, 0.04), Inf),
	border_color=NA,
	cluster_cols = F,
	cluster_rows = F,
	annotation_row = annot_row,
	annotation_col = annot_col, 
	show_colnames=F, 
	show_rownames=F,
	gaps_col = cumsum(as.vector(table(annot_col$DEGs)[c("WTHigh_PD1KO", "WTHigh_Both", "WTHigh_TOXKO")])),
	gaps_row = cumsum(as.vector(table(annot_row$Cluster)[c("WT_1","TOX_KO", "PD1_KO")])),
	labels_row = F,
	fontsize=15))
dev.off()

pdf("./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_Pos_Heatmap_Cellmean.pdf", width=10, height=10)
# png("tmp.png", width=1000, height=1000)
print(pheatmap(t(Expression.Cellmean),
	color=colorRampPalette(c("green", "black", "red"))(102),
	legend_breaks =  c(-2, -1, 0, 1, 2),
	breaks =  c(-Inf, seq(-2, 2, 0.04), Inf),
	border_color=NA,
	cluster_cols = F,
	cluster_rows = F,
	#annotation_row = annot_row,
	annotation_row = annot_col, 
	show_colnames=T, 
	show_rownames=F,
	gaps_row = cumsum(as.vector(table(annot_col$DEGs)[c("WTHigh_PD1KO", "WTHigh_Both", "WTHigh_TOXKO")])),
	#gaps_row = cumsum(as.vector(table(annot_row$Cluster)[c("WT_1","TOX_KO", "PD1_KO")])),
	labels_row = T,
	fontsize=15))
dev.off()


DEGs.WTvsPD1KO.Pos <- row.names(DEGs.WTvsPD1KO)[DEGs.WTvsPD1KO$avg_log2FC > 0 & DEGs.WTvsPD1KO$p_val_adj < 0.0001]
DEGs.WTvsPD1KO.Neg <- row.names(DEGs.WTvsPD1KO)[DEGs.WTvsPD1KO$avg_log2FC < 0 & DEGs.WTvsPD1KO$p_val_adj < 0.0001]
DEGs.WTvsTOXKO.Pos <- row.names(DEGs.WTvsTOXKO)[DEGs.WTvsTOXKO$avg_log2FC > 0 & DEGs.WTvsTOXKO$p_val_adj < 0.0001]
DEGs.WTvsTOXKO.Neg <- row.names(DEGs.WTvsTOXKO)[DEGs.WTvsTOXKO$avg_log2FC < 0 & DEGs.WTvsTOXKO$p_val_adj < 0.0001]

WTHigh <- list(WTHigh_PD1KO = DEGs.WTvsPD1KO.Pos, WTHigh_TOXKO = DEGs.WTvsTOXKO.Pos)

genes.intersect <- intersect(WTHigh[[1]], WTHigh[[2]])
genes <- c(setdiff(WTHigh[[1]], genes.intersect), genes.intersect, setdiff(WTHigh[[2]], genes.intersect))

Expression.cut <- Expression[,genes]
Expression.cut <- Expression.cut[order(Treg.intd$Cluster.int),]

annot_col <- as.data.frame(c(rep("WTHigh_PD1KO",length(setdiff(WTHigh[[1]], genes.intersect))), rep("WTHigh_Both", length(genes.intersect)), rep("WTHigh_TOXKO", length(setdiff(WTHigh[[2]], genes.intersect)))))
row.names(annot_col) <- genes
colnames(annot_col) <- "DEGs"
annot_col$DEGs <- as.factor(annot_col$DEGs)
annot_col$DEGs <- factor(annot_col$DEGs, levels=c("WTHigh_PD1KO", "WTHigh_Both", "WTHigh_TOXKO"))

annot_row <- as.data.frame(sort(Treg.intd$Cluster.int))
colnames(annot_row) <- "Cluster"

Expression.cut <- Expression.cut[annot_row$Cluster %in% c("WT_1", "TOX_KO", "PD1_KO"),]
annot_row <- as.data.frame(annot_row[annot_row$Cluster %in% c("WT_1", "TOX_KO", "PD1_KO"),])
colnames(annot_row) <- "Cluster"
row.names(annot_row) <- row.names(Expression.cut)
annot_row$Cluster <- as.factor(as.character(annot_row$Cluster))
annot_row$Cluster <- factor(annot_row$Cluster, levels=c("WT_1","TOX_KO", "PD1_KO"))

Expression.Cellmean <- data.frame()
for(group in levels(annot_row$Cluster)){
	Expression.Cellmean <- rbind(Expression.Cellmean, colMeans(Expression.cut[annot_row$Cluster == group,]))
}
colnames(Expression.Cellmean) <- colnames(Expression.cut)
row.names(Expression.Cellmean) <- levels(annot_row$Cluster)

pdf("./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_Pos_Heatmap_Cellmean_p.0001.pdf", width=10, height=10)
# png("tmp.png", width=1000, height=1000)
print(pheatmap(t(Expression.Cellmean),
	color=colorRampPalette(c("green", "black", "red"))(102),
	legend_breaks =  c(-2, -1, 0, 1, 2),
	breaks =  c(-Inf, seq(-2, 2, 0.04), Inf),
	border_color=NA,
	cluster_cols = F,
	cluster_rows = F,
	#annotation_row = annot_row,
	annotation_row = annot_col, 
	show_colnames=T, 
	show_rownames=F,
	gaps_row = cumsum(as.vector(table(annot_col$DEGs)[c("WTHigh_PD1KO", "WTHigh_Both", "WTHigh_TOXKO")])),
	#gaps_row = cumsum(as.vector(table(annot_row$Cluster)[c("WT_1","TOX_KO", "PD1_KO")])),
	labels_row = T,
	fontsize=15))
dev.off()




DEGs.intersect.Pos <- intersect(DEGs.WTvsPD1KO.Pos,DEGs.WTvsTOXKO.Pos)
DEGs.intersect.Neg <- intersect(DEGs.WTvsPD1KO.Neg,DEGs.WTvsTOXKO.Neg)
# write.table(t(DEGs.intersect.Pos), "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_PosGene.tsv", quote = F, sep = '\t', col.names=NA)
# write.table(t(DEGs.intersect.Neg), "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_NegGene.tsv", quote = F, sep = '\t', col.names=NA)

write.table(DEGs.intersect.Pos, "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_PosGene_433.tsv", quote = F, sep = '\t', col.names=NA)
write.table(setdiff(DEGs.WTvsPD1KO.Pos, DEGs.intersect.Pos), "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_PosGene_562.tsv", quote = F, sep = '\t', col.names=NA)
write.table(setdiff(DEGs.WTvsTOXKO.Pos, DEGs.intersect.Pos), "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_PosGene_230.tsv", quote = F, sep = '\t', col.names=NA)


colnames(DEGs.WTvsPD1KO) <- paste0("PD1_", colnames(DEGs.WTvsPD1KO))
colnames(DEGs.WTvsTOXKO) <- paste0("TOX_", colnames(DEGs.WTvsTOXKO))
DEGs.df <- cbind(DEGs.WTvsPD1KO, DEGs.WTvsTOXKO[row.names(DEGs.WTvsPD1KO),])
DEGs.df$gene <- row.names(DEGs.df)
DEGs.df$group <- "NotSig"
DEGs.df$group[DEGs.df$PD1_p_val_adj < 0.05] <- "PD1Sig"
DEGs.df$group[DEGs.df$TOX_p_val_adj < 0.05] <- "TOXSig"
DEGs.df$group[row.names(DEGs.df) %in% c(DEGs.intersect.Pos, DEGs.intersect.Neg)] <- "DoubleSig"
DEGs.df$group <- factor(DEGs.df$group, levels=c("NotSig", "PD1Sig", "TOXSig", "DoubleSig"))
write.table(DEGs.df, "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO.tsv", quote = F, sep = '\t', col.names=NA)

write.table(DEGs.df[DEGs.df$group == "TOXSig" & DEGs.df$TOX_avg_log2FC > 0,], "./../Results/Treg_Int/Intd_DEG.Intd2000_TOXSig_WTHigh.tsv", sep='\t', quote=F, col.names=NA)
write.table(DEGs.df[DEGs.df$group == "PD1Sig" & DEGs.df$PD1_avg_log2FC > 0,], "./../Results/Treg_Int/Intd_DEG.Intd2000_PD1Sig_WTHigh.tsv", sep='\t', quote=F, col.names=NA)

pdf("./../Results/Treg_Int/Intd_DEG.Intd2000_PD1KO.TOXKO_Scatterplot.pdf")
ggplot(DEGs.df) +
	geom_point(aes(x = PD1_avg_log2FC, y = TOX_avg_log2FC, colour = group)) +
	#ggrepel::geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = ifelse(genelabels == T, as.character(gene),"")), max.overlaps=Inf, force=20, segment.color="grey") +
	xlab("PD-1 KO <---> WT") +
	ylab("TOX KO <---> WT") +
	theme(
		plot.title = element_text(size = rel(1.5), hjust = 0.5),
		axis.title = element_text(size = rel(1.25))) +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=0) +
	scale_color_manual(values=c("gray", "blue", "green", "red"))
dev.off()

# GSVA
Treg.PD1 <- Treg.intd[,Treg.intd@active.ident %in% c("WT_1", "PD1_KO")]
Treg.TOX <- Treg.intd[,Treg.intd@active.ident %in% c("WT_1", "TOX_KO")]
gsva.PD1 <- gsva(as.matrix(Treg.PD1@assays$integrated@data), GS.MSigDB, parallel.sz = 28)
gsva.TOX <- gsva(as.matrix(Treg.TOX@assays$integrated@data), GS.MSigDB, parallel.sz = 28)
saveRDS(gsva.PD1, "./../Results/Treg_Int/Intd_Whole_GSVA_PD1_MSigDB.norm.rds")
saveRDS(gsva.TOX, "./../Results/Treg_Int/Intd_Whole_GSVA_TOX_MSigDB.norm.rds")
gsva.PD1 <- readRDS("./../Results/Treg_Int/Intd_Whole_GSVA_PD1_MSigDB.norm.rds")
gsva.TOX <- readRDS("./../Results/Treg_Int/Intd_Whole_GSVA_TOX_MSigDB.norm.rds")

volcano.df <- data.frame()
for(term in row.names(gsva.PD1)){
	term1.df <- as.data.frame(cbind(GSVA = gsva.PD1[term, ], cluster = as.character(Treg.PD1@active.ident)))
	term1.df$GSVA <- as.numeric(as.character(term1.df$GSVA))
	term2.df <- as.data.frame(cbind(GSVA = gsva.TOX[term, ], cluster = as.character(Treg.TOX@active.ident)))
	term2.df$GSVA <- as.numeric(as.character(term2.df$GSVA))

	term1.fc <- mean(term1.df$GSVA[term1.df$cluster == "WT_1"]) - mean(term1.df$GSVA[term1.df$cluster == "PD1_KO"])
	term1.pval <- wilcox.test(term1.df$GSVA[term1.df$cluster == "WT_1"], term1.df$GSVA[term1.df$cluster == "PD1_KO"])$p.value
	term2.fc <- mean(term2.df$GSVA[term2.df$cluster == "WT_1"]) - mean(term2.df$GSVA[term2.df$cluster == "TOX_KO"])
	term2.pval <- wilcox.test(term2.df$GSVA[term2.df$cluster == "WT_1"], term2.df$GSVA[term2.df$cluster == "TOX_KO"])$p.value

	volcano.df <- rbind(volcano.df, c(term, term1.fc, term1.pval, term2.fc, term2.pval))
}

colnames(volcano.df) <- c("Term", "PD1_Diff", "PD1_pvalue", "TOX_Diff", "TOX_pvalue")
volcano.df$PD1_Diff <- as.numeric(volcano.df$PD1_Diff)
volcano.df$PD1_pvalue <- as.numeric(volcano.df$PD1_pvalue)
volcano.df$TOX_Diff <- as.numeric(volcano.df$TOX_Diff)
volcano.df$TOX_pvalue <- as.numeric(volcano.df$TOX_pvalue)
volcano.df$Term <- str_split_fixed(volcano.df$Term, "[_]", 2)[,-1]
volcano.df$PD1_adjp <- p.adjust(volcano.df$PD1_pvalue)
volcano.df$TOX_adjp <- p.adjust(volcano.df$TOX_pvalue)

volcano.df$group <- "NotSig"
volcano.df$group[volcano.df$PD1_adjp < 0.05 & volcano.df$PD1_Diff > 0] <- "PD1_WTHigh"
volcano.df$group[volcano.df$TOX_adjp < 0.05 & volcano.df$TOX_Diff > 0] <- "TOX_WTHigh"
volcano.df$group[volcano.df$PD1_adjp < 0.05 & volcano.df$PD1_Diff > 0 & volcano.df$TOX_adjp < 0.05 & volcano.df$TOX_Diff > 0] <- "Double_WTHigh"
volcano.df$group <- factor(volcano.df$group, levels=c("NotSig", "PD1_WTHigh", "TOX_WTHigh", "Double_WTHigh"))

volcano.df$group[volcano.df$PD1_Diff < 0 & volcano.df$TOX_Diff < 0] <- "NotSig"
volcano.df$group[(volcano.df$PD1_Diff < 0 | volcano.df$TOX_Diff < 0) & volcano.df$group == "DoubleSig"] <- "NotSig"

pdf("./../Results/Treg_Int/Intd_GSVA.Intd2000_PD1KO.TOXKO_MSigDB.norm_Volcano.pdf", width=11, height=11)
for(sz in seq(4, 5, 0.1)){
	print(ggplot(volcano.df, aes(x=PD1_Diff, y=TOX_Diff, color=group)) +
		geom_point() +
		#xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
		geom_text_repel(aes(x = PD1_Diff, y = TOX_Diff, label = ifelse(group != "NotSig", as.character(Term),"")), size=sz, max.overlaps=Inf, force=20, segment.color="grey") +
		geom_hline(yintercept=0) +
		geom_vline(xintercept=0) +
		scale_color_manual(values=c("gray", "blue", "green", "red")) +
		theme_bw())
}
dev.off()

pdf("./../Results/Treg_Int/Intd_GSVA.Intd2000_PD1KO.TOXKO_MSigDB.norm_Volcano_nofont1.pdf", width=11, height=11)
ggplot(volcano.df, aes(x=PD1_Diff, y=TOX_Diff, color=group)) +
	geom_point() +
	#xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
	geom_text_repel(aes(x = PD1_Diff, y = TOX_Diff, label = ifelse(group != "NotSig", " ","")), max.overlaps=Inf, force=20, segment.color="grey") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=0) +
	scale_color_manual(values=c("gray", "blue", "green", "red")) +
	theme_bw() +
  theme(text = element_text(size=0))
dev.off()

pdf("./../Results/Treg_Int/Intd_GSVA.Intd2000_PD1KO.TOXKO_MSigDB.norm_Volcano_nofont2.pdf", width=11, height=11)
ggplot(volcano.df, aes(x=PD1_Diff, y=TOX_Diff, color=group)) +
	geom_point() +
	#xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
	geom_text_repel(aes(x = PD1_Diff, y = TOX_Diff, label = ifelse(group != "NotSig", "","")), max.overlaps=Inf, force=20, segment.color="grey") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=0) +
	scale_color_manual(values=c("gray", "blue", "green", "red")) +
	theme_bw() +
  theme(text = element_text(size=0))
dev.off()

pdf("./../Results/Treg_Int/Intd_GSVA.Intd2000_PD1KO.TOXKO_MSigDB.norm_Barplot.pdf", width=11, height=11)
for(case in  c("PD1Sig", "TOXSig")){
	volcano.df.cut <- subset(volcano.df, group == case)
	volcano.df.cut[,"logP"] <- -log(volcano.df.cut[,paste0(gsub("Sig", "", case), "_adjp")])
	volcano.df.cut <- volcano.df.cut[order(volcano.df.cut[,paste0(gsub("Sig", "", case), "_Diff")]),]
	print(ggplot(volcano.df.cut, aes_string(x="Term", y=paste0(gsub("Sig", "", case), "_Diff"), fill="logP")) + 
		geom_bar(stat="identity") + 
		scale_x_discrete(limits=volcano.df.cut$Term) + 
		coord_flip() + 
		theme(axis.text.y = element_text(size = 15),axis.title.x = element_text(size = 15)) + 
		ylab(paste0(gsub("Sig", "", case), " KO <---> WT")) + 
		scale_fill_continuous(name="-log(adjP)") + xlab("")
	)
}
dev.off()


# Randomize DEGs































DEGs$threshold <- DEGs$p_val_adj < 0.05 & DEGs$avg_log2FC > 0.25
DEGs.Sig <- DEGs[DEGs$threshold == T,]
DEGs.table <- table(DEGs.Sig$gene)
DEGs.excl <- names(DEGs.table)[DEGs.table == 1]
DEGs.exclusive <- DEGs.Sig[DEGs.Sig$gene %in% DEGs.excl,]
write.table(DEGs.exclusive, "./../Results/NewTreg/Seurat/Int_DEG.Exclusive_Cluster.tsv", quote = F, sep = '\t', col.names=NA)


pdf("./../Results/NewTreg/Seurat/Int_Dotplot_Cluster.pdf", width=21)
DotPlot(Int, features = feat, dot.scale = 8) + RotatedAxis()
dev.off()
pdf("./../Results/NewTreg/Seurat/Int_Dotplot.FACS_Cluster.pdf", width=21)
DotPlot(Int, features = c("Pdcd1", "Mki67", "Il2ra", "Cd69", "Cd44", "Havcr2", "Lag3", "Tigit"), dot.scale = 16) + RotatedAxis()
dev.off()

Int <- readRDS("./../Results/NewTreg/Seurat/Int_SeuratObj.rds")
Int$clonotype_id[Int$orig.ident == "Treg_1"] <- gsub("Iso", "Treg_1", Int$clonotype_id[Int$orig.ident == "Treg_1"])
Int$clonotype_id[Int$orig.ident == "Treg_2"] <- paste0("Treg_2_", Int$clonotype_id[Int$orig.ident == "Treg_2"])
Int$clonotype_id[Int$orig.ident == "Treg_3"] <- paste0("Treg_3_", Int$clonotype_id[Int$orig.ident == "Treg_3"])

# Divide clonotypes into multi-clonotype and single-clonotype. Then assign each cell to each type.
div_clonotype <- function(seurat_obj){
	# Extracting clonotype information
	clonotable <- table(seurat_obj$clonotype_id)
	clonotable <- clonotable[clonotable >= 1]
	multiclone <- names(clonotable[clonotable > 1])
	singleclone <- names(clonotable[clonotable == 1])

	# Making object
	clonotype <- rep("NA", ncol(seurat_obj))
	clonotype[which(seurat_obj$clonotype_id %in% multiclone)] <- "multi-clonotype"
	clonotype[which(seurat_obj$clonotype_id %in% singleclone)] <- "single-clonotype"

	# Calculating the number of each clones
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotable[as.character(seurat_obj$clonotype_id)], col.name = "cloneNum")
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotype, col.name = "clonotype")
	return(seurat_obj)
}

Int <- div_clonotype(Int)
saveRDS(Int, "./../Results/NewTreg/Seurat/Int_SeuratObj.rds")


# Visualizing chemokines
Chemokines <- read.delim("./../Data/Signature/FunctionalGeneset/Chemokines", sep='\t', header=F)[,1]
ChemokineReceptors <- read.delim("./../Data/Signature/FunctionalGeneset/ChemokineReceptors", header=F, sep='\t')[,1]

Chemokine.list <- list(Chemokines, ChemokineReceptors)
names(Chemokine.list) <- c("Chemokines", "ChemokineReceptors")

groups <- levels(Int$Cluster)
combs <- expand.grid(groups, groups)
combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
colnames(my.comps) <- NULL
rownames(my.comps) <- NULL
my.comps <- as.list(my.comps)

p <- list()
for(geneset in names(Chemokine.list)){
  genes <- Chemokine.list[[geneset]]
  
  i <- 1
  q <- list()
  for(gene in genes){
    if(!(gene %in% row.names(Int@assays$RNA@data))){next}
    VlnObj <- Int[row.names(Int) == gene,]
    VlnMat <- VlnObj@assays$RNA@data

		Vln.df <- as.data.frame(t(as.matrix(VlnMat)))
		Vln.df$Cluster <- Int$Cluster
		
		q[[i]] <- ggplot(Vln.df, aes_string(x="Cluster", y=gene)) + geom_violin() + geom_boxplot() +
		stat_compare_means(comparisons=my.comps, label = "p.signif", method='wilcox.test') +
		ggtitle(gene) +
		theme_bw() +
		theme(axis.title.x = element_blank(), axis.title.y = element_blank())

		i <- i + 1
  }
	p[[geneset]] <- q
}
pdf(paste0("./../Results/NewTreg/Seurat/Chemokines_Chemokines.pdf"), width=28, height=16)
print(grid.arrange(grobs = p[["Chemokines"]], ncol=7, left=grid::textGrob(gene)))
dev.off()
pdf(paste0("./../Results/NewTreg/Seurat/Chemokines_ChemokineReceptors.pdf"), width=24, height=16)
print(grid.arrange(grobs = p[["ChemokineReceptors"]], ncol=6, left=grid::textGrob(gene)))
dev.off()

for(geneset in names(Chemokine.list)){
  genes <- Chemokine.list[[geneset]]
  pdf(paste0("./../Results/NewTreg/Chmokines_", geneset, "_VolcanoPlot.pdf"), width=14)
  for(ct in levels(ChemokineCells$State)){
    p <- list()
     # Comparison between BLIA vs Rest
    VlnObj <- ChemokineCells[row.names(ChemokineCells) %in% genes, ChemokineCells$State == ct]
    VlnMat <- VlnObj@assays$RNA@data
    min.VlnMat <- min(VlnMat[VlnMat != 0])

    TNBCSubtype <- VlnObj$TNBCSubtype
    levels(TNBCSubtype) <- c("BLIA", "Other", "Other", "Other")

    DEG.df <- data.frame()
    for(gene in row.names(VlnMat)){
      wilcox.gene <- wilcox.test(VlnMat[gene,TNBCSubtype == "BLIA"], VlnMat[gene,TNBCSubtype == "Other"])
      logFC.gene <- log2((mean(VlnMat[gene,TNBCSubtype == "BLIA"]) + min.VlnMat) / (mean(VlnMat[gene,TNBCSubtype == "Other"]) + min.VlnMat))
      DEG.df <- rbind(DEG.df, cbind(gene, p.value=wilcox.gene$p.value, logFC.gene))
    }
    DEG.df$adj.p.value <- p.adjust(DEG.df$p.value)
    DEG.df$logFC.gene <- as.numeric(as.character(DEG.df$logFC.gene))
    p[[1]] <- ggplot(DEG.df, aes(x=logFC.gene, y=-log(adj.p.value))) + geom_point() +
      ggrepel::geom_text_repel(aes(x = logFC.gene, y = -log(adj.p.value), label = gene), max.overlaps=Inf, force=20, segment.color="grey") +
      ggtitle(paste0("Other <--> BLIA"))

    # Comparison between BLIA vs BLIS
    VlnObj <- VlnObj[,VlnObj$TNBCSubtype %in% c("BLIA", "BLIS")]
    VlnMat <- VlnObj@assays$RNA@data
    min.VlnMat <- min(VlnMat[VlnMat != 0])

    TNBCSubtype <- VlnObj$TNBCSubtype

    DEG.df <- data.frame()
    for(gene in row.names(VlnMat)){
      wilcox.gene <- wilcox.test(VlnMat[gene,TNBCSubtype == "BLIA"], VlnMat[gene,TNBCSubtype == "BLIS"])
      logFC.gene <- log2((mean(VlnMat[gene,TNBCSubtype == "BLIA"]) + min.VlnMat) / (mean(VlnMat[gene,TNBCSubtype == "BLIS"]) + min.VlnMat))
      DEG.df <- rbind(DEG.df, cbind(gene, p.value=wilcox.gene$p.value, logFC.gene))
    }
    DEG.df$adj.p.value <- p.adjust(DEG.df$p.value)
    DEG.df$logFC.gene <- as.numeric(as.character(DEG.df$logFC.gene))
    p[[2]] <- ggplot(DEG.df, aes(x=logFC.gene, y=-log(adj.p.value))) + geom_point() +
      ggrepel::geom_text_repel(aes(x = logFC.gene, y = -log(adj.p.value), label = gene), max.overlaps=Inf, force=20, segment.color="grey") +
      ggtitle(paste0("BLIS <--> BLIA"))
    
    print(grid.arrange(grobs = p, ncol=2, top = grid::textGrob(ct)))
  }
  dev.off()
}






pdf("./../Results/NewTreg/Seurat/Int_UMAP_Clonotype.pdf")
print(DimPlot(Int, reduction = "umap", group.by = "clonotype", pt.size = 1, cols = c("grey", "red"), order = c("multi-clonotype")))
dev.off()


clonality.df <- data.frame()
for(ct in levels(Int@active.ident)){
	for(sp in unique(Int$orig.ident)){
		tmp <- Int[,Int$orig.ident == sp]
		H <- asbio::alpha.div(table(tmp$clonotype_id[tmp@active.ident == ct]),index='shan')
		Hmax <- log(length(unique(tmp$clonotype_id[tmp@active.ident == ct])))
		clonality.df <- rbind(clonality.df, c(sp, ct, 1 - H/Hmax))
	}
}
colnames(clonality.df) <- c("Sample", "Celltype", "Clonality")
clonality.df$Celltype <- as.factor(clonality.df$Celltype)
clonality.df$Clonality <- as.numeric(clonality.df$Clonality)

pdf("./../Results/NewTreg/Seurat/Int_Clonality_Cluster.pdf")
	groups <- levels(clonality.df$Celltype)
	combs <- expand.grid(groups, groups)
	combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

	my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
	colnames(my.comps) <- NULL
	rownames(my.comps) <- NULL
	my.comps <- as.list(my.comps)

	print(
		ggplot(clonality.df, aes(x = Celltype, y = Clonality)) +
		geom_boxplot() +
		stat_compare_means(comparisons = my.comps, label="p.format",method='wilcox.test', method.args = list(alternative = "less")) +
		xlab("Cell types") +
		ylab("Clonality") +
		#scale_x_discrete(limits = c("Iso PD-1+ Treg", "aPD-1 PD-1+ Treg", "Iso PD-1- Treg", "aPD-1 PD-1- Treg")) +
		theme_cowplot()
	)
	print(
		ggplot(clonality.df, aes(x = Celltype, y = Clonality)) +
		geom_boxplot() +
		stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
		xlab("Clusters") +
		ylab("Clonality") +
		#scale_x_discrete(limits = c("Iso PD-1+ Treg", "aPD-1 PD-1+ Treg", "Iso PD-1- Treg", "aPD-1 PD-1- Treg")) +
		theme_cowplot()
	)
dev.off()

Gini.df <- data.frame()
for(ct in levels(Int@active.ident)){
	for(sp in unique(Int$orig.ident)){
		tmp <- Int[,Int$orig.ident == sp]
		Gini.df <- rbind(Gini.df, c(sp, ct, ineq::ineq(table(tmp$clonotype_id[tmp@active.ident == ct]), type='Gini')))
	}
}
colnames(Gini.df) <- c("Sample", "Celltype", "Gini")
Gini.df$Celltype <- as.factor(Gini.df$Celltype)
Gini.df$Gini <- as.numeric(Gini.df$Gini)

pdf("./../Results/NewTreg/Seurat/Int_GiniIndex_Cluster.pdf")
	groups <- levels(Gini.df$Celltype)
	combs <- expand.grid(groups, groups)
	combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

	my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
	colnames(my.comps) <- NULL
	rownames(my.comps) <- NULL
	my.comps <- as.list(my.comps)

	print(
		ggplot(Gini.df, aes(x = Celltype, y = Gini)) +
			geom_boxplot() +
			stat_compare_means(comparisons = my.comps, label="p.format",method='wilcox.test', method.args = list(alternative = "less")) +
			xlab("Cell types") +
			ylab("Gini Index") +
			#scale_x_discrete(limits = c("Iso PD-1+ Treg", "aPD-1 PD-1+ Treg", "Iso PD-1- Treg", "aPD-1 PD-1- Treg")) +
			theme_cowplot()
	)
	print(
		ggplot(Gini.df, aes(x = Celltype, y = Gini)) +
			geom_boxplot() +
			stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
			xlab("Cell types") +
			ylab("Gini Index") +
			theme_cowplot()
	)	
dev.off()


#######################################################################################3
# Trajectory Analysis
library(monocle)
library(reshape2)
library(Seurat)

# Reading the Seurat file.
SeuratObj <- Int

# Extracting information from Seurat object
Seur.exp <- as.matrix(SeuratObj@assays$RNA@counts)
Seur.metadata <- SeuratObj@meta.data
Seur.gene <- cbind(gene_short_name = row.names(SeuratObj@assays$RNA@meta.features), SeuratObj@assays$RNA@meta.features)

pd <- new("AnnotatedDataFrame", data = Seur.metadata)
fd <- new("AnnotatedDataFrame", data = Seur.gene)

# Making monocle data object
cds <- newCellDataSet(Seur.exp, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# Calculating additional data
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))
#cds <- cds[,pData(cds)$Total_mRNAs < 1e6]
#upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) + 2*sd(log10(pData(cds)$Total_mRNAs)))
#lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) - 2*sd(log10(pData(cds)$Total_mRNAs)))

pdf("./../Results/NewTreg/Monocle2/Int_Dist_RNA.pdf")
qplot(Total_mRNAs, data = pData(cds), geom = "density")
dev.off()

# I will not filter out cells.
#cds <- cds[,pData(cds)$Total_mRNAs > lower_bound & pData(cds)$Total_mRNAs < upper_bound]

# Finding expressed genes
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# Checking whether data follows standard distribution.
# Log transformation
L <- log(exprs(cds[expressed_genes,]))
# Standardization
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

pdf("./../Results/NewTreg/Monocle2/Int_Dist_RNA.log.pdf")
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(UMI)") + ylab("Density")
dev.off()


cds$Cluster <- SeuratObj@active.ident[colnames(cds)]

clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = '~Cluster + orig.ident', reducedModelFormulaStr = '~orig.ident', cores = 28)
saveRDS(clustering_DEG_genes, "./../Results/NewTreg/Monocle2/Int_Cluster_DEG.rds")
saveRDS(cds, "./../Results/NewTreg/Monocle2/Int_MonocleObj.rds")

feat.T <- c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Satb1", "Tnf", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

for (n in seq(100, 2000, 10)){
	message(paste0("Calculating gene number : ", n))
	#fl <- plyr::round_any(n, 100, f=floor)
	#ordering_genes <- clustering_DEG_genes_name[1:n]
  ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:n]
	cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
	cds <- reduceDimension(cds, method = 'DDRTree')
	cds <- orderCells(cds)

	traj <- t(cds@reducedDimS)
	colnames(traj) <- paste0("Trajectory_", 1:2)
	row.names(traj) <- colnames(cds)
	# We will now store this as a custom dimensional reduction called 'monocle2'
	SeuratObj[["trajectory"]] <- CreateDimReducObject(embeddings = as.matrix(traj), key = "Trajectory_", assay = DefaultAssay(SeuratObj))

	png(paste0("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_Feature.T_", n, ".png"), width=5600, height=4000, res=144)
	print(FeaturePlot(SeuratObj, features = feat.T, reduction = "trajectory", ncol = 7))
	dev.off()
}

n <- 740
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:n]
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
saveRDS(cds, "./../Results/NewTreg/Monocle2/Int_MonocleObj.rds")

pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_densityPeak_Cluster.pdf")
p <- plot_cell_trajectory(cds, color_by = "Cluster")
print(p)
dev.off()

pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_densityPeak_State.pdf")
p <- plot_cell_trajectory(cds, color_by = "State")
print(p)
dev.off()

traj <- t(cds@reducedDimS)
colnames(traj) <- paste0("Trajectory_", 1:2)
row.names(traj) <- colnames(cds)
# We will now store this as a custom dimensional reduction called 'monocle2'
Int[["trajectory"]] <- CreateDimReducObject(embeddings = as.matrix(traj), key = "Trajectory_", assay = DefaultAssay(Int))

pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_Sample.pdf", width = 8)
print(DimPlot(Int, reduction = "trajectory", group.by="orig.ident", label=F))
dev.off()


png("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_Feature.TF.png", width=4800, height=4800, res=144)
print(FeaturePlot(Int, features = feat.T, reduction = "trajectory", ncol = 6))
dev.off()

Int$State <- cds$State
Int$Cluster <- cds$Cluster
Int@active.ident <- Int$State
DEGs <- FindAllMarkers(Int, logfc.threshold = 0, min.pct = 0)

PlottingVolcano(DEGs, "./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_DEG.pdf", c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1","Lef1","Tbx1","Eomes","Prdm1","Bcl6","Id2","Id3","Gata3","Nfil3","Batf","Bcl11b","Stat3","Stat4","Stat5a","Irf4","Runx3","Bach2","Nr4a1","Foxo1","Zeb1","Zfp683","Klf2","Cebpb","Rela","Relb","Hopx","Rbpj","Nfkb2","Ikzf2","Myb","Junb","Runx1"))

cds <- orderCells(cds, root_state = 1) # You should change the number of state. Using plot_cell_trajectory(HSMM_myo, color_by = "State") command.

# Plotting both branches simultaneously
TF_genes <- row.names(subset(fData(cds), gene_short_name %in% feat.T))
pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_GBP.pdf", width = 13)
plot_genes_branched_pseudotime(cds[TF_genes,], branch_point = 1, color_by = "State", ncol = 4)
dev.off()




pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_Clonotype.pdf")
print(DimPlot(Int, reduction = "trajectory", group.by="clonotype", label=F))
dev.off()


clonality.df <- data.frame()
for(ct in levels(Int@active.ident)){
	for(sp in unique(Int$orig.ident)){
		tmp <- Int[,Int$orig.ident == sp]
		H <- asbio::alpha.div(table(tmp$clonotype_id[tmp@active.ident == ct]),index='shan')
		Hmax <- log(length(unique(tmp$clonotype_id[tmp@active.ident == ct])))
		clonality.df <- rbind(clonality.df, c(sp, ct, 1 - H/Hmax))
	}
}
colnames(clonality.df) <- c("Sample", "Celltype", "Clonality")
clonality.df$Celltype <- as.factor(clonality.df$Celltype)
clonality.df$Clonality <- as.numeric(clonality.df$Clonality)

pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_Clonality.pdf")
	groups <- levels(clonality.df$Celltype)
	combs <- expand.grid(groups, groups)
	combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

	my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
	colnames(my.comps) <- NULL
	rownames(my.comps) <- NULL
	my.comps <- as.list(my.comps)

	print(
		ggplot(clonality.df, aes(x = Celltype, y = Clonality)) +
		geom_boxplot() +
		stat_compare_means(comparisons = my.comps, label="p.format",method='wilcox.test', method.args = list(alternative = "less")) +
		xlab("Cell types") +
		ylab("Clonality") +
		#scale_x_discrete(limits = c("Iso PD-1+ Treg", "aPD-1 PD-1+ Treg", "Iso PD-1- Treg", "aPD-1 PD-1- Treg")) +
		theme_cowplot()
	)
	print(
		ggplot(clonality.df, aes(x = Celltype, y = Clonality)) +
		geom_boxplot() +
		stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
		xlab("Clusters") +
		ylab("Clonality") +
		#scale_x_discrete(limits = c("Iso PD-1+ Treg", "aPD-1 PD-1+ Treg", "Iso PD-1- Treg", "aPD-1 PD-1- Treg")) +
		theme_cowplot()
	)
dev.off()

Gini.df <- data.frame()
for(ct in levels(Int@active.ident)){
	for(sp in unique(Int$orig.ident)){
		tmp <- Int[,Int$orig.ident == sp]
		Gini.df <- rbind(Gini.df, c(sp, ct, ineq::ineq(table(tmp$clonotype_id[tmp@active.ident == ct]), type='Gini')))
	}
}
colnames(Gini.df) <- c("Sample", "Celltype", "Gini")
Gini.df$Celltype <- as.factor(Gini.df$Celltype)
Gini.df$Gini <- as.numeric(Gini.df$Gini)

pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_GiniIndex.pdf")
	groups <- levels(Gini.df$Celltype)
	combs <- expand.grid(groups, groups)
	combs <- combs[!duplicated(t(apply(combs, 1, sort))) & apply(combs, 1, function(x){x[1] != x[2]}),]

	my.comps <- as.data.frame(t(combs), stringsAsFactors=FALSE)
	colnames(my.comps) <- NULL
	rownames(my.comps) <- NULL
	my.comps <- as.list(my.comps)

	print(
		ggplot(Gini.df, aes(x = Celltype, y = Gini)) +
			geom_boxplot() +
			stat_compare_means(comparisons = my.comps, label="p.format",method='wilcox.test', method.args = list(alternative = "less")) +
			xlab("Cell types") +
			ylab("Gini Index") +
			#scale_x_discrete(limits = c("Iso PD-1+ Treg", "aPD-1 PD-1+ Treg", "Iso PD-1- Treg", "aPD-1 PD-1- Treg")) +
			theme_cowplot()
	)
	print(
		ggplot(Gini.df, aes(x = Celltype, y = Gini)) +
			geom_boxplot() +
			stat_compare_means(comparisons = my.comps, label="p.signif",method='wilcox.test', method.args = list(alternative = "less")) +
			xlab("Cell types") +
			ylab("Gini Index") +
			theme_cowplot()
	)	
dev.off()


MosaicData <- Int
MosaicData$State <- as.factor(as.character(MosaicData$State))
MetaInfo <- as.factor(as.character(MosaicData$Cluster))
ClusterInfo <- MosaicData$State
MosaicTable <- table(MetaInfo, ClusterInfo)

Mosaic <- vcd::structable(ClusterInfo ~ MetaInfo)
#Mosaic <- Mosaic[c(3, 14, 8, 11, 2, 12, 4, 9, 6, 5, 1, 13, 10, 7),]
pdf("./../Results/NewTreg/Monocle2/Int_Cluster_Trajectory.DDRTree_MosaicPlot.pdf", width=8)
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

saveRDS(Int, "./../Results/NewTreg/Seurat/Int_SeuratObj.rds")

# Compare between state 2 and state 3
State23 <- Int[,Int$State %in% c(2,3)]
State23@active.ident <- State23$State
DEGs <- FindAllMarkers(State23, logfc.threshold = 0, min.pct = 0)


PlottingVolcano(DEGs, "./../Results/NewTreg/Monocle2/State23_State_DEG.pdf", c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1","Lef1","Tbx1","Eomes","Prdm1","Bcl6","Id2","Id3","Gata3","Nfil3","Batf","Bcl11b","Stat3","Stat4","Stat5a","Irf4","Runx3","Bach2","Nr4a1","Foxo1","Zeb1","Zfp683","Klf2","Cebpb","Rela","Relb","Hopx","Rbpj","Nfkb2","Ikzf2","Myb","Junb","Runx1"))


# GSVA!!!
library(GSVA)
library(GSVAdata)
library(Seurat)
library(reshape2)
library(ggplot2)
library("RColorBrewer")
State23 <- ScaleData(State23, features = row.names(State23))

GS.MSigDB <- getGmt("./../Data/Signature/MSigDB/Hallmark.symbols_mouse.gmt")

gsva.MSigDB <- gsva(as.matrix(State23@assays$RNA@scale.data), GS.MSigDB, parallel.sz = 28)

gsva.tmp <- gsva.MSigDB
p <- list()
i <- 1
volcano.df <- data.frame()
for(term in row.names(gsva.tmp)){
	term.df <- as.data.frame(cbind(GSVA = gsva.tmp[term, ], celltype = as.character(State23$State)))
	term.df$GSVA <- as.numeric(as.character(term.df$GSVA))

  my_comparisons <- list(c("2", "3"))

	p[[i]] <- ggplot(term.df, aes(x = celltype, y = GSVA)) +
	  geom_violin(adjust = 3) +
	  geom_boxplot() +
	  ggtitle(term) +
	  xlab("State") +
	  ylab("GSVA score") +
    stat_compare_means(comparisons=my_comparisons, label = "p.signif")

	i <- i + 1
  volcano.df <- rbind(volcano.df, c(term, mean(term.df$GSVA[term.df$celltype == "2"]) - mean(term.df$GSVA[term.df$celltype == "3"]), wilcox.test(term.df$GSVA[term.df$celltype == "2"], term.df$GSVA[term.df$celltype == "3"])$p.value))
}
pdf("./../Results/NewTreg/Monocle2/State23_GSVA_MSigDB.scale_Violin.pdf", width = 36, height = 24)
grid.arrange(grobs = p, ncol = 9)
dev.off()

colnames(volcano.df) <- c("Term", "Diff", "pvalue")
volcano.df$Diff <- as.numeric(volcano.df$Diff)
volcano.df$pvalue <- as.numeric(volcano.df$pvalue)
volcano.df$Term <- str_split_fixed(volcano.df$Term, "[_]", 2)[,-1]

pdf("./../Results/NewTreg/Monocle2/State23_GSVA_MSigDB.scale_Volcano.pdf", width=11, height=11)
ggplot(volcano.df, aes(x=Diff, y=-log10(pvalue))) +
	geom_point() +
  xlim(c(-max(abs(volcano.df$Diff)), max(abs(volcano.df$Diff)))) +
	geom_text_repel(aes(x = Diff, y = -log10(pvalue), label = ifelse(pvalue < 0.05, as.character(Term),"")), max.overlaps=Inf, force=20, segment.color="grey") +
  geom_hline(yintercept=-log10(0.05))
dev.off()







pdf("tmp.pdf")
ggplot(tmp.df, aes(x = Celltype, y = Exp_Lgmn)) +
	  geom_violin(adjust = 3) +
	  geom_boxplot() +
	  xlab("Celltype") +
	  ylab("Expression of Lgmn") +
		theme_cowplot() +
		scale_x_discrete(limits=c("Iso_PD-1_WT_Treg", "Iso_PD-1_KO_Treg"))
dev.off()



#Plotting each branch separately.
cds.mem <- cds[,cds$State != 2]
cds.exh <- cds[,cds$State == 1 | cds$State == 2]

to_be_tested <- row.names(subset(fData(cds), gene_short_name %in% c("HAVCR2", "PDCD1", "KLRG1", "CD244", "CTLA4", "LAG3", "CD160", "TOX", "TIGIT", "ENTPD1", "ITGAE", "SELL", "IL7R")))
#mem_subset <- mono.mel.mem[to_be_tested,]
mem_subset <- cds.mem[to_be_tested,]
exh_subset <- cds.exh[to_be_tested,]

#diff_test_mem <- differentialGeneTest(mem_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_mem <- differentialGeneTest(mem_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_exh <- differentialGeneTest(exh_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")

#pdf(paste0("./../Results/Monocle/Pseudotime_State.", ngene, "_mem.pdf"), height = 30)
#p <- plot_genes_in_pseudotime(mem_subset, color_by = "State")
#print(p)
#dev.off()

pdf("./../Results/Monocle2/CD8.Seurat_Trajectory.Memory_Gene.pdf", height = 30)
p <- plot_genes_in_pseudotime(mem_subset, color_by = "State")
print(p)
dev.off()


pdf("./../Results/Monocle2/CD8.Seurat_Trajectory.Exhaustion_Gene.pdf", height = 30)
p <- plot_genes_in_pseudotime(exh_subset, color_by = "State")
print(p)
dev.off()

exh <- cds[,cds$State == 2]
exh.df <- pData(exh)

pdf("./../Results/Monocle2/CD8.Seurat_Trajectory.Exhaustion_Pseudotime_TNBCSubtype.pdf")
ggplot(exh.df, aes(x=TNBCSubtype, y=Pseudotime)) +
  geom_violin() +
  geom_boxplot(width=0.1)
dev.off()

mem <- cds[,cds$State != 1 & cds$State != 2]
mem.df <- pData(mem)

pdf("./../Results/Monocle2/CD8.Seurat_Trajectory.Memory_Pseudotime_TNBCSubtype.pdf")
ggplot(mem.df, aes(x=TNBCSubtype, y=Pseudotime)) +
  geom_violin() +
  geom_boxplot(width=0.1)
dev.off()

exh.df$Branch <- "Exhausted"
mem.df$Branch <- "Memory"

Branch.df <- rbind(mem.df, exh.df)










# Reading the Seurat file.
Treg <- Int[,Int$Celltype == "Treg_Norm"]

SeuratObj <- Treg

# Extracting information from Seurat object
Seur.exp <- as.matrix(SeuratObj@assays$RNA@counts)
Seur.metadata <- SeuratObj@meta.data
Seur.gene <- cbind(gene_short_name = row.names(SeuratObj@assays$RNA@meta.features), SeuratObj@assays$RNA@meta.features)

pd <- new("AnnotatedDataFrame", data = Seur.metadata)
fd <- new("AnnotatedDataFrame", data = Seur.gene)

# Making monocle data object
cds <- newCellDataSet(Seur.exp, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# Calculating additional data
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))
#cds <- cds[,pData(cds)$Total_mRNAs < 1e6]
#upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) + 2*sd(log10(pData(cds)$Total_mRNAs)))
#lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) - 2*sd(log10(pData(cds)$Total_mRNAs)))

pdf("./../Results/NewTreg/Monocle2/Treg_Dist_RNA.pdf")
qplot(Total_mRNAs, data = pData(cds), geom = "density")
dev.off()

# I will not filter out cells.
#cds <- cds[,pData(cds)$Total_mRNAs > lower_bound & pData(cds)$Total_mRNAs < upper_bound]

# Finding expressed genes
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# Checking whether data follows standard distribution.
# Log transformation
L <- log(exprs(cds[expressed_genes,]))
# Standardization
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

pdf("./../Results/NewTreg/Monocle2/Treg_Dist_RNA.log.pdf")
qplot(value, geom = "density", data = melted_dens_df) + stat_function(fun = dnorm, size = 0.5, color = 'red') + xlab("Standardized log(UMI)") + ylab("Density")
dev.off()


cds$Cluster <- SeuratObj@active.ident[colnames(cds)]

clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = '~Cluster + orig.ident', reducedModelFormulaStr = '~orig.ident', cores = 28)
saveRDS(clustering_DEG_genes, "./../Results/NewTreg/Monocle2/Treg_Cluster_DEG.rds")
saveRDS(cds, "./../Results/NewTreg/Monocle2/Treg_MonocleObj.rds")

for (n in seq(100, 2000, 100)){
	message(paste0("Calculating gene number : ", n))
	#ordering_genes <- clustering_DEG_genes_name[1:n]
  ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:n]
	cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
	cds <- reduceDimension(cds, method = 'DDRTree')
	cds <- orderCells(cds)

	pdf(paste0("./../Results/NewTreg/Monocle2/Treg_Cluster_Trajectory.DDRTree_densityPeak_Cluster_", n, ".pdf"))
	p <- plot_cell_trajectory(cds, color_by = "Cluster")
	print(p)
	dev.off()

	pdf(paste0("./../Results/NewTreg/Monocle2/Treg_Cluster_Trajectory.DDRTree_densityPeak_State_", n, ".pdf"))
	p <- plot_cell_trajectory(cds, color_by = "State")
	print(p)
	dev.off()

	traj <- t(cds@reducedDimS)
	colnames(traj) <- paste0("Trajectory_", 1:2)
	row.names(traj) <- colnames(cds)
	# We will now store this as a custom dimensional reduction called 'monocle2'
	SeuratObj[["trajectory"]] <- CreateDimReducObject(embeddings = as.matrix(traj), key = "Trajectory_", assay = DefaultAssay(SeuratObj))

	pdf(paste0("./../Results/NewTreg/Monocle2/Treg_Cluster_Trajectory.DDRTree_Sample_", n, ".pdf"), width = 8)
	print(DimPlot(SeuratObj, reduction = "trajectory", group.by="orig.ident", label=F))
	dev.off()

	feat.T <- c("Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Satb1", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Cd36", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")
	png(paste0("./../Results/NewTreg/Monocle2/Treg_Cluster_Trajectory.DDRTree_Feature.T_", n, ".png"), width=4800, height=4800, res=144)
	print(FeaturePlot(SeuratObj, features = feat.T, reduction = "trajectory", ncol = 6))
	dev.off()
}

n <- 1200
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:n]
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)










################################################################################
# 2-1-4 : Finding multiclonotype overlap...
Norm <- list()
for(l in 1:3){
	Norm[[l]] <- Int.list[[l]][,(Int.list[[l]]$Celltype == "Treg_Norm" | Int.list[[l]]$Celltype == "Mki67High" | Int.list[[l]]$Celltype == "Tconv_PD1High") & Int.list[[l]]$orig.ident == "Iso"]
}
saveRDS(Norm, "./../Results/NewTreg/Seurat/Norm.list_SeuratObj.rds")

for(l in 1:3){
	Norm[[l]] <- NormalizeData(object = Norm[[l]], verbose = FALSE)
	Norm[[l]] <- FindVariableFeatures(object = Norm[[l]], selection.method = "vst", verbose = FALSE, nfeatures = 2000)
	Norm[[l]]$orig.ident <- paste0("Norm_", l)
}

features <- SelectIntegrationFeatures(Norm)

#Using only 30 dimensions for integration
Norm.anchors <- FindIntegrationAnchors(Norm, dims = 1:30, anchor.features = features)
Norm.Int <- IntegrateData(anchorset = Norm.anchors, dims = 1:30)
Norm <- Norm.Int

DefaultAssay(Norm) <- "integrated"

# Run the standard workflow for visualization and clustering
#Int <- FindVariableFeatures(object = Int, selection.method = "vst", verbose = FALSE, nfeatures = 2000)

Norm <- ScaleData(Norm, verbose = FALSE)
Norm <- RunPCA(Norm, npcs = 100, verbose = FALSE)

Norm <- FindNeighbors(Norm, dims = 1:50)
Norm <- FindClusters(Norm, resolution = seq(0.1, 1, 0.1))

Norm <- RunUMAP(Norm, reduction = "pca", dims = 1:50)
# Changing the default assays into RNA!!!!!!!!!!
DefaultAssay(Norm) <- "RNA"

#Norm <- readRDS("./../Results/NewTreg/Seurat/Norm_SeuratObj.rds")
Norm$clonotype_id[Norm$orig.ident == "Norm_1"] <- gsub("Iso", "Norm_1", Norm$clonotype_id[Norm$orig.ident == "Norm_1"])
Norm$clonotype_id[Norm$orig.ident == "Norm_2"] <- paste0("Norm_2_", Norm$clonotype_id[Norm$orig.ident == "Norm_2"])
Norm$clonotype_id[Norm$orig.ident == "Norm_3"] <- paste0("Norm_3_", Norm$clonotype_id[Norm$orig.ident == "Norm_3"])
saveRDS(Norm, "./../Results/NewTreg/Seurat/Norm_SeuratObj.rds")


# Divide clonotypes into multi-clonotype and single-clonotype. Then assign each cell to each type.
div_clonotype <- function(seurat_obj){
	# Extracting clonotype information
	clonotable <- table(seurat_obj$clonotype_id)
	clonotable <- clonotable[clonotable >= 1]
	multiclone <- names(clonotable[clonotable > 1])
	singleclone <- names(clonotable[clonotable == 1])

	# Making object
	clonotype <- rep("NA", ncol(seurat_obj))
	clonotype[which(seurat_obj$clonotype_id %in% multiclone)] <- "multi-clonotype"
	clonotype[which(seurat_obj$clonotype_id %in% singleclone)] <- "single-clonotype"

	# Calculating the number of each clones
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotable[as.character(seurat_obj$clonotype_id)], col.name = "cloneNum")
	seurat_obj <- AddMetaData(seurat_obj, metadata = clonotype, col.name = "clonotype")
	return(seurat_obj)
}

Norm <- div_clonotype(Norm)
saveRDS(Norm, "./../Results/NewTreg/Seurat/Norm_SeuratObj.rds")


cut <- "
Int.Treg.multi.Ids <- unique(Int$clonotype_id[Int@active.ident == 2 & Int$clonotype == "multi-clonotype"])
Int <- AddMetaData(Int, metadata = Int$clonotype_id %in% Int.Treg.multi.Ids, col.name = "Treg_multi_over")
tmp <- as.factor(Int$Treg_multi_over)
Int$Treg_multi_over <- tmp

# Plotting bar plot
Int.df <- as.data.frame(cbind(as.character(Int@active.ident[Int$clonotype_id %in% Int.Treg.multi.Ids]), Int$orig.ident[Int$clonotype_id %in% Int.Treg.multi.Ids]))
colnames(Int.df) <- c("Cluster", "Group")

pdf("./../Results/Seurat/Treg_Bar_Treg.multi.over.pdf")
ggplot(Int.df, aes(Group, fill = Cluster)) + geom_bar(position = 'fill') + scale_x_discrete(limits = c("Iso", "aPD1"))
dev.off()
"

# Plotting overlap between clusters.
Norm$orig.ident <- as.factor(as.character(Norm$orig.ident))

Prop <- data.frame()
for(sp in levels(Norm$orig.ident)){
	SeuratObj <- Norm[,Norm$orig.ident == sp & Norm$clonotype == "multi-clonotype"]
	SeuratObj$Celltype <- as.factor(as.character(SeuratObj$Celltype))

	# Making data frame for all cell venn diagram
	SeuratObj.split <- split(SeuratObj$clonotype_id, SeuratObj$Celltype)

	venn.diagram(
	x = SeuratObj.split,
		main = "Multi clonotype ovelap",
		category.names = levels(SeuratObj$Celltype),
		filename = paste0("./../Results/NewTreg/Seurat/CloneOver_", sp, "_whole_VennDiagram.png"),
		imagetype="png",
		height = 800,
		width = 800,
		resolution = 200,
		compression = "lzw",
		cat.pos = 0,

		output=F
	)

	for(cell.1 in levels(SeuratObj$Celltype)){
		for(cell.2 in levels(SeuratObj$Celltype)){
			if(cell.1 >= cell.2){next}

			SeuratObj.cell.1 <- subset(SeuratObj, Celltype == cell.1)
			SeuratObj.cell.2 <- subset(SeuratObj, Celltype == cell.2)

			if(length(intersect(SeuratObj.cell.1$clonotype_id, SeuratObj.cell.2$clonotype_id)) == 0){ next;}

			clonotype.list <- list(SeuratObj.cell.1$clonotype_id, SeuratObj.cell.2$clonotype_id)
			clonotype.inter <- intersect(SeuratObj.cell.1$clonotype_id, SeuratObj.cell.2$clonotype_id)

			Prop <- rbind(Prop, c(sp, cell.1, cell.2, length(clonotype.inter) / length(unique(SeuratObj.cell.1$clonotype_id)) * 100, length(clonotype.inter) / length(unique(SeuratObj.cell.2$clonotype_id)) * 100))

			venn.diagram(
				x = clonotype.list,
				main = paste0("Multi clonotype ovelap\nbetween ", cell.1, " and ", cell.2, " from sample ", sp),
				category.names = c(cell.1, cell.2),
				filename = paste0("./../Results/NewTreg/Seurat/CloneOver_", sp, "_", cell.1, "n", cell.2, "_VennDiagram.png"),
				imagetype="png",
				height = 800,
				width = 800,
				resolution = 200,
				compression = "lzw",
				cat.pos = 0,

				output=F
			)
		}
	}
}
colnames(Prop) <- c("Sample", "Cell1", "Cell2", "Cell1_Prop", "Cell2_Prop")

Prop.cut <- Prop[Prop$Cell1 == "Mki67High" | Prop$Cell2 == "Mki67High",]
Prop.df <- melt(Prop.cut, c("Sample", "Cell1", "Cell2"))
Prop.df$value <- as.numeric(Prop.df$value)
levels(Prop.df$variable) <- c("Mki67High_Prop", "Other_Prop")

pdf(paste0("./../Results/NewTreg/Seurat/CloneOver_Boxplot.pdf"))
ggplot(Prop.df, aes(x=Cell2, y=value)) +
	facet_wrap(~variable) +
	geom_boxplot()
dev.off()



	# Plotting UMAP
	pdf("./../Results/Seurat/Int_UMAP_Treg.multi.over.pdf")
	DimPlot(Int, reduction = "umap", group.by = "Treg_multi_over", pt.size = 1, cols = c("grey", "red"), order = c("Overlapping"))
	dev.off()


