# =============================================================================
# Reanalysis of public CD4 T cells (GSE152022)
# Manuscript: "TOX enforces the immunosuppressive program of tumor-infiltrating
#              regulatory T cells" (Park et al., Nature Immunology)
# Figures: Fig 1c-f
# Environment: scRNA env (R 4.0.5)   (exact versions in ../envs/ and README.md)
# Provenance: original project script "Analyzing_scRNAseq_GSE152022.R" — analysis logic UNCHANGED; only this
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
Data1.raw <- Read10X("./../Data/PublicData/GSE152022/GSM4598898/")
Data2.raw <- Read10X("./../Data/PublicData/GSE152022/GSM4598899/")

Data1 <- CreateSeuratObject(counts = Data1.raw, project = "Tumor", min.cells = 3, min.features = 200)
Data2 <- CreateSeuratObject(counts = Data2.raw, project = "Normal", min.cells = 3, min.features = 200)

Data1[["percent.mt"]] <- PercentageFeatureSet(Data1, pattern = "^mt-")
Data2[["percent.mt"]] <- PercentageFeatureSet(Data2, pattern = "^mt-")


plot1 <- FeatureScatter(Data1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Data1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Data2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Data2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


Data1 <- subset(Data1, subset = nFeature_RNA > 200 & percent.mt < 5)
Data2 <- subset(Data2, subset = nFeature_RNA > 200 & percent.mt < 5)

Data.list <- list(Data1, Data2)
Data.list <- lapply(X = Data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Data.list)

Data.anchors <- FindIntegrationAnchors(object.list = Data.list, anchor.features = features)

Data.intd <- IntegrateData(anchorset = Data.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Data.intd) <- "integrated"

# Run the standard workflow for visualization and clustering
Data.intd <- ScaleData(Data.intd, verbose = FALSE)
Data.intd <- RunPCA(Data.intd, npcs = 30, verbose = FALSE)
Data.intd <- RunUMAP(Data.intd, reduction = "pca", dims = 1:30)
Data.intd <- FindNeighbors(Data.intd, reduction = "pca", dims = 1:30)
Data.intd <- FindClusters(Data.intd, resolution = seq(0.1, 1, 0.1))

DefaultAssay(Data.intd) <- "RNA"

feat <- c("Cd4", "Cd8", "Tox", "Pdcd1", "Foxp3", "Havcr2", "Tigit", "Lag3", "Ctla4", "Anxa5", "Mki67", "Il10", "Tgfb1", "Tgfb2", "Tcf7", "Tbx21", "Prdm1", "Ifng", "Tnf", "Il2", "Il2ra", "Eomes", "Gata3", "Nfil3", "Batf", "Bcl11b", "Lef1", "Id2", "Ezh2", "Gzmb", "Tnfrsf4", "Entpd1", "Ifit3", "Cd5", "Infg", "Cd44", "Cd69", "Tnfrsf9", "Ccr8", "Il1rl1")

png(paste0("./../Results/GSE152022/Intd.2000_UMAP_Feature.png"), width = 6400, height = 4000, res = 144)
print(FeaturePlot(Data.intd, reduction = "umap", features = feat, ncol = 8))
dev.off()

pdf(paste0("./../Results/GSE152022/Intd.2000_Sample.pdf"))
print(DimPlot(Data.intd, group.by = "orig.ident"))
dev.off()

pdf("./../Results/GSE152022/Intd.2000_UMAP_Clusters.pdf")
for(j in seq(0.1, 1, 0.1)){
	print(DimPlot(Data.intd, reduction = "umap", label = T, pt.size = 0.5, group.by = paste0("integrated_snn_res.", j)))
}
dev.off()

Data.intd@active.ident <- as.factor(Data.intd$integrated_snn_res.0.4)

pdf("./../Results/GSE152022/Intd.2000_UMAP_Cluster.pdf")
print(DimPlot(Data.intd,label = T))
dev.off()

saveRDS(Data.intd, "./../Results/GSE152022/Intd.2000_SeuratObj.rds")

png(paste0("./../Results/GSE152022/Intd.2000_VlnPlot_Feature.png"), width = 6400, height = 4000, res = 144)
VlnPlot(Data.intd, features = feat, ncol=8)
dev.off()


# For paper
feat <- c( "Foxp3", "Il2ra", "Pdcd1", "Tox")
png(paste0("./../Results/GSE152022/Intd.2000_VlnPlot_Feature.Int.png"), width = 3200, height = 800, res = 144)
VlnPlot(Data.intd, features = feat, ncol=4)
dev.off()

feat <- c("Cd4","Foxp3","Gata3","Nfil3","Batf","Lef1","Tcf7","Id2","Il2ra","Tox","Pdcd1","Havcr2","Tigit","Ctla4","Lag3","Tnfrsf4","Tnfrsf9","Tnfrsf18","Entpd1","Il1rl1","Ccr8","Cd44","Il10","Ccr7","Il7r","Sell")

pdf("./../Results/GSE152022/Intd.2000_Dotplot_Feature.pdf", width=14, height=4)
DotPlot(Data.intd, features = feat) + RotatedAxis()
dev.off()







