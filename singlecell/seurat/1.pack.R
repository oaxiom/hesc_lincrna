library(Seurat)
library(ggplot2)

setwd('/Users/andrew/Projects/linc_te/hesc_lincrna/singlecell/seurat')

data1 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.C11Solo.out/"), project="primed_c11_ipsc")
data2 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.WIBR3PriSolo.out"), project="primed_wibr3_esc")
data3 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam1Solo.out"), project="primed_wtc_ipsc_rp1")
data4 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam3Solo.out"), project="primed_wtc_ipsc_rp3")
data5 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam4Solo.out"), project="primed_wtc_ipsc_rp4")
data6 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam5Solo.out"), project="primed_wtc_ipsc_rp5")

# Filter and integrate:
data1 <- subset(x=data1, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)
data2 <- subset(x=data2, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)
data3 <- subset(x=data3, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)
data4 <- subset(x=data4, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)
data5 <- subset(x=data5, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)
data6 <- subset(x=data6, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)

sams = list(data1, data2, data3, data4, data5, data6)

sams <- lapply(sams, NormalizeData)

anchors <- FindIntegrationAnchors(object.list=sams, dims=1:20)
all_data = IntegrateData(anchorset=anchors, dims=1:20)

DefaultAssay(all_data) <- "integrated"

all_data <- ScaleData(all_data, verbose=FALSE)
all_data <- RunPCA(all_data, npcs=100, verbose = FALSE)
all_data <- RunUMAP(all_data, reduction = "pca", dims = 1:50)
all_data <- FindNeighbors(all_data, reduction = "pca", dims = 1:50)
all_data <- FindClusters(all_data, resolution = 0.5)

# Quick QC:
pdf("QC.pdf", width=18, height=6)
ElbowPlot(all_data, ndims = 100)
VlnPlot(object=all_data, features=c("nCount_RNA"), ncol=1, pt.size=0, group.by='orig.ident')
VlnPlot(object=all_data, features=c("nFeature_RNA"), ncol=1, pt.size=0, group.by='orig.ident')
dev.off()

pdf("learning.pdf", width=6, height=6)
DimPlot(all_data, pt.size=1, reduction="umap", group.by='orig.ident')
DimPlot(all_data, pt.size=1, reduction="umap", label=T)
DimHeatmap(all_data, dims=1:9, cells=500, balanced=TRUE)
DimHeatmap(all_data, dims=41:50, cells=500, balanced=TRUE)
FeaturePlot(esc, pt.size=1, features = c("POU5F1", "SOX2", "NANOG", "LIN28A"), min.cutoff = "q4")
#PC1:
FeaturePlot(esc, pt.size=1, features = c("MALAT1", "MTATP6P1", "SNHG14", "XACT"), min.cutoff = "q4")
FeaturePlot(esc, pt.size=1, features = c("RPL18A", "RPS2", "RPL15", "RPS17"), min.cutoff = "q4")
#PC2:
FeaturePlot(esc, pt.size=1, features = c("LINC01356", "ESRG", "POU5F1", "TDGF1"), min.cutoff = "q4")
FeaturePlot(esc, pt.size=1, features = c("LINC01158", "CKS2", "NR2F2", "PTN"), min.cutoff = "q4")
#PC3:
FeaturePlot(esc, pt.size=1, features = c("ENO1", "PGK1", "ALDOA", "PKM"), min.cutoff = "q4")
FeaturePlot(esc, pt.size=1, features = c("SNHG19", "CKS2", "SNHG7", "AC106864.1"), min.cutoff = "q9")
FeaturePlot(esc, pt.size=1, features = c("HNRNPU", "HNRNPK", "BRCA1", "BRCA2"), min.cutoff = "q9")

DimHeatmap(esc, dims=1:9, cells=500, balanced=TRUE)
dev.off()

#VlnPlot(object=esc, features=c("nFeature_RNA", "nCount_RNA"), ncol=2,  group.by='orig.ident')
save(esc, file="human_escs.Rda")

# And export to loom:
hescs <- as.loom(esc, filename="human_escs.loom")
hescs
hescs$close_all()
