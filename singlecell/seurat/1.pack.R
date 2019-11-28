library(Seurat)

setwd('/Users/andrew/Projects/linc_te/hesc_lincrna/singlecell/seurat')

data1 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.C11Solo.out/"), project="primed_c11_ipsc")
data2 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.WIBR3PriSolo.out"), project="primed_wibr3_esc")
data3 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam1Solo.out"), project="primed_wtc_ipsc#1")
data4 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam3Solo.out"), project="primed_wtc_ipsc#3")
data5 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam4Solo.out"), project="primed_wtc_ipsc#4")
data6 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam5Solo.out"), project="primed_wtc_ipsc#5")

sams = list(data1, data2, data3, data4, data5, data6)

sams <- lapply(sams, subset, subset=nFeature_RNA > 2000 & nFeature_RNA < 5000)
sams <- lapply(sams, NormalizeData)

sams <- FindIntegrationAnchors(object.list=sams, anchor.features = 4000, dims=1:20)
all_data = IntegrateData(anchorset=sams, dims=1:20)

DefaultAssay(all_data) <- "integrated"
all_data <- FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 4000)
all_data <- ScaleData(all_data)

# Quick QC:
# get number of cells kept:
cell.num <- table(all_data$orig.ident)

pdf("QC.pdf", width=10, height=6)
op <- par(mar = c(5,14, 4,2) + 0.1)
barplot(cell.num,horiz=T, main='Cell numbers', space=0.4, las=2)
par(op)
VlnPlot(object=all_data, features=c("nCount_RNA"), ncol=1, pt.size=0, group.by='orig.ident')
VlnPlot(object=all_data, features=c("nFeature_RNA"), ncol=1, pt.size=0, group.by='orig.ident')
dev.off()

# Learning
all_data <- RunPCA(all_data, npcs=100)

pdf("jackstraw.pdf", width=10, height=6)
all_data <- JackStraw(all_data, num.replicate = 100)
all_data <- ScoreJackStraw(all_data, dims = 1:20)
JackStrawPlot(all_data, dims = 1:20)
ElbowPlot(all_data, ndims = 100)
dev.off()

# This data still has some pretty bad batch effects, a lot of which are low down the PCS, so I use a lot less than the 
# ElbowPlot and JackStraw suggest: This matches the expected variance in this dataset which should be not really so high
# And matches closer to tSNE seen in the Powell dataset alone.
all_data <- RunUMAP(all_data, reduction="pca", dims=1:20)
all_data <- RunTSNE(all_data, reduction='pca', dims=1:20)
all_data <- FindNeighbors(all_data, reduction="pca", dims=1:20)
all_data <- FindClusters(all_data, resolution=0.2)

pdf('pca.pdf', width=12, height=12)
DimPlot(all_data, dims=1:2, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=3:4, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=5:6, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=7:8, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=9:10, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=11:12, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=13:14, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=15:16, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=17:18, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimPlot(all_data, dims=19:20, pt.size=0.4, reduction="pca", group.by='orig.ident')
DimHeatmap(all_data, dims=1:9, cells=500, balanced=TRUE)
DimHeatmap(all_data, dims=10:19, cells=500, balanced=TRUE)
DimHeatmap(all_data, dims=20:29, cells=500, balanced=TRUE)
DimHeatmap(all_data, dims=30:39, cells=500, balanced=TRUE)
DimHeatmap(all_data, dims=39:40, cells=500, balanced=TRUE)
dev.off()

pdf("learning.pdf", width=12, height=12)
all_data <- RunTSNE(all_data, reduction='pca', dims=1:5)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:5)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:10)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:10)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:15)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:15)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:20)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:20)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:25)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:25)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:30)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:30)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:35)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:35)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:40)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:40)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:45)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:45)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
all_data <- RunTSNE(all_data, reduction='pca', dims=1:50)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:50)
DimPlot(all_data, pt.size=0.5, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.5, reduction="umap", group.by='orig.ident')
dev.off()

all_data <- RunTSNE(all_data, reduction='pca', dims=1:20)
all_data <- RunUMAP(all_data, reduction="pca", dims=1:20)

pdf("learning-final.pdf", width=13, height=12)
DimPlot(all_data, pt.size=1, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=1, reduction="umap", group.by='orig.ident')
DimPlot(all_data, pt.size=1, reduction="tsne")
DimPlot(all_data, pt.size=1, reduction="umap")
dev.off()

pdf('pc_genes.pdf', width=12, height=12)
FeaturePlot(all_data, pt.size=0.1, features = c("POU5F1", "SOX2", "NANOG", "LIN28A"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
FeaturePlot(all_data, pt.size=0.1, features = c("L1TD1", "HSP90AB1", "IGFBP2", "FTH1"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
#PC1:
FeaturePlot(all_data, pt.size=0.1, features = c("MALAT1", "MTATP6P1", "SNHG14", "XACT"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
FeaturePlot(all_data, pt.size=0.1, features = c("RPL18A", "RPS2", "RPL15", "RPS17"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
#PC2:
FeaturePlot(all_data, pt.size=0.1, features = c("LINC01356", "ESRG", "POU5F1", "TDGF1"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
FeaturePlot(all_data, pt.size=0.1, features = c("LINC01158", "CKS2", "NR2F2", "PTN"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
#PC3:
FeaturePlot(all_data, pt.size=0.1, features = c("ENO1", "PGK1", "ALDOA", "PKM"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
FeaturePlot(all_data, pt.size=0.1, features = c("SNHG19", "CKS2", "SNHG7", "AC106864.1"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))
FeaturePlot(all_data, pt.size=0.1, features = c("HNRNPU", "HNRNPK", "BRCA1", "BRCA2"), min.cutoff="0", max.cutoff='1.5', cols=c('purple', "lightgrey", "orange"))

# cluster-specific markers:
# Clust1:
FeaturePlot(all_data, pt.size=0.1, features = c("PIF1", "UBE2S", "RPL39", "RPS29"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust2:
FeaturePlot(all_data, pt.size=0.1, features = c("ENO1", "LDHB", "PPIA", "PRDX1"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust2:
FeaturePlot(all_data, pt.size=0.1, features = c("STMN1", "ANXA5", "EPCAM", "ILF2"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust3:
FeaturePlot(all_data, pt.size=0.1, features = c("HSP90AA1", "NPM1", "LDHA", "HDAC2"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust4:
FeaturePlot(all_data, pt.size=0.1, features = c("ID3", "PMAIP1", "SNHG7", "AL118516.1"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust5:
FeaturePlot(all_data, pt.size=0.1, features = c("BUB3", "AURKB", "CDK1", "CDC20"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust6:
FeaturePlot(all_data, pt.size=0.1, features = c("IGFBP2", "TOMM7", "SNHG6", "FUS"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust7:
FeaturePlot(all_data, pt.size=0.1, features = c("SDCCAG8", "LINC01356", "GAL", "SDCCAG8"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
# Clust8:
FeaturePlot(all_data, pt.size=0.1, features = c("PDIA3", "CENPF", "SEMA6A", "LEFTY1"), min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
dev.off()

viol_and_umap <- function(x) {
  print(x)
  plot1 = VlnPlot(object=all_data, features=x, y.min=0, y.max=3, pt.size=0)
  plot2 = FeaturePlot(object=all_data, features=x, pt.size=0.02, min.cutoff="0", reduction="umap", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
  plot3 = FeaturePlot(object=all_data, features=x, pt.size=0.02, min.cutoff="0", reduction="tsne", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
  CombinePlots(plots = list(plot1, plot2, plot3), ncol=3)
}

# Patterns from the Powell paper:
genes = list('LEFTY2', 'OTX2', 'GDF3', 'DNMT3B', 'DPPA5', 'UTF1', 'NODAL', 'DPPA2', 'KLF4', 'SDC2', 'ZFP42', 'HESX1',
             'EGFR', 'PDGFA', 'ESR2', 'ENO1', 'NMNAT1', 'NETO2')
pdf('powell.pdf', width=13, height=3)
par(mfrow=c(1,1))
lapply(genes, viol_and_umap)
dev.off()

# Naive Primed
genes = list('ESRRB',  'TFCP2L1', 'ZFP42', 'MT1H', 'DPPA3', 'DPPA4', 'DPPA5', 'ZNF486', 'CR1L', 'DNMT3L', 'HNRNPU', 'HNRNPK')
pdf('naive.pdf', width=13, height=3)
par(mfrow=c(1,1))
lapply(genes, viol_and_umap)
dev.off()

# Marker genes; 
all_data.markers <- FindAllMarkers(all_data, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25)
genes = all_data.markers['gene'][,1]
pdf('violins.pdf', width=13, height=3)
lapply(genes, viol_and_umap)
dev.off()
write.table(all_data.markers, file='marker_genes.tsv', sep='\t')

save(all_data, file="human_escs.Rda")

