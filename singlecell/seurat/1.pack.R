library(Seurat)

setwd('/Users/andrew/Projects/linc_te/hesc_lincrna/singlecell/seurat')

data1 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.C11Solo.out/"), project="primed_c11_ipsc")
data2 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.WIBR3PriSolo.out"), project="primed_wibr3_esc")
data3 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam1Solo.out"), project="primed_wtc_ipsc_rp1")
data4 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam3Solo.out"), project="primed_wtc_ipsc_rp3")
data6 = CreateSeuratObject(counts=Read10X(data.dir="../rawdata/ss.hesc_posam5Solo.out"), project="primed_wtc_ipsc_rp5")

sams = list(data1, data2, data3, data4, data6)

sams <- lapply(sams, subset, subset=nFeature_RNA > 1000 & nFeature_RNA < 5000)
sams <- lapply(sams, NormalizeData)

sams <- FindIntegrationAnchors(object.list=sams, dims=1:20)
all_data = IntegrateData(anchorset=sams, dims=1:20)

DefaultAssay(all_data) <- "integrated"

all_data <- FindVariableFeatures(all_data, selection.method = "vst", nfeatures = 2000)

all_data <- ScaleData(all_data)

# Quick QC:
pdf("QC.pdf", width=10, height=6)
VlnPlot(object=all_data, features=c("nCount_RNA"), ncol=1, pt.size=0, group.by='orig.ident')
VlnPlot(object=all_data, features=c("nFeature_RNA"), ncol=1, pt.size=0, group.by='orig.ident')
dev.off()

# Learning
all_data <- RunPCA(all_data, npcs=100)

pdf("jackstraw.pdf", width=10, height=6)
all_data <- JackStraw(all_data, num.replicate = 100)
all_data <- ScoreJackStraw(all_data, dims = 1:30)
JackStrawPlot(all_data, dims = 1:20)
ElbowPlot(all_data, ndims = 100)
dev.off()

all_data <- RunUMAP(all_data, reduction="pca", dims=1:20)
all_data <- RunTSNE(all_data, reduction='pca', dims=1:20)
all_data <- FindNeighbors(all_data, reduction="pca", dims=1:20)
all_data <- FindClusters(all_data, resolution=0.2)

pdf("learning.pdf", width=12, height=12)
DimPlot(all_data, pt.size=0.4, reduction="tsne", group.by='orig.ident')
DimPlot(all_data, pt.size=0.4, reduction="tsne", label=T)
DimPlot(all_data, pt.size=0.4, reduction="umap", group.by='orig.ident')
DimPlot(all_data, pt.size=0.4, reduction="umap", label=T)
DimHeatmap(all_data, dims=1:9, cells=500, balanced=TRUE)
DimHeatmap(all_data, dims=32:40, cells=500, balanced=TRUE)
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
  plot1 = VlnPlot(object=all_data, features=x, y.min=0, y.max=3, pt.size=0, group.by='orig.ident')
  plot2 = FeaturePlot(object=all_data, features=x, pt.size=0.02, min.cutoff="0", max.cutoff='1.5',cols=c('purple', "lightgrey", "orange"))
  CombinePlots(plots = list(plot1, plot2))
}

# Marker genes; 
all_data.markers <- FindAllMarkers(all_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
genes = all_data.markers['gene'][,1]
pdf('violins.pdf', width=8, height=3)
lapply(genes, viol_and_umap)
dev.off()
write.table(all_data.markers, file='marker_genes.tsv', sep='\t')


# Patterns from the Powell paper:
genes = list('LEFTY2', 'OTX2', 'GDF3', 'DNMT3B', 'DPPA5', 'UTF1', 'NODAL', 'DPPA2', 'KLF4', 'SDC2', 'ZFP42', 'HESX1',
             'EGFR', 'PDGFA', 'ESR2', 'ENO1', 'NMNAT1', 'NETO2')
pdf('powell.pdf', width=8, height=3)
par(mfrow=c(1,1))
lapply(genes, viol_and_umap)
dev.off()

# Naive Primed
genes = list('ESRRB',  'TFCP2L1', 'ZFP42', 'MT1H', 'DPPA3', 'DPPA4', 'DPPA5', 'ZNF486', 'CR1L', 'DNMT3L', 'HNRNPU', 'HNRNPK')
pdf('naive.pdf', width=8, height=3)
par(mfrow=c(1,1))
lapply(genes, viol_and_umap)
dev.off()

save(all_data, file="human_escs.Rda")

# And export to loom: All these conversions are currently broken:
#hescs <- as.loom(all_data, filename="human_escs.loom")
#hescs
#hescs$close_all()
#pfile <- Convert(from=all_data, to="loom", filename="human_escs.loom", display.progress=T)
#pbmc.sce <- as.SingleCellExperiment(pbmc)

# Can load in scanpy:
DefaultAssay(all_data) <- 'RNA'

cell.order <- colnames(all_data)
feature.order <- rownames(all_data)

# Save metadata:
meta.data <- all_data[[]][cell.order, ]
meta.data <- meta.data[c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'seurat_clusters')]
write.table(meta.data, file='metadata.tsv', sep='\t')

# Save matrix:
counts <- as.matrix(t(as.data.frame(GetAssayData(all_data, assay="RNA", slot='counts')))) # [feature.order, cell.order]))

keep_genes <- colSums(counts) > 200  # Filter out empty genes:
keep_cells <- rowSums(counts) > 1000
norm_table = counts[keep_genes]
norm_table = norm_table[keep_cells]
keep_genes = feature.order[keep_genes]
keep_cells = cell.order[keep_cells]

table(keep_cells)
table(keep_genes)

write.table(norm_table, file='matrix.tsv', sep='\t', col.names=F, row.names=F)
write.table(keep_genes, file='gene_names.tsv', sep='\t')
write.table(keep_cells, file='cell_names.tsv', sep='\t')

# Get the metadata of the samples:


