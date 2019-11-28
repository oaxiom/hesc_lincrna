library(Seurat)
library(data.table)

load(file="human_escs.Rda")

# Can load in scanpy:
DefaultAssay(all_data) <- 'integrated'
cell.order <- colnames(all_data)
feature.order <- rownames(all_data)

DefaultAssay(all_data) <- 'RNA'
feature.order.all <- rownames(all_data)

# Save metadata:
meta.data <- all_data[[]][cell.order, ]
meta.data <- meta.data[c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'seurat_clusters')]
fwrite(data.table(data.frame("barcode"=rownames(meta.data),meta.data)), file='export/metadata.tsv',col.names=T, sep='\t')

# Save matrix:
DefaultAssay(all_data) <- 'integrated'
counts <- data.frame(t(GetAssayData(all_data, assay='integrated', slot='scale.data'))) # [feature.order, cell.order]))
raw_counts <- data.frame(t(as.matrix(GetAssayData(all_data, assay="RNA", slot='data'))))

norm_table = counts

# Make sure in correcto order?
#norm_table = counts[,keep_genes]
#norm_table = norm_table[keep_cells,]
#keep_genes = feature.order[keep_genes]
#keep_cells = cell.order[keep_cells]
#table(keep_cells)
#table(keep_genes)

fwrite(data.table(raw_counts), file='export/raw_matrix.tsv', sep='\t', col.names=F, row.names=F)
fwrite(data.table(feature.order.all), file='export/raw_gene_names.tsv', sep='\t', col.names=F)
fwrite(data.table(norm_table), file='export/matrix.tsv', sep='\t', col.names=F, row.names=F)
fwrite(data.table(feature.order), file='export/gene_names.tsv', sep='\t', col.names=F)
fwrite(data.table(cell.order), file='export/cell_names.tsv', sep='\t', col.names=F)



