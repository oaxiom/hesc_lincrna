library(Seurat)
library(data.table)

load(file="human_escs.Rda")

# Can load in scanpy:
DefaultAssay(all_data) <- 'integrated'

cell.order <- colnames(all_data)
feature.order <- rownames(all_data)

# Save metadata:
meta.data <- all_data[[]][cell.order, ]
meta.data <- meta.data[c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'seurat_clusters')]
fwrite(data.table(data.frame("barcode"=rownames(meta.data),meta.data)), file='export/metadata.tsv',col.names=T, sep='\t')

# Save matrix:
counts <- data.frame(t(GetAssayData(all_data, slot='scale.data'))) # [feature.order, cell.order]))

norm_table = counts

norm_table = counts[,keep_genes]
norm_table = norm_table[keep_cells,]
keep_genes = feature.order[keep_genes]
keep_cells = cell.order[keep_cells]

table(keep_cells)
table(keep_genes)

fwrite(data.table(norm_table), file='export/matrix.tsv', sep='\t', col.names=F, row.names=F)
fwrite(data.table(keep_genes), file='export/gene_names.tsv', sep='\t', col.names=F)
fwrite(data.table(keep_cells), file='export/cell_names.tsv', sep='\t', col.names=F)



