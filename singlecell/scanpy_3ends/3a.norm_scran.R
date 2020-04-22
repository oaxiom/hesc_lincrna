
library('data.table')
library('scran')

setwd('/Users/andrew/Projects/linc_te/hesc_lincrna/singlecell/scanpy_3ends')

data_mat = data.matrix(fread('dense_array.tsv', sep='\t'))
clusters <- quickCluster(data_mat)
size_factors = calculateSumFactors(data_mat, clusters=clusters, min.mean=0.1)
fwrite(data.frame(size_factors), file='size_factors.csv')
