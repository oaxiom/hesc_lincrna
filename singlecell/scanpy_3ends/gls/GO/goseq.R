
library(goseq)
library(geneLenDataBase)
library(org.Hs.eg.db)
library(GO.db)
#library(TxDb.Mmusculus.UCSC.mm10.ensGene)

setwd("/Users/andrew/Projects/linc_te/hesc_lincrna/singlecell/scanpy_3ends/gls/GO")

filenames = list.files(path=".", pattern="*.tsv$")

#all = read.delim("all_genes.tsv")$ensg
all = read.delim("../../../../gencode/hg38_ensembl_v95_ensg.tsv")$ensg

for (filename in filenames) {
  print(filename)
  
  data = read.delim(filename, row.names=1) # id for deseq
  ensg = rownames(data)
  
  vec = as.integer(all%in%ensg)
  names(vec) = all
  
  pdf(paste("barplots.", filename, ".pdf", sep=""), width=6, height=1.8)
  par(mar=c(4, 14, 3, 3)) 
  
  #pwf = nullp(vec, "mm10", "ensGene")#, bias.data=rep(100, length(vec)))
  # make a pwf, with no bias data
  pwf = data.frame(DEgenes=vec, bias.data=rep(1, length(all)), pwf=rep(1, length(all)), row.names=all)                
  GO.sig = goseq(pwf, "hg38", "ensGene", method='Wallenius', use_genes_without_cat=FALSE) # Hypergeometric
  trimmedGO = GO.sig[GO.sig$numInCat > 20 & GO.sig$numInCat < 300,]# & GO.sig$numDEInCat>=4,]
  
  #trimmedGO$over_represented_pvalue = p.adjust(trimmedGO$over_represented_pvalue, method="BH")+1e-40
  enriched = trimmedGO[trimmedGO$over_represented_pvalue <0.01,] # q-value
  #enriched[, "Ontology"] = Ontology(enriched$category)
  
  # Build a table for a barplot:
  
  for (ont in c("CC", "BP", "MF")) {
    this_ont = enriched[enriched$ontology == ont,]
    if (dim(this_ont)[1] > 1) {
      tab = data.frame(p=this_ont$over_represented_pvalue, numInCat=this_ont$numInCat, row.names=paste(this_ont$category, " ", this_ont$term, sep=""))
      write.table(tab, paste("tabs/tab", ont, '_', filename, sep=""), sep="\t", col.names=NA)
    }
    # Do the same for the barcharts
    this_ont = enriched[enriched$ontology == ont,]
    if (dim(this_ont)[1] > 1) { # It is a bug that the barcharts cannot deal with only a single sig entry... Not clear how to fix
      tab = c("GO", "q")
      for (r in 1:min(100, dim(this_ont)[1])) {
        pp = paste(Term(this_ont[r,]$category), " (", this_ont[r,]$category, ")", sep="")
        tab = rbind(tab, c(pp, this_ont[r,]$over_represented_pvalue)) 
      }
      
      colnames(tab) = tab[1,]
      tab = (tab[-1,])
      tab = data.frame(q=as.numeric(tab[,2]), row.names=tab[,1])
      tab = -log10(tab)
      tab = tab[order(tab[,1], decreasing=T),,drop=F]
      
      par(las=2, mar=c(2,27,1,0.5))#, mgp=c(3,1,0)) # perpendicular and margins
      barplot(as.numeric(rev(tab[1:20,1])), names.arg=rev(rownames(tab)[1:20]), horiz=T, main=ont, cex.axis=0.7, 
              space=0.4, cex=0.5)
    }
  }
  dev.off()
}

