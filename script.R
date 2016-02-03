rm(list=ls())
# Load required packages
require(ggdendro)
require(pheatmap)
# Load custom functions
setwd('~/github/ngs_r/')
source('RNASeq.R')
source('hclust.R')
# Read in RNA seq, filter and log transform
tpmFile = '/home/adam/sc_gene.tpm'
expFile = '/home/adam/sc_gene.exp'
filterTPM <- filterRNACounts(tpm=tpmFile, exp=expFile, minratio=0.25)
hc <- count2hclust(filterTPM)
# Create groups
groups <- list()
groups$'celltype' <- gsub('^S\\d\\.(P\\d).*?$','\\1',colnames(filterTPM))
names(groups$'celltype') <- colnames(filterTPM)
groups$'replicate' <- gsub('^(S\\d).*?$', '\\1', colnames(filterTPM))
names(groups$'replicate') <- colnames(filterTPM)
# Open 
pdf('/home/adam/test.pdf', onefile=T, height=7, width=10, paper='a4r')
pheatmap(
  log2(filterTPM+1),
  clustering_distance_rows='euclidean',
  clustering_distance_cols='euclidean',
  clustering_method='complete',
  show_rownames=F,
  show_colnames=F
)
plotGroupDendro(hc, groups)
dev.off()
pheatmap(log2Filter)

require(Rtsne)
tsneData <- Rtsne(
  X = t(log2(filterTPM+1)),
  dims=2,
  pca=F,
  perplexity=10,
  theta=0.1
)


