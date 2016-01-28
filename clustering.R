rm(list=ls())
# Load packages
require(ggplot2)
require(pheatmap)
require(fpc)
require(ggplot2)
# Set values for analysis
tpmFile = '/farm/scratch/rs-bio-lif/rabino01/bonet/singleCell576Analysis/sc_gene.tpm'
expFile = '/farm/scratch/rs-bio-lif/rabino01/bonet/singleCell576Analysis/sc_gene.exp'
minAligned = 4e5
minExpressed = 2e3
minTPM = 10
minTPMCount = 5
kRange <- 2:10
# Read in data
tpmMatrix = read.table(tpmFile)
expMatrix = read.table(expFile)
# Extract accpeted samples
acceptedSamples <- colnames(expMatrix)[
  colSums(expMatrix) >= minAligned &
  apply(expMatrix, 2, function(z) {sum(z > 0)}) >= minExpressed
]
# Extract accepted genes
acceptedGenes <- row.names(tpmMatrix)[
  apply(
    tpmMatrix[,acceptedSamples],
    1,
    function(z) {
      sum(z >= minTPM)
    }
  ) >= minTPMCount
]
# Extract tpm counts
filteredTPM <- tpmMatrix[acceptedGenes,acceptedSamples]
# Generate dissimilarity matrix
processDist <- function(tm, outFile) {
  sampleNames = colnames(tm)
  tmDist <- dist(t(tm))
  # Find jaccard simmilarity values for range of k
  outJaccard <- data.frame()
  outCluster <- data.frame(row.names = sampleNames)
  outCluster$celltype <- gsub('^S\\d\\.(P\\d)\\.\\w\\d+','\\1', sampleNames)
  for (k in kRange) {
    clb <- clusterboot(
      tmDist,
      bootmethod = 'boot',
      B = 100,
      distances = T,
      multipleboot = F,
      clustermethod = kmeansCBI,
      k = k,
      seed = 1979
    )
    jaccard <- rowMeans(clb$bootresult)
    jaccardDF <- data.frame(
      k = rep(k, k),
      jaccard = jaccard
    )
    outJaccard <- rbind(outJaccard, jaccardDF)
    outCluster[,paste0('k',k)] <- as.factor(clb$partition)
  }
  outJaccard$k <- as.factor(outJaccard$k)
  # Plot heatmap
  pdf(outFile, onefile = T, paper='a4r', width=9, height=7)
  pheatmap(tmDist)
  # Plot jaccard similarities
  p <- ggplot(outJaccard, aes(x=k, y = jaccard)) +
    geom_point(shape = 20, size = 8) +
    geom_hline(aes(yintercept = c(0.5,0.75)), col = c('red','green')) +
    ylim(0,1) + 
    stat_summary(fun.y=mean, geom="point", color="blue", shape = 20, size = 8) +
    theme(axis.text=element_text(size=16),
      axis.title=element_text(size=18,face="bold"))
  print(p)
  # Plot contribution of each cell type to cluster
  for (k in kRange) {
    clID <- paste0('k',k)
    p <- ggplot(outCluster, aes_string(clID, fill='celltype')) +
      geom_bar(aes(fill=celltype)) +
      theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) 
    print(p)
  }
  dev.off()
}
# Create plots
processDist(filteredTPM,'/farm/home/rabino01/bonet/tpmCluster.pdf')
log2filterTPM <- log2(filteredTPM + 1)
processDist(log2filterTPM, '/farm/home/rabino01/bonet/log2tpmCluster.pdf')
log2tpmVar <- apply(log2filterTPM, 1, var)
mostVar <- names(log2tpmVar[order(-log2tpmVar)])[1:1000]
processDist(log2filterTPM[mostVar,], '/farm/home/rabino01/bonet/log2tpmVarCluster.pdf')
