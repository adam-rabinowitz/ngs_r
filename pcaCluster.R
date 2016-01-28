# Empty workspace
rm(list=ls())
# Load packages
require(ggplot2)
require(pheatmap)
require(fpc)
require(ggplot2)
require(Rtsne)
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
# Log data
log2TPM <- t(log2(filteredTPM + 1))
# Create function to cluster pca analysis
kmeans2DGroups <- function(X, outFile, krange=2:10, groups=NULL, seed=1) {
  # Check for dimensionality and rename columns
  if (ncol(X) != 2) {
    stop('Input data should have two columns')
  }
  X <- as.data.frame(X)
  colnames(X) <- c('x','y')
  # Process groups
  if (is.null(groups)) {
    X$group <- 1
  } else {
    X$group <- groups
  }
  X$group <- as.factor(X$group)
  # Extract and check row names
  sampleNames = row.names(X)
  if (is.null(sampleNames)) {
    sampleNames = 1:nrow(X)
  }
  # Create output objects
  outJaccard <- data.frame()
  outCluster <- data.frame(
    group = X$group,
    row.names = sampleNames
  )
  # Calcualte euclidean distance
  eDist <- dist(X)
  # Perform a range of kmeans clustering
  for (k in kRange) {
    # Perform boostrap of kmeans
    clb <- clusterboot(
      eDist,
      bootmethod = 'boot',
      B = 10,
      distances = T,
      multipleboot = F,
      clustermethod = kmeansCBI,
      k = k,
      seed = seed
    )
    # Add jaccard values to output
    jaccard <- rowMeans(clb$bootresult)
    jaccardDF <- data.frame(
      k = rep(k, k),
      jaccard = jaccard
    )
    outJaccard <- rbind(outJaccard, jaccardDF)
    # Add cluster data to output
    outCluster[,paste0('k',k)] <- as.factor(clb$partition)
  }
  outJaccard$k <- as.factor(outJaccard$k)
  # Open output file
  pdf(outFile, onefile = T, paper='a4r', width=8, height=7)
  # Create theme for plots
  t <- theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16),
    legend.title=element_text(size=16)
  )
  # Plot clustering
  p <- ggplot(X, aes(x=x, y=y, color=group)) + geom_point() + 
    t +
    labs(x='Dimension 1', y='Dimension 2')
  print(p)
  # Plot heatmap
  pheatmap(eDist)
  # Plot jaccard similarities
  p <- ggplot(outJaccard, aes(x=k, y = jaccard)) +
    geom_hline(aes(yintercept = c(0.5,0.75)), col = c('red','green')) +
    geom_point(shape = 20, size = 8) +
    t +
    ylim(0,1) + 
    stat_summary(fun.y=mean, geom="point", color="blue", shape = 20, size = 8)  
  print(p)
  # Plot contribution of each cell type to cluster
  for (k in kRange) {
    clID <- paste0('k',k)
    p <- ggplot(outCluster, aes_string(clID, fill='group')) +
      geom_bar(aes(fill=group)) +
      t
    print(p)
  }
  dev.off()
}
# Perform tsne
tsneData <- Rtsne(
  X = log2TPM,
  dims = 2,
  pca = FALSE,
  theta = 0.1,
  perplexity = 10
)
kmeans2DGroups()
# Perform pca analysis
pcaData <- prcomp(log2TPM)
# 

