require('ggplot2')
require('Rtsne')

# check clustereval

##############################################################################
## filter RNA-Seq counts
##############################################################################
filterRNACounts <- function(tpm, exp, reads=400000, genes=2000, mintpm=10,
  minratio=0.1, minsample=5) {
  # Read in data
  message('Reading TPM count file')
  tpmMatrix = read.table(tpm)
  message('Reading expected count file')
  expMatrix = read.table(exp)
  # Extract samples with sufficient aligned reads and expressed genes
  acceptedSamples <- colnames(expMatrix)[
    colSums(expMatrix) >= reads &
      apply(expMatrix, 2, function(z) {sum(z > 0)}) >= genes
    ]
  message(paste(length(acceptedSamples), 'of', ncol(tpmMatrix),
      'sample accepted'))
  # Calculate number of samples in which a gene must be expresed
  sampleNo <- ceiling(max(minsample, length(acceptedSamples) * minratio))
  message(paste('Finding Genes expressed in >=', sampleNo, 'samples'))
  # Extract accepted genes
  acceptedGenes <- row.names(tpmMatrix)[
    apply(
      tpmMatrix[,acceptedSamples], 1, function(z) {sum(z >= mintpm)}
    ) >= sampleNo
  ]
  message(paste(length(acceptedGenes), 'of', nrow(tpmMatrix), 'accepted'))
  # Extract and return filtered counts
  filteredTPM <- tpmMatrix[acceptedGenes,acceptedSamples]
  return(filteredTPM)
}

##############################################################################
## Perform PCA
##############################################################################
performPCA <- function(tpm, log=T) {
  # Manipulate data
  if (log) {
    tpm <- t(log2(tpm+1))
  } else {
    tpm <- t(tpm)
  }
  # Perform PCA and return data
  pca <- prcomp(tpm)
  return(pca)
}

##############################################################################
## Generate euclidean distance matrix from pca data
##############################################################################
pcaDist <- function(pca, pc=c(1,2)) {
  pcData <- pca$x[,pc]
  dist <- dist(pcData, method='euclidean')
  return(dist)
}

##############################################################################
## Check concordane between pca data and group data
##############################################################################
checkPCAGroups <- function(pca, groups) {
  # Check group data is character class
  if (!is.character(groups)) {
    stop('"groups" must be a character vector') 
  }
  # Check group variable name length
  gname = sort(names(groups))
  pname = sort(row.names(pca$x))
  # Check group variable length
  if (!identical(gname,pname)) {
    stop('Group names and sample names are not identical')
  }
  # Reorder group, turn to factor and return
  groups <- groups[row.names(pca$x)]
  groups <- factor(groups, levels = unique(groups))
  return(groups)
}

##############################################################################
## Function to plot pca data
##############################################################################
plotPCA <- function(pca, outFile, groups=NULL) {
  # Extract pca data as data frame
  pcDF <- as.data.frame(pca$x)
  # Check concordance between PC data and groups if supplied
  if (!is.null(groups)) {
    groups <- checkPCAGroups(pca, groups)
  }
  # Extract variance data
  varRatio <- (pca$sdev^2) / sum(pca$sdev^2)
  varDF <- data.frame(
    'PC' = factor(paste0('PC', 1:10), levels=paste0('PC', 1:10)),
    'VAR' = varRatio[1:10]
  )
  # Create theme
  thm <- theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16),
    legend.title=element_text(size=16)
  )
  # Open PDF file
  pdf(file = outFile, onefile=T, paper='a4r', width=8, height=7)
  # Create Variance plot
  p <- ggplot(varDF, aes(x = PC, y=VAR)) +
    geom_bar(stat='identity') +
    xlab('Principal Component') +
    ylab('Proportion Of Variance') +
    ylim(0,max(varDF$VAR)) +
    thm
  print(p)
  # Create PC1 vs PC2 plot
  p <- ggplot(pcDF, aes(x = PC1, y = PC2)) +
    thm
  if (is.factor(groups)) {
    p <- p + geom_point(aes(color=groups))
  } else {
    p <- p + geom_point()
  }
  print(p)
  # Create density plot for first ten principal components
  for (PCn in paste0('PC',1:10)) {
    p <- ggplot(pcDF, aes_string(x = PCn)) +
      thm +
      ylab('Density')
    if (is.factor(groups)) {
      p <- p + geom_density(aes(color=groups))
    } else {
      p <- p + geom_density()
    }
    geom_density(aes(color=groups))
    print(p)
  }
  dev.off()
}

##############################################################################
## Cluster PCA data
##############################################################################
clustPCA <- function(pca, pc=c(1:2), groups=NULL) {
  # Extract pca data as data frame
  pcDF <- as.data.frame(pca$x)
  # Check concordance between PC data and groups if supplied
  if (!is.null(groups)) {
    groups <- checkGroups(pca, groups)
  }
  # Create distance matrix
  pcaDist <- dist(pcDF[,pc])
  # 
}

##############################################################################
## Plot tsne of multiple principal components
##############################################################################
tsnePlot <- function(pca, pc=c(1,2), perplexity = 10, groups=NULL) {
  # Check concordance between PC data and groups if supplied
  if (!is.null(groups)) {
    groups <- checkGroups(pca, groups)
  }
  # Perform tsne
  pcad <- pcaDist(pca, pc=pc)
  tsneData <- Rtsne(pcad, theta=0, perplexity=perplexity,  pca=F)
  plot(tsneData$Y)
}
tsnePlot(pcaData)


p <- pcaDist(pcaData)
plotPCA(pcaData,outFile,groups=groups)