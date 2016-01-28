require('ggplot2')
##############################################################################
## Generate distance matrix from pca data
##############################################################################
pcaDist <- function(pca, method='euclidean' , pc=c(1,2)) {
  pcData <- pca$x[,pc]
  dist <- dist(pcData, method=method)
  return(dist)
}

##############################################################################
## Check concordane between pca data and group data
##############################################################################
checkGroups <- function(pca, groups) {
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
plotPCA <- function(pca, outFile, groups=NA) {
  # Extract pca data as data frame
  pcDF <- as.data.frame(pca$x)
  # Check concordance between PC data and groups if supplied
  groups <- checkGroups(pca, groups)
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


p <- pcaDist(pcaData)
plotPCA(pcaData,outFile,groups=groups)