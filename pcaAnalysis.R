require('ggplot2')
require('Rtsne')

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
## Generate distance matrix from pca data
##############################################################################
pcaDist <- function(pca, dmethod='euclidean' pc=c(1,2)) {
  pcData <- pca$x[,pc]
  dist <- dist(pcData, method=dmethod)
  return(dist)
}

##############################################################################
## Function to plot pca data
##############################################################################
plotPCA <- function(pca, groups=NULL) {
  # Check concordance between PC data and groups if supplied
  print(names(groups))
  if (!is.null(groups)) {
    groups <- checkGroups(row.names(pca$x), groups)
  }
  print(names(groups))
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
    legend.title=element_text(size=16),
    title=element_text(size=18,face='bold')
  )
  # Create Variance plot
  p <- ggplot(varDF, aes(x = PC, y=VAR)) +
    geom_bar(stat='identity') +
    xlab('Principal Component') +
    ylab('Proportion Of Variance') +
    ylim(0,max(varDF$VAR)) +
    thm
  print(p)
  # Create PC1 vs PC2 plot
  if (is.list(groups)) {
    for (i in 1:length(groups)) {
      # Extract group data
      groupName = names(groups)[i]
      groupFactor = groups[[i]]
      # Assemble plot data
      plotDF <- as.data.frame(pca$x)[,1:2]
      plotDF[,groupName] <- groupFactor
      # Create plot
      fig <- ggplot(plotDF, aes_string(x = 'PC1', y = 'PC2', color=groupName)) +
        geom_point() +
        thm +
        ggtitle('PCA Analysis - PC1 vs PC2')
      # Print plot and remove
      print(fig)
      rm(fig)
    }
  } else {
    # Assemble plot data
    plotDF <- as.data.frame(pca$x)[,1:2]
    p <- ggplot(plotDFDF, aes(x = PC1, y = PC2)) +
      thm +
      geom_point() +
      ggtitle('PCA Analysis - PC1 vs PC2')
    print(p)
  }
  # Create density plot for top principal components
  for (PCn in paste0('PC',1:10)) {
    # Create plots with group data
    if (is.list(groups)) {
      for (i in 1:length(groups)) {
        # Extract group data
        groupName = names(groups)[i]
        groupFactor = groups[[i]]
        # Create plot data
        plotDF <- as.data.frame(pca$x[,PCn])
        colnames(plotDF) <- PCn
        plotDF[,groupName] = groupFactor
        # Create plot
        fig <- ggplot(plotDF, aes_string(x = PCn, color=groupName)) +
          thm +
          geom_density() +
          ggtitle(paste(PCn, 'Distribution')) +
          ylab('Density')
        print(fig)
        rm(fig)
      }
    # Create plots without groups data
    } else {
      for (i in 1:length(groups)) {
        # Create plot
        p <- ggplot(pcDF, aes_string(x = PCn)) +
          thm +
          geom_density() +
          ggtitle(paste(PCn, 'Distribution')) +
          ylab('Density')
        print(p)
      }
    }
  }
}

pdf('/home/adam/test.pdf',onefile=T,paper='a4r')
plotPCA(pca, groups)
dev.off()

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