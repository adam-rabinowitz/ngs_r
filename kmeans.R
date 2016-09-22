require(fpc)
source('/farm/home/rabino01/github/ngs_r/ggplot/plotThemes.R')

##############################################################################
## Extract clusters and their stability scores for different values of K
##############################################################################
# Create method to help determine kmeans clustering
kmeansKChoice <- function(X, krange=2:10, seed=1) {
  # Check for dimensionality and rename columns
  if (ncol(X) != 2) {
    stop('Input data should have two columns')
  }
  X <- as.data.frame(X)
  colnames(X) <- c('x','y')
  # Extract and check row names
  sampleNames = row.names(X)
  if (is.null(sampleNames)) {
    sampleNames = 1:nrow(X)
  }
  # Create output objects
  outJaccard <- data.frame()
  outCluster <- data.frame(
    row.names = sampleNames
  )
  # Calcualte euclidean distance
  eDist <- dist(X)
  # Perform a range of kmeans clustering
  for (k in krange) {
    # Perform boostrap of kmeans
    clb <- clusterboot(
      eDist,
      bootmethod = 'boot',
      B = 1000,
      distances = T,
      multipleboot = F,
      clustermethod = kmeansCBI,
      k = k,
      seed = seed
    )
    # Add jaccard values to output
    jaccardDF <- data.frame(
      k = rep(k, k),
      cluster = 1:k,
      jaccard = clb$bootmean
    )
    outJaccard <- rbind(outJaccard, jaccardDF)
    # Add cluster data to output
    outCluster[,paste0('k',k)] <- as.factor(clb$partition)
  }
  outJaccard$k <- as.factor(outJaccard$k)
  outList <- list(
    'coordinates' = X,
    'jaccard' = outJaccard,
    'clusters'=outCluster
  )
  return(outList)
}

##############################################################################
## Plot kmeans data and cluster stability
##############################################################################
kmeansChoicePlot <- function(X, outpdf) {
  # Extract plot data
  plotData <- cbind(X$coordinates, X$clusters)
  # Open pdf and create unclustered plot
  pdf(file = outpdf, paper = 'a4r', width = 8, height =7, onefile = T)
  p <- ggplot(
    plotData,
    aes_string('x', 'y')
  ) +
    geom_point() +
    pdfA4SquareTheme
  print(p)
  # Loop through clusters
  for (clust in colnames(X$clusters)) {
    print(clust)
    # extract value of k
    kValue <- as.numeric(substr(clust,2,10))
    jValues <- subset(X$jaccard, k == kValue)
    jValues <- jValues[order(jValues$cluster),]
    jValues$jaccard <- round(jValues$jaccard, 2)
    # Extract jaccard cluster values
    p <- ggplot(
      plotData,
      aes_string('x', 'y', color = clust)
    ) +
    geom_point() +
    ggtitle(clust) +
    scale_color_discrete(
      name = 'Cluster\nJaccard',
      labels = paste(
        jValues$cluster,
        sprintf('%.2f', jValues$jaccard),
        sep = ' - '
      )
    ) +
    pdfA4SquareTheme
    print(p)
  }
  dev.off()
}

##############################################################################
## Plot how supplied groupings are places in clusters
##############################################################################
kmeansGroupPlot <- function(X, groups, outpdf) {
  # Extract plot data
  plotData <- cbind(X$coordinates, X$clusters)
  # Add groups to plot data
  if (!identical(names(groups), row.names(plotData))) {
    stop('group names do not match kmeans data')
  }
  plotData$group <- as.factor(groups)
  # Open pdf
  pdf(file = outpdf, paper = 'a4r', width = 8, height =7, onefile = T)
  # Create scatter plot
  p <- ggplot(
    plotData,
    aes_string('x', 'y', color = 'group')
  ) +
    geom_point() +
    ggtitle('Group Distribution') +
    pdfA4SquareTheme
  print(p)
  # Loop through clusters
  print('clust')
  for (clust in colnames(X$clusters)) {
    # Create plot
    p <- ggplot(
      plotData,
      aes_string(clust, fill = 'group')
    ) +
    geom_bar() +
    ggtitle(clust) +
    pdfA4SquareTheme
    print(p)
  }
  dev.off()
}



