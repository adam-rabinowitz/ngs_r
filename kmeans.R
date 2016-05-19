require(fpc)
##############################################################################
## Plot kmeans clsutering
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
      B = 100,
      distances = T,
      multipleboot = F,
      clustermethod = kmeansCBI,
      k = k,
      seed = seed
    )
    # Add jaccard values to output
    jaccardDF <- data.frame(
      k = rep(paste0('k',k), k),
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

kmeansKChoicePlot <- function(X, outpdf) {
  # Extract plot data
  plotData <- cbind(X$coordinates, X$clusters)
  # Open pdf
  pdf(file = outpdf, paper = 'a4r', width = 8, height =7, onefile = T)
  # Loop through clusters
  for (clust in colnames(X$clusters)) {
    print(clust)
    # extract value of k
    jValues <- subset(X$jaccard, k == clust)
    jValues <- jValues[order(jValues$cluster),]
    jValues$jaccard <- round(jValues$jaccard, 2)
    # Extract jaccard cluster values
    kplot <- ggplot(
      plotData,
      aes_string('x', 'y', color = clust)
    ) +
    geom_point() +
    ggtitle(clust) +
    scale_color_discrete(
      name = 'Cluster\nJaccard',
      labels = paste(
        jValues$cluster,
        jValues$jaccard,
        sep = ' - '
      )
    )
    print(kplot)
  }
  dev.off()
}


