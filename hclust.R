require(snow)
require(pheatmap)
require(dendextend)

##############################################################################
## Calculate distcance of count matrix
##############################################################################
count2Dist <- function(tpm, method='euclidean', log=T) {
  if (isTRUE(log)) {
    tpm = t(log2(tpm + 1))
  } else {
    tpm = t(tpm)
  }
  d <- dist(tpm, method=method)
  return(d)
}

##############################################################################
## Perform hierachical clustering on count matrix
##############################################################################
count2hclust <- function(tpm, dmethod='euclidean', cmethod='complete', log=T){
  # Calculate distance
  d <- count2Dist(tpm, method=dmethod, log=log)
  # Find clusters and return
  hc <- hclust(d, method=cmethod)
  return(hc)
}

##############################################################################
## Plot dendro with groups
##############################################################################
plotDendro <- function(hc, groups=NULL) {
  # Calculate y-limits for plot
  ylim <- c(0, max(hc$height))
  # convert hclust object to dendrogram
  dend <- as.dendrogram(hc, ylim=ylim, hang=-1)
  # Create plots for groups
  if (!is.null(groups)) {
    # Check groups
    groups <- checkGroups(labels(dend), groups)
    # Create colours for plots 
    colours <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65,
      h.start = 0, direction = 1)
    for (i in 1:length(groups)) {
      # Extract group data
      groupName <- names(groups)[i]
      groupFactor <- groups[[i]]
      # Convert groups to numeric and generate colours
      groupLevels <- as.numeric(groupFactor)
      colVector = colours(max(groupLevels))
      # Format dendrogram
      groupColours <- colVector[groupLevels]
      dend <- branches_color(dend, clusters=groupLevels)
      # Create plot
      plot(
        dend,
        xlab=paste('Distance:', hc$dist.method),
        sub=paste('Cluster:', hc$method),
        leaflab='none',
        main='Hierarchical Clustering'
      )
      legend(
        'topright',
        legend = levels(groupFactor),
        col=colVector,
        lty=1,
        lwd=3,
        title=groupName,
        bg='white'
      )
    }
  } else {
    # Create plot
    plot(
      dend,
      xlab=paste('Distance:', hc$dist.method),
      sub=paste('Cluster:', hc$method),
      leaflab='none',
      main='Hierarchical Clustering'
    )
  }
}

##############################################################################
## Perform pvclust
##############################################################################
clusterConfidence <- function(tpm, dmethod='euclidean', cmethod='complete',
  log=T, nboot=100, threads=2) {
  # Convert to log if required
  if (log) {
    tpm = log2(tpm +1)
  }
  # Create cluster
  cl <- makeCluster(threads, type='SOCK')
  # Perform pvclust
  pv.tpm <- pvclust(
    tpm,
    method.hclust = cmethod,
    method.dist = dmethod,
    nboot = nboot, 
    parallel = cl
  )
  # stop cluster
  stopCluster(cl)
  # Return results
  return(pv.tpm)
}

##############################################################################
## Plot pvclust
##############################################################################
plotPVClust <- function(pvcData, alpha=0.95) {
  plot(
    pvcData$hclust,
    labels = F,
    xlab = paste('Distance:', pvcData$hclust$dist.method),
    sub = paste('Clustering:', pvcData$hclust$method),
    main = 'Raw Clusters'
  )
  plot(
    pvcData,
    labels = F,
    xlab = paste('Distance:', pvcData$hclust$dist.method),
    sub = paste('Clustering:', pvcData$hclust$method),
    main = 'Cluster P-Value'
  )
  plot(
    pvcData$hclust,
    labels = F,
    xlab = paste('Distance:', pvcData$hclust$dist.method),
    sub = paste('Clustering:', pvcData$hclust$method),
    main = paste('Cluster P-Value >=', alpha)
  )
  pvrect(pvcData, alpha=alpha)
}

##############################################################################
## heatmap
##############################################################################
drawHeat <- function(tpm, log=T, dmethod='euclidean', cmethod='complete') {
  # Log counts if requested
  if (isTRUE(log)) {
    tpm <- log2(tpm + 1)
  }
  pheatmap(
    tpm,
    clustering_distance_rows = dmethod,
    clustering_distance_cols = dmethod,
    clustering_method = cmethod,
    show_rownames = F,
    show_colnames = F
  )
}

##############################################################################
## Complete cluster analysis
##############################################################################
heatClusterCombo <- function(tpm, log=T, dmethod='euclidean',
  cmethod='complete', groups=NULL) {
  # Draw heatmap
  drawHeat(tpm, log=log, dmethod=dmethod, cmethod=cmethod)
  # Perform hierachical clustering and plot
  hc <- count2hclust(tpm, log=log, dmethod=dmethod, cmethod=cmethod)
  plotDendro(hc, groups=groups)
}

##############################################################################
## Cut tree
##############################################################################
extractClusters <- function(tpm, k, minsize=1, log=T, dmethod='euclidean',
  cmethod='complete') {
  # Perform hierachical clustering and plot
  hc <- count2hclust(tpm, log=log, dmethod=dmethod, cmethod=cmethod)
  # Store the required and max value for k
  requiredK <- k
  # Find clusters
  while (TRUE) {
    # Check cluster number
    clusters <- cutree(hc, k=k)
    clustTable <- table(clusters)
    clustCount <- sum(clustTable >= minsize)
    # Find clusters of sufficient size
    if (clustCount == requiredK) {
      # Extract clsters of the required size
      acceptedClusters <- names(clustTable[clustTable >= minsize])
      accpetedClusters <- as.numeric(acceptedClusters)
      clusters <- clusters[clusters %in% acceptedClusters]
      # Factorise clusters
      clusterFactor <- factor(clusters, levels=unique(clusters))
      # Reorganise factor and return
      returnVector <- as.numeric(clusterFactor)
      names(returnVector) <- names(clusters)
      return(returnVector)
    # Else increase k size
    } else {
      k <- k + 1
    }
    if (k > ncol(tpm)) {
      warning('No clusters of required size identified')
      return(NULL)
    }
  }
}



