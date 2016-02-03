require(snow)

##############################################################################
## Calculate distcance of count matrix
##############################################################################
count2Dist <- function(tpm, method='euclidean', log=T) {
  if (log) {
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
count2hclust <- function(tpm, dmethod='euclidean', cmethod='complete', log=T) {
  # Calculate distance
  d <- count2Dist(tpm, method=dmethod, log=log)
  # Find clusters and return
  hc <- hclust(d, method=cmethod)
  return(hc)
}

##############################################################################
## Plot dendro with groups
##############################################################################
plotGroupDendro <- function(hc, groups) {
  # convert hclust object to dendrogram and order  groups
  dend <- as.dendrogram(hc)
  groups <- checkGroups(labels(dend), groups)
  # Calculate y-limits for plot
  ylim <- c(min(hc$height), max(hc$height))
  # Create colours for plots 
  colours <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65,
    h.start = 0, direction = 1)
  # Create plots for each groups
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
    dend <- hang.dendrogram(dend, 0.1)
    # Create plot
    plot(
      dend,
      ylim=ylim,
      xlab=paste('Distance:', hc$dist.method),
      sub=paste('Cluster:', hc$method),
      bty='l',
      leaflab='none'
    )
    legend(
      'topright',
      legend = levels(groupFactor),
      col=colVector,
      lty=1,
      lwd=3,
      title=groupName
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
## Complete cluster analysis
##############################################################################


