##############################################################################
## Plot kmeans clsutering
##############################################################################
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