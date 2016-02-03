# Load required packages
require(ggplot2)
require(ggdendro)
require(amap)
require(pheatmap)
require
# Create function to draw dendrogram
filterTPM <- filterRNACounts(tpm=tpmFile, exp=expFile, minratio=0.2)
log2Filter <- log2(filterTPM + 1)
dlog2 <- dist(t(log2Filter))
hc <- hclust(dlog2)
plot(hc, label=F)



# Calculate groups
groups <- gsub('^S\\d\\.(P\\d).*?$','\\1',colnames(filterTPM))
names(groups) <- colnames(filterTPM)
# Create distance object
pcaData <- performPCA(filterTPM, log=T)
plotPCA(pcaData, outFile = '/home/rabino01/test.pdf', groups=groups)


##############################################################################
## Create distance
##############################################################################
calcDist <- function(tpm, dmethod='euclidean', log=T) {
  if (log) {
    tpm = t(log2(tpm + 1))
  } else {
    tpm = t(tpm)
  }
  d <- Dist(tpm, method=dmethod)
}

##############################################################################
## Check groups for tpm matrices
##############################################################################
checkGroups <- function(tpm, groups) {
  # Check group data is character class
  if (!is.character(groups)) {
    stop('"groups" must be a character vector') 
  }
  # Check group variable name length
  gname = sort(names(groups))
  pname = sort(colnames(tpm))
  # Check group variable length
  if (!identical(gname,pname)) {
    stop('Group names and sample names are not identical')
  }
  # Reorder group, turn to factor and return
  groups <- groups[colnames(tpm)]
  groups <- factor(groups, levels = unique(groups))
  return(groups)
} 

##############################################################################
## Calculate robustness of clusters for varying k
##############################################################################
plotHClustK <- function(tpm, outfile, groups=NULL, dmethod='euclidean',
  krange=2:10) {
  # Check groups if supplied
  if (!is.null(groups)) {
    groups <- checkGroups(tpm, groups)
  }
  # Calculate distance
  tpmDist <- calcDist(tpm, dmethod=dmethod)
  # Create variables to store data
  outJaccard <- data.frame()
  outCluster <- data.frame(row.names = colnames(tpm))
  # Loop through k and extract data
  for (k in krange) {
    clb <- clusterboot(
      tpmDist,
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
  pdf(outfile, onefile = T, paper='a4r', width=9, height=7)
  pheatmap(tpmDist)
  # Create theme
  t <- theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16),
    legend.title=element_text(size=16)
  ) 
  # Plot jaccard similarities
  p <- ggplot(outJaccard, aes(x=k, y = jaccard)) +
    geom_hline(aes(yintercept = c(0.5,0.75)), col = c('red','green')) +
    geom_point(shape = 20, size = 8) +
    ylim(0,1) + 
    stat_summary(fun.y=mean, geom="point", color="blue", shape = 20, size = 8) +
    t
  print(p)
  # Plot contribution of each cell type to cluster
  for (k in krange) {
    # Create plot
    clID <- paste0('k',k)
    p <- ggplot(outCluster, aes_string(clID)) +
      t
    # Add group data
    if (!is.null(groups)) {
      p <- p + geom_bar(aes(fill=groups))
    } else {
      p <- p + geom_bar()
    }
    # Print
    print(p)
  }
  # Close file
  dev.off()
}



plotDendro <- function(d) {
  hc <- hclust(d)
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc, type = "rectangle")
  p <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0.2, 0))
  p
}