require(Rtsne)
require(ggplot2)

##############################################################################
## Perform TSNE
##############################################################################
tsneReduce <- function(tpm, log=T, perplexity=10, theta=0.1, seed=1979,
  dims=2) {
  # Transpose and (optionally) log transform data
  if (log) {
    tpm <- t(log2(tpm + 1))
  } else {
    tpm <- t(tpm)
  }
  # Perform tsne and return data
  set.seed(seed)
  tData <- Rtsne(X=tpm, dims=dims, pca=F, perplexity=perplexity, theta=theta)
  # Add sample names to output and return
  tData$samples <- row.names(tpm)
  return(tData)
}

##############################################################################
## Plot TSNE
##############################################################################
tsnePlot <- function(tsne, groups=NULL) {
  # Create theme
  thm <- theme(
    axis.text=element_text(size=16),
    axis.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16),
    legend.title=element_text(size=16),
    title=element_text(size=18,face='bold')
  )
  # Create plots with groups
  if(!is.null(groups)){
    # Check groups
    groups <- checkGroups(tsne$samples, groups)
    # Create plot for each group
    for (i in 1:length(groups)) {
      # Extract group data
      groupName = names(groups)[i]
      groupFactor = groups[[i]]
      # Create dataframe for plots
      plotDF <- as.data.frame(tsne$Y)
      colnames(plotDF) <- paste0('Dim', 1:ncol(plotDF))
      plotDF[,groupName] = groupFactor
      # Create 2D plot if possible
      if (ncol(plotDF) >= 3) {
        fig <- ggplot(plotDF, aes_string(x='Dim1', y='Dim2', color=groupName)) + 
          geom_point() +
          ggtitle('tSNE - Dim1 vs Dim2') +
          thm
        print(fig)
      }
      # Modify data for 1D plot
      plotDFalt <- plotDF
      plotDF[,groupName] <- 'All'
      plotDF <- rbind(plotDF, plotDFalt)
      # Create plot for each dimension
      for (i in 1:(ncol(plotDF)-1)) {
        dim <- paste0('Dim', i)
        fig <- ggplot(plotDF, aes_string(x=dim, color=groupName)) + 
          geom_density() +
          ggtitle(paste('tSNE -', dim, 'Distribution')) +
          thm
        print(fig)
      }
    }
  } else {
    # Create dataframe for plots
    plotDF <- as.data.frame(tsne$Y)
    colnames(plotDF) <- paste0('Dim', 1:ncol(plotDF))
    # Create 2D plot if possible
    if (ncol(plotDF) >= 2) {
      fig <- ggplot(plotDF, aes_string(x='Dim1', y='Dim2')) + 
        geom_point() +
        ggtitle('tSNE - Dim1 vs Dim2') +
        thm
      print(fig)
    }
    # Create plot for each dimension
    for (i in 1:ncol(plotDF)) {
      dim <- paste0('Dim', i)
      fig <- ggplot(plotDF, aes_string(x=dim)) + 
        geom_density() +
        ggtitle(paste('tSNE -', dim, 'Distribution')) +
        thm
      print(fig)
    }
  }
}
