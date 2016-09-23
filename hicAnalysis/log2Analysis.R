require('ggplot2')
require('ggdendro')
require('reshape2')

###############################################################################
## Generate plot for log2 clustering
###############################################################################
log2.clustering.plot <- function(log2.data, sample.data) {
  # Perform clustetring
  trim.data <- log2.data[,6:ncol(log2.data)]
  log2.hclust <- hclust(dist(t(trim.data)))
  # Extract dendrogram data
  dd <- as.dendrogram(log2.hclust)
  ddata <- dendro_data(dd)
  # Extract condition for each sample
  sample.match <- match(
    as.character(ddata_x$labels$label),
    as.character(sample.data$sample)
  )
  conditions <- as.character(sample.data$condition)[sample.match]
  # Create dendrogram
  plot <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    theme(
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      panel.grid=element_blank()
    ) +
    ggtitle('Clustering of euclidean distance between log2 vectors') +
    geom_text(
      data=label(ddata_x),
      aes(label=label, x=x, y=0, colour=conditions, vjust=1.5, xpd=T)
    )
  # Return plot
  return(plot)
}

###############################################################################
## Generate plot for variance distribution
###############################################################################
log2.density.plot <- function(log2.data, sample.data) {
  # Perform clustetring
  trim.data <- as.matrix(log2.data[,6:ncol(log2.data)])
  trim.data <- na.omit(trim.data)
  # Extract data for each condition
  condition.list <- split(sample.data$sample, sample.data$condition)
  log2.list <- lapply(condition.list, function(z) {
    c(trim.data[,as.character(z)])
  })
  # Create plot
  plot.data <- do.call(cbind, log2.list)
  plot.data <- melt(plot.data)
  colnames(plot.data) <- c('_','Condition','Log2')
  # Create plot and return
  plot <- ggplot(plot.data, aes(x=Log2, colour=Condition)) +
    geom_density() +
    scale_x_continuous(limits=c(-1,1)) +
    ggtitle('Distribution of log2 values')
  return(plot)
}

###############################################################################
## Functions to calculate log2 ratio
###############################################################################
# Calculates mean log2 ratio of pairs
log2.pair.calc <- function(
  upstream.prob, downstream.prob 
) {
  log2.values <- log2(upstream.prob / downstream.prob)
  log2.mean <- mean(log2.values)
  return(log2.mean)
}
# Calculates log2 ratio of sum of probabilities
log2.sum.calc <- function(
  upstream.prob, downstream.prob 
) {
  return(log2(sum(upstream.prob) / sum(downstream.prob)))
}

###############################################################################
## Function to calculate log2 data
###############################################################################
calc.log2.direction <- function(
  norm.matrix, log2.function, max.window.no, min.window.no
) {
  # Remove zero values
  norm.matrix[norm.matrix == 0] <- NA
  # Calculate bin number and generate output variable
  number.bins <- ncol(norm.matrix)
  log2.out <- rep(NA, number.bins)
  # Loop through bins
  for (index in 1:number.bins) {
    # Calculate window size for statistic
    actual.window.no <- min(index - 1, number.bins - index - 1, max.window.no)
    if (actual.window.no < min.window.no) {
      next
    }
    # Extract upstream and downstream probabilities
    upstream.indices <- (index - actual.window.no) : (index - 1)
    upstream.prob <- rev(norm.matrix[upstream.indices, index])
    downstream.indices <- (index + 1) : (index + actual.window.no)
    downstream.prob <- norm.matrix[downstream.indices, index]
    # Calculate accpetable indices
    acceptable.indices <- !(is.na(upstream.prob) | is.na(downstream.prob))
    if (sum(acceptable.indices) < min.window.no) {
      next
    }
    # Extract accpetable probabilities and calculate log2 value
    upstream.prob <- upstream.prob[acceptable.indices]
    downstream.prob <- downstream.prob[acceptable.indices]
    log2.ratio <- log2.function(upstream.prob, downstream.prob)
    log2.out[index] <- log2.ratio
  }
  # Create and return output data frame
  out.df <- as.data.frame(
    do.call(rbind, strsplit(colnames(norm.data), '[:-]')),
    stringsAsFactors = F
  )
  colnames(out.df) <- c('chr', 'start', 'end')
  out.df$start <- as.numeric(out.df$start)
  out.df$end <- as.numeric(out.df$end)
  out.df$log2 <- log2.out
  row.names(out.df) <- colnames(norm.matrix)
  return(out.df)
}

###############################################################################
## Function to test whether two vectors have different variance
###############################################################################
deltaVarianceTest <- function(v1, v2, iterations=10000, seed=1234) {
  # Calculate reference delta variation
  deltaVar <- var(v1) - var(v2)
  # Concatenate values and make output vector
  jointValues <- c(v1, v2)
  deltaVarResample <- numeric(length=iterations)
  # Resample joint values and calculate delta variance distribution
  set.seed(seed)
  for (i in 1:iterations) {
    reorderValues <- sample(jointValues)
    var1 <- var(head(reorderValues, length(v1)))
    var2 <- var(tail(reorderValues, length(v2)))
    deltaVarResample[i] <- var1 - var2
  }
  # Calculate p.value
  extreme <- abs(deltaVarResample) >= abs(deltaVar)
  pValue <- max((sum(extreme) / iterations), (1 / iterations)) 
  # Create plot
  plot(density(deltaVarResample))
  abline(v=deltaVar)
  output = list()
  output$'var1' = var(v1)
  output$'var2' = var(v2)
  output$'pvalue' = pValue
  return(output)
}

###############################################################################
## Compare variance across clusters
###############################################################################


###############################################################################
## Function to extract sample names grouped by condition
###############################################################################
samplesByCondition <- function(sampleData) {
  outList = list()
  for (condition in unique(sampleData$condition)) {
    samples <- sampleData$sample[sampleData$condition == condition]
    outList[[condition]] = as.character(samples)
  }
  return(outList)
}

###############################################################################
## Function to extract vectors of log2 values for supplied conditions
###############################################################################
extractConditionLog2 <- function(log2Data, sampleData) {
  sampleConditions <- samplesByCondition(sampleData)
  output <- list()
  for (condition in names(sampleConditions)) {
    samples <- sampleConditions[[condition]]
    vector <- c(as.matrix(log2Data[,samples]))
    output[[condition]] <- vector
  }
  output <- data.frame(do.call(cbind, output))
  return(output)
}

###############################################################################
## Function to extract tad boundaries
###############################################################################
splitLog2RegionSample <- function(log2Data) {
  samples <- colnames(log2Data)[6:ncol(log2Data)]
  regionList <- split(log2Data, log2Data$region)
  regionSampleList <- list()
  for (region in names(regionList)) {
    regionData <- regionList[[region]]
    sampleList <- list()
    for (sample in samples) {
      sampleIndex <- match(sample, colnames(regionData))
      sampleDF <- regionData[,c(1:5, sampleIndex)]
      colnames(sampleDF)[6] <- 'log2'
      sampleList[[sample]] <- sampleDF
    }
    regionSampleList[[region]] <- sampleList
  }
  return(regionSampleList)
}









