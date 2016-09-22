require(tsne)




data <- read.table('~/test.df', header=T, sep='\t', check.names=F)
log2Data <- na.omit(data[,6:ncol(data)])
sampleData <- read.table('~/yasuMatrices/sample_data.txt', header=T, sep='\t')


log2Clustering <- function(log2Data) {
  trimData <- log2Data[,6:ncol(log2Data)]
  log2Dist <- dist(t(trimData))
  plot(hclust(log2Dist))
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
  pValue <- sum(extreme) / iterations
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
conditionLog2 <- extractConditionLog2(log2Data, sampleData)
deltaVarianceTest(conditionLog2$Mitosis.AuB, conditionLog2$Mitosis)

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

###############################################################################
## Function to calculate jaccard similarity scores across samples
###############################################################################
calculate.jaccard <- function(value.list) {
  output <- matrix(
    ncol=length(value.list),
    nrow=length(value.list),
    dimnames=list(
      names(value.list),
      names(value.list)  
    )
  )
  output['NGS-8179','NGS-8179'] <- 1
  combinations <- expand.grid(names(value.list), names(value.list))
  for (i in 1:nrow(combinations)) {
    sample1 <- combinations[i,1]
    sample2 <- combinations[i,2]
    values1 <- value.list[[sample1]]
    values2 <- value.list[[sample2]]
    intersectNo <- length(intersect(values1, values2))
    totalNo <- length(values1) + length(values2)
    jaccard <- intersectNo / (totalNo - intersectNo)
    output[sample1, sample2] <- jaccard
  }
  return(output)
}

###############################################################################
## Function to calculate mean jaccard similarity scores across conditions
###############################################################################
mean.jaccard.similarity <- function(jaccard.matrix, sample.data) {
  # Check row and column names and remove diagonal
  if (!all.equal(colnames(jaccard.matrix), row.names(jaccard.matrix))) {
    stop('Row and column names must match')
  }
  diag(jaccard.matrix) <- NA
  # Extract conditions and generate output matrix
  conditions <- split(sample.data$sample, sample.data$condition)
  output <- matrix(
    ncol=length(conditions),
    nrow=length(conditions),
    dimnames=list(
      names(conditions),
      names(conditions)  
    )
  )
  # Fill output matrix and return
  combinations <- expand.grid(names(conditions), names(conditions))
  for (i in 1:nrow(combinations)) {
    condition1 <- combinations[i,1]
    condition2 <- combinations[i,2]
    samples1 <- as.character(conditions[[condition1]])
    samples2 <- as.character(conditions[[condition2]])
    subMatrix <- jaccard.matrix[samples1, samples2]
    print(samples1)
    print(samples2)
    print(subMatrix)
    subMatrixMean <- mean(c(subMatrix), na.rm=T)
    print(subMatrixMean)
    output[condition1,condition2] <- subMatrixMean
  }
  return(output)
}
