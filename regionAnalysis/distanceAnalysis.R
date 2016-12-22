require('GenomicRanges')
require('ggplot2')
###############################################################################
## Convert bed to grange
###############################################################################
# Create fuction to read in bed data
bed2Grange <- function(file, chr=NULL) {
  # Read in data, name columns and ajust start
  bedData <- read.table(file, header=F, sep='\t')
  colnames(bedData) <-c('chr', 'start', 'end')
  bedData$start <- bedData$start + 1
  # Convert to genomic range without chromosome data
  if (is.null(chr)) {
    bedGrange <- GRanges(
      seqnames = as.character(bedData$chr),
      ranges = IRanges(bedData$start, bedData$end)
    )
  # Convert to genomic range with chromosome data
  } else {
    bedData <- bedData[bedData$chr %in% chr,]
    bedGrange <- GRanges(
      seqnames = as.character(bedData$chr),
      ranges = IRanges(bedData$start, bedData$end),
      seqinfo = Seqinfo(seqnames=chr)
    )
  }
  # Return genomic range
  return(bedGrange)
}

###############################################################################
## Extract cumulative distances between genomic ranges objects
###############################################################################
calculateDistanceCumFreq <- function(query, subject) {
  dtn <- distanceToNearest(query, subject)
  d <- mcols(dtn)$distance
  d <- d[order(d)]
  drle <- rle(d)
  outdf <- data.frame(
    'distance' = drle$values,
    'cumfreq' = cumsum(drle$lengths) / sum(drle$lengths)
  )
  return(outdf)
}

###############################################################################
## Extract cumulative distance between genomic range object lists
###############################################################################
calculateDistanceCumFreqList <- function(queryList, subjectList) {
  # Extract all potential combinations
  combinations <- expand.grid(names(queryList), names(subjectList))
  # Calculate distance for each combination and add query and subject name
  listData <- apply(combinations, 1, function(z) {
    names <- as.character(z)
    distanceDF <- calculateDistanceCumFreq(
      queryList[[names[1]]], subjectList[[names[2]]])
    distanceDF$query <- names[1]
    distanceDF$subject <- names[2]
    return(distanceDF)
  })
  # Create and return output dataframe
  outData <- do.call(rbind, listData)
  outData <- outData[,c('query', 'subject', 'distance', 'cumfreq')]
  outData$query <- factor(outData$query, levels=names(queryList))
  outData$subject <- factor(outData$subject, levels=names(subjectList))
  return(outData)
}

###############################################################################
## Create cumulative distance plot
###############################################################################
createDistanceCumFreqPlot <- function(plotData, maxDist = 100000) {
  plot <- ggplot(
      plotData,
      aes(x=distance, y=cumfreq, group=query, colour=query)
    ) +
    geom_line() +
    facet_wrap(~subject) +
    scale_x_continuous(limits=c(0, maxDist)) +
    scale_y_continuous(limits=c(0, 1)) +
    labs(x='Distance', y='Cumulative Frequency')
  return(plot)
}

###############################################################################
## Calculate p-value for overlap btween 2 query samples and reference
###############################################################################
pvalue_overlap <- function(
  query1, query2, subject, maxgap=0, minoverlap=1
) {
  # Perform overlaps
  q1Length <- length(query1)
  q1Overlaps <- length(queryHits(findOverlaps(
    query1, subject, maxgap=maxgap, minoverlap=minoverlap
  )))
  q2Length <- length(query2)
  q2Overlaps <-length(queryHits(findOverlaps(
    query2, subject, maxgap=maxgap, minoverlap=minoverlap
  )))
  # Generate query ratio difference
  q1Ratio <- q1Overlaps / q1Length
  q2Ratio <- q2Overlaps / q2Length
  deltaRatio <- q1Ratio - q2Ratio
  # Calculate pvalue
  countMatrix <- matrix(c(
    q1Overlaps, q1Length - q1Overlaps,
    q2Overlaps, q2Length - q2Overlaps
  ), ncol=2)
  pvalue <- fisher.test(countMatrix)$p.value
  print(countMatrix)
  return(list(
    'maxgap'=maxgap,
    'minoverlap'=minoverlap,
    'query1ratio'=q1Ratio,
    'query2ratio'=q2Ratio,
    'deltaratio'=deltaRatio,
    'pvalue'=pvalue
  ))
}

###############################################################################
## Calculate p-value for overlap between paired multiple query samples and 
## multiple reference samples
###############################################################################
multiple_pvalue_overlap <- function(
    queryList, subjectList, maxgap=0, minoverlap=1
) {
  # Create output comparison dataframe
  comparisons <- combn(names(queryList), 2)
  outputDF <- data.frame(
    'query1' = rep(comparisons[1,], length(subjectList)),
    'query2' = rep(comparisons[2,], length(subjectList)),
    'subject' = rep(names(subjectList), each=ncol(comparisons)),
    stringsAsFactors = F
  )
  # Add query and subject length
  outputDF$query1_length <- sapply(queryList, length)[outputDF$query1]
  outputDF$query2_length <- sapply(queryList, length)[outputDF$query2]
  outputDF$subject_length <- sapply(subjectList, length)[outputDF$subject]
  # Add blank columns
  outputDF$query1_ratio <- NA
  outputDF$query2_ratio <- NA
  outputDF$pvalue <- NA
  # Populate blank columns
  for (i in 1:nrow(outputDF)) {
    query1 <- outputDF[i, 'query1']
    query2 <- outputDF[i, 'query2']
    subject <- outputDF[i, 'subject']
    data <- pvalue_overlap(
      queryList[[query1]], queryList[[query2]], subjectList[[subject]], 
      maxgap=maxgap, minoverlap=minoverlap
    )
    outputDF[i, 'query1_ratio'] <- data$query1ratio
    outputDF[i, 'query2_ratio'] <- data$query2ratio
    outputDF[i, 'pvalue'] <- data$pvalue
  }
  return(outputDF)
}

###############################################################################
## Calculate p-value for overlap
###############################################################################
pvalue_median_distance <- function(
  query, reference, subject, iterations=10000, seed=1234
) {
  # Extract distance
  qDist <- mcols(distanceToNearest(query, subject))$distance
  rDist <- mcols(distanceToNearest(reference, subject))$distance
  combinedDist <- c(qDist, rDist)
  # Generate reference ratio difference
  qMedian <- median(qDist)
  rMedian <- median(rDist)
  deltaMedian <- qMedian - rMedian
  # Create variables for iterations
  sampleDeltaMedian <- numeric(length=iterations)
  qLength <- length(qDist)
  rLength <- length(rDist)
  # Perform iterations to calculate difference between resampled medians
  set.seed(seed)
  for (i in 1:iterations) {
    # Resample test and control samples
    resample <- sample(combinedDist)
    qMedianResample <- median(head(resample, qLength))
    rMedianResample <- median(tail(resample, rLength))
    # Calculate and store ratio
    deltaMedianResample <- qMedianResample - rMedianResample
    sampleDeltaMedian[i] <- deltaMedianResample
  }
  # Calculate pvalue
  count <- sum(abs(sampleDeltaMedian) >= abs(deltaMedian))
  pvalue <- max(3 / iterations, sum(count) / iterations)
  return(list(
    'iterations'=iterations,
    'querymedian'=qMedian,
    'referencemedian'=rMedian,
    'deltamedian'=deltaMedian,
    'resample'=sampleDeltaMedian,
    'pvalue'=pvalue
  ))
}

###############################################################################
## Function to create genomic range object from diffreps output file
###############################################################################
diffreps2grange <- function(file, log2fc, chr=NULL) {
  # Read in data
  data <- read.table(file, header=T, sep='\t')
  if (log2fc == '+') {
    data <- subset(data, data$log2FC > 0)
  } else if (log2fc == '-') {
    data <- subset(data, data$log2FC < 0)
  } else if (log2fc != '*') {
    stop("log2fc must be one of '+', '-', or '*'" )
  }
  # Convert to genomic range without chromosome data
  if (is.null(chr)) {
    gr <- GRanges(
      seqnames = as.character(data$Chrom),
      ranges = IRanges(data$Start, data$End),
      treatment = data$Treatment.avg,
      control = data$Control.avg,
      log2fc = data$log2FC,
      pval = data$pval,
      padj = data$padj
    )
    # Convert to genomic range with chromosome data
  } else {
    data <- data[data$Chrom %in% chr,]
    gr <- GRanges(
      seqnames = as.character(data$Chrom),
      ranges = IRanges(data$Start, data$End),
      seqinfo = Seqinfo(seqnames=chr),
      treatment = data$Treatment.avg,
      control = data$Control.avg,
      log2fc = data$log2FC,
      pval = data$pval,
      padj = data$padj
    )
  }
  # Return genomic range
  return(gr)
}


