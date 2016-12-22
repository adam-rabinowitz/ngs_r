###############################################################################
## Function to create DESeq2 dataset
###############################################################################
createDESeqDataSet <- function(counts, samples, batch=F) {
  # Read sample data file
  sampleData <- read.table(samples, sep = '\t', header = T,
    row.names = 1, check.names = F)
  # Convert condition and replicate to factors
  sampleData$condition <- factor(sampleData$condition,
    levels=unique(sampleData$condition))
  sampleData$replicate <- factor(as.character(sampleData$replicate),
    levels=unique(as.character(sampleData$replicate)))
  # Read expected count file
  countData <- read.table(counts, sep = '\t', row.names = 1,
    header = T, check.names = 'F')
  # Check all samples in count data and resort
  inter <- intersect(row.names(sampleData), colnames(countData))
  if (length(inter) != length(row.names(sampleData))) {
    stop('Not all samples in count data')
  }
  countData <- countData[,row.names(sampleData)]
  # Round expected counts
  countData <- round(countData)
  # Create DESeq2 dataset and perform analysis
  if (batch) {
    DESeqData <- DESeqDataSetFromMatrix(
      countData = countData,
      colData = sampleData,
      design = ~ replicate + condition
    )
  } else {
    DESeqData <- DESeqDataSetFromMatrix(
      countData = countData,
      colData = sampleData,
      design = ~ condition
    )
  }
  # Return deseq object
  return(DESeqData)
}
  
###############################################################################
## Create MA plot from DESeqDataSet
###############################################################################
maPlotDESeq <- function(DESeqDataSet, padj) {
  # Check input object
  if (class(DESeqDataSet)[1] != 'DESeqDataSet') {
    stop('DESeqDataSet object must be supplied')
  }
  # Extract results as dataframe
  resultsData <- as.data.frame(results(DESeqDataSet))
  resultsData <- subset(resultsData, is.na(resultsData$log2FoldChange) == F)
  # Find y axis limits
  ylim = max(abs(range(resultsData$log2FoldChange)))
  # Create colour for plot
  resultsData$sig <- ifelse(resultsData$padj < padj, '+', '-') 
  # Create plots
  p <- ggplot(resultsData, aes(x = log2(baseMean), y = log2FoldChange,
    colour = sig)) +
    geom_hline(yintercept=0) +
    geom_point() +
    ylim(-ylim,ylim) +
    xlab('Log2 Base Mean') +
    ylab('Log2 Fold Change') +
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=18,face="bold"),
      legend.text=element_text(size=16),
      legend.title=element_text(size=16)
    ) 
  return(p)
}

###############################################################################
## Create Volcano plot from DESeqDataSet
###############################################################################
volcanoPlotDESeq <- function(DESeqDataSet, padj) {
  # Check input object
  if (class(DESeqDataSet)[1] != 'DESeqDataSet') {
    stop('DESeqDataSet object must be supplied')
  }
  # Extract results as dataframe
  resultsData <- as.data.frame(results(DESeqDataSet))
  resultsData <- subset(resultsData, is.na(resultsData$padj) == F)
  # Find y axis limits
  xlim = max(abs(range(resultsData$log2FoldChange)))
  # Create colour for plot
  resultsData$sig <- ifelse(resultsData$padj < padj, '+', '-') 
  # Create plots
  p <- ggplot(resultsData, aes(x = log2FoldChange, y = -log10(padj),
    colour = sig)) +
    geom_point() +
    xlim(-xlim,xlim) +
    xlab('Log2 Base Mean') +
    ylab('-Log10 Adjusted P-Value') +
    theme(
      axis.text=element_text(size=16),
      axis.title=element_text(size=18,face="bold"),
      legend.text=element_text(size=16),
      legend.title=element_text(size=16)
    ) 
  return(p)
}

###############################################################################
## Create hierachical cluster
###############################################################################
rlogClusteringDESeq <- function(DESeqDataSet) {
  # Check input object
  if (class(DESeqDataSet)[1] != 'DESeqDataSet') {
    stop('DESeqDataSet object must be supplied')
  }
  # Extract rlog counts and perform clustering
  rlogData <- rlog(DESeqDataSet, blind=T)
  rlogCounts <- assay(rlogData)
  # Perform clustring and create dendrogram
  hc <- hclust(dist(t(rlogCounts)))
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc, type = "rectangle")
  # Create plot
  p <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() +
    geom_text(data = label(ddata), aes(x = x, y = y, label = label),
      hjust=0) +
    scale_y_reverse(expand = c(0.2, 0))+
    theme(
      panel.grid.minor.y=element_blank(),
      panel.grid.major.y=element_blank()
    )
  # Return plot and rlog data
  return(p)
}

###############################################################################
## Save results from DESEq data object to file
###############################################################################
extractDataDESeq  <- function(DESeqDataSet, prefix, altname=NA) {
  # Extract file for gene name conversion
  if (!is.na(altname)) {
    nameData <- read.table(altname, sep = '\t', header = F,
      stringsAsFactors=F)
    nameData <- split(nameData[,2], nameData[,1])
    nameData <- sapply(nameData, paste0, collapse=' ')
  }
  # Extract results and sort by p-value
  resultData <- as.data.frame(results(DESeqDataSet))
  resultData <- resultData[order(resultData$pvalue,
    resultData$log2FoldChange),]
  # Add gene names to dataframe and reorder
  resultData$gene <- row.names(resultData)
  if (!is.na(altname)) {
    resultData$altname <- nameData[row.names(resultData)]
    colNo <- ncol(resultData)
    resultData <- resultData[,c(colNo - 1, colNo, 1:(colNo-2))]
  } else {
    colNo <- ncol(resultData)
    resultData <- resultData[,c(colNo, 1:(colNo-1))]
  }
  # Extract and sort counts
  countData <- as.data.frame(counts(DESeqDataSet, normalized=T))
  countData <- countData[resultData$gene,]
  # Add genenames to datframe and reorder
  countData$gene <- row.names(countData)
  if (!is.na(altname)) {
    countData$altname <- nameData[row.names(resultData)]
    colNo <- ncol(countData)
    countData <- countData[,c(colNo - 1, colNo, 1:(colNo-2))]
  } else {
    colNo <- ncol(countData)
    countData <- countData[,c(colNo, 1:(colNo-1))]
  }
  # Save rdata
  write.table(resultData, file = paste0(prefix,'.results.txt'), sep = '\t',
    quote = F, col.names = T, row.names = F)
  write.table(countData, file = paste0(prefix,'.counts.txt'), sep = '\t',
    quote = F, col.names = T, row.names = F)
}

###############################################################################
## Create database from DESeq2 dataset
###############################################################################


