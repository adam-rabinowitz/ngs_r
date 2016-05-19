# Clear workspace
rm(list=ls())
# Create command line oprions
suppressPackageStartupMessages(require('optparse'))
option_list = list(
  make_option(c("-E", "--expData"), action="store", default=NA,
              type='character', help="expected counts file"),
  make_option(c("-S", "--sampleData"), action="store", default=NA,
              type="character", help="sample data"),
  make_option(c("-O", "--outPrefix"), action="store", default=NA,
              type='character', help='output prefix'),
  make_option(c("-B", "--modelBatch"), action="store_true", default=F,
              type='logical', help='add batch effects to model [FALSE]'),
  make_option(c("-N", "--altName"), action="store", default=NA,
              type='character', help='alternative gene name file [NA]'),
  make_option(c("-T", "--bioType"), action="store", default=NA,
              type='character', help='gene biotype file [NA]'),
  make_option(c("-P", "--padj"), action="store", default=0.05,
              type='numeric', help='adjusted p-value [0.05]')
)
# Import arguments
opt = parse_args(OptionParser(option_list=option_list))
# Check required arguments have been supplied
if (is.na(opt$expData)) {
  stop('Expected count file must be supplied')
}
if (is.na(opt$sampleData)) {
  stop('Sample file must be supplied')
}
if (is.na(opt$outPrefix)) {
  stop('Output prefix must be supplied')
}
# Load additional required modules
require('DESeq2')
require('ggplot2')
require('ggdendro')
# Read sample data file
sampleData <- read.table(opt$sampleData, sep = '\t', header = T,
  row.names = 1, check.names = F)
sampleData$condition <- factor(sampleData$condition,
  levels=unique(sampleData$condition))
sampleData$replicate <- factor(sampleData$replicate,
  levels=unique(sampleData$replicate))
# Read expected count file
expData <- read.table(opt$expData, sep = '\t', row.names = 1,
  header = T, check.names = 'F')
# Check column names in expected count file and resort
inter <- intersect(row.names(sampleData), colnames(expData))
if (length(inter) != length(row.names(sampleData))) {
    stop('Not all samples in count data')
}
expData <- expData[,row.names(sampleData)]
# Round expected counts
expData <- round(expData)
# Extract file for gene name conversion
if (!is.na(opt$altName)) {
  nameData <- read.table(opt$altName, sep = '\t', header = T,
    stringsAsFactors=F)
  nameData <- split(nameData[,2], nameData[,1])
  nameData <- sapply(nameData, paste0, collapse=' ')
}
# Extract file for biotype conversion
if (!is.na(opt$bioType)) {
  bioData <- read.table(opt$bioType, sep = '\t', header = T,
    stringsAsFactors=F)
  bioData <- split(bioData[,2], bioData[,1])
  bioData <- sapply(bioData, paste0, collapse=' ')
}
# Create DESeq2 dataset
if (opt$modelBatch) {
  DESeqData <- DESeqDataSetFromMatrix(
    countData = expData,
    colData = sampleData,
    design = ~ replicate + condition
  )
} else {
  DESeqData <- DESeqDataSetFromMatrix(
    countData = expData,
    colData = sampleData,
    design = ~ condition
  )
}
# Perform DESeq2 analysis
DESeqAnalysis <- DESeq(DESeqData)
# Create MA plot
maPlot <- function(DESeqObject, padj = opt$padj) {
  # Extract results as dataframe
  resultsData <- as.data.frame(results(DESeqObject))
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
# Create volcano plot
volcanoPlot <- function(DESeqObject, padj = opt$padj) {
  # Extract results as dataframe
  resultsData <- as.data.frame(results(DESeqObject))
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
# Create hierachical cluster
cluster <- function(DESeqObject) {
  countData <- counts(DESeqObject, normalized = T)
  countData <- countData[apply(countData, 1, max) > 0,]
  hc <- hclust(dist(t(log2(countData+1))))
  dhc <- as.dendrogram(hc)
  ddata <- dendro_data(dhc, type = "rectangle")
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
  return(p)
}
# Function to store results
extractResults  <- function(DESeqObject, prefix = opt$outPrefix) {
  # Extract results and sort by p-value
  resultData <- as.data.frame(results(DESeqObject))
  resultData <- resultData[order(resultData$pvalue),]
  # Add biotype and reorder
  if (!is.na(opt$bioType)) {
    resultData$biotype <- bioData[row.names(resultData)]
    resultData <- resultData[,c(ncol(resultData),1:(ncol(resultData)-1))]
  }
  # Add alternative gene names and reorder
  if (!is.na(opt$altName)) {
    resultData$altname <- nameData[row.names(resultData)]
    resultData <- resultData[,c(ncol(resultData),1:(ncol(resultData)-1))]
  }
  # Add gene names to dataframe and reorder
  resultData$gene <- row.names(resultData)
  resultData <- resultData[,c(ncol(resultData),1:(ncol(resultData)-1))]
  # Extract and sort counts
  countData <- as.data.frame(counts(DESeqObject, normalized=T))
  countData <- countData[resultData$gene,]
  countData$gene <- row.names(countData)
  countData <- countData[,c(ncol(countData),1:(ncol(countData)-1))]
  # Save results data
  write.table(resultData, file = paste0(prefix,'.results'), sep = '\t',
    quote = F, col.names = T, row.names = F)
  # Save count data
  write.table(countData, file = paste0(prefix,'.counts'), sep = '\t',
    quote = F, col.names = T, row.names = F)
}
# Store results and count data
extractResults(DESeqAnalysis)
# Store plots
pdf(file = paste0(opt$outPrefix,'.pdf'), paper = 'a4r', width = 9,
  height = 7, onefile = T)
print(cluster(DESeqAnalysis))
print(maPlot(DESeqAnalysis))
print(volcanoPlot(DESeqAnalysis))
dev.off()
