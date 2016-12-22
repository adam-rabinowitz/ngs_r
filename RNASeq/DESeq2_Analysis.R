# Clear workspace
rm(list=ls())
# Create command line options and parse
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
source('~/github/ngs_r/RNASeq/DESeq2_Functions.R')
DESeqData <- createDESeqDataSet <- function(
  opt$expData,
  opt$sampleData,
  opt$modelBatch
)
DESeqAnalysis <- DESeq(DESeqData)
# Store results
extractDataDESeq(DESeqAnalysis, opt$outPrefix, opt$altName)
# Create plots
clusterPlot <- rlogClusteringDESeq(DESeqAnalysis)
maPlot <- maPlotDESeq(DESeqAnalysis, opt$padj)
volcanoPlot <- volcanoPlotDESeq(DESeqAnalysis, opt$padj)
pdf(file = paste0(opt$outPrefix,'.pdf'), paper = 'a4r', width = 9,
  height = 7, onefile = T)
print(clusterPlot)
print(maPlot)
print(volcanoPlot)
dev.off()
