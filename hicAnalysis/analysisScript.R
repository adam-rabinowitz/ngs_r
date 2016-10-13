rm(list=ls())
suppressPackageStartupMessages(require('optparse'))
option_list = list(
  make_option(c('-L', '--log2Data'), action='store', default=NA,
    type='character', help="log2 data file"),
  make_option(c('-S', '--sampleData'), action='store', default=NA,
    type="character", help='file listing condition for each sample'),
  make_option(c('-W', '--halfWindowSize'), action='store', default=NA,
    type="character", help='size of half window for tad identification'),
  make_option(c("-M", "--minShift"), action="store", default=NA,
    type='character', help='minimum shift for tad filtering'),
  make_option(c("-O", "--outPrefix"), action="store", default=NA,
    type='character', help='Full path to prefix of output files')
)
# Import arguments
opt <- parse_args(OptionParser(option_list=option_list))
opt$log2Data <- '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/log2.df.txt'
opt$sampleData <- '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/2kbBinData/auroraBInterphaseMitosis/sample_data.txt'
opt$halfWindowSize <- 5
opt$minShift <- 0.2
opt$outPrefix <- '/farm/home/rabino01/Yasu/auroraB'
# Check required arguments have been supplied
if (is.na(opt$log2Data)) {
  stop('Log2 data file must be supplied')
}
if (is.na(opt$sampleData)) {
  stop('Sample file must be supplied')
}
if (is.na(opt$outPrefix)) {
  stop('Output prefix must be supplied')
}
# Source required functions
source('~/github/ngs_r/hicAnalysis/log2Analysis.R')
source('~/github/ngs_r/hicAnalysis/tadAnalysis.R')
source('~/github/ngs_r/hicAnalysis/similarityAnalysis.R')
# Read in log2 data and check all samples are present in sample data
log2.data <- read.table(opt$log2Data, header=T, sep='\t', check.names=F)
sample.data <- read.table(opt$sampleData, header=T, sep='\t', check.names=F)
sample.names <- colnames(log2.data)[6:ncol(log2.data)]
if (anyNA(match(c(sample.names), sample.data$sample))) {
  stop('Sample data not present for all samples in log2 table')
}
# Create plot for clustering
clusterPlot <- log2.clustering.plot(log2.data, sample.data)
varPlot <- log2.density.plot(log2.data, sample.data)
# Calculate pvalue for variation in absolute log2
varPvalue <- createVarPvalueMatrix(log2.data, sample.data, iteration=10000)
# Split data and identify tads
region.sample.log2 <- splitLog2RegionSample(log2.data)
region.sample.tad <- lapply(region.sample.log2, function(r) {
  lapply(r, function(s) {
    find.filter.tad.window(
      s, half.window=opt$halfWindowSize, min.shift=opt$minShift, na.rm=T
    )
  })
})
# Calculate overlap of tads
tad.centres <- lapply(region.sample.tad, function(r) {
  lapply(r, function(s) {
    s$centre
  })
})
# Calculate mean jaccard scores for each value
similarity.matrices <- lapply(tad.centres, function(z) {
  sample.similarity <- calculate.similarity(z, diff=2000)
  condition.similarity <- condition.similarity(sample.similarity, sample.data)
  return(condition.similarity)
})
mean.matrices <- Reduce('+', similarity.matrices) / length(similarity.matrices)

tad.plot <- plot.tads(region.sample.tad[[1]], sample.data)
  
  








