# Clear workspace
rm(list=ls())
# Create command line oprions
suppressPackageStartupMessages(require('optparse'))
option_list = list(
  make_option(c("-N", "--normData"), action="store", default=NA,
              type='character', help="normalised counts file"),
  make_option(c("-O", "--outFile"), action="store", default=NA,
              type="character", help="output pdf file")
)
# Import arguments
opt = parse_args(OptionParser(option_list=option_list))
# Read in counts
normCounts = read.table(opt$normData, header=T, sep='\t', check.names=F)
# Create prob plot
source('~/github/ngs_r/hicAnalysis/tadAnalysis.R')
normPlot <- create_prob_plot(
    normCounts,
    max_prob = 1e-1,
    min_prob = 1e-6
)
pdf(opt$outFile, paper='a4')
print(normPlot)
dev.off()
