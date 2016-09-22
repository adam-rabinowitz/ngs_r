##############################################################################
## Intialize analysis 
##############################################################################
# Empty workspace and load modules
rm(list= ls())
require('scran')
require('ggplot2')
require('Rtsne')
require('scde')
require('DESeq2')
require('limma')
source('/farm/home/rabino01/github/ngs_r/ggplot/plotThemes.R')
source('/farm/home/rabino01/github/ngs_r/kmeans.R')
# Set analysis parameters
acceptObsv <- c('1')
acceptDates <- c(110216, 150216, 140316)
minExpressed <- 2000
minAligned <- 250000
maxRRNA <- 0.1
minERCC <- 0.01
maxERCC <- 0.5
minRatioNorm <- 0.5
minRatioCluster <- 0.25
minRatioSCDE <- 0.25
# Paths
sdPath <- '/farm/home/rabino01/Calado/completeSingleCellData.txt'
tpmPath <- '/farm/scratch/rs-bio-lif/rabino01/Calado/singleCell/GRCm38.84_combined_gene.tpm'
expPath <- '/farm/scratch/rs-bio-lif/rabino01/Calado/singleCell/GRCm38.84_combined_gene.exp'
qcPath <- '/farm/scratch/rs-bio-lif/rabino01/Calado/singleCell/GRCm38.84_combined.qc'
erccPath <- '/farm/home/rabino01/ERCC/ercc.conc.txt'
biotypePath <- '/farm/scratch/rs-bio-lif/rabino01/Calado/genomeData/GRCm38.84/Mus_musculus.GRCm38.84.primary_assembly_exon.biotype'
outPrefix <- '/farm/home/rabino01/Calado/GRCm38.84'
geneLists <- list(
  'Myc' = 'ENSMUSG00000022346',
  'Xbp1' = 'ENSMUSG00000020484',
  'Irf4' = 'ENSMUSG00000021356',
  'Prdm1'= 'ENSMUSG00000038151',
  'Polh' = 'ENSMUSG00000023953',
  'Ccnd2' = 'ENSMUSG00000000184',
  'Cd83' = 'ENSMUSG00000015396',
  'Nfkbia' = 'ENSMUSG00000021025',
  'Cd86' = 'ENSMUSG00000022901',
  'Cd40' = 'ENSMUSG00000017652',
  'eGFP' = 'eGFP',
  'Gapdh' = 'ENSMUSG00000057666',
  'Actb' = 'ENSMUSG00000029580',
  'Hprt' = 'ENSMUSG00000025630'
)

##############################################################################
## Extract acceptable samples
##############################################################################
# Read in sample data and find single cell samples
sd <- read.table(sdPath, header = T, sep = '\t', row.names = 1,
  stringsAsFactors = F)
# Read in expected data
exp <- read.table(expPath, header = T, sep = '\t', row.names = 1,
  stringsAsFactors = F, check.names = F)
# Read in biptype data amd extract rRNA genes
biotype <- read.table(biotypePath)
rRNAgenes <- subset(biotype[,1], biotype[,2] == 'rRNA')
# Extract gene types
geneTypes <- factor( c( ENSM="MM", ERCC="ERCC", eGFP="eGFP")[
  substr( rownames(exp), 1, 4 ) ] )
mouseGenes <- row.names(exp)[geneTypes == 'MM']
erccGenes <- row.names(exp)[geneTypes == 'ERCC']

###############################################################################
## Generate qc data
###############################################################################
qc <- data.frame(row.names = colnames(exp))
# Extract gene count
qc$expressed_genes <- apply(exp[mouseGenes,], 2, function(z) {sum(z > 0)})
# Extract gene read count
qc$read_count <- apply(exp[mouseGenes,], 2, sum)
# Extract gene read count
qc$gene_read_count <- apply(exp[mouseGenes,], 2, sum)
# Extract errc read count
qc$ercc_read_count <- apply(exp[erccGenes,], 2, sum)
# Extract ercc ratio
qc$ercc_ratio <- qc$ercc_read_count / 
  (qc$ercc_read_count + qc$gene_read_count)
# Extract rRNA read count
qc$rrna_read_count <- apply(exp[rRNAgenes,], 2, sum)
# Extract rRNA ration
qc$rrna_ratio <- qc$rrna_read_count / qc$gene_read_count

################################################################################
## Extract passed samples and normalize counts
################################################################################
# Accepted genes
passedSamples <- intersect(
  subset(
    row.names(qc),
    qc$expressed_genes >= minExpressed &
    qc$gene_read_count >= minAligned &
    qc$rrna_ratio <= 0.05 &
    qc$ercc_ratio >= 0.1 &
    qc$ercc_ratio <= 0.5
  ),
  subset(
    row.names(sd),
    sd$obsv %in% acceptObsv
  )
)

###############################################################################
## Find most variable genes
###############################################################################
# Extract size factors using ERCC controls
expERCC <- exp[geneTypes == 'ERCC', passedSamples]
erccSizeFactors <- estimateSizeFactorsForMatrix(expERCC[rowMeans(expERCC) > 100,])
# Normalize, using ERCC factors, and log counts
normExp <- t(t(exp[,passedSamples]) /erccSizeFactors)
log2NormExp <- log2(normExp + 1)
# Create sce object
sce <- newSCESet(exprsData=log2NormExp)
isSpike(sce) <- geneTypes == 'ERCC'
# Perform fitting
fit <- trendVar(sce, use.spikes=T, assay = 'exprs', trend = 'loess')
decomp <- decomposeVar(sce, fit)
# Extract variant genes
varGenes <- subset(row.names(decomp), decomp$bio > 0 )
widelyExpressed <- row.names(log2NormExp)[
  apply(log2NormExp, 1, function(z) {
    (sum(z > 0) / length(z)) >= 0.5
  })
]
corGenes <- intersect(widelyExpressed, setdiff(varGenes, rRNAgenes))
log2CorGenes <- log2NormExp[corGenes,]
# Plot results
plot(
  decomp$mean,
  decomp$total,
  xlab="Mean log-expression",
  ylab="Variance",
  pch=19,
  cex=0.2
)
o <- order(decomp$mean)
lines(decomp$mean[o], decomp$tech[o], col="red", lwd=2)
points(fit$mean, fit$var, col="red", pch=19, cex=.5)

corGeneBiotype <- biotype$V2[match(corGenes, biotype$V1)]


###############################################################################
## Generate correlation matrix
###############################################################################
exprCorGenes <- apply(log2CorGenes, 2, function(z) {sum(z > 0)})
design <- model.matrix( ~ exprCorGenes )
fit <- lmFit(log2CorGenes, design )
resid <- residuals( fit, log2CorGenes )
corMatrix <- cor(log2CorGenes)

###############################################################################
## Perform PCA
###############################################################################
pca <- prcomp(t(log2CorGenes))
varContribution <- (pca$sdev^2) / sum(pca$sdev^2)
barplot(varContribution[1:10])
pcaPlotData <- cbind(
  pca$x[,1:10],
  qc[
    row.names(pca$x),
    c('rrna_ratio', 'gene_read_count', 'expressed_genes')
    ]
)
ggplot(pcaPlotData, aes(x = PC1, y = PC2, colour = rrna_ratio)) +
  geom_point() +
  scale_colour_gradientn(colours=rainbow(4)) +
  pdfA4SquareTheme
testDF <- data.frame(
  'pc' = rep(c('PC1','PC2','PC3','PC4','PC5'), each = nrow(pcaPlotData)),
  'value' = c(pca$x[,'PC1'], pca$x[,'PC2'], pca$x[,'PC3'], pca$x[,'PC4'], pca$x[,'PC5']),
  'exprs' = rep(pcaPlotData$expressed_genes, 5),
  'rrna' = rep(pcaPlotData$rrna_ratio, 5)
)
p<-ggplot(testDF, aes(x=pc, y=value, color=rrna)) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_colour_gradientn(colours=rainbow(4))
p


##############################################################################
## Perform tsne analyis and create plot data
##############################################################################
# Perform tsne
set.seed(1)
tsneData <- Rtsne(t(log2CorGenes), dims = 2, theta = 0.0, pca = FALSE)
# Create tsne plot data
tsnePlotData <- cbind(
  sd[
    passedSamples,
    c('date', 'plate', 'well', 'obsv')
  ],
  qc[
    passedSamples,
    c('rrna_ratio', 'gene_read_count', 'expressed_genes')
  ]
)
tsnePlotData$expressed_cor_genes <- exprCorGenes
tsnePlotData[,c('x','y')] <- tsneData$Y
# Set date as factor
tsnePlotData$date <- factor(
  as.character(tsnePlotData$date),
  levels = as.character(acceptDates)
)


##############################################################################
## Create plots for experimental variables
##############################################################################
# Open pdf
pdf(paste0(outPrefix,'.correlates.pdf'), paper = 'a4r', onefile = T,
    width = 10, height = 7)
# Create plot based on sample date
ggplot(tsnePlotData, aes(x = x, y = y, colour = date)) +
  geom_point() +
  pdfA4SquareTheme
# Plot observation
ggplot(tsnePlotData, aes(x = x, y = y, colour = obsv)) +
  geom_point() +
  pdfA4SquareTheme
# Plot aligned counts
ggplot(tsnePlotData, aes(x=x, y=y, colour=gene_read_count)) +
  geom_point() + 
  scale_colour_gradientn(colours=rainbow(4)) +
  pdfA4SquareTheme
# Plot rRNA rate
ggplot(tsnePlotData, aes(x=x, y=y, colour=rrna_ratio)) +
  geom_point() + 
  scale_colour_gradientn(colours=rainbow(4)) +
  pdfA4SquareTheme
# Plot expressed genes
ggplot(tsnePlotData, aes(x=x, y=y, colour=expressed_genes)) +
  geom_point() + 
  scale_colour_gradientn(colours=rainbow(4)) +
  pdfA4SquareTheme
ggplot(tsnePlotData, aes(x=x, y=y, colour=expressed_cor_genes)) +
  geom_point() + 
  scale_colour_gradientn(colours=rainbow(4)) +
  pdfA4SquareTheme
dev.off()

##############################################################################
## Perform kmeans clustering
##############################################################################
# Extract cluster stability values for different k
kData <- kmeansKChoice(tsnePlotData[,c('x','y')])
# Plot cluster stability values
kmeansChoicePlot(kData, paste0(outPrefix, '.kmeans.stability.pdf'))

##############################################################################
## Plot data for k of 4
##############################################################################
tsnePlotData$cluster <- kData$clusters$k4
ggplot(tsnePlotData, aes(x=rRNA_rate, y=expressed_genes,
  colour=cluster)) +
  geom_point() + 
  pdfA4SquareTheme
ggplot(tsnePlotData, aes(x=rRNA_rate, colour=cluster)) +
  geom_density()
t.test(
  subset(
    tsnePlotData$rRNA_rate,
    tsnePlotData$cluster == 1 
  ),
  subset(
    tsnePlotData$rRNA_rate,
    tsnePlotData$cluster == 2 
  )
)

##############################################################################
## Perform scde
##############################################################################
# Generate pair comparisons
cmp <- data.frame(
  condition = c('control', 'test'),
  cluster = c(1,2)
)
# Extract counts
controlSamples <- row.names(
  tsnePlotData[tsnePlotData$cluster %in% cmp$cluster[cmp$condition == 'control'],]
)
testSamples <- row.names(
  tsnePlotData[tsnePlotData$cluster %in% cmp$cluster[cmp$condition == 'test'],]
)
allSamples <- c(controlSamples, testSamples)
# Find acceptable genes
genePassSCDE <- row.names(exp[,allSamples])[
  apply(exp[,allSamples], 1, function(z) {
    sum(z > 0) / length(z)
  }) >= minRatioSCDE
]
# Extract counts
expInt <- round(exp[genePassSCDE,allSamples])
expInt <- apply(
  expInt, 2, function(z) {
    storage.mode(z) <- 'integer'; z
  }
)
o.ifm <- scde.error.models(
  counts = expInt,
  n.cores = 3
)


