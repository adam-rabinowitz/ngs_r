require('scde')
require('WriteXLS')

##############################################################################
## Function to extract clusters from a named character vector
##############################################################################
extractClusters <- function(clusters) {
  clusterValues <- unique(clusters)
  clusterValues <- clusterValues[order(clusterValues)]
  comparisons <- lapply(
    clusterValues,
    function(z) {
      factors <- factor(
        clusters == z,
        levels=c(TRUE,FALSE),
        labels=c(paste0('cl.',z), paste0('cl.!',z))
      )
      return(factors)
    }
  )
  if (length(comparisons) == 2) {
    comparisons <- comparisons[1]
  }
}

###############################################################################
## Perform scde
###############################################################################
performScde <- function(counts, clusters, outFile, cores = 6) {
  # Generate comparisons
  comparisons <- extractClusters(clusters)
  # Convert counts to integers
  counts <- apply(
    counts, 2, function(z) {
      storage.mode(z) <- 'integer'; z}
  )
  # Perform scde
  scdeResults <- lapply(
    comparisons, function(z) {
      # Create error models
      o.ifm <- scde.error.models(
        counts = counts, groups = z, n.cores = cores,
        threshold.segmentation = TRUE, save.crossfit.plots = FALSE,
        save.model.plots = FALSE, verbose = 1)
      # Create priors
      o.prior <- scde.expression.prior(
        models = o.ifm, counts = counts, length.out = 400, show.plot = F)
      # Perform differential expression
      ediff <- scde.expression.difference(
        models = o.ifm, counts = counts, prior = o.prior, groups = z,
        n.randomizations = 100, n.cores = cores, verbose = 1,
        return.posteriors = T)
      # Add raw and adjusted p-values
      ediff$results$cZ <- NULL
      ediff$results$pval <- 2*pnorm(-abs(ediff$results$Z))
      ediff$results$padj <- p.adjust(ediff$results$pval, method = 'BH')
      ediff$results <- ediff$results[order(ediff$results$padj,
        ediff$results$pval),]
      return(ediff$results)
    }
  )
  WriteXLS(
    'scdeResults',
    ExcelFileName = outFile,
    row.names=T
  )
}

# load example dataset
data(es.mef.small)
# factor determining cell types
clusters <- gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small))
names(clusters) <- colnames(es.mef.small)
# clean up the dataset
counts <- es.mef.small[rowMeans(es.mef.small) > 1,]
performScde(counts, clusters, '/farm/home/rabino01/Calado/scdeTest.xls', 6)


