##############################################################################
## filter RNA-Seq counts
##############################################################################
filterRNACounts <- function(tpm, exp, reads=400000, genes=2000, mintpm=10,
  minratio=0.1, minsample=5) {
  # Read in data if files supplied
  if (is.character(tpm)) {
    message('Reading TPM count file')
    tpm = read.table(tpm, sep='\t', header=T, row.names=1)
  }
  if (is.character(exp)) {
    message('Reading expected count file')
    exp = read.table(exp, sep='\t', header=T, row.names=1)
  }
  # Extract samples with sufficient aligned reads and expressed genes
  acceptedSamples <- colnames(exp)[
    colSums(exp) >= reads &
      apply(exp, 2, function(z) {sum(z > 0)}) >= genes
    ]
  message(paste(length(acceptedSamples), 'of', ncol(tpm),
    'samples accepted'))
  # Calculate number of samples in which a gene must be expresed
  sampleNo <- ceiling(max(minsample, length(acceptedSamples) * minratio))
  message(paste('Finding Genes expressed in >=', sampleNo, 'samples'))
  # Extract accepted genes
  acceptedGenes <- row.names(tpm)[
    apply(
      tpm[,acceptedSamples], 1, function(z) {sum(z >= mintpm)}
    ) >= sampleNo
    ]
  message(paste(length(acceptedGenes), 'of', nrow(tpm),
    'genes accepted'))
  # Extract and return filtered counts as matrix
  filteredTPM <- tpm[acceptedGenes,acceptedSamples]
  filteredTPM <- as.matrix(filteredTPM)
  return(filteredTPM)
}

##############################################################################
## Check concordane between samples and group list
##############################################################################
checkGroups <- function(samples, groups) {
  # Sort sample names
  samplesSorted <- sort(samples)
  # Loop through group list
  lapply(groups, function(g) {
    # Check groups are character vectors
    if (!is.character(g) & !is.factor(g) & !is.numeric(g)) {
      stop('groups must be a list of character vectors') 
    }
    # Check group names are identical to sample names
    gname = sort(names(g))
    if (!identical(gname, samplesSorted)) {
      stop('Group names and sample names are not identical')
    }
    # Sort and factor groups
    g <- g[samples]
    g <- factor(g, levels = unique(g))
    return(g)
  })
}
