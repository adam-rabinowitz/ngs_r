##############################################################################
## filter RNA-Seq counts
##############################################################################
filterTPMCounts <- function(tpm, exp, reads=400000, genes=2000, mintpm=10,
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
## Filter Exp counts
###############################################################################
filterExpCounts <- function(tpm, exp, reads=400000, genes=2000, mintpm=10,
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
  filteredExp <- exp[acceptedGenes,acceptedSamples]
  filteredExp <- as.matrix(filteredExp)
  return(filteredExp)
}

##############################################################################
## Create gene translation vecotr
##############################################################################
# Extract file for gene name conversion
createTran <- function(tranfile) {
  # Read in tran data
  tranData <- read.table(tranfile, sep = '\t', header = T,
    stringsAsFactors=F)
  # Split alternative names by source name
  tranData <- split(tranData[,2], tranData[,1])
  # Create character vecor and return
  tranData <- sapply(tranData, paste0, collapse=' ')
  return(tranData)
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

###############################################################################
## Function to perform scde analysis
###############################################################################
performSCDE <- function(counts, groups, altname=NULL, cores=4) {
  # Load required module
  require('scde')
  # Check groups
  groups <- checkGroups(colnames(counts), groups)
  # Check that each group has two factors
  factorCount <- sapply(groups, function(z) {length(levels(z))})
  countTable <- table(factorCount)
  if (countTable['2'] != sum(countTable)) {
    stop('groups must have exactly two factors')
  }
  # Round counts
  counts <- round(counts,digits=0)
  storage.mode(counts) <- 'integer'
  # Perform scde
  results <- lapply(groups, function(z){
    # Build models
    scde.fitted.model <- scde.error.models(
      counts=counts,
      n.cores=cores,
      save.model.plots=F
    )
    scde.prior <- scde.expression.prior(
      models=scde.fitted.model,
      counts=counts
    )
    # Calculate differential expression
    ediff <- scde.expression.difference(
      scde.fitted.model,
      counts,
      scde.prior,
      groups = groups[[1]],
      n.cores = 4
    )
    # Calculate p-values
    ediff$pvalue <- 2*pnorm(abs(ediff$Z),lower.tail=F) # 2-tailed p-value
    ediff$padj <- p.adjust(ediff$pvalue,'BH')
    ediff <- ediff[order(ediff$padj),]
    # Add gene names and adjust column order
    ediff$gene <- row.names(ediff)
    if (!is.null(altname)) {
      ediff$altname <- altname[ediff$gene]
      ediff <- ediff[,c(9:10,1:8)]
    } else {
      ediff <- ediff[,c(9,1:8)]
    }
    # rename rows and return
    row.names(ediff) <- 1:nrow(ediff)
    return(ediff)
  })
  return(results)
}






