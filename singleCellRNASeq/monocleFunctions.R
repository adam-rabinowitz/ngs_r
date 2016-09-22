###############################################################################
## Function to extract filter gene data for tow specified states
###############################################################################
extractGeneStateData <- function(cds, state1, state2) {
  # Check states
  states <- as.numeric(unique(pData(cds)$State))
  states <- states[order(states)]
  if (! state1 %in% states) {
    stop('state1 not a valid state')
  }
  if (! state2 %in% states) {
    stop('state2 not a valid state')
  }
  # Extract samples
  samples <- row.names(pData(cds))[
    pData(cds)$State %in% c(state1, state2)
    ]
  # Extract filter genes
  filterGenes <- fData(cds)[
    fData(cds)$use_for_ordering == T,
    'gene_short_name',
    drop = F
    ]
  # Extract gene expression data and melt
  plotData <- exprs(cds)[
    row.names(filterGenes),
    samples
    ]
  plotData <- melt(plotData)
  colnames(plotData) <- c('name', 'sample', 'value')
  # Add state, short-name and pseudotime
  plotData$state <- pData(cds)[
    as.character(plotData$sample), 'State']
  plotData$short_name <- filterGenes[
    as.character(plotData$name), 'gene_short_name']
  plotData$time <- pData(cds)[
    as.character(plotData$sample), 'Pseudotime']
  # Add colour data
  hues = seq(15, 375, length = length(states) + 1)
  cols = hcl(h = hues, l = 65, c = 100)[1:length(states)]
  plotData$col <- cols[as.numeric(plotData$state)]
  # Return data
  return(plotData)
}

###############################################################################
## Function to plot gene-state against time
###############################################################################
plotGeneStateTime <- function(cds, state1, state2) {
  # Extract plot data
  plotData <- extractGeneStateData(cds, state1, state2)
  # Extract colours
  states <- unique(plotData$state)
  cols <- sapply(states, function(z) {
    plotData$col[match(z, plotData$state)]
  })
  names(cols) <- states
  # Create and return plot
  p <- ggplot(plotData, aes(x = time, y = log10(round(value) + 1))) +
    geom_point(aes(colour = state)) +
    facet_wrap(~ short_name, scales = 'free_y') +
    scale_color_manual(values = cols) +
    scale_y_continuous(breaks=c(0,1,2,3,4,5,6)) +
    stat_smooth(method = 'glm', method.args = list(family = 'quasipoisson'),
      formula = y ~ ns(x, 2), se = F) +
    labs(title = paste('State', state1, 'to State', state2),
      y = 'Log10 Expression', x = 'Pseudotime')
  return(p)
}

###############################################################################
## Function to plot gene expression against time
###############################################################################
plotAllGeneStateTime <- function(cds, fileName) {
  # Extract states
  states <- unique(pData(cds)$State)
  if (is.null(states)) {
    stop('States must be set')
  }
  states <- states[order(states)]
  # Generate list of comparisons
  comparison <- list()
  for (i in 2:length(states)) {
    name <- paste0(1,'vs',i)
    comparison[[name]] <- c(1,i)
  }
  # Create list of plots
  plots = list(span = plot_spanning_tree(cds))
  for (i in 1:length(comparison)) {
    name = names(comparison)[i]
    state1 = comparison[[i]][1]
    state2 = comparison[[i]][2]
    plots[[name]] = plotGeneStateTime(cds, state1, state2)
  }
  # Draw plots to file
  pdf(fileName, onefile = T, paper = 'a4r', width = 10, height = 7)
  lapply(plots, print)
  dev.off()
}

###############################################################################
## Function to plot gene-state
###############################################################################
plotGeneState <- function(cds, state1, state2) {
  # Extract plot data and add log of values
  plotData <- extractGeneStateData(cds, state1, state2)
  plotData$logvalue <- log2(plotData$value + 1)
  # Extract colours
  states <- unique(plotData$state)
  cols <- sapply(states, function(z) {
    plotData$col[match(z, plotData$state)]
  })
  names(cols) <- states
  # Create and return plot
  p <- ggplot(plotData, aes(x=state, y = logvalue, fill=state, color = state)) +
    geom_jitter(width = 0.6, height = 0) +
    facet_wrap(~ short_name) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(breaks=(seq(0,20,2))) +
    stat_summary(fun.y = mean, geom = 'point', col = 'black') +
    labs(title = paste('State', state1, 'to State', state2),
      y = 'Log2 Expression', x = 'State') +
    theme(legend.position="none", panel.grid.minor = element_blank())
  return(p)  
}

###############################################################################
## Function to plot gene expression against time
###############################################################################
plotAllGeneState <- function(cds, fileName) {
  # Extract states
  states <- unique(pData(cds)$State)
  if (is.null(states)) {
    stop('States must be set')
  }
  states <- states[order(states)]
  # Generate list of comparisons
  comparison <- combn(states, 2)
  colnames(comparison) <- apply(comparison, 2, function(z) {
    paste(z, collapse='vs')})
  # Create list of plots
  plots = list(span = plot_spanning_tree(cds))
  for (i in 1:ncol(comparison)) {
    name = colnames(comparison)[i]
    state1 = as.numeric(comparison[1,i])
    state2 = as.numeric(comparison[2,i])
    plots[[name]] = plotGeneState(cds, state1, state2)
  }
  # Draw plots to file
  pdf(fileName, onefile = T, paper = 'a4r', width = 10, height = 7)
  lapply(plots, print)
  dev.off()
}

###############################################################################
## Function to plot gene expression against state
###############################################################################
plotAllGeneStateData <- function(cds, fileName) {
  # Extract states
  states <- unique(pData(cds)$State)
  if (is.null(states)) {
    stop('States must be set')
  }
  states <- states[order(states)]
  # Generate list of comparisons
  comparison <- list()
  for (i in 2:length(states)) {
    name <- paste0(1,'vs',i)
    comparison[[name]] <- c(1,i)
  }
  # Create list of plots
  plots = list(span = plot_spanning_tree(cds))
  for (i in 1:ncol(comparison)) {
    print(i)
    name = names(comparison)[i]
    state1 = comparison[[i]][1]
    state2 = comparison[[i]][2]
    plots[[name]] = plotGeneStateData(cds, state1, state2)
  }
  # Draw plots to file
  pdf(fileName, onefile = T, paper = 'a4r', width = 10, height = 7)
  lapply(plots, print)
  dev.off()
}

###############################################################################
## Generate cds objects with ordering filters
###############################################################################
createFilterCds <- function(cellData, featureData, phenoData, filterGeneList,
  minCellExp, numberPaths, seed = 1066) {
  # Create cds object and count gene expression
  cds <- newCellDataSet(cellData = cellData, phenoData = phenoData,
    featureData = featureData)
  cds <- detectGenes(cds)
  expressed_genes <- row.names(
    subset(fData(cds), num_cells_expressed >= minCellExp))
  # Find structure in data using filter genes
  cdsGoi <- lapply(filterGeneList, function(z) {
    filterGenes <- intersect(z, expressed_genes)
    cdsGoi <- setOrderingFilter(cds, filterGenes)
    set.seed(seed)
    cdsGoi <- reduceDimension(cdsGoi, use_irlba=FALSE)
    cdsGoi <- orderCells(cdsGoi, num_paths=numberPaths, reverse=TRUE)
  })
  # Return data
  return(cdsGoi)
}

###############################################################################
## Extract states
###############################################################################
extractCdsStates <- function(cdsList) {
  output = data.frame(row.names = row.names(pData(cdsList[[1]])))
  stateDF <- for (i in 1:length(cdsGoi)) {
    df <- pData(cdsGoi[[i]])
    name <- names(cdsGoi)[i]
    output[[name]] <- df[row.names(output),'State']
  }
  return(output)
}

###############################################################################
## Create state overlap
###############################################################################
plotStateOverlap  <- function(stateDF) {
  if (ncol(stateDF) != 2) {
    stop('input dataframe should have two columns')
  }
  plotDF <- melt(table(stateDF[,1], stateDF[,2]))
  plotDF <- plotDF[plotDF$value > 0,]
  p <- ggplot(plotDF, aes(x = Var1, y = Var2, size = value)) +
    geom_point() +
    labs(x = colnames(stateDF)[1], y = colnames(stateDF)[2])
  return(p)
}