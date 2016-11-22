# Load packages
require(ggplot2)
require(reshape2)
require(IRanges)
require(zoo)

###############################################################################
## Window shift detection of tad
###############################################################################
find.tad.window <- function(log2.data, half.window = 3, na.rm = T) {
  # Create empty data frame for output
  output <- data.frame()
  # Loop through windows
  for(index in 1:(nrow(log2.data) - half.window)) {
    # Extract data for upstream bins
    up.indices <- index : (index - 1 + half.window)
    up.data <- log2.data[up.indices,]
    if (all(is.na(up.data$log2))) {
      next
    }
    # Extract max log2 value for upstream bins
    max.log2 <- max(up.data$log2, na.rm = na.rm)
    if (is.na(max.log2)) {
      next
    }
    # Extract data for max log2 value
    max.indices <- which(up.data$log2 == max.log2)
    up.data.max <- up.data[tail(max.indices, 1),]
    # Extract data and min log2 for downstream bins
    down.indices <- (index + half.window) : (index + (2 * half.window) - 1)
    down.data <- log2.data[down.indices,]
    if (all(is.na(down.data$log2))) {
      next
    }
    # Extract min log2 value for downstream bins
    min.log2 <- min(down.data$log2, na.rm = na.rm)
    if (is.na(min.log2)) {
      next
    }
    # Extract data for min log2 value
    min.indices <- which(down.data$log2 == min.log2)
    down.data.min <- down.data[head(min.indices, 1),]
    # Calculate switch values and skip if negative
    log2.switch <- up.data.max$log2 - down.data.min$log2
    if (log2.switch < 0) {
      next
    }
    # Create and append dataframe for output
    tad.window <- data.frame(
      'chr' = up.data.max$chr,
      'start' = up.data.max$start,
      'end' = down.data.min$end,
      'centre' = (up.data.max$start + down.data.min$end) / 2,
      'max' = up.data.max$log2,
      'min' = down.data.min$log2,
      'shift' = up.data.max$log2 - down.data.min$log2
    )
    output <- rbind(output, tad.window)
  }
  # Filter and order output
  output <- output[order(-output$shift),]
  return(output)
}

###############################################################################
## Funcion to filter tads identified through window search
###############################################################################
filter.tad.window <- function(tad.data, min.shift = 0) {
  # Order and filter data by shift
  tad.data <- tad.data[order(-tad.data$shift),]
  tad.data <- subset(tad.data, tad.data$shift >= min.shift)
  # Build iranges object and find overlaps
  tad.ranges <- IRanges(tad.data$start, tad.data$end)
  tad.overlaps <- findOverlaps(tad.ranges, tad.ranges)
  # Find intervals overlapping higher shift intervals
  lower.shift.overlaps <- subjectHits(tad.overlaps) > queryHits(tad.overlaps)
  remove.rows <- subjectHits(tad.overlaps)[lower.shift.overlaps]
  remove.rows <- unique(remove.rows)
  remove.rows <- remove.rows[order(remove.rows)]
  # Remove unwanted intervals and return data
  tad.data <- tad.data[-remove.rows,]
  rownames(tad.data) <- 1:nrow(tad.data)
  return(tad.data)
}

###############################################################################
## Function to find and filter tads
###############################################################################
find.filter.tad.window <- function(
  log2.data, half.window = 3, min.shift=0, na.rm = T
) {
  tad.data <- find.tad.window(log2.data, half.window, na.rm)
  tad.filtered <- filter.tad.window(tad.data, min.shift)
  return(tad.filtered)
}


###############################################################################
## Function to create plot
###############################################################################
create_prob_plot <- function(
    norm_data, tad_data = NA, max_prob = 1e-1, min_prob = 1e-6
  ) {
  # Check min and max prob values
  if (log10(max_prob) %% 1 > 0) {
    stop('incorrect max_prob value')
  }
  if (log10(min_prob) %% 1 > 0) {
    stop('incorrect min_prob value')
  }
  plot_breaks <- 10^(seq(log10(min_prob), log10(max_prob)))
  plot_labels <- paste0('1e', log10(plot_breaks))
  plot_labels[1] <- paste0('<=', plot_labels[1])
  plot_labels[length(plot_labels)] <- paste0('>=', plot_labels[length(plot_labels)])
  # Copy plot data and rename rows and columns
  plot_data <- norm_data
  name_list <- strsplit(colnames(plot_data), '[:-]')
  centre <- sapply(name_list, function(z) {mean(as.numeric(z[2:3]))})
  colnames(plot_data) <- centre
  rownames(plot_data) <- centre
  # Manipulate plot data
  plot_data[is.na(plot_data)] <- min_prob
  plot_data[plot_data < min_prob] <- min_prob
  plot_data[plot_data > max_prob] <- max_prob
  plot_data <- as.matrix(plot_data)
  plot_data <- melt(plot_data)
  # Create plot
  p <- ggplot(plot_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_raster() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient(
      name = "prob",
      limits = c(min_prob, max_prob),
      trans = "log",
      breaks = plot_breaks,
      labels = plot_labels,
      low = 'black',
      high = 'cyan2'
    ) +
    coord_fixed(ratio = 1)
  if (!is.na(tad_data)) {
    p <- p + geom_point(
      data = tad_data,
      aes(x = centre, y = centre),
      inherit.aes = F,
      fill = 'red',
      colour = 'red'
    )
  }
  return(p)
}

###############################################################################
## Function to find and annotate tad boundaries
###############################################################################
annotate_tad_boundary <- function(log2_data) {
  # Check all intervals are on same chromosome
  chrom <- unique(log2_data$chrom)
  if (length(chrom) != 1) {
    stop('log2 values must be from same chromosome')
  }
  # Using sign of log2 identify last base in tad
  sign_change <- diff(sign(log2_ratio$log2))
  sign_change[is.na(sign_change)] <- 0
  last_in_tad <- sign_change == -2
  # Split dataframe into tads
  tad_factor <- factor(cumsum(c(F,last_in_tad)))
  tad_list <- split(log2_data, tad_factor)
  # Extract positions of tad
  position <- sapply(tad_list, function(z) {tail(z$end,1)})
  position <- head(position, -1)
  # Extract log2 values proximal to tad
  prox_up <- sapply(tad_list, function(z) {tail(z$log2,1)})
  prox_up <- head(prox_up, -1)
  prox_down <- sapply(tad_list, function(z) {head(z$log2,1)})
  prox_down <- tail(prox_down, -1)
  prox_switch <- prox_up - prox_down
  # Calculate log2 values across tad
  max_up <- sapply(tad_list, function(z) {max(z$log2)})
  max_up <- head(max_up, -1)
  max_down <- sapply(tad_list, function(z) {min(z$log2)})
  max_down <- tail(max_down, -1)
  max_switch <- max_up - max_down
  # Create output dataframe
  out_df <- data.frame(
    'chrom' = chrom,
    'pos' = position,
    'prox_up' = prox_up,
    'prox_down' = prox_down,
    'prox_switch' = prox_switch,
    'max_up' = max_up,
    'max_down' = max_down,
    'max_switch' = max_switch
  )
  return(out_df)
}

###############################################################################
## Extract diagonals
###############################################################################
extract_diagonals <- function(data_matrix) {
  # Generate factors for diagonals and split matrix
  diag_factors <- col(data_matrix) - row(data_matrix) + 1
  diag_data <- split(data_matrix, diag_factors)
  # Select, reorder and rename diagonals
  diag_data <- rev(diag_data[1:ncol(data_matrix)])
  names(diag_data) <- 1:ncol(data_matrix)
  return(diag_data)
}

###############################################################################
## Function to calculate values for diagonals
###############################################################################
calc_diagonals <- function(
  data_matrix, diagonals, calc_function, rm_zero = T
) {
  # Remove zeros if required
  if (rm_zero) {
    data_matrix[data_matrix == 0] <- NA
  }
  # Extract diagonals and rename them
  diag <- extract_diagonals(data_matrix)
  bin_names <- colnames(data_matrix)
  diag <- lapply(diag, function(z) {
    names(z) <- head(bin_names, length(z))
    return(z)
  })
  select_diag <- diag[diagonals]
  print(select_diag)
  diag_df <- data.frame(
    bin = unlist(sapply(select_diag, names)), 
    value = unlist(select_diag)
  )
  data_list <- do
  output <- tapply(diag_df$value, diag_df$bin, calc_function)
  return(output)
}

###############################################################################
## Extract squares of diagonal
###############################################################################
extract_squares <- function(
  data_matrix, size, offset = 1, rm_zero_bins = F, rm_zero_values = T
) {
  # Process zero bins
  if (rm_zero_bins) {
    low_bins <- which(colSums(data_matrix) < 0.9)
    data_matrix <- data_matrix[-low_indices, -low_indices]
  }
  # Process zero values
  if (rm_zero_values) {
    data_matrix[data_matrix == 0] <- NA
  }
  # Check object is square
  dimensions <- dim(dhalf_window = 3, na.rm = Tata_matrix)
  if (dimensions[1] != dimensions[2]) {
    stop('matrix not square')
  }
  # Calculate corners of sub-matrices
  col_starts <- ((2*offset) + size):((dimensions[1] + 1) - size)
  row_starts <- 1:length(col_starts)
  col_ends <- col_starts + (size - 1)
  row_ends <- row_starts + (size - 1)
  # Extract and return matrices
  output <- list()
  for (i in 1:length(col_starts)) {
    output[[i]] <- data_matrix[
      row_starts[i]:row_ends[i],
      col_starts[i]:col_ends[i]
    ]
  }
  # Rename output and return data
  bin_names <- colnames(data_matrix)
  names(output) <- head(
    tail(
      bin_names,
      -(size + offset - 1)
    ), 
    length(output)
  )
  return(output)
}

###############################################################################
## Calculate metrics for offset
###############################################################################
extract_square_metrics <- function(matrix_list, span = 0.75, half_window = 15) {
  # Create dataframe summarising matrices
  split_names <- strsplit(names(matrix_list), ('-|:'))
  centre <- sapply(split_names, function(z) {mean(as.numeric(z[2:3]))})
  matrix_df <- data.frame(
    row.names = 1:length(matrix_list),
    bin_no = 1:length(matrix_list),
    name = names(matrix_list),
    centre = centre,
    median = sapply(matrix_list, function(z) {median(z, na.rm = T)})
  )
  # Fit loess and find difference
  matrix_df$loess <- loess(median ~ centre, matrix_df, span = span)$fitted
  matrix_df$diff <- matrix_df$median - matrix_df$loess
  # Find slope upstream
  slope <- rollapply(matrix_df$median, half_window, function(z) {
    df <- data.frame(
      x = 1:half_window,
      y = z
    )
    fit <- lm(y ~ x, df)
    slope <- sign(fit$coefficients['x'])
    return(slope)
  })
  matrix_df$up_slope <- c(rep(NA, half_window - 1), slope)
  matrix_df$down_slope <- c(slope, rep(NA, half_window - 1))
  # Return data
  return(matrix_df)
}

###############################################################################
## Function to filter minima 
###############################################################################
filter_square_metrics <- function(minima, half_window = 16) {
  # Remove non-minima bins
  filter_minima <- square_metrics[
    which((square_metrics$up_slope - square_metrics$down_slope) == -2),
  ]
  # Extract minium value for each minima
  filter_minima <- split(
    filter_minima,
    c(0, cumsum(diff(filter_minima$bin_no) > 1))
  )
  filter_minima <- lapply(filter_minima, function(z) {
    print(z)
    return(z[which.min(z$median),])
  })
  filter_minima <- do.call(rbind, filter_minima)
  # Remove minima above trend, sort and return
  filter_minima <- subset(filter_minima, filter_minima$diff < 0)
  row.names(filter_minima) <- 1:nrow(filter_minima)
  return(filter_minima)
}

###############################################################################
## Calculate loess fit and standard deviation of probability to distance
###############################################################################
extract_distance_loess <- function(
  data.matrix, zero.values = F, zero.dist = F, span = 0.01
) {
  # Create distance matrix
  split.names <- strsplit(colnames(data.matrix), ('-|:'))
  centre <- sapply(split.names, function(z) {mean(as.numeric(z[2:3]))})
  dist.matrix <- abs(outer(centre, centre, '-'))
  # Create and process distance vs proablity dataframe
  dist.df <- data.frame(
    dist = dist.matrix[upper.tri(dist.matrix, diag = zero.dist)],
    prob = data.matrix[upper.tri(data.matrix, diag = zero.dist)]
  )
  dist.df <- dist.df[order(dist.df$dist),]
  if (!zero.values) {
    dist.df <- dist.df[dist.df$prob > 0,]
  }
  # Perform loess and extract results into dataframe
  loess.data <- loess(prob ~ dist, dist.df, span = span)
  loess.df <- data.frame(
    dist = dist.df$dist,
    fit = loess.data$fitted,
    residual = loess.data$residuals
  )
  # Extract fit and standard deviation from loess data and return
  out.df <- data.frame(
    row.names = 1:length(unique(loess.df$dist)),
    dist = tapply(loess.df$dist, loess.df$dist, mean),
    no = tapply(loess.df$dist, loess.df$dist, length),
    fit = tapply(loess.df$fit, loess.df$dist, mean),
    sd = tapply(loess.df$residual, loess.df$dist, function(z) {
      sqrt(sum(z^2) / length(z))
    })
  )
  return(out.df)
}

###############################################################################
## Function to create matrix of z-scores
###############################################################################
extract_z_scores <- function(data.matrix, span = 0.05) {
  # Create distance matrix
  split.names <- strsplit(colnames(data.matrix), ('-|:'))
  centre <- sapply(split.names, function(z) {mean(as.numeric(z[2:3]))})
  dist.matrix <- abs(outer(centre, centre, '-'))
  # Extract loess data and create fit and sd matrices
  loess.data <- extract_distance_loess(norm_data, span=span, zero.dist=T,
    zero.values=F)
  fit.matrix <- matrix(
    loess.data$fit[match(dist.matrix, loess.data$dist)],
    ncol = ncol(dist.matrix)
  )
  sd.matrix <- matrix(
    loess.data$sd[match(dist.matrix, loess.data$dist)],
    ncol = ncol(dist.matrix)
  )
  # Process data matrix and calculate z.matrix
  data.matrix[data.matrix == 0] <- NA
  z.matrix = (data.matrix - fit.matrix) / sd.matrix
  return(z.matrix)
}

################################################################################
## Function to create z-score plot
################################################################################
create_zscore_plot <- function(
  z.data, max.z = 3, min.z = -3
) {
  # Check min and max prob values
  if (max.z %% 1 > 0) {
    stop('incorrect max.z value')
  }
  if (min.z %% 1 > 0) {
    stop('incorrect min.z value')
  }
  plot.breaks <- seq(min.z, max.z, 1)
  plot.labels <- as.character(plot.breaks)
  plot.labels[1] <- paste0('<=', plot.labels[1])
  plot.labels[length(plot.labels)] <- paste0('>=', plot.labels[length(plot.labels)])
  # Copy plot data and rename rows and columns
  plot.data <- z.data
  name.list <- strsplit(colnames(plot.data), '[:-]')
  centre <- sapply(name.list, function(z) {mean(as.numeric(z[2:3]))})
  colnames(plot.data) <- centre
  rownames(plot.data) <- centre
  # Manipulate plot data
  plot.data[plot.data < min.z] <- min.z
  plot.data[plot.data > max.z] <- max.z
  plot.data[is.na(plot.data)] <- 0
  plot.data <- melt(plot.data)
  # Create plot
  p <- ggplot(plot.data, aes(x = Var1, y = Var2, fill = value)) +
    geom_raster() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_gradient2(
      name = "prob",
      limits = c(min.z, max.z),
      breaks = plot.breaks,
      labels = plot.labels,
      low = 'blue',
      high = 'red',
      mid = 'white'
    ) +
    coord_fixed(ratio = 1)
  return(p)
}

###############################################################################
## Use bootstrap to calculate reporducibility of tad regions
###############################################################################
tad_overlap_pvalue <- function(
  query.bins, ref.bins, permutations = 100
) {
  ref.overlap <- sum(!is.na(match(query.bins, ref.bins)))
  combined.bins <- c(query.bins, ref.bins)
  overlap.dist <- vector(length=permutations)
  for (i in 1:permutations) {
    sampled.bins <- sample(
      combined.bins, 
      size = length(combined.bins),
      replace = F
    )
    query.sample <- head(sampled.bins, length(query.bins))
    ref.sample <- tail(sampled.bins, length(ref.bins))
    overlap <- sum(!is.na(match(query.sample, ref.sample)))
    overlap.dist[i] <- overlap 
  }
  return(density(overlap.dist))
}

###############################################################################
## Function to plot tads
###############################################################################
plot.tads <- function(tad.list, sample.data) {
  # Sort tad.list bu condition
  sample.data <- sample.data[order(sample.data$condition),]
  tad.list <- tad.list[sample.data$sample]
  plot.data <- data.frame()
  for (sample in names(tad.list)) {
    sample.tad <- tad.list[[sample]]
    sample.index <- match(sample, sample.data$sample)
    sample.df <- data.frame(
      sample = sample,
      condition = sample.data$condition[sample.index],
      centre = sample.tad$centre,
      shift = sample.tad$shift
    )
    plot.data <- rbind(plot.data, sample.df)
  }
  plot <- ggplot(
    plot.data, 
    aes(x=centre, y=sample, colour=condition)
  ) + geom_point()
  return(plot)
}
