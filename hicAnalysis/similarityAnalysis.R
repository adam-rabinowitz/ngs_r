###############################################################################
## Function to calculate overlap score with a supplied allowed difference
###############################################################################
similarity.score <- function(Var1, Var2, diff=0) {
  grid <- expand.grid(Var1,Var2)
  grid$diff <- grid$Var1 - grid$Var2
  grid <- grid[abs(grid$diff)<=diff,]
  Var1Pass <- unique(grid$Var1)
  Var2Pass <- unique(grid$Var2)
  totalPass <- length(Var1Pass) + length(Var2Pass)
  ratio <- totalPass / (length(Var1) + length(Var2))
  return(ratio)
}

###############################################################################
## Function to calculate similarity scores across samples
###############################################################################
calculate.similarity <- function(value.list, diff=0) {
  # Create output matrix and set diagonal to 1
  output <- matrix(
    ncol=length(value.list),
    nrow=length(value.list),
    dimnames=list(
      names(value.list),
      names(value.list)  
    )
  )
  diag(output) <- 1
  # Calculate similarity scores for all other elements of the matrix
  combinations <- t(combn(names(value.list), 2))
  for (i in 1:nrow(combinations)) {
    # Extract samples and their corresponding values
    sample1 <- combinations[i,1]
    sample2 <- combinations[i,2]
    values1 <- value.list[[sample1]]
    values2 <- value.list[[sample2]]
    # Calculate similarity score and add to matrix
    similarity <- similarity.score(values1, values2, diff=diff)
    output[sample1, sample2] <- similarity
    output[sample2, sample1] <- similarity
  }
  # Return matrix
  return(output)
}

###############################################################################
## Function to calculate mean similarity scores across conditions
###############################################################################
condition.similarity <- function(similarity.matrix, sample.data) {
  # Check row and column names and remove diagonal
  if (!all.equal(colnames(similarity.matrix), row.names(similarity.matrix))) {
    stop('Row and column names must match')
  }
  diag(similarity.matrix) <- NA
  # Extract conditions and generate output matrix
  conditions <- split(sample.data$sample, sample.data$condition)
  output <- matrix(
    ncol=length(conditions),
    nrow=length(conditions),
    dimnames=list(
      names(conditions),
      names(conditions)  
    )
  )
  # Fill output matrix and return
  combinations <- expand.grid(names(conditions), names(conditions))
  for (i in 1:nrow(combinations)) {
    condition1 <- combinations[i,1]
    condition2 <- combinations[i,2]
    samples1 <- as.character(conditions[[condition1]])
    samples2 <- as.character(conditions[[condition2]])
    subMatrix <- similarity.matrix[samples1, samples2]
    subMatrixMean <- mean(c(subMatrix), na.rm=T)
    output[condition1,condition2] <- subMatrixMean
  }
  return(output)
}


