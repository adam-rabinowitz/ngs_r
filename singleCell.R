require(scran)
require(matrixStats)
scranSizeFactors <- function(exprMatrix, minRatio = 0.75) {
  # Convert dataframe
  if (class(exprMatrix) == 'data.frame') {
    exprMatrix <- as.matrix(exprMatrix)
  }
  # Find genes expressed in the minimum ratio of cells
  sampleRatio <- apply(
    exprMatrix, 1, function (z) {
      sum(z > 0) / length(z)
    }
  )
  passed <- which(sampleRatio >= minRatio)
  # Generate sizes
  sizes = seq(20, ncol(exprMatrix) / 2, 20)
  # Calculate size factors
  message(paste('Estimate size factors from', length(passed), 'genes'))
  sumFactors <- computeSumFactors(
    exprMatrix[passed,],
    sizes = sizes
  )
  return(sumFactors)
}
sf <- scranSizeFactors(exp[,passedSamples], minRatio = 0.75)

findVariantGenes <- function(exprMatrix, sizeFactors, minMean) {
  # Normalize expression
  if (ncol(exprMatrix) != length(sizeFactors)) {
    sys.exit('Requires single size factor for each sample')
  }
  exprMatrix <- t(t(exprMatrix) / sizeFactors)
  # Extract mean and variance
  geneMean <- rowMeans(exprMatrix)
  geneVar <- rowVars(exprMatrix)
  cv2 <- geneVar / geneMean^2
  # Perform linear model
  useForFit <- geneMean >= minMean
  fit <- glmgam.fit(
    cbind(
      a0 = 1,
      a1tilde = 1/geneMean[useForFit]
    ), 
    cv2[useForFit]
  )
  # Extract values
  xi <- mean( 1 / sizeFactors)
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"] - xi )
  print(a1)
  # Create plot
  plot(geneMean, cv2, log='xy', col='red', pch=19, cex=0.2)
  xg <- 10^seq( -2, 6, length.out=1000 )
  df <- ncol(exprMatrix) - 1
  lines( xg, (xi+a1)/xg + a0, col='black', lwd=3 )
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .975, df ) / df, 
    col="black", lwd=2, lty="dashed" )
  lines( xg, ( (xi+a1)/xg + a0  ) * qchisq( .025, df ) / df, 
    col="black", lwd=2, lty="dashed" )  
  
  #plot( NULL, xaxt="n", yaxt="n", log="xy", xlim = c( 1e-1, 3e5 ),
  #  ylim = c( .005, 8 ), xlab = "average normalized read count",
  #  ylab = "squared coefficient of variation (CV^2)")
  #axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
  #  expression(10^4), expression(10^5) ) )
  #axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10", "100" ), las=2 )
  #abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  ## Add the data points
  #points(geneMean, cv2, pch=20, cex=.2, col='red' )
}
findVariantGenes(exp[notErcc,passedSamples], sf, 200)

normExpr <- t(t(exp[,passedSamples]) / sf)
logNormExpr <- log2(normExpr[rowMeans(normExpr) > 10,] + 1)
rMean <- apply(logNormExpr, 1, mean)
rVar <- apply(logNormExpr, 1, var)
df <- data.frame(
  'mean' = rMean,
  'var' = rVar,
  'disp' = rVar / rMean
)
ggplot(df, aes(mean, var)) +
  geom_point() +
  geom_smooth()
sf <- scranSizeFactors(exp[,passedSamples], minRatio = 0.75)




