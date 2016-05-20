require(ggplot2)
pdfA4SquareTheme <- theme(
  aspect.ratio = 1,
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14, face = 'bold'),
  title = element_text(size = 15, face = 'bold')
)