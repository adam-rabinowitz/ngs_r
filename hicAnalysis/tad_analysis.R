# Set parameters
rm(list = ls())
input_file <- '/farm/scratch/rs-bio-lif/rabino01/Yasu/Yasu_HiC/intraArmData/NGS-8178.2000.I_ArmR.500.normMatrix.gz'
output_plot = '/farm/home/rabino01/Yasu/prob_plot.tiff'
source('/farm/home/rabino01/Yasu/tad_functions.R')
# Read in data
norm_data <- as.matrix(read.table(input_file, header = T, check.names = F))
# Extract squares
square_list <- extract_squares(norm_data, size=50, offset=5, rm_zero=T)
# Calculate square metrics
square_metrics <- extract_square_metrics(square_list, span=0.75, half_window=10)
# Find tads
tad_data <- filter_square_metrics(square_metrics, half_window=10)
# Plot tad_data
prob_plot <- create_prob_plot(norm_data, tad_data, max_prob=1e-2, min_prob=1e-4)
ggsave('/farm/home/rabino01/Yasu/test.tiff', prob_plot)
plot(square_metrics$centre, square_metrics$median, type = 'l')
points(square_metrics$centre, square_metrics$loess, col = 'red', type = 'l')
abline(v = tad_data$centre, col = 'green')
abline(v=4719566)

test_plot <- subset(
  square_metrics,
  square_metrics$centre >= 4720000 &
  square_metrics$centre <= 4750000  
)
plot(test_plot$centre, test_plot$media, type = 'l')
abline(v = test_plot$centre[test_plot$min_diff])

log2_direction <- calc_log2_direction(
  norm_data, log2_pair_calc, max_window_no=100, min_window_no=2 )
tad_window <- find_tad_window(log2_direction, half_window = 5)
filtered_tad <- filter_tad(tad_window, min_shift = 0.2)
create_prob_plot(norm_data, filtered_tad, max_prob=1e-2, min_prob=1e-4)

z.scores <- extract_z_scores(norm_data, span = 0.05)
create_zscore_plot(z.scores)


