################################################################################
# Script Information
# Author: Eisnecker, Philipp
# Date: 2025/02/18
# Purpose: Visualization of the dependent variable 'Vol' (wood volume) distribution
#          comparing original and multiple transformed versions side-by-side.
#          Bars colored by relative distance from the median (|x - median| / median)
#          using a custom Wes Anderson palette.
################################################################################

# --- Load libraries -----------------------------------------------------------
library(ggplot2)
library(dplyr)
library(e1071)
library(wesanderson)
library(patchwork)
library(MASS)

# --- Define file path and read data -------------------------------------------
method <- "very_strict"

file_path <- paste0("S:/Projekte/VorratAktuell/Arbeitspakete/AP3_Eingangsdaten/Parameterselektion/all/",
                    method, "/stp_all_", method, ".csv")
df <- read.csv(file_path)

vol_df <- data.frame(Vol = df$Vol)

# --- Determine optimal Box-Cox lambda -----------------------------------------
boxcox_model <- boxcox(lm(Vol ~ 1, data = vol_df), lambda = seq(-2, 2, 0.1), plotit = FALSE)
lambda_opt <- boxcox_model$x[which.max(boxcox_model$y)]
cat("Optimal lambda for Box-Cox:", lambda_opt, "\n")

# --- Create transformations ---------------------------------------------------
df$Vol_sqrt <- sqrt(df$Vol)
df$Vol_log <- ifelse(df$Vol > 0, log(df$Vol), NA)
df$Vol_boxcox <- if(lambda_opt == 0) log(df$Vol) else (df$Vol^lambda_opt - 1) / lambda_opt
df$Vol_cbrt <- sign(df$Vol) * abs(df$Vol)^(1/3)  # Cube root
df$Vol_scaled <- scale(df$Vol)

# --- Function to prepare histogram data ----------------------------------------
prepare_hist_data <- function(data) {
  data <- data[!is.na(data)]
  median_val <- median(data)
  hist_data <- hist(data, breaks = 50, plot = FALSE)
  hist_df <- data.frame(mid = hist_data$mids, counts = hist_data$counts)
  hist_df$rel_dist_median <- abs(hist_df$mid - median_val) / max(abs(hist_df$mid - median_val))
  return(list(hist_df = hist_df, median_val = median_val))
}

# --- Prepare list of datasets to plot ------------------------------------------
data_list <- list(
  Original = df$Vol,
  SquareRoot = df$Vol_sqrt,
  Logarithm = df$Vol_log,
  BoxCox = df$Vol_boxcox,
  CubeRoot = df$Vol_cbrt,
  Scaled = as.numeric(df$Vol_scaled)
)

# --- Create histogram data for all ---------------------------------------------
hist_list <- lapply(data_list, prepare_hist_data)

# --- Define common max relative distance for color scaling ---------------------
max_rel_dist <- max(sapply(hist_list, function(x) max(x$hist_df$rel_dist_median, na.rm=TRUE)))

# --- Define color palette ------------------------------------------------------
a_palette <- c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF",
               "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF")

# --- Function to create histogram plot -----------------------------------------
create_histogram_plot <- function(hist_data, median_val, skew_val, n_points, title, max_rel_dist, palette) {
  ggplot(hist_data, aes(x = mid, y = counts, fill = rel_dist_median)) +
    geom_col(color = "black") +
    scale_fill_gradientn(colors = palette,
                         limits = c(0, max_rel_dist),
                         values = scales::rescale(c(0, 0.05, 0.1, 0.3, 0.6, 1))) +
    geom_vline(xintercept = median_val, color = "darkred", linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = paste0("N=", n_points, " | Skewness: ", round(skew_val, 2)),
      x = "Volume (transformed)",
      y = "Frequency",
      fill = "Relative Distance\nto Median",
      caption = "Dark red dashed line = Median"
    )
}

# --- Generate plots ------------------------------------------------------------
n_points <- nrow(df)
plots <- list()
for(name in names(data_list)) {
  data_vec <- data_list[[name]]
  data_vec <- data_vec[!is.na(data_vec)]
  skew_val <- skewness(data_vec)
  hist_data <- hist_list[[name]]$hist_df
  median_val <- hist_list[[name]]$median_val
  title <- paste0(name, ifelse(name == "BoxCox", paste0(" (lambda=", round(lambda_opt, 2), ")"), ""))
  plots[[name]] <- create_histogram_plot(hist_data, median_val, skew_val, length(data_vec), title, max_rel_dist, a_palette)
}

# --- Arrange and display plots in grid ----------------------------------------
final_plot <- (plots$Original + plots$SquareRoot + plots$Logarithm) /
  (plots$BoxCox + plots$CubeRoot + plots$Scaled) + plot_layout(guides = 'collect')
print(final_plot)

# --- Additional plot: Original vs best transformation (Box-Cox) ----------------
best_plots <- (plots$Original + plots$BoxCox) + plot_layout(guides = 'collect')
print(best_plots)
