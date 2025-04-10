rm(list = ls())

library(MASS)
library(Matrix)     # for bdiag
library(ggplot2)
library(dplyr)
library(glue)

work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'
source('/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/R/convergence.R')
source(glue("{work_directory}R/anchored_lasso.R"))



# ---------------------------
# Settings
# ---------------------------
n_sim <- 10
n_control <- 100
n_treatment <- 100
p <- 50
block_size <- 5
n_blocks <- p / block_size

# Make block covariance (AR(1) inside each block)
make_block <- function(size, rho = 0.6) {
  mat <- outer(1:size, 1:size, function(i, j) rho^abs(i - j))
  return(mat)
}

block_list <- replicate(n_blocks, make_block(block_size), simplify = FALSE)
Sigma_block_diag <- as.matrix(bdiag(block_list))  # Full p x p covariance

# Means
mu_control <- rep(0, p)
# mu_tr <- c(rep(0.4, 5), rep(0, p - 5))  # Only first block has signal
mu_tr <- c(c(10,rep(0,block_size - 1)), c(10,rep(0,block_size - 1)),rep(0, p - 2*block_size))  # Only first block has signal

# Group vector: repeat 1, 2, ..., for each block
group <- rep(1:n_blocks, each = block_size)

# ---------------------------
# Storage for Results
# ---------------------------
sim_results <- data.frame(test_statistic = numeric(n_sim),
                          p_value = numeric(n_sim))
all_active_feature_counts <- list()

set.seed(123)
for (i in 1:n_sim) {
  message(glue("repeat is {i}"))
  
  control <- as.data.frame(mvrnorm(n_control, mu_control, Sigma_block_diag))
  treatment <- as.data.frame(mvrnorm(n_treatment, mu_tr, Sigma_block_diag))
  
  res <- mean_comparison_anchor(
    control = control,
    treatment = treatment,
    pca_method = "dense_pca",
    classifier_method = "group_lasso",
    lambda_type = "lambda.1se",
    n_folds = 10,
    group = group,
    verbose = FALSE
  )
  # res <- mean_comparison_anchor(
  #   control = control,
  #   treatment = treatment,
  #   pca_method = "dense_pca",
  #   classifier_method = "lasso",
  #   lambda_type = "lambda.1se",
  #   n_folds = 10,
  #   verbose = FALSE
  # )
  
  sim_results[i, ] <- c(res$test_statistic, res$p_value)
  all_active_feature_counts[[i]] <- collect_active_features(res)
}

# =========================
# Visualization
# =========================

# Histogram of test statistics
ggplot(sim_results, aes(x = test_statistic)) +
  geom_histogram(color = "black", bins = 15) +
  labs(title = "Histogram of Test Statistics",
       x = "Test Statistic", y = "Count")

# Q-Q plot with y = x reference line
ggplot(sim_results, aes(sample = test_statistic)) +
  stat_qq(distribution = qnorm) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Q-Q Plot of Test Statistics (y = x)",
       x = "Theoretical Quantiles", y = "Sample Quantiles")

# Histogram of p-values
ggplot(sim_results, aes(x = p_value)) +
  geom_histogram(color = "black", bins = 15) +
  labs(title = "Histogram of P-values",
       x = "P-value", y = "Count")

# =========================
# Active Feature Frequency Plot
# =========================

feature_freq <- table(unlist(all_active_feature_counts))
feature_df <- as.data.frame(feature_freq)
colnames(feature_df) <- c("Feature", "Count")
feature_df$Feature <- factor(feature_df$Feature, levels = feature_df$Feature[order(-feature_df$Count)])

# Mark true active features
true_active_names <- paste0("V", c(1,6))
feature_df$is_true_active <- feature_df$Feature %in% true_active_names

ggplot(feature_df, aes(x = Feature, y = Count, fill = is_true_active)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "tomato")) +
  labs(
    title = "Feature Selection Frequency Across Simulations",
    x = "Feature",
    y = "Selection Count",
    fill = "Truly Active"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
