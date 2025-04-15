rm(list = ls())
set.seed(125)
library(MASS)
library(Matrix)     # for bdiag
library(ggplot2)
library(dplyr)
library(glue)
library(tidyr)

work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'
source('/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/R/convergence.R')
source(glue("{work_directory}R/anchored_lasso.R"))

# ---------------------------
# Settings
# ---------------------------
n_sim <- 30
n_control <- 100
n_treatment <- 100
p <- 100
block_size <- 5
n_blocks <- p / block_size

make_block <- function(size, rho = 0.6) {
  outer(1:size, 1:size, function(i, j) rho^abs(i - j))
}

# Generate variance multipliers for each block (e.g., from Uniform[0.5, 2])
block_variances <- runif(n_blocks, min = 0.5, max = 5)

# Create blocks with different variance levels
block_list <- lapply(1:n_blocks, function(b) {
  block_var <- block_variances[b]
  block_var * make_block(block_size)
})

Sigma_block_diag <- as.matrix(bdiag(block_list))
baseline_value <- 100
mu_control <- rep(baseline_value, p)
mu_tr <- c(c(baseline_value+0.7, rep(baseline_value, block_size - 1)), 
           c(baseline_value+0.7, rep(baseline_value, block_size - 1)), 
           rep(baseline_value, p - 2 * block_size))

group <- rep(1:n_blocks, each = block_size)

# ---------------------------
# Storage
# ---------------------------
sim_results <- data.frame(test_statistic = numeric(n_sim), p_value = numeric(n_sim))
all_active_feature_counts <- list()
all_projections <- list()


for (i in 1:n_sim) {
  message(glue("repeat is {i}"))
  
  control <- as.data.frame(mvrnorm(n_control, mu_control, Sigma_block_diag))
  treatment <- as.data.frame(mvrnorm(n_treatment, mu_tr, Sigma_block_diag))
  
  res <- mean_comparison_anchor(
    control = control,
    treatment = treatment,
    pca_method = "dense_pca",
    classifier_method = "lasso",
    lambda_type = "lambda.1se",
    n_folds = 5,
    verbose = FALSE,
    standardize_feature = TRUE
  )
  
  sim_results[i, ] <- c(res$test_statistic, res$p_value)
  all_active_feature_counts[[i]] <- collect_active_features(res)
  
  for (j in 1:length(res$fold_data)) {
    all_projections[[length(all_projections) + 1]] <- data.frame(
      Feature = paste0("V", 1:p),
      proj_direction = res$fold_data[[j]]$proj_direction,
      classifier_coef = res$fold_data[[j]]$classifier_coef,
      leanding_pc = res$fold_data[[j]]$leanding_pc,
      Fold = j,
      Sim = i
    )
  }
}

# =========================
# Combine & Plot Projection Components
# =========================
projection_df <- bind_rows(all_projections)

# Ensure proper ordering of features (V1, V2, ..., Vp)
long_proj_df <- pivot_longer(projection_df,
                             cols = c("proj_direction", "classifier_coef", "leanding_pc"),
                             names_to = "Type", values_to = "Value")
long_proj_df$Feature <- factor(long_proj_df$Feature,
                               levels = paste0("V", 1:p),
                               ordered = TRUE)
# Subset to first 3 simulations only
subset_proj_df <- long_proj_df %>% filter(Sim <= 3)

ggplot(subset_proj_df, aes(x = Feature, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Sim ~ Fold, labeller = label_both) +
  labs(title = "Comparison of Projection Components (Sim 1–3)",
       y = "Component Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# =========================
# QQ Plot for Test Statistics
# =========================
sim_results <- sim_results %>%
  arrange(test_statistic) %>%
  mutate(expected = qnorm(ppoints(n_sim)))

ggplot(sim_results, aes(x = expected, y = test_statistic)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 2, color = "darkred") +
  geom_abline(slope = 1, intercept = 0)+
  labs(title = "QQ Plot of Test Statistics",
       x = "Expected N(0,1) Quantiles",
       y = "Observed Test Statistics") +
  theme_minimal()
