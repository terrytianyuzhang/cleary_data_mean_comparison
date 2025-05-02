rm(list = ls())
# Load required packages
library(MASS)       # For mvrnorm
library(PMA)        # For sparse PCA
library(irlba)      # For dense PCA
library(Matrix)     # For norm
library(glue)

work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'
source('/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/R/convergence.R')
source(glue("{work_directory}R/anchored_lasso.R"))

# ------- Simulate low-rank signal in mean difference -------
simulate_pca_shift_data <- function(n_control = 100, n_treatment = 100, 
                                    p = 50, signal_strength = 1.5, 
                                    shift_rank = 1, rho = 0.6) {
  # Create block-correlated covariance matrix
  Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
  
  # Generate a low-rank mean shift (in the first few directions)
  U <- matrix(rnorm(p * shift_rank), p, shift_rank)
  U <- apply(U, 2, function(u) u / sqrt(sum(u^2)))  # orthonormalize
  mu_diff <- signal_strength * U[,1]
  
  # Generate samples
  control <- MASS::mvrnorm(n_control, mu = rep(0, p), Sigma = Sigma)
  treatment <- MASS::mvrnorm(n_treatment, mu = mu_diff, Sigma = Sigma)
  
  # Assign column names
  col_names <- paste0("gene_", seq_len(p))
  colnames(control) <- col_names
  colnames(treatment) <- col_names
  
  list(control = control, treatment = treatment, mu_diff = mu_diff)
}


# Simulate
sim_data <- simulate_pca_shift_data()

# Dense PCA test
result_dense <- mean_comparison_anchor(
  control = sim_data$control,
  treatment = sim_data$treatment,
  pca_method = "dense_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  n_folds = 5,
  verbose = TRUE
)

# Sparse PCA test
result_sparse <- mean_comparison_anchor(
  control = sim_data$control,
  treatment = sim_data$treatment,
  pca_method = "sparse_pca",
  classifier_method = "lasso",
  lambda_type = "lambda.min",
  n_folds = 5,
  verbose = TRUE
)

cat("\n--- Dense PCA ---\n")
print(result_dense$p_value)

cat("\n--- Sparse PCA ---\n")
print(result_sparse$p_value)

plot(result_dense$fold_data[[1]]$leading_pc)
plot(result_sparse$fold_data[[1]]$leading_pc)
