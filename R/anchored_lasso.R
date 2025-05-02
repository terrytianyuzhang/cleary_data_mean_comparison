
# ---------------------------
# Fold Functions
# ---------------------------

process_fold_mean_diff <- function(fold_index, control, treatment,
                                   control_split_index, tr_split_index,
                                   pca_method, classifier_method, lambda_type, group, verbose) {
  if (verbose) message(paste0("Processing fold ", fold_index))
  
  # Extract training and testing datasets
  control_train <- control[-control_split_index[[fold_index]], ]
  control_test <- control[control_split_index[[fold_index]], ]
  
  tr_train <- treatment[-tr_split_index[[fold_index]], ]
  tr_test <- treatment[tr_split_index[[fold_index]], ]
  
  # Fit Lasso model for control vs treatment
  classifier_coef <- tryCatch({
    fit_lasso(control_train, tr_train, lambda_type, classifier_method, group)
  }, error = function(e) {
    message("LASSO failed for first training set: ", e$message)
    return(NULL)
  })
  
  # Estimate Principal Component and Projection Direction
  leading_pc <- estimate_leading_pc(control_train, pca_method)
  
  # ===================================================
  # Compute adjusted projection direction
  # ===================================================
  n_effect <- 2 * min(nrow(control_train), nrow(tr_train))  # Effective training size
  a_n <- n_effect^(1/3)  # Scaling factor based on effective sample size
  
  # Adjust projection direction using Lasso coefficients and normalize
  proj_direction <- (leading_pc + a_n * classifier_coef)  # Adjustment
  proj_direction <- proj_direction / norm(proj_direction, type = "2")  # Normalize
  
  # ===================================================
  # Compute test statistic, scores, and variance
  # ===================================================
  # Extract relevant gene data from control and treatment groups
  control_matrix <- as.matrix(control_test)
  tr_matrix <- as.matrix(tr_test)
  
  # Compute scores for control and treatment groups
  control_score <- control_matrix %*% proj_direction
  tr_score <- tr_matrix %*% proj_direction
  
  # Compute the test statistic as the difference of means
  T_stat <- mean(control_score) - mean(tr_score)
  
  # Get sample sizes
  n_x <- nrow(control_test)
  n_z <- nrow(tr_test)
  
  # Compute the variance of the test statistic
  T_variance <- (var(control_score) * (n_x - 1) / (n_x^2)) + (var(tr_score) * (n_z - 1) / (n_z^2))
  
  # Take notes for each split
  return(list(
    statistic = T_stat,  # T_stat is the mean instead of the one with a standard normal distribution
    variance = T_variance,
    control_score = control_score,
    tr_score = tr_score,
    proj_direction = proj_direction,
    classifier_coef = classifier_coef,
    leading_pc = leading_pc,
    group = group
  ))
}

combine_folds_mean_diff <- function(fold_data, verbose = FALSE) {
  
  n_folds <- length(fold_data)
  numerator_test_statistic<-0
  denominator_variance<-0
  
  # Identify successful folds (non-degenerate cases)
  valid_folds <- 1:n_folds
  first_valid_fold <- valid_folds[1]  # Select the first non-degenerate case as the baseline
  
  # ===================================================
  # compute test statistics
  # ===================================================
  
  for (i in valid_folds) {
    denominator_variance <- denominator_variance + fold_data[[i]]$variance
    projection_sign_match <- sign(crossprod(fold_data[[first_valid_fold]]$proj_direction, fold_data[[i]]$proj_direction))
    
    if (projection_sign_match == 0) {
      if (verbose) message("The projection directions are orthogonal")
      projection_sign_match <- 1
    }
    
    numerator_test_statistic <- numerator_test_statistic + projection_sign_match * fold_data[[i]]$statistic
  }
  
  numerator_test_statistic <- numerator_test_statistic / n_folds
  standard_error <- sqrt(denominator_variance / (n_folds^2))
  test_statistic <- numerator_test_statistic / standard_error
  p_value <- 2 * pnorm(-abs(test_statistic))
  
  return(list(
    p_value = p_value,
    test_statistic = test_statistic,
    fold_data = fold_data
  ))
}

# ---------------------------
# Main Function
# ---------------------------

mean_comparison_anchor <- function(
    control, treatment,
    pca_method = c("dense_pca", "sparse_pca"),
    classifier_method = c("lasso", "group_lasso"),
    lambda_type = 'lambda.1se',
    n_folds = 10,
    group = NULL,
    standardize_feature = TRUE,
    verbose = TRUE
) {
  # ============================================
  # Data Preprocessing: Validation and Conversion
  # ============================================
  control <- validate_and_convert_data(control, "control")
  treatment <- validate_and_convert_data(treatment, "treatment")
  
  check_non_null_and_identical_colnames(list(control, treatment))
  
  if(standardize_feature){
    # Normalize and split
    normalized_list <- normalize_and_split(control, treatment)
    
    # Access results
    control <- normalized_list$df1
    treatment <- normalized_list$df2
  }
  
  pca_method <- match.arg(pca_method) #match.arg takes the first argument of pca_method as the default value if not specified
  classifier_method <- match.arg(classifier_method)
  
  if (!is.null(group) && classifier_method == 'lasso'){
    message("the grouping vector is not NULL but the method is normal LASSO, set classifier_method as group_lasso in mean_comparison_anchor()")
    
  }
  
  if (!is.null(group) && (!is.vector(group) || length(group) != ncol(control))) {
    stop("Error: `group` must be NULL or a vector of the same length as the number of columns in `control`.")
  }
  
  # ============================================
  # Split Datasets into Folds
  # ============================================
  split_indices <- lapply(list(control, treatment), function(data) {
    check_data_for_folds(data, n_folds)
    index_spliter(1:nrow(data), n_folds)
  })
  
  control_split_index <- split_indices[[1]]
  tr_split_index <- split_indices[[2]]
  
  fold_data <- vector("list", n_folds)
  
  # ============================================
  # Process data for each fold
  # ============================================  
  for(i in 1:n_folds){
    fold_data[[i]] <- process_fold_mean_diff(i, control, treatment,
                                             control_split_index, tr_split_index, 
                                             pca_method, classifier_method, lambda_type, group, verbose)
  }
  
  # ===================================================
  # Now combine the folds
  # ===================================================
  return(combine_folds_mean_diff(fold_data, verbose))
  
}

# ---------------------------
# Post-hoc Analysis Function
# ---------------------------

compute_predictive_contributions <- function(result, group, group_threshold = 5) {
  # ============================================
  # Step 1: Stack classifier coefficients
  # ============================================
  n_folds <- length(result$fold_data)
  beta_matrix <- do.call(rbind, lapply(1:n_folds, function(i) result$fold_data[[i]]$classifier_coef))
  colnames(beta_matrix) <- names(result$fold_data[[1]]$classifier_coef)
  
  # ============================================
  # Step 2: Average beta coefficients across folds
  # ============================================
  avg_beta <- colMeans(beta_matrix)
  
  # ============================================
  # Step 3: Identify active groups based on threshold
  # ============================================
  active <- collect_active_features_proj(result, group = group, group_threshold = group_threshold)
  active_groups <- active$active_groups
  
  # ============================================
  # Step 4: Compute group-wise contributions
  # ============================================
  contribution_scores <- numeric(length(active_groups) + 1)
  names(contribution_scores) <- c(active_groups, "other")
  
  for (group_index in seq_along(active_groups)) {
    sub_classifier <- avg_beta[group == active_groups[group_index]]
    contribution_scores[group_index] <- sum(sub_classifier^2) / sum(avg_beta^2) 
  }
  
  contribution_scores["other"] <- sum(avg_beta^2) - sum(contribution_scores)
  
  # ============================================
  # Step 5: Convert to data frame
  # ============================================
  contribution_df <- data.frame(
    group = names(contribution_scores),
    score = as.numeric(contribution_scores)
  )
  
  return(contribution_df)
}


collect_active_features_proj <- function(test_result, voting_method = c("majority_voting"), 
                                    group = NULL, group_threshold = 1) {
  fold_data <- test_result$fold_data
  voting_method <- match.arg(voting_method)
  n_folds <- length(fold_data)
  active_features_list <- vector("list", n_folds)
  
  # Collect non-zero features for each fold
  for (i in 1:n_folds) {
    beta <- result$fold_data[[i]]$proj_direction
    names(beta) <- names(result$fold_data[[i]]$classifier_coef)
    
    non_zero_features <- names(beta[abs(beta) > 1e-10])
    active_features_list[[i]] <- non_zero_features
  }
  
  # Flatten and count
  all_active_features <- unlist(active_features_list)
  feature_counts <- table(all_active_features)
  
  # Apply majority voting
  if (voting_method == 'majority_voting') {
    active_features <- names(feature_counts[feature_counts > n_folds / 2])
  }
  
  # Group handling
  if (!is.null(group)) {
    if (is.null(names(group))) {
      names(group) <- names(fold_data[[1]]$classifier_coef)
    }
    
    group_nonzero_counts <- table(group[active_features])
    active_groups <- names(group_nonzero_counts[group_nonzero_counts >= group_threshold])
    
    return(list(
      active_features = active_features,
      active_groups = active_groups
    ))
  }
  
  return(active_features)
}
