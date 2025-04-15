rm(list = ls())
set.seed(123)

library(glue)
library(data.table)
library(grpreg)
library(HMC)
library(foreach)
library(doParallel)

# Setup directories
data_directory <- '/raid6/Tianyu/convergence_risk_gene/try_Cleary_data/'
work_directory <- '/raid6/Tianyu/SingleCell/cleary_data_mean_comparison/'
batch_name <- '11_all_groups'

source(glue('{data_directory}R/convergence.R'))
source(glue("{work_directory}R/anchored_lasso.R"))

# Create necessary folders
dir.create(glue('{work_directory}data/intermediate/{batch_name}/'), recursive = TRUE, showWarnings = FALSE)
dir.create(glue('{work_directory}log/'), recursive = TRUE, showWarnings = FALSE)

# Load and preprocess data
control_subsetting_sample_size <- 500

residual_subset <- readRDS(glue('{data_directory}data/intermediate_data/residual_matrix_all_in_paper.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

clustering_file <- glue("{data_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

# Get treatment names
all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != 'non-targeting']

# Preprocessing function
preprocess_one_setting <- function(residual_subset, treatment_name, clustering){
  control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
  treatment <- residual_subset[Guides_collapsed_by_gene == treatment_name, -"Guides_collapsed_by_gene"]
  
  matched_genes <- clustering$gene_name[clustering$gene_name %in% colnames(control)]
  control <- control[, ..matched_genes]
  treatment <- treatment[, ..matched_genes]
  
  if (nrow(control) < control_subsetting_sample_size) {
    warning(glue("Only {nrow(control)} non-targeting samples available; using all of them."))
  } else {
    control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
  }
  
  return(list(control = control, treatment = treatment))
}

# Register parallel backend
n_cores <- 10
cl <- makeCluster(n_cores)
registerDoParallel(cl)

log_file <- glue('{work_directory}log/logs.txt')

# Parallel processing
foreach(treatment_name = all_treatment_names, .packages = c("data.table", "glue", "grpreg", "HMC")) %dopar% {
  tryCatch({
    message(glue("Processing {treatment_name}"))
    log_message <- glue("{Sys.time()} - Processing {treatment_name}\n")
    write(log_message, file = log_file, append = TRUE)
    
    process_data <- preprocess_one_setting(residual_subset, treatment_name, clustering)
    
    test_result <- mean_comparison_anchor(
      control = process_data$control,
      treatment = process_data$treatment,
      pca_method = "sparse_pca",
      classifier_method = "group_lasso",
      lambda_type = "lambda.min",
      n_folds = 5,
      group = clustering$cluster_index,
      verbose = FALSE
    )
    
    output_filename <- glue('{work_directory}data/intermediate/{batch_name}/{treatment_name}.rds')
    saveRDS(test_result, file = output_filename)
  }, error = function(e) {
    warning(glue("Failed to process {treatment_name}: {e$message}"))
  })
}

# Stop cluster
stopCluster(cl)
