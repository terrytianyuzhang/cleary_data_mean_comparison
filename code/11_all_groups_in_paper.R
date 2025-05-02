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
batch_name <- '11_all_groups_sparse'

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
n_cores <- 25
cl <- makeCluster(n_cores)
registerDoParallel(cl)

log_file <- glue('{work_directory}log/logs.txt')
log_dir <- glue('{work_directory}log/')
log_files <- list.files(log_dir, full.names = TRUE)
file.remove(log_files)

foreach(treatment_name = all_treatment_names,
        .packages = c("data.table", "glue", "grpreg", "HMC")) %dopar% {
          tryCatch({
            # Define log file path
            treatment_log_file <- paste0(work_directory, 'log/', treatment_name, '.log')
            
            # Open file connections
            log_con_out <- file(treatment_log_file, open = "wt")
            log_con_msg <- file(treatment_log_file, open = "at")
            
            # Ensure cleanup on exit
            on.exit({
              sink(NULL, type = "message")
              sink(NULL)
              close(log_con_out)
              close(log_con_msg)
            }, add = TRUE)
            
            # Redirect output
            sink(log_con_out, split = TRUE)
            sink(log_con_msg, type = "message", append = TRUE)
            
            cat(glue("Processing {treatment_name}\n"))
            cat(glue("{Sys.time()} - Starting...\n"))
            
            process_data <- preprocess_one_setting(residual_subset, treatment_name, clustering)
            
            test_result <- mean_comparison_anchor(
              control = process_data$control,
              treatment = process_data$treatment,
              pca_method = "sparse_pca",
              classifier_method = "group_lasso",
              lambda_type = "lambda.min",
              n_folds = 5,
              group = clustering$cluster_index,
              verbose = TRUE
            )
            
            output_filename <- paste0(work_directory, 'data/intermediate/', batch_name,'/', treatment_name, '.rds')
            saveRDS(test_result, file = output_filename)
            
            cat(glue("{Sys.time()} - Finished successfully.\n"))
            
          }, error = function(e) {
            # Capture error in a central log (optional)
            central_log <- glue('{work_directory}log/errors.log')
            write(glue("[{Sys.time()}] Error for {treatment_name}: {e$message}"), 
                  file = central_log, append = TRUE)
          })
        }

# Stop cluster
stopCluster(cl)
