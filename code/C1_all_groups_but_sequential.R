rm(list = ls())
set.seed(123)

library(glue)
library(data.table)
library(grpreg)
library(HMC)
library(foreach)
library(doParallel)

# Setup directories
work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'
data_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
# data_directory <- '/raid6/Tianyu/convergence_risk_gene/try_Cleary_data/'
# work_directory <- '/raid6/Tianyu/SingleCell/cleary_data_mean_comparison/'

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

for (treatment_name in all_treatment_names) {
    cat(glue("Processing {treatment_name}\n"))
    
    process_data <- preprocess_one_setting(residual_subset, treatment_name, clustering)
    
    test_result <- mean_comparison_anchor(
      control = process_data$control[,1:50],
      treatment = process_data$treatment[, 1:50],
      pca_method = "sparse_pca",
      classifier_method = "group_lasso",
      lambda_type = "lambda.min",
      n_folds = 5,
      group = clustering$cluster_index[1:50],
      verbose = TRUE
    )
    
    output_filename <- glue('{work_directory}data/intermediate/{batch_name}/{treatment_name}.rds')
    # saveRDS(test_result, file = output_filename)
    
    cat(glue("{Sys.time()} - Finished successfully.\n"))
    
  
}

