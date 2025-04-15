rm(list = ls())
set.seed(123)
library(glue)
library(data.table)
library(grpreg)
library(HMC)
data_directory <- '/raid6/Tianyu/convergence_risk_gene/try_Cleary_data/'
work_directory <- '/raid6/Tianyu/SingleCell/cleary_data_mean_comparison/'
batch_name <- '01_small_scale'
# data_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
# work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'
source(glue('{data_directory}R/convergence.R'))
source(glue("{work_directory}R/anchored_lasso.R"))

control_subsetting_sample_size <- 500

residual_subset <- readRDS(glue('{data_directory}data/intermediate_data/residual_matrix_small.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

####
clustering_file <- glue("{data_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

dir.create(glue('{work_directory}data/intermediate/{batch_name}/'), recursive = TRUE, showWarnings = FALSE)
dir.create(glue('{work_directory}log/'), recursive = TRUE, showWarnings = FALSE)

#####
all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != 'non-targeting']
table(residual_subset$Guides_collapsed_by_gene)


preprocess_one_setting <- function(residual_subset, treatment_name, clustering){
  
  control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
  treatment <- residual_subset[Guides_collapsed_by_gene == treatment_name, -"Guides_collapsed_by_gene"]

  control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]
  treatment <- treatment[, match(clustering$gene_name, colnames(treatment)), with = FALSE]

  control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
  return(list(control = control,
              treatment = treatment))
}

log_file <- glue('{work_directory}log/logs.txt')

for(treatment_name in all_treatment_names){

  message(glue("Processing {treatment_name}"))
  log_message <- glue("{Sys.time()} - Processing {treatment_name}\n")
  write(log_message, file = log_file, append = TRUE)
  
  process_data <- preprocess_one_setting(residual_subset, treatment_name, clustering)
  
  test_result <- mean_comparison_anchor(control = process_data$control,
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
}

# summary_results_dt <- rbindlist(summary_results)
# fwrite(summary_results_dt, glue('{work_directory}data/intermediate/01_small_scale/all_results_summary.csv'))
# message("Processing completed. Summary saved.")

