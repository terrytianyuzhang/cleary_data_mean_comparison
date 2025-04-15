rm(list = ls())

library(data.table)
library(glue)
library(ggplot2)
library(tidyr)
library(dplyr)

data_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
# work_directory <- '/raid6/Tianyu/SingleCell/cleary_data_mean_comparison/'
source(glue('{data_directory}R/convergence.R'))
work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'

batch_name <- '01_small_scale'
# Define target folder
result_folder <- glue('{work_directory}data/intermediate/{batch_name}/')

# Get all .rds files in the folder
result_files <- list.files(result_folder, pattern = "\\.rds$", full.names = TRUE)
results_list <- lapply(result_files, readRDS)
names(results_list) <- gsub("\\.rds$", "", basename(result_files))

clustering_file <- glue("{data_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)


# Initialize storage
summary_list <- vector("list", length(results_list))
# Loop over treatments
for (treatment_name in names(results_list)) {
  result <- results_list[[treatment_name]]
  active <- collect_active_features(result, 
                                    group = clustering$cluster_index,
                                    group_threshold = 5)
  
  # Convert active features to comma-separated string
  feature_string <- if (is.list(active) && !is.null(active$active_features)) {
    paste(active$active_features, collapse = ", ")
  } else {
    paste(active, collapse = ", ")
  }
  
  # Convert active groups to comma-separated string
  group_string <- if (is.list(active) && !is.null(active$active_groups)) {
    paste(active$active_groups, collapse = ", ")
  } else {
    NA_character_
  }
  
  # Store result
  summary_list[[treatment_name]] <- data.table(
    treatment = treatment_name,
    p_value = result$p_value,
    test_statistic = result$test_statistic,
    active_features = feature_string,
    active_groups = group_string
  )
}

# Combine into one data.table
summary_dt <- rbindlist(summary_list, use.names = TRUE)

# View
print(summary_dt)



# Assume summary_dt is your current result table
summary_dt$p_value <- summary_dt$p_value.V1
summary_dt$test_statistic <- summary_dt$test_statistic.V1

# Step 1: expand the active group list
summary_long <- summary_dt[, .(
  group = unlist(strsplit(active_groups, ",\\s*")),
  p_value = p_value,
  treatment = treatment
), by = treatment]

# Step 2: convert group to integer (optional)
summary_long[, group := as.integer(group)]

# Step 3: for visualization, use -log10(p) as color intensity
summary_long[, log10_p := -log10(p_value)]

# Step 4: spread to wide format
heatmap_data <- summary_long %>%
  pivot_wider(names_from = treatment, values_from = log10_p, values_fill = 0)

# Step 5: convert back to long format for ggplot
heatmap_long <- heatmap_data %>%
  pivot_longer(cols = -group, names_to = "treatment", values_to = "log10_p")

# Step 6: plot
ggplot(heatmap_long, aes(x = treatment, y = factor(group), fill = log10_p)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = expression(-log[10](p))) +
  labs(x = "Treatment", y = "Active Group", title = "Heatmap of Active Groups by Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

