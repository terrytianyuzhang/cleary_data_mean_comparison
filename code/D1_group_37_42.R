rm(list = ls())
set.seed(123)

library(glue)
library(data.table)
library(grpreg)
library(HMC)
library(foreach)
library(doParallel)

# Setup directories
data_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'
source(glue('{data_directory}R/convergence.R'))
source(glue("{work_directory}R/anchored_lasso.R"))

# Load and preprocess data
control_subsetting_sample_size <- 500

residual_subset <- readRDS(glue('{data_directory}data/intermediate_data/residual_matrix_all_in_paper.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

clustering_file <- glue("{data_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

control_residual <- residual_subset[Guides_collapsed_by_gene == 'non-targeting',]
# Clusters to analyze
target_clusters <- c(37, 42, 31,15,1,10)

# Initialize result list
cluster_stats_list <- list()

# Loop through target clusters
for (cluster_id in target_clusters) {
  selected_genes <- clustering$gene_name[clustering$cluster_index == cluster_id]
  
  # Subset control residual for selected genes
  group_data <- control_residual[, ..selected_genes]
  
  # Compute summary statistics
  mean_vec <- colMeans(group_data)
  sd_vec <- apply(group_data, 2, sd)
  
  cluster_stats_list[[as.character(cluster_id)]] <- data.table(
    gene = names(mean_vec),
    mean = mean_vec,
    sd = sd_vec,
    cluster = factor(cluster_id)
  )
}

# Combine all clusters
cluster_summary_dt <- rbindlist(cluster_stats_list)

# Plot mean vs. sd
library(ggplot2)
ggplot(cluster_summary_dt, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.7, color = "#a65552") +
  facet_wrap(~ cluster) +        # axes are fixed by default
  coord_fixed() +                # optional: enforces equal unit length for x & y
  labs(
    title = "Mean vs. SD of Genes in Selected Clusters",
    x = "Mean",
    y = "Standard Deviation"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14)
  )


