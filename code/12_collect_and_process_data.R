rm(list = ls())

library(data.table)
library(glue)
library(ggplot2)
library(tidyr)
library(dplyr)

data_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'

source(glue('{data_directory}R/convergence.R'))
source(glue('{work_directory}R/anchored_lasso.R'))

batch_name <- '11_all_groups_sparse'
result_folder <- glue('{work_directory}data/intermediate/{batch_name}/')

result_files <- list.files(result_folder, pattern = "\\.rds$", full.names = TRUE)
results_list <- lapply(result_files, readRDS)
names(results_list) <- gsub("\\.rds$", "", basename(result_files))

clustering_file <- glue("{data_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

# ============================================
# Step 1: Summarize test results
# ============================================
summary_list <- vector("list", length(results_list))
group_threshold <- 1
for (treatment_name in names(results_list)) {
  result <- results_list[[treatment_name]]
  active <- collect_active_features_proj(
    result,
    group = clustering$cluster_index,
    group_threshold = group_threshold
  )
  
  feature_string <- if (is.list(active) && !is.null(active$active_features)) {
    paste(active$active_features, collapse = ", ")
  } else {
    paste(active, collapse = ", ")
  }
  
  group_string <- if (is.list(active) && !is.null(active$active_groups)) {
    paste(active$active_groups, collapse = ", ")
  } else {
    NA_character_
  }
  
  summary_list[[treatment_name]] <- data.table(
    treatment = treatment_name,
    p_value = result$p_value,
    test_statistic = result$test_statistic,
    active_groups = group_string
  )
}

summary_dt <- rbindlist(summary_list, use.names = TRUE)
print(summary_dt)

# Clean p-values and test statistics
summary_dt$p_value <- summary_dt$p_value.V1
summary_dt <- summary_dt[p_value < 0.05, ]
summary_dt$test_statistic <- summary_dt$test_statistic.V1

# ============================================
# Step 2: Expand active group list
# ============================================
summary_long <- summary_dt[, .(
  group = unlist(strsplit(active_groups, ",\\s*")),
  p_value = p_value
), by = treatment]

summary_long[, group := as.integer(group)]
summary_long <- summary_long[!is.na(group), ]

# ============================================
# Step 3: Compute group-level contribution scores
# ============================================
score_list <- vector("list", length(results_list))

for (treatment_name in names(results_list)) {
  result <- results_list[[treatment_name]]
  score_df <- compute_predictive_contributions(
    result,
    clustering$cluster_index,
    group_threshold = group_threshold
  )
  score_df$treatment <- treatment_name
  score_list[[treatment_name]] <- score_df
}

score_long <- rbindlist(score_list, use.names = TRUE)

# Ensure consistent types
score_long[, group := as.character(group)]
summary_long[, group := as.character(group)]

# ============================================
# Step 4: Merge contributions with p-values
# ============================================
summary_with_score <- merge(summary_long, score_long, by = c("treatment", "group"), all.x = TRUE)
print(summary_with_score)

# ============================================
# Step 5: Heatmap visualization
# ============================================
summary_with_score[, group := factor(group)]

ggplot(summary_with_score, aes(x = treatment, y = group, fill = score)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white", 
                      name = "Predictive\nScore") +
  labs(x = "Treatment", y = "Group", title = "Predictive Score Heatmap") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# ============================================
# Step 6: Matrix + base heatmap
# ============================================
score_matrix <- reshape2::acast(
  summary_with_score,
  group ~ treatment,
  value.var = "score"
)

score_matrix[is.na(score_matrix)] <- 0
score_matrix[is.infinite(score_matrix)] <- 0

heatmap(
  score_matrix,
  col = colorRampPalette(c("white", "red"))(100),
  scale = "none",
  na.rm = FALSE,
  margins = c(8, 8),
  main = "Clustering Heatmap of Predictive Scores"
)

# Step 1: Get clustered row and column orders from base::heatmap
hm <- heatmap(
  score_matrix,
  col = colorRampPalette(c("white", "red"))(100),
  scale = "none",
  na.rm = FALSE,
  margins = c(8, 8),
  main = "Clustering Heatmap of Predictive Scores"
)

row_order <- rownames(score_matrix)[hm$rowInd]
col_order <- colnames(score_matrix)[hm$colInd]

# Step 2: Convert matrix to long format
score_long_for_plot <- as.data.table(as.table(score_matrix))
colnames(score_long_for_plot) <- c("group", "treatment", "score")

# Step 3: Apply clustering order as factor levels
score_long_for_plot[, group := factor(group, levels = row_order)]
score_long_for_plot[, treatment := factor(treatment, levels = col_order)]

score_long_for_plot <- score_long_for_plot[group %in% row_order[1:16], ]
score_long_for_plot <- score_long_for_plot[treatment %in% col_order[1:25], ]
# Step 4: Plot using ggplot2
ggplot(score_long_for_plot, aes(x = treatment, y = group, fill = score)) +
  geom_tile(color = "#555555") +
  scale_fill_gradient2(high = "#a65552", name = "Load. Conc.") +
  labs(x = "Perturbation", y = "Group", title = "Loading Concentration for Significant Perturbations") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),   # 👈 Resize y-axis labels
    axis.title = element_text(size = 16),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

dir.create(file.path(work_directory, "report/", batch_name), showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(work_directory, "report/", batch_name, "/loading_contribution.pdf"),
       width = 8, height = 5)


