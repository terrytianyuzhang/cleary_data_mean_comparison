library(ggplot2)
library(ggpubr)

rm(list = ls())

work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'

plot_dir <- paste0(work_directory, '/report/31_GO_plot/')
GO_plt <- readRDS(paste0(plot_dir, 'GO_plot.rds'))

batch_name <- '11_all_groups_sparse'
heatmap_plt <- readRDS(paste0(work_directory, "report/", batch_name, "/loading_contribution.rds"))

# Modify GO_plt: add bigger side margins
GO_plt <- GO_plt + theme(
  plot.margin = margin(t = 5, r = 80, b = 5, l = 80)  # Top, Right, Bottom, Left
)

# Now arrange
combined_plot <- ggarrange(heatmap_plt, GO_plt, 
                           ncol = 1, nrow = 2, 
                           labels = c('A', 'B'))

# View
combined_plot
plot_dir <- paste0(work_directory, '/report/41_paper_plot/')
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
ggsave(paste0(plot_dir, "loading_con_and_GO.pdf"),
       width = 11, height = 9)
