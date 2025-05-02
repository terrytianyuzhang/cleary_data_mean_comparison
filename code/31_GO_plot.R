rm(list = ls())

library(data.table)
library(ggplot2)

work_directory <- '/Users/tianyuzhang/Documents/SingleCell/code_paper/cleary_data_mean_comparison/'

all_GO_result <- data.table()
for(module_index in c(15,31,37)){
  GO_file <- paste0(work_directory, 'data/final/module_GO/module_', module_index,'.csv')
  GO_result <- fread(GO_file)
  GO_result[, module_index := module_index]
  all_GO_result <- rbind(all_GO_result, GO_result)
}

GO_to_present <- c(###module 15
                   "signaling receptor binding",
                   "growth factor activity",
                   "cell-cell signaling",
                   "response to virus", ###module 31
                   "regulation of response to biotic stimulus",
                   "regulation of innate immune response",
                   "positive regulation of response to biotic stimulus",
                   "cytokine-mediated signaling pathway",
                   "response to chemokine", ##module 37
                   "leukocyte migration",
                   "cellular homeostasis",
                   "G protein-coupled receptor signaling pathway"
                   )
GO_result_plot <- all_GO_result[Description %in% GO_to_present, .(Description, module_index, p.adjust, ID)]
GO_result_plot <- GO_result_plot[p.adjust < 0.05, ]
GO_result_plot[, module_index := as.factor(module_index)]
GO_result_plot[, Description := factor(Description, levels = unique(Description))]

library(ggplot2)

# Shorten the Description labels
GO_result_plot[, Short_Description := fcase(
  Description == "signaling receptor binding", "Signaling receptor binding",
  Description == "growth factor activity", "Growth factor activity",
  Description == "cell-cell signaling", "Cell-cell signaling",
  Description == "response to virus", "Response to virus",
  Description == "regulation of response to biotic stimulus", "Response to biotic stimulus",
  Description == "regulation of innate immune response", "Innate immune response",
  Description == "positive regulation of response to biotic stimulus", "Response to biotic stimulus",
  Description == "cytokine-mediated signaling pathway", "Cytokine-mediated pathway",
  Description == "response to chemokine", "Chemokine response",
  Description == "leukocyte migration", "Leukocyte migration",
  Description == "cellular homeostasis", "Cell homeostasis",
  Description == "G protein-coupled receptor signaling pathway", "GPCR signaling"
)]

# Make sure the new Short_Description factor follows original order
GO_result_plot[, Short_Description := factor(Short_Description, levels = unique(Short_Description))]

p1 <- ggplot(GO_result_plot) +
  geom_point(aes(x = Short_Description, y = module_index, size = -log10(p.adjust)), color = "steelblue") +
  scale_size_continuous(name = expression(-log[10]~"(GO p-value)")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),  # <-- bigger x-axis font
    axis.text.y = element_text(size = 11),                                    # optional: make y-axis font nice too
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "right"
  ) +
  labs(
    y = "Module Index"
    # title = "GO Terms by Module"
  )


plot_dir <- paste0(work_directory, '/report/31_GO_plot/')
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Save the plot
ggsave(filename = paste0(plot_dir, 'GO_plot.pdf'), 
       plot = last_plot(),  # or you can put your ggplot object here
       width = 8, height = 4.8)
saveRDS(p1, paste0(plot_dir, 'GO_plot.rds'))

# Select the columns
table_data <- unique(GO_result_plot[, .(Short_Description, ID)])

# Create a vector to store each line
latex_lines <- c()

# Add beginning of table
latex_lines <- c(latex_lines, "\\begin{table}[ht]")
latex_lines <- c(latex_lines, "\\centering")
latex_lines <- c(latex_lines, "\\begin{tabular}{ll}")
latex_lines <- c(latex_lines, "\\hline")
latex_lines <- c(latex_lines, "\\textbf{Short Description} & \\textbf{GO ID} \\\\")
latex_lines <- c(latex_lines, "\\hline")

# Add table rows
for (i in 1:nrow(table_data)) {
  latex_lines <- c(latex_lines, paste0(table_data$Short_Description[i], " & ", table_data$ID[i], " \\\\"))
}

# End table
latex_lines <- c(latex_lines, "\\hline")
latex_lines <- c(latex_lines, "\\end{tabular}")
latex_lines <- c(latex_lines, "\\caption{Short descriptions and corresponding GO IDs.}")
latex_lines <- c(latex_lines, "\\label{tab:short_desc_go_ids}")
latex_lines <- c(latex_lines, "\\end{table}")

# Now output all at once
cat(paste(latex_lines, collapse = "\n"))
