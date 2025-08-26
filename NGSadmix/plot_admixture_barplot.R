#================================================================
# R script to plot NGSadmix results for multiple runs
# Author: Your Name
# Date: 2025-08-26
#================================================================

#------------------------
# Set working directory
#------------------------
# Replace with your own path or use a relative path
setwd("/path/to/NGSadmix/output/qopt/")

#------------------------
# Load population assignment file
#------------------------
pop <- read.table("pop_list.txt", header = FALSE)
head(pop)

#------------------------
# List of NGSadmix .qopt files to plot
#------------------------
q_files <- c(...)

# Read all Q matrices into a list
q_list <- lapply(q_files, read.table)

#------------------------
# Order individuals by population
#------------------------
ord <- order(pop[, 1])

#------------------------
# Function to plot admixture barplot
#------------------------
plot_admixture <- function(q_matrix, pop_order, pop_labels = pop[,1], colors = 2:10) {
  barplot(
    t(q_matrix)[, pop_order],
    col = colors,
    space = 0,
    border = NA,
    xlab = "Individuals",
    ylab = "Admixture proportions"
  )
  
  # Add population labels below the barplot
  text(
    tapply(1:length(pop_order), pop_labels[pop_order], mean),
    -0.05,
    unique(pop_labels[pop_order]),
    xpd = TRUE
  )
  
  # Add vertical lines separating populations
  abline(
    v = cumsum(sapply(unique(pop_labels[pop_order]), function(x) sum(pop_labels[pop_order] == x))),
    col = 1,
    lwd = 1.2
  )
}

#------------------------
# Plot all Q matrices
#------------------------
for (q in q_list) {
  plot_admixture(q, ord)
}
