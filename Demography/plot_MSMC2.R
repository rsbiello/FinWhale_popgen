#================================================================
# Plot MSMC2 Results with Historical Periods Highlighted
# Author: Roberto Biello
# Date: 2025
#================================================================

library(ggplot2)

#------------------------
# User-configurable variables
#------------------------
setwd("/path/to/MSMC2/output/")  # replace with your path

mu <- 1.38e-8    # mutation rate per site per generation
gen <- 25.9      # generation time in years

# Sample groups and corresponding files
samples <- list(
  POP1 = c("ind1.msmc.final.txt"),
  POP2 = c("ind2.msmc.final.txt"),
  POP3 = c("ind3.msmc.final.txt"),
  POP4 = c("ind4.msmc.final.txt"),
  POP5 = c("ind5.msmc.final.txt")
)

# Colors per group
group_colors <- c(POP1="cyan", POP2="yellow", POP3="red", POP4="blue", POP5="darkgreen")

# Historical periods (years ago)
periods <- data.frame(
  start = c(3e6, 1.2e6, 6e5, 5.3e5, 4e5, 3.3e5, 2.4e5, 1.3e5, 11700),
  end   = c(2.5e6, 7e5, 5.7e5, 5.1e5, 3.7e5, 3e5, 2.1e5, 1.15e5, 1.1e4),
  color = c(rgb(0.7,0.2,0.3,0.3), rgb(0.2,0.7,0.2,0.3), rep(rgb(0.7,0.5,0.1,0.3),5), rgb(0.7,0.7,0.2,0.3)),
  label = c("PPT","MPT","MIS15e","MIS13e","MIS11e","MIS9e","MIS7e","Eemian","Holocene")
)

#------------------------
# Function to plot a single MSMC line
#------------------------
plot_msmc <- function(file, color) {
  data <- read.table(file, header=TRUE)
  x <- data$left_time_boundary / mu * gen
  y <- 1 / data$lambda / mu
  lines(x, y, type="s", col=color)
}

#------------------------
# Initialize plot
#------------------------
plot(NA, NA, log="x",
     xlim=c(1.1e4, 1e7), ylim=c(0, 2.5e5),
     xlab="Years ago", ylab="Effective population size",
     main="MSMC2 Inference")

#------------------------
# Add shaded historical periods
#------------------------
for(i in 1:nrow(periods)) {
  rect(periods$start[i], 0, periods$end[i], 2.5e5,
       col=periods$color[i], border=NA)
}

#------------------------
# Add MSMC lines for each sample group
#------------------------
for(group in names(samples)) {
  for(file in samples[[group]]) {
    plot_msmc(file, group_colors[group])
  }
}

#------------------------
# Add legend
#------------------------
legend("topright", legend=names(samples),
       col=group_colors, lty=1, cex=0.8)
