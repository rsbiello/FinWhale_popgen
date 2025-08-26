library(ggplot2)
library(tidyverse)

# Set working directory
setwd("dir")

# Load table
data <- read.table("froh", header = TRUE)

# Ensure factors are ordered (optional: order by ROH or by population)
data$ind <- factor(data$ind, levels = data$ind[order(data$ROH)])

# Bar Plot
plot <- ggplot(data, aes(x = ind, y = ROH, fill = pop)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(x = "Individuals", y = "FROH") +
  coord_flip() +                               # flip axes
  coord_cartesian(ylim = c(0, 0.5)) +          # <- adjust to your scale (maybe 0–0.5 instead of 0–50?)
  theme(
    axis.text.y = element_text(size = 8),      # smaller labels if many inds
    legend.position = "top"
  )

print(plot)
