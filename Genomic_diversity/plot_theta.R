#------------------------
# Load libraries
#------------------------
library(ggplot2)
library(viridis)       # optional: color palettes
library(hrbrthemes)    # optional: themes

#------------------------
# Set working directory
#------------------------
setwd("/path/to/analyses/theta/")  # update path as needed

#------------------------
# Load data
#------------------------
theta <- read.table("theta.txt", header = TRUE)
head(theta)

#------------------------
# Create dot plot
#------------------------
dotplot <- ggplot(theta, aes(x = pop, y = theta, color = pos, shape = pos)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_viridis(discrete = TRUE, option = "D") + # optional: nicer colors
  ylim(0.0010, 0.0028) +
  xlab("Populations") +
  ylab("Theta W") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))

#------------------------
# Display plot
#------------------------
print(dotplot)
