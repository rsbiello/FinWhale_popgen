# Load necessary libraries
library(ggplot2)
library(reshape2)
library(readr)

# Set working directory
setwd("/dir/")

# Load the data
file_path <- "count_geno.GERP_0.CDS"  # Replace with your file path
df <- read_table(file_path)

# Melt the data for ggplot
df_melted <- melt(df, id.vars = "POP", 
                  measure.vars = c("HETEROZYGOUS_SITES", "HOMOZYGOUS_DERIVED_SITES"),
                  variable.name = "Genotype_Type", 
                  value.name = "Proportion")

head(df_melted)

# Specify the order for the POP factor levels
df_melted$POP <- factor(df_melted$POP, levels = c("NPA", "SOC", "SVA", "ICE", "MED"))

# Custom colors for each population
pop_colors <- c(
  "POP1" = "#21E000",
  "POP2" = "#FFF803",
  "POP3" = "#3C4400",
  "POP4" = "#3131FF",
  "POP5" = "#FF0037"
)

# Create the boxplot with vertical faceting by genotype type
p_grid <- ggplot(df_melted, aes(x = POP, y = Proportion, fill = POP)) +
  geom_boxplot() +
  facet_wrap(~Genotype_Type, scales = "free_y", nrow = 2, strip.position = "right") +
  scale_fill_manual(values = pop_colors) +
  theme_bw(base_size = 14) +
  labs(y = "Proportion of Derived Genotype", x = "Population") +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

# Display plot
p_grid
