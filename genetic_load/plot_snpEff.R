# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory
setwd("...")

# Load the data
df <- read.table("count_snpeff.txt", header = TRUE)
head(df)

# -----------------------------
# 1) Reshape ratio data
# -----------------------------
df_long_ratios <- df %>%
  pivot_longer(
    cols = c("HeterozygousRatio", "HomoDerivedRatio"),
    names_to = "RatioType",      
    values_to = "RatioValue"     
  )

# Set factor levels for populations
df_long_ratios$Population <- factor(df_long_ratios$Population, levels = c("NPA", "SOC", "SVA", "ICE", "MED"))

# Custom population colors
my_pop_colors <- c(
  "POP1" = "#3131FF",
  "POP2" = "#3C4400",
  "POP3" = "#FF0037",
  "POP4" = "#21E000",
  "POP5" = "#FFF803"
)

# Boxplot for ratios (RatioType x Category)
p_ratios <- ggplot(df_long_ratios, aes(x = Population, y = RatioValue, fill = Population)) +
  geom_boxplot() +
  facet_grid(RatioType ~ Category, scales = "free_y") +
  scale_fill_manual(values = my_pop_colors) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "Ratio", title = "Heterozygous and Homozygous Derived Ratios by Population")
p_ratios

# -----------------------------
# 2) Reshape counts data
# -----------------------------
df_long_counts <- df %>%
  pivot_longer(
    cols = c(Heterozygous, HomoDerived),
    names_to = "Genotype_Type",
    values_to = "Genotype_Count"
  )

# Boxplot for counts (Genotype_Type x Category)
p_counts <- ggplot(df_long_counts, aes(x = Population, y = Genotype_Count, fill = Population)) +
  geom_boxplot() +
  facet_grid(Category ~ Genotype_Type, scales = "free_y") +
  scale_fill_manual(values = my_pop_colors) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(y = "Genotype Counts", title = "Heterozygous and Homozygous Derived Counts by Population")
p_counts
