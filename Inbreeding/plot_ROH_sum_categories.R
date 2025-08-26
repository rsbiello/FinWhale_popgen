###############################################
# ROH Analysis & Visualization
# Author: [Roberto Biello]
# Date:   [2024]
# ---------------------------------------------
# Description:
#   - Classifies ROH lengths into categories
#   - Summarises counts and sums per individual & population
#   - Produces stacked barplots and scatterplots
# ---------------------------------------------
###############################################

# --- Load libraries ---
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)   # for comma axis labels

# --- Set working directory & input ---
workdir <- "path/to/roh/results"
setwd(workdir)

input_file <- "all_samples.roh_length"
output_summary <- "classified_summary.txt"
output_filtered <- "filtered_classified_summary.txt"

# --- Load data ---
# Expected format: Pop, ID, Value (ROH length in bp)
data <- read.table(input_file, header = FALSE,
                   col.names = c("Pop", "ID", "Value"))

head(data)

# --- Classify ROH length into bins ---
data <- data %>%
  mutate(Range = case_when(
    Value >= 1e6   & Value <= 2.5e6  ~ "D_1-2.5Mb",
    Value >  2.5e6 & Value <= 5e6    ~ "C_2.5-5Mb",
    Value >  5e6   & Value <= 1e7    ~ "B_5-10Mb",
    Value >  1e7                     ~ "A_>10Mb",
    TRUE ~ "NA"
  ))

# --- Summarise counts & sums ---
summary_data <- data %>%
  group_by(ID, Pop, Range) %>%
  summarise(
    Count = n(),
    Sum   = sum(Value),
    .groups = "drop"
  )

# Save full summary
write.table(summary_data, output_summary,
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Filter out "NA" bins ---
filtered_data <- summary_data %>%
  filter(Range != "NA")

write.table(filtered_data, output_filtered,
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Plot: ROH counts per individual ---
ggplot(filtered_data, aes(x = ID, y = Count, fill = Range)) +
  geom_col(position = "stack") +
  labs(title = "ROH counts by length class",
       x = "Individual", y = "Count") +
  facet_grid(~ Pop, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- Plot: ROH total lengths per individual ---
ggplot(filtered_data, aes(x = ID, y = Sum, fill = Range)) +
  geom_col(position = "stack") +
  labs(title = "Total ROH length by class",
       x = "Individual", y = "Total ROH length (bp)") +
  facet_grid(~ Pop, scales = "free_x", space = "free") +
  scale_y_continuous(labels = comma) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- Reshape to wide format ---
summary_data_wide <- filtered_data %>%
  pivot_wider(names_from = Range, values_from = c(Sum, Count),
              names_glue = "{Range}_{.value}") %>%
  rowwise() %>%
  mutate(
    Total_Sum   = sum(c_across(ends_with("_Sum")), na.rm = TRUE),
    Total_Count = sum(c_across(ends_with("_Count")), na.rm = TRUE)
  ) %>%
  ungroup()

head(summary_data_wide)

# --- Scatterplot: Total ROH burden ---
ggplot(summary_data_wide, aes(x = Total_Sum, y = Total_Count, color = Pop)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "Total ROH burden per individual",
       x = "Total ROH length (bp)", y = "Total ROH count") +
  scale_x_continuous(labels = comma) +
  theme_bw()

# --- Scatterplot faceted by population ---
ggplot(summary_data_wide, aes(x = Total_Sum, y = Total_Count, color = Pop)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Total ROH burden by population",
       x = "Total ROH length (bp)", y = "Total ROH count") +
  scale_x_continuous(labels = comma) +
  facet_wrap(~ Pop, scales = "free") +
  theme_bw()
