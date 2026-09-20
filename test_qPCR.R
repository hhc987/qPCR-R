library(tidyverse)
library(janitor)
library(ggplot2)
library(RSelenium)

raw_data <- Class_qPCR_Data %>% 
  clean_names()

# Extract no. from well str and map wells to their respective grps
cleaned_groups <- raw_data %>%
  mutate(
    well_column = as.numeric(str_extract(well, "\\d+")),
    group_num = case_when(
      well_column %in% c(1, 2) ~ "Group 1",
      well_column %in% c(3, 4) ~ "Group 2",
      well_column %in% c(5, 6) ~ "Group 3",
      well_column %in% c(7, 8) ~ "Group 4",
      well_column %in% c(9, 10) ~ "Group 5",
      well_column %in% c(11, 12) ~ "Group 6",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(cq))

### THE CORRECTIONS

# AVERAGING THE REPLICATES IN EACH GRP
group_averages <- cleaned_groups %>%
  group_by(group_num, target, sample) %>%
  summarise(avg_cq = mean(cq), .groups = "drop")

# Assign experiments and clean replicate numbers (1, 2, 3) to grps
biological_data <- group_averages %>%
  mutate(
    experiment = case_when(
      group_num %in% c("Group 1", "Group 2", "Group 3") ~ "TNF-a",
      group_num %in% c("Group 4", "Group 5", "Group 6") ~ "IL-6"
    ),
    # match grps to replicate, here we use each grp ave. as a replicate for calc.
    replicate = case_when(
      group_num %in% c("Group 1", "Group 4") ~ 1,
      group_num %in% c("Group 2", "Group 5") ~ 2,
      group_num %in% c("Group 3", "Group 6") ~ 3
    )
  )

# Separate GAPDH and Targets
gapdh_df <- biological_data %>%
  filter(target == "GAPDH") %>%
  select(experiment, replicate, sample, gapdh_cq = avg_cq)

targets_df <- biological_data %>%
  filter(target != "GAPDH") %>%
  select(experiment, target, replicate, sample, target_cq = avg_cq)

# qPCR calc.
qpcr_results <- targets_df %>%
  # join by experiment AND replicate AND sample
  inner_join(gapdh_df, by = c("experiment", "replicate", "sample")) %>%
  
  # Calculate dCt per replicate
  mutate(delta_ct = target_cq - gapdh_cq)

# Format into matrix
matrix <- qpcr_results %>%
  pivot_wider(
    id_cols = c(experiment, replicate),
    names_from = sample,
    values_from = c(gapdh_cq, target_cq, delta_ct)
  ) %>%
  mutate(
    delta_delta_ct = delta_ct_Infected - delta_ct_Uninfected,
    fold_change = 2^(-delta_delta_ct)
  )

print(matrix)
view(matrix)

fold_change_table <- matrix %>%
  group_by(experiment) %>%
  summarise(
    Infected = mean(fold_change),
    Uninfected = 1,
    .groups = "drop"
  ) %>%
  # Flip the columns into rows
  pivot_longer(
    cols = c(Uninfected, Infected),
    names_to = "Sample",
    values_to = "Relative Fold Change"
  ) %>%
  select("Target Gene" = experiment, Sample, "Relative Fold Change")

print(as.data.frame(fold_change_table))
View(as.data.frame(fold_change_table))

# graph

ggplot(fold_change_table, aes(x = `Target Gene` , y = `Relative Fold Change`, fill = Sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(
    aes(label = round(`Relative Fold Change`, 2)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  # Custom y label for superscript adapted from stack overflow, labs doesn't work?
  ylab(expression(paste("Relative Fold Change ( ", 2, phantom()^{-Delta*Delta*Cq}, ")"))) +
  labs(
    title = "Relative Fold Change of Target Genes",
    subtitle = "Comparing Fold Change of Infected and Uninfected samples",
    x = "Target Gene"
  ) +
  scale_fill_manual(values = c("Infected" = "#a00000", "Uninfected" = "blue4")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey50"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )



