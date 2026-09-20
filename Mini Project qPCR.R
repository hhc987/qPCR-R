library(tidyverse)
library(janitor)

Raw_data <- Class_qPCR_Data |> 
  clean_names()

# Extract no. from well str and map wells to their respective grps
grps <- Raw_data |> 
  mutate(
    well = as.numeric(str_extract(well, "\\d+")),
    grp_no = case_when(
      well %in% c(1, 2) ~ "Group 1",
      well %in% c(3, 4) ~ "Group 2",
      well %in% c(5, 6) ~ "Group 3",
      well %in% c(7, 8) ~ "Group 4",
      well %in% c(9, 10) ~ "Group 5",
      well %in% c(11, 12) ~ "Group 6",
      TRUE ~ NA_character_
    )
  ) |> 
  filter(!is.na(cq))

### THE CORRECTIONS

# AVERAGING THE REPLICATES IN EACH GRP
grp_ave <- grps |> 
  group_by(grp_no, target, sample) |> 
  summarise(ave_cq = mean(cq), .groups = "drop")

# Assign experiments and clean replicate numbers (1, 2, 3) to grps
formatted_data <- grp_ave |> 
  mutate(
    expt = case_when(
      grp_no %in% c("Group 1", "Group 2", "Group 3") ~ "TNF-a",
      grp_no %in% c("Group 4", "Group 5", "Group 6") ~ "IL-6"
    ),
    # match grps to replicate, here we use each grp ave. as a replicate for calc.
    replicate = case_when(
      grp_no %in% c("Group 1", "Group 4") ~ 1,
      grp_no %in% c("Group 2", "Group 5") ~ 2,
      grp_no %in% c("Group 3", "Group 6") ~ 3
    )
  )

# Separate GAPDH and Targets
gapdh_df <- formatted_data |> 
  filter(target == "GAPDH") |> 
  select(expt, target, replicate, sample, gapdh_cq = ave_cq)
targets_df <- formatted_data |> 
  filter(target != "GAPDH") |> 
  select(expt, target, replicate, sample, target_cq = ave_cq)

# calc.
qPCR_calc <- targets_df |> 
  inner_join(gapdh_df, by = c("expt", "replicate", "sample")) |> 
  mutate(dct = target_cq - gapdh_cq)

# Format into matrix
matrix <- qPCR_calc |> 
  pivot_wider(
    id_cols = c(expt, replicate),
    names_from = sample,
    values_from = c(gapdh_cq, target_cq, dct)
  ) |> 
  mutate(
    ddct = dct_Infected - dct_Uninfected,
    fold_change = 2^(-ddct)
  )

print(matrix)
# view(matrix)

fold_change_table <- matrix |> 
  group_by(expt) |> 
  summarise(
    Infected = mean(fold_change),
    Uninfected = 1,
    .groups = "drop"
  ) |> 
  # Flip the columns into rows
  pivot_longer(
    cols = c(Uninfected, Infected),
    names_to = "Sample",
    values_to = "Relative Fold Change"
  ) |> 
  select("Target Gene" = expt, Sample, "Relative Fold Change")

print(as.data.frame(fold_change_table))
# View(as.data.frame(fold_change_table))

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
  # Custom y label for superscript adapted from stack overflow
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



