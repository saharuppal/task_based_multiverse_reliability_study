library(dplyr)
library(tidyr)
library(flextable)
library(readr)
library(officer)
# -------------------------
# Read CSV
# -------------------------
voxel_data <- read_csv("all_mean_zstat_voxelwise.csv")

# -------------------------
# Filter same tasks/copes
# -------------------------
filtered_data <- voxel_data %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope3") |
      (task == "oasis" & cope %in% c("cope4", "cope5"))
  ) %>%
  mutate(
    task_cope = paste(task, cope, sep = "_"),
    pipeline  = gsub("_", "-", pipeline),
    mean_disp = sprintf("%.2f (%.2f)", mean_z, sd_z),
    vox_disp  = format(n_voxels, big.mark = ",")
  )

# -------------------------
# Numeric voxel table (for bolding)
# -------------------------
voxel_numeric <- filtered_data %>%
  select(pipeline, task_cope, n_voxels) %>%
  pivot_wider(names_from = task_cope, values_from = n_voxels)

mean_numeric <- filtered_data %>%
  select(pipeline, task_cope, mean_z) %>%
  pivot_wider(names_from = task_cope, values_from = mean_z) %>%
  rename(!!!setNames(task_map, names(task_map)))

# -------------------------
# Two separate wide tables: mean and voxels
# -------------------------
mean_wide <- filtered_data %>%
  select(pipeline, task_cope, mean_disp) %>%
  pivot_wider(names_from = task_cope, values_from = mean_disp)

vox_wide <- filtered_data %>%
  select(pipeline, task_cope, vox_disp) %>%
  pivot_wider(names_from = task_cope, values_from = vox_disp)

# -------------------------
# Rename columns
# -------------------------
task_map <- c(
  Pipeline          = "pipeline",
  GoNoGo            = "gonogo_cope3",
  MID               = "mid_cope3",
  `N-Back`          = "nback_cope3",
  `OASIS (pos>neu)` = "oasis_cope4",
  `OASIS (neg>neu)` = "oasis_cope5"
)

rename_df <- function(df) {
  df %>% rename(!!!setNames(task_map, names(task_map)))
}

mean_wide    <- rename_df(mean_wide)
vox_wide     <- rename_df(vox_wide)
voxel_numeric <- rename_df(voxel_numeric)

task_cols <- setdiff(names(mean_wide), "Pipeline")

# -------------------------
# Interleave mean and voxel columns: mean_GoNoGo, vox_GoNoGo, mean_MID, ...
# -------------------------
combined <- mean_wide %>% select(Pipeline)

for (col in task_cols) {
  combined[[paste0("mean_", col)]] <- mean_wide[[col]]
  combined[[paste0("vox_",  col)]] <- vox_wide[[col]]
}

# After creating combined, replace "NA (NA)" with "NA"
combined <- combined %>%
  mutate(across(everything(), ~gsub("NA \\(NA\\)", "NA", .)))


# -------------------------
# Build flextable
# -------------------------
ft <- flextable(combined) %>%
  theme_booktabs()

# ---- Header row 1: span task name across mean + vox pair
all_cols   <- colnames(combined)
mean_cols  <- grep("^mean_", all_cols, value = TRUE)
vox_cols   <- grep("^vox_",  all_cols, value = TRUE)

# Add a second header row with sub-labels
ft <- set_header_labels(ft,
                        values = setNames(
                          c("Pipeline", rep(c("Mean Z (SD)", "Voxels"), length(task_cols))),
                          all_cols
                        )
)

# Add spanning top header for each task
ft <- add_header_row(ft,
                     values = c("Pipeline", task_cols),
                     colwidths = c(1, rep(2, length(task_cols)))
)

# ---- Formatting
all_data_cols <- c(mean_cols, vox_cols)

ft <- align(ft,  j = all_data_cols, align = "center", part = "all")
ft <- align(ft,  j = "Pipeline",    align = "left",   part = "all")
ft <- valign(ft, j = all_data_cols, valign = "center")
ft <- line_spacing(ft, space = 1.1)
ft <- autofit(ft)

# After autofit(ft), add right border after each task's vox column
for (col in task_cols) {
  vox_col_name <- paste0("vox_", col)
  ft <- border(ft, 
               j = vox_col_name,
               border.right = fp_border(color = "white", width = 1.5),
               part = "all")
}

# -------------------------
# Bold highest voxels and highest mean Z separately
# -------------------------
for (col in task_cols) {
  
  vox_col_name  <- paste0("vox_", col)
  mean_col_name <- paste0("mean_", col)
  
  # highest voxel count
  max_vox   <- max(voxel_numeric[[col]], na.rm = TRUE)
  vox_rows  <- which(voxel_numeric[[col]] == max_vox)
  
  # highest mean Z
  max_mean  <- max(mean_numeric[[col]], na.rm = TRUE)
  mean_rows <- which(mean_numeric[[col]] == max_mean)
  
  # bold voxels
  ft <- bold(ft, i = vox_rows, j = vox_col_name, bold = TRUE)
  
  # bold mean Z
  ft <- bold(ft, i = mean_rows, j = mean_col_name, bold = TRUE)
}

# -------------------------
# Bold highest n_voxels per task column
# -------------------------
# After autofit(ft)
for (col in task_cols) {
  vox_col_name <- paste0("vox_", col)
  ft <- padding(ft, j = vox_col_name, padding.right = 20, part = "all")
}

ft


# -------------------------# -------------------------# -------------------------# -------------------------# -------------------------
# -------------------------# -------------------------# -------------------------# -------------------------# -------------------------

library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

# -------------------------
# Read data
# -------------------------
voxel_data <- read_csv("all_mean_zstat_voxelwise.csv") %>%
  mutate(pipeline = gsub("_", "-", pipeline))

# Preserve table order
pipeline_order <- voxel_data %>%
  distinct(pipeline) %>%
  pull(pipeline)

# -------------------------
# Function for single-bar tasks
# -------------------------
make_single_task_plot <- function(data, title_text, bar_colour) {
  
  mean_val <- mean(data$n_voxels, na.rm = TRUE)
  
  ggplot(data, aes(x = factor(pipeline, levels = pipeline_order),
                   y = n_voxels)) +
    geom_col(fill = bar_colour, width = 0.7) +
    geom_hline(yintercept = mean_val,
               linetype = "dashed",
               linewidth = 0.7) +
    labs(
      x = "Pipeline",
      y = "Number of voxels",
      title = title_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
}

# -------------------------
# Task datasets
# -------------------------
gonogo_df <- voxel_data %>%
  filter(task == "gonogo", cope == "cope3")

mid_df <- voxel_data %>%
  filter(task == "mid", cope == "cope3")

nback_df <- voxel_data %>%
  filter(task == "nback", cope == "cope3")

# OASIS grouped dataset
oasis_df <- voxel_data %>%
  filter(task == "oasis", cope %in% c("cope4", "cope5")) %>%
  mutate(
    cope_label = case_when(
      cope == "cope4" ~ "pos>neu",
      cope == "cope5" ~ "neg>neu"
    ),
    pipeline = factor(pipeline, levels = pipeline_order)
  )

# -------------------------
# Build plots
# -------------------------
p1 <- make_single_task_plot(gonogo_df, "(A) Go/No-Go", "aquamarine4")
p2 <- make_single_task_plot(mid_df, "(B) MID", "darkslateblue")
p3 <- make_single_task_plot(nback_df, "(C) N-back", "darkgreen")

# OASIS grouped bars
mean_oasis <- mean(oasis_df$n_voxels, na.rm = TRUE)

p4 <- ggplot(oasis_df,
             aes(x = pipeline,
                 y = n_voxels,
                 fill = cope_label)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.5   # thinner bars
  ) +
  geom_hline(yintercept = mean_oasis,
             linetype = "dashed",
             linewidth = 0.7) +
  scale_fill_manual(
    name = "Contrast",
    values = c(
      "pos>neu" = "#C040c0",
      "neg>neu" = "#804080"
    )
  ) +
  labs(
    x = "Pipeline",
    y = "Number of voxels",
    title = "(D) OASIS"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# -------------------------
# Combine panel evenly
# -------------------------
panel_plot <- (p1 | p2) / (p3 | p4)

panel_plot




# -------------------------# -------------------------# -------------------------# -------------------------# -------------------------
# -------------------------# -------------------------# -------------------------# -------------------------# -------------------------
install.packages("ggtext")
library(dplyr)
library(ggtext)
library(ggplot2)
library(readr)
library(patchwork)

# -------------------------
# Read data
# -------------------------
voxel_data <- read_csv("all_mean_zstat_voxelwise.csv") %>%
  mutate(pipeline = gsub("_", "-", pipeline))

# Preserve table order
pipeline_order <- voxel_data %>%
  distinct(pipeline) %>%
  pull(pipeline)

# -------------------------
# Function for single-bar tasks
# -------------------------
make_single_task_plot <- function(data, title_text, bar_colour) {
  
  mean_val <- mean(data$n_voxels, na.rm = TRUE)
  
  ggplot(data, aes(x = factor(pipeline, levels = pipeline_order),
                   y = n_voxels)) +
    geom_col(fill = bar_colour, width = 0.7) +
    geom_hline(yintercept = mean_val,
               linetype = "dashed",
               linewidth = 0.7) +
    labs(
      x = "Pipeline",
      y = "Number of voxels",
      title = title_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_markdown(face = "bold")  # key change
    )
}

# -------------------------
# Task datasets
# -------------------------
gonogo_df <- voxel_data %>%
  filter(task == "gonogo", cope == "cope3")

mid_df <- voxel_data %>%
  filter(task == "mid", cope == "cope3")

nback_df <- voxel_data %>%
  filter(task == "nback", cope == "cope3")

# OASIS grouped dataset
oasis_df <- voxel_data %>%
  filter(task == "oasis", cope %in% c("cope4", "cope5")) %>%
  mutate(
    cope_label = case_when(
      cope == "cope4" ~ "pos>neu",
      cope == "cope5" ~ "neg>neu"
    ),
    pipeline = factor(pipeline, levels = pipeline_order)
  )

# -------------------------
# Build plots
# -------------------------
p1 <- make_single_task_plot(
  gonogo_df,
  "(A) <span style='font-size:18pt'>Go/No-Go</span>",
  "aquamarine4"
)

p2 <- make_single_task_plot(
  mid_df,
  "(B) <span style='font-size:18pt'>MID</span>",
  "darkslateblue"
)

p3 <- make_single_task_plot(
  nback_df,
  "(C) <span style='font-size:18pt'>N-back</span>",
  "darkgreen"
)

# OASIS grouped bars
mean_oasis <- mean(oasis_df$n_voxels, na.rm = TRUE)

p4 <- ggplot(oasis_df,
             aes(x = pipeline,
                 y = n_voxels,
                 fill = cope_label)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.5   # thinner bars
  ) +
  geom_hline(yintercept = mean_oasis,
             linetype = "dashed",
             linewidth = 0.7) +
  scale_fill_manual(
    name = "Contrast",
    values = c(
      "pos>neu" = "#C040c0",
      "neg>neu" = "#804080"
    )
  ) +
  labs(
    x = "Pipeline",
    y = "Number of voxels",
    title = "(D) <span style='font-size:18pt'>OASIS</span>"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_markdown(face = "bold")
  )

# -------------------------
# Combine panel evenly
# -------------------------
panel_plot <- (p1 | p2) / (p3 | p4)

panel_plot
