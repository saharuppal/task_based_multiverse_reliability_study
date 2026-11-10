library(dplyr)
library(tidyr)
library(flextable)
library(readr)

# Read ICC CSV
all_icc_values <- read_csv("all_icc_masked_values.csv")

# Filter required copes
filtered_data <- all_icc_values %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope-3") |
      (task == "oasis" & cope %in% c("cope-4", "cope-5"))
  ) %>%
  mutate(
    task_cope = paste(task, cope, sep = "_"),
    pipeline = gsub("_", "-", pipeline)
  )

# ---- numeric ICC table (for max detection)
icc_numeric <- filtered_data %>%
  select(pipeline, task_cope, icc) %>%
  pivot_wider(names_from = task_cope, values_from = icc)

# ---- formatted display table
table_data <- filtered_data %>%
  mutate(
    display = sprintf("%.3f\n(%.3f, %.3f)", icc, ci_lower, ci_upper)
  ) %>%
  select(pipeline, task_cope, display) %>%
  pivot_wider(names_from = task_cope, values_from = display)

# ---- Rename columns FIRST
icc_numeric <- icc_numeric %>%
  rename(
    Pipeline = pipeline,
    GoNoGo = `gonogo_cope-3`,
    MID = `mid_cope-3`,
    `N-Back` = `nback_cope-3`,
    `OASIS (pos>neu)` = `oasis_cope-4`,
    `OASIS (neg>neu)` = `oasis_cope-5`
  )

table_data <- table_data %>%
  rename(
    Pipeline = pipeline,
    GoNoGo = `gonogo_cope-3`,
    MID = `mid_cope-3`,
    `N-Back` = `nback_cope-3`,
    `OASIS (pos>neu)` = `oasis_cope-4`,
    `OASIS (neg>neu)` = `oasis_cope-5`
  )

# ---- Create flextable
ft <- flextable(table_data) %>%
  autofit() %>%
  theme_booktabs()

cols <- setdiff(colnames(table_data), "Pipeline")

# Center text
ft <- align(ft, j = cols, align = "center", part = "all")
ft <- valign(ft, j = cols, valign = "center")
ft <- line_spacing(ft, space = 1.2)

# ---- Bold highest ICC per column
for (col in cols) {
  max_val <- max(icc_numeric[[col]], na.rm = TRUE)
  
  ft <- bold(
    ft,
    i = icc_numeric[[col]] == max_val,
    j = col,
    bold = TRUE
  )
}

ft




library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

# -------------------------
# Read ICC data
# -------------------------
icc_data <- read_csv("all_icc_masked_values.csv") %>%
  mutate(pipeline = gsub("_", "-", pipeline))

# -------------------------
# Keep copes of interest
# -------------------------
filtered_data <- icc_data %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope-3") |
      (task == "oasis" & cope %in% c("cope-4", "cope-5"))
  ) %>%
  mutate(
    cope_label = case_when(
      task == "gonogo" ~ "GoNoGo",
      task == "mid" ~ "MID",
      task == "nback" ~ "N-Back",
      task == "oasis" & cope == "cope-4" ~ "OASIS (pos>neu)",
      task == "oasis" & cope == "cope-5" ~ "OASIS (neg>neu)"
    )
  )

# -------------------------
# Convert to matrix
# -------------------------
heatmap_df <- filtered_data %>%
  select(pipeline, cope_label, icc) %>%
  pivot_wider(names_from = cope_label, values_from = icc)

heatmap_mat <- heatmap_df %>%
  column_to_rownames("pipeline") %>%
  as.matrix()

# -------------------------
# Column order requested
# -------------------------
heatmap_mat <- heatmap_mat[, c(
  "N-Back",
  "OASIS (pos>neu)",
  "OASIS (neg>neu)",
  "MID",
  "GoNoGo"
)]

# -------------------------
# Colour scaling to data
# -------------------------
min_icc <- min(heatmap_mat, na.rm = TRUE)
max_icc <- max(heatmap_mat, na.rm = TRUE)

breaks <- seq(min_icc, max_icc, length.out = 100)

heat_colors <- colorRampPalette(c(
  "brown1",
  "#FEE08B",
  "#1A9850"
))(99)

# -------------------------
# Plot heatmap
# -------------------------
# pheatmap(
# heatmap_mat,
# cluster_rows = TRUE,
# cluster_cols = FALSE,
# clustering_distance_rows = "correlation",
# clustering_method = "average",
# display_numbers = TRUE,
# number_format = "%.2f",
# number_color = "black",   # ← black font
# color = heat_colors,
# breaks = breaks,
# angle_col = 45,
# border_color = NA
# )

# ── Adjustable font sizes & angles ────────────────
cell_number_size  <- 14
axis_label_size   <- 12
col_label_angle   <- 45    # try 315 for the other direction, or 90 for vertical
# ──────────────────────────────────────────────────

pipeline_order <- c(
  "6HMP",
  "24HMP",
  "6HMP_WM_CSF",
  "24HMP_WM_CSF",
  "6HMP_aCompCor",
  "24HMP_aCompCor",
  "6HMP_spike_reg",
  "24HMP_spike_reg",
  "6HMP_WM_CSF_spike_reg",
  "24HMP_WM_CSF_spike_reg",
  "6HMP_aCompCor_spike_reg",
  "24HMP_aCompCor_spike_reg"
  
)

pipeline_order <- gsub("_", "-", pipeline_order)
heatmap_mat <- heatmap_mat[pipeline_order, , drop = FALSE]

pheatmap(
  heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_row = 12,
  cellheight = 35,
  color = heat_colors,
  breaks = breaks,
  angle_col = col_label_angle,
  border_color = NA,
  fontsize = axis_label_size,
  fontsize_number = cell_number_size
)

setdiff(pipeline_order, rownames(heatmap_mat))   # in your order, missing from matrix
setdiff(rownames(heatmap_mat), pipeline_order)   # in matrix, missing from your order
dim(heatmap_mat)


library(dplyr)
library(tidyr)
library(flextable)
library(readr)

# Read ICC CSV
all_icc_values <- read_csv("all_icc_masked_values.csv")

# Filter required copes
filtered_data <- all_icc_values %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope-3") |
      (task == "oasis" & cope %in% c("cope-4", "cope-5"))
  ) %>%
  mutate(
    task_cope = paste(task, cope, sep = "_"),
    pipeline = gsub("_", "-", pipeline)
  )

# ---- numeric ICC table (for max detection)
icc_numeric <- filtered_data %>%
  select(pipeline, task_cope, icc) %>%
  pivot_wider(names_from = task_cope, values_from = icc)

# ---- formatted display table
table_data <- filtered_data %>%
  mutate(
    display = sprintf("%.3f\n(%.3f, %.3f)", icc, ci_lower, ci_upper)
  ) %>%
  select(pipeline, task_cope, display) %>%
  pivot_wider(names_from = task_cope, values_from = display)

# ---- Rename columns FIRST
icc_numeric <- icc_numeric %>%
  rename(
    Pipeline = pipeline,
    GoNoGo = `gonogo_cope-3`,
    MID = `mid_cope-3`,
    `N-Back` = `nback_cope-3`,
    `OASIS (pos>neu)` = `oasis_cope-4`,
    `OASIS (neg>neu)` = `oasis_cope-5`
  )

table_data <- table_data %>%
  rename(
    Pipeline = pipeline,
    GoNoGo = `gonogo_cope-3`,
    MID = `mid_cope-3`,
    `N-Back` = `nback_cope-3`,
    `OASIS (pos>neu)` = `oasis_cope-4`,
    `OASIS (neg>neu)` = `oasis_cope-5`
  )

# ---- Create flextable
ft <- flextable(table_data) %>%
  autofit() %>%
  theme_booktabs()

cols <- setdiff(colnames(table_data), "Pipeline")

# Center text
ft <- align(ft, j = cols, align = "center", part = "all")
ft <- valign(ft, j = cols, valign = "center")
ft <- line_spacing(ft, space = 1.2)

# ---- Bold highest ICC per column
for (col in cols) {
  max_val <- max(icc_numeric[[col]], na.rm = TRUE)
  
  ft <- bold(
    ft,
    i = icc_numeric[[col]] == max_val,
    j = col,
    bold = TRUE
  )
}

ft





