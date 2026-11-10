library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

rm(list = ls())

# ── Load subject-level data ───────────────────────────────────────────────────
raw <- read_csv("n_suprathreshold.csv")

# ── Step 1: average sessions per subject ──────────────────────────────────────
per_subject <- raw %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope3") |
      (task == "oasis" & cope %in% c("cope4", "cope5"))
  ) %>%
  group_by(subject, task, pipeline, cope) %>%
  summarise(n_suprathreshold = mean(n_suprathreshold, na.rm = TRUE), .groups = "drop")

# ── Step 2: mean across subjects per pipeline x task x cope ───────────────────
summary_data <- per_subject %>%
  group_by(task, pipeline, cope) %>%
  summarise(
    mean_n_suprathreshold = mean(n_suprathreshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pipeline = gsub("_", "-", pipeline))

pipeline_order <- c(
  "6HMP",
  "24HMP",
  "6HMP-WM-CSF",
  "24HMP-WM-CSF",
  "6HMP-aCompCor",
  "24HMP-aCompCor",
  "6HMP-spike-reg",
  "24HMP-spike-reg",
  "6HMP-WM-CSF-spike-reg",
  "24HMP-WM-CSF-spike-reg",
  "6HMP-aCompCor-spike-reg",
  "24HMP-aCompCor-spike-reg"
)

# ── Plot function ─────────────────────────────────────────────────────────────
# Draws a dashed line at the group mean across all pipelines for that panel
make_single_task_plot <- function(data, title_text, bar_colour) {
  panel_mean <- mean(data$mean_n_suprathreshold, na.rm = TRUE)
  
  ggplot(data, aes(x = factor(pipeline, levels = pipeline_order),
                   y = mean_n_suprathreshold)) +
    geom_col(fill = bar_colour, width = 0.7) +
    geom_hline(yintercept = panel_mean, linetype = "dashed", linewidth = 0.5) +
    labs(x = "Pipeline", y = "Number of voxels", title = title_text) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.title      = element_text(face = "bold")
    )
}

# ── Task datasets ─────────────────────────────────────────────────────────────
gonogo_df <- summary_data %>% filter(task == "gonogo", cope == "cope3")
mid_df    <- summary_data %>% filter(task == "mid",    cope == "cope3")
nback_df  <- summary_data %>% filter(task == "nback",  cope == "cope3")

oasis_df <- summary_data %>%
  filter(task == "oasis", cope %in% c("cope4", "cope5")) %>%
  mutate(
    cope_label = case_when(
      cope == "cope4" ~ "pos>neu",
      cope == "cope5" ~ "neg>neu"
    ),
    cope_label = factor(cope_label, levels = c("neg>neu", "pos>neu")),
    pipeline   = factor(pipeline, levels = pipeline_order)
  )

# ── Build plots ───────────────────────────────────────────────────────────────
p1 <- make_single_task_plot(gonogo_df, "(A) Go/No-Go", "aquamarine4")
p2 <- make_single_task_plot(mid_df,    "(B) MID",      "darkslateblue")
p3 <- make_single_task_plot(nback_df,  "(C) N-back",   "darkgreen")

oasis_mean <- mean(oasis_df$mean_n_suprathreshold, na.rm = TRUE)

p4 <- ggplot(oasis_df,
             aes(x = pipeline, y = mean_n_suprathreshold, fill = cope_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.5) +
  geom_hline(yintercept = oasis_mean, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(
    name   = "Contrast",
    values = c("neg>neu" = "#804080", "pos>neu" = "#C040c0")
  ) +
  labs(x = "Pipeline", y = "Number of voxels", title = "(D) OASIS") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold")
  )

# ── Panel ─────────────────────────────────────────────────────────────────────
panel_plot <- (p1 | p2) / (p3 | p4)
panel_plot


raw %>%
  filter(pipeline == "6HMP", task == "gonogo", cope == "cope3") %>%
  arrange(subject, session) %>%
  print(n = Inf)


# 1. What's actually in the raw file for each pipeline, this task/cope
raw %>%
  filter(pipeline %in% c("6HMP", "24HMP"), task == "gonogo", cope == "cope3") %>%
  arrange(pipeline, subject, session) %>%
  print(n = Inf)

# 2. Same, after session-averaging into per_subject
per_subject %>%
  filter(pipeline %in% c("6HMP", "24HMP"), task == "gonogo", cope == "cope3") %>%
  arrange(pipeline, subject) %>%
  print(n = Inf)

# 3. Final pipeline-level means going into the plot
summary_data %>%
  filter(pipeline %in% c("6HMP", "24HMP"), task == "gonogo", cope == "cope3")

raw %>%
  filter(pipeline == "6HMP", task == "gonogo", cope == "cope3") %>%
  arrange(subject, session) %>%
  print(n = Inf)


library(dplyr)
library(readr)

raw <- read_csv("n_suprathreshold.csv")

# ── Step 1: average sessions within subject ────────────────────────────────
per_subject <- raw %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope3") |
      (task == "oasis" & cope %in% c("cope4", "cope5"))
  ) %>%
  group_by(subject, task, pipeline, cope) %>%
  summarise(n_suprathreshold = mean(n_suprathreshold, na.rm = TRUE), .groups = "drop")

# ── Step 2: average across subjects, per pipeline ──────────────────────────
calculated_averages <- per_subject %>%
  group_by(task, pipeline, cope) %>%
  summarise(
    n_subjects            = n(),
    mean_n_suprathreshold = mean(n_suprathreshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(task, cope, pipeline)

print(calculated_averages, n = Inf)


# Number of suprathreshold voxels per pipeline, four-panel figure
# (A) Go/No-Go, (B) MID, (C) N-back, (D) OASIS (pos>neu / neg>neu)
# Recovered from chat of 6 March 2026. Input produced by collate_mean_zstat_voxelwise.py
# Expected HPC location of input:
# /rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/
#   mean_zstat_voxelwise_summaries/all_mean_zstat_voxelwise.csv

# ── Adjustable font sizes ──────────────────────────
base_font_size   <- 12   # axis titles, legend
title_font_size  <- 13   # bold panel titles
axis_text_size   <- 10   # x and y tick labels
# ──────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

voxel_data <- read_csv("~/Desktop/all_mean_zstat_voxelwise_HPC.csv") %>%
  mutate(pipeline = gsub("_", "-", pipeline))

# Preserve table order
pipeline_order <- voxel_data %>%
  distinct(pipeline) %>%
  pull(pipeline)

make_single_task_plot <- function(data, title_text, bar_colour) {
  mean_val <- mean(data$n_voxels, na.rm = TRUE)
  
  ggplot(data, aes(x = factor(pipeline, levels = pipeline_order),
                   y = n_voxels)) +
    geom_col(fill = bar_colour, width = 0.7) +
    geom_hline(yintercept = mean_val, linetype = "dashed", linewidth = 0.7) +
    labs(x = "Pipeline", y = "Number of voxels", title = title_text) +
    theme_minimal(base_size = base_font_size) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y     = element_text(size = axis_text_size),
      plot.title      = element_text(face = "bold", size = title_font_size)
    )
}

gonogo_df <- voxel_data %>% filter(task == "gonogo", cope == "cope3")
mid_df    <- voxel_data %>% filter(task == "mid",    cope == "cope3")
nback_df  <- voxel_data %>% filter(task == "nback",  cope == "cope3")

oasis_df <- voxel_data %>%
  filter(task == "oasis", cope %in% c("cope4", "cope5")) %>%
  mutate(
    cope_label = case_when(
      cope == "cope4" ~ "pos>neu",
      cope == "cope5" ~ "neg>neu"
    ),
    pipeline = factor(pipeline, levels = pipeline_order)
  )

p1 <- make_single_task_plot(gonogo_df, "(A) Go/No-Go", "aquamarine4")
p2 <- make_single_task_plot(mid_df,    "(B) MID",      "darkslateblue")
p3 <- make_single_task_plot(nback_df,  "(C) N-back",   "darkgreen")

mean_oasis <- mean(oasis_df$n_voxels, na.rm = TRUE)

p4 <- ggplot(oasis_df, aes(x = pipeline, y = n_voxels, fill = cope_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.5) +
  geom_hline(yintercept = mean_oasis, linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(
    name = "Contrast",
    values = c("pos>neu" = "#C040c0", "neg>neu" = "#804080")
  ) +
  labs(x = "Pipeline", y = "Number of voxels", title = "(D) OASIS") +
  theme_minimal(base_size = base_font_size) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    plot.title  = element_text(face = "bold", size = title_font_size)
  )

panel_plot <- (p1 | p2) / (p3 | p4)
panel_plot

# ggsave("n_voxels_by_pipeline.png", panel_plot, width = 12, height = 9, dpi = 300)




# Number of suprathreshold voxels per pipeline, four-panel figure
# Left column: (A) N-back, (B) MID; right column: (C) Go/No-Go, (D) OASIS
# Input: all_mean_zstat_voxelwise.csv (whole-brain version, dated 2 Mar 2026) copied from
# /rds/general/project/task_reliability_study/live/MultiverseReliability/group_levels/
#   mean_zstat_voxelwise_summaries/all_mean_zstat_voxelwise.csv

# ── Adjustable font sizes ──────────────────────────
base_font_size   <- 12   # axis titles, legend
title_font_size  <- 13   # bold panel titles
axis_text_size   <- 10   # x and y tick labels
# ──────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

voxel_data <- read_csv("~/Desktop/all_mean_zstat_voxelwise_HPC.csv") %>%
  mutate(pipeline = gsub("_", "-", pipeline))

# Fixed pipeline order (underscores converted to match the data)
pipeline_order <- gsub("_", "-", c(
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
))

make_single_task_plot <- function(data, title_text, bar_colour) {
  mean_val <- mean(data$n_voxels, na.rm = TRUE)
  
  ggplot(data, aes(x = factor(pipeline, levels = pipeline_order),
                   y = n_voxels)) +
    geom_col(fill = bar_colour, width = 0.7) +
    geom_hline(yintercept = mean_val, linetype = "dashed", linewidth = 0.7) +
    scale_x_discrete(drop = FALSE) +
    labs(x = "Pipeline", y = "Number of voxels", title = title_text) +
    theme_minimal(base_size = base_font_size) +
    theme(
      panel.grid      = element_blank(),
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y     = element_text(size = axis_text_size),
      plot.title      = element_text(face = "bold", size = title_font_size)
    )
}

nback_df  <- voxel_data %>% filter(task == "nback",  cope == "cope3")
mid_df    <- voxel_data %>% filter(task == "mid",    cope == "cope3")
gonogo_df <- voxel_data %>% filter(task == "gonogo", cope == "cope3")

oasis_df <- voxel_data %>%
  filter(task == "oasis", cope %in% c("cope4", "cope5")) %>%
  mutate(
    cope_label = case_when(
      cope == "cope4" ~ "pos>neu",
      cope == "cope5" ~ "neg>neu"
    ),
    pipeline = factor(pipeline, levels = pipeline_order)
  )

p_nback  <- make_single_task_plot(nback_df,  "(A) N-back",   "darkgreen")
p_mid    <- make_single_task_plot(mid_df,    "(B) MID",      "darkslateblue")
p_gonogo <- make_single_task_plot(gonogo_df, "(C) Go/No-Go", "aquamarine4")

mean_oasis <- mean(oasis_df$n_voxels, na.rm = TRUE)

p_oasis <- ggplot(oasis_df, aes(x = pipeline, y = n_voxels, fill = cope_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.5) +
  geom_hline(yintercept = mean_oasis, linetype = "dashed", linewidth = 0.7) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    name = "Contrast",
    values = c("pos>neu" = "#C040c0", "neg>neu" = "#804080")
  ) +
  labs(x = "Pipeline", y = "Number of voxels", title = "(D) OASIS") +
  theme_minimal(base_size = base_font_size) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_size),
    axis.text.y = element_text(size = axis_text_size),
    plot.title  = element_text(face = "bold", size = title_font_size)
  )

# Left column: N-back over MID; right column: Go/No-Go over OASIS
panel_plot <- (p_nback | p_gonogo) / (p_mid | p_oasis)
panel_plot

# ggsave("~/Desktop/n_voxels_by_pipeline.png", panel_plot, width = 12, height = 9, dpi = 300)
