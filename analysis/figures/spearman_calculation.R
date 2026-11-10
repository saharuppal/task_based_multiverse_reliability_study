# Install packages (run once)
install.packages(c("dplyr", "tidyr", "flextable", "readr"))

library(dplyr)
library(tidyr)
library(flextable)
library(readr)

# Read CSV
all_spearman_values <- read_csv("all_spearman_values.csv")

# Filter required copes
filtered_data <- all_spearman_values %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope3") |
      (task == "oasis" & cope %in% c("cope4", "cope5"))
  ) %>%
  mutate(
    task_cope = paste(task, cope, sep = "_"),
    pipeline = gsub("_", "-", pipeline)
  )

# Reshape
table_data <- filtered_data %>%
  select(pipeline, task_cope, spearman) %>%
  pivot_wider(names_from = task_cope, values_from = spearman)

# Rename safely (by column name)
table_data <- table_data %>%
  rename(
    Pipeline = pipeline,
    GoNoGo = gonogo_cope3,
    MID = mid_cope3,
    `N-Back` = nback_cope3,
    `OASIS (pos>neu)` = oasis_cope4,
    `OASIS (neg>neu)` = oasis_cope5
  )

# Optional: round values
numeric_cols <- setdiff(colnames(table_data), "Pipeline")
table_data[numeric_cols] <- round(table_data[numeric_cols], 3)

# Create flextable
ft <- flextable(table_data) %>%
  autofit() %>%
  theme_booktabs()

ft <- align(ft,  j = cols, align = "center", part = "all")
ft <- align(ft,  j = "Pipeline", align = "left", part = "all")
ft <- valign(ft, j = cols, valign = "center")

# ---- Bold highest per column EXCEPT GoNoGo ----

cols_to_bold <- numeric_cols

for (col in cols_to_bold) {
  max_val <- max(table_data[[col]], na.rm = TRUE)
  
  ft <- bold(
    ft,
    i = table_data[[col]] == max_val,
    j = col,
    bold = TRUE
  )
}

ft




# ── Adjustable font sizes ──────────────────────────
body_font_size   <- 16    # cell text (Spearman values)
header_font_size <- 18   # column headers
# ──────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(flextable)
library(readr)

all_spearman_values <- read_csv("all_spearman_values.csv")

filtered_data <- all_spearman_values %>%
  filter(
    (task %in% c("nback", "gonogo", "mid") & cope == "cope3") |
      (task == "oasis" & cope %in% c("cope4", "cope5"))
  ) %>%
  mutate(
    task_cope = paste(task, cope, sep = "_"),
    pipeline = gsub("_", "-", pipeline)
  )

table_data <- filtered_data %>%
  select(pipeline, task_cope, spearman) %>%
  pivot_wider(names_from = task_cope, values_from = spearman)

table_data <- table_data %>%
  rename(
    Pipeline = pipeline,
    GoNoGo = gonogo_cope3,
    MID = mid_cope3,
    "N\u2011Back" = nback_cope3,
    "OASIS (pos>neu)" = oasis_cope4,
    "OASIS (neg>neu)" = oasis_cope5
  )

numeric_cols <- setdiff(colnames(table_data), "Pipeline")
table_data[numeric_cols] <- round(table_data[numeric_cols], 3)

cols <- numeric_cols   # fix: cols was undefined in original

ft <- flextable(table_data) %>%
  autofit() %>%
  theme_booktabs()

ft <- align(ft,  j = cols, align = "center", part = "all")
ft <- align(ft,  j = "Pipeline", align = "left", part = "all")
ft <- valign(ft, j = cols, valign = "center")
ft <- fontsize(ft, size = body_font_size,   part = "body")
ft <- fontsize(ft, size = header_font_size, part = "header")

for (col in numeric_cols) {
  max_val <- max(table_data[[col]], na.rm = TRUE)
  ft <- bold(ft, i = table_data[[col]] == max_val, j = col, bold = TRUE)
}

ft

