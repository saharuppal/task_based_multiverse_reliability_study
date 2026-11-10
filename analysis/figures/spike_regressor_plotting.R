library(dplyr)
library(flextable)
library(officer)

# ── Read in the data ──────────────────────────────────────────────────────
df <- read.csv("final_spike_regressor_counts.csv", stringsAsFactors = FALSE)

# Quick check of structure
str(df)
# Expect columns like: subject, session, task, n_motion_outliers, file_found

# ── Optional: drop the file_found column and any missing rows for the tables ─
df_clean <- df %>%
  filter(file_found == "yes") %>%
  select(subject, session, task, n_motion_outliers)

# ── Split into one dataframe per task ────────────────────────────────────
tasks <- c("gonogo", "oasis", "nback", "mid")

task_tables <- lapply(tasks, function(t) {
  df_clean %>%
    filter(task == t) %>%
    select(-task) %>%
    arrange(subject, session)
})
names(task_tables) <- tasks

# ── Build a flextable for one task's data ────────────────────────────────
make_flextable <- function(data, title_text) {
  ft <- flextable(data) %>%
    set_header_labels(
      subject = "Subject",
      session = "Session",
      n_motion_outliers = "N Motion Outliers"
    ) %>%
    theme_vanilla() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all")
  return(ft)
}

# ── Create the Word document ─────────────────────────────────────────────
doc <- read_docx()

task_titles <- c(
  gonogo = "Go/No-Go",
  oasis  = "OASIS",
  nback  = "N-back",
  mid    = "MID"
)

for (t in tasks) {
  doc <- doc %>%
    body_add_par(paste("Table:", task_titles[[t]], "- Spike Regressor Counts"),
                 style = "heading 1") %>%
    body_add_flextable(make_flextable(task_tables[[t]])) %>%
    body_add_par("", style = "Normal")  # spacer
  
  # Add a page break between tables except after the last one
  if (t != tasks[length(tasks)]) {
    doc <- doc %>% body_add_break()
  }
}

print(doc, target = "spike_regressor_tables_by_task.docx")