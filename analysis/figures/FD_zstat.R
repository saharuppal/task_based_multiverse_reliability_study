# =============================================================================
# Head-motion (FD) vs activation outcome analysis
# Spearman correlations per pipeline + linear mixed models
# Sessions analysed separately throughout
# =============================================================================



library(tidyverse)
library(lme4)
library(lmerTest)   # p-values for lmer via Satterthwaite
library(broom.mixed)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Load and merge
# -----------------------------------------------------------------------------

fd        <- read_csv("mean_fd.csv")           # subject, session, task, mean_fd
mean_z    <- read_csv("mean_z.csv")            # subject, session, task, pipeline, cope, mean_z
top1      <- read_csv("mean_top1pct_z.csv")      # subject, session, task, pipeline, cope, mean_top1pct_z
n_supra   <- read_csv("n_suprathreshold.csv")  # subject, session, task, pipeline, cope, n_suprathreshold

# FD is per subject/session/task — join to the outcome frames
df <- mean_z |>
  left_join(top1,    by = c("subject","session","task","pipeline","cope")) |>
  left_join(n_supra, by = c("subject","session","task","pipeline","cope")) |>
  left_join(fd,      by = c("subject","session","task"))

# Enforce factor types
df <- df |>
  mutate(
    pipeline = factor(pipeline),
    cope     = factor(cope),
    task     = factor(task),
    session  = factor(session),
    subject  = factor(subject)
  ) |>
  mutate(pipeline = relevel(pipeline, ref = "6HMP"))

# Convenience: list of outcome variables
outcomes <- c("mean_z", "mean_top1pct_z", "n_suprathreshold")

# Split once — reuse throughout
df_s1 <- filter(df, session == "ses-1")
df_s2 <- filter(df, session == "ses-2")

# -----------------------------------------------------------------------------
# 2. Spearman correlations: FD vs each outcome, per pipeline
#    Run within each task × cope combination, then pool or inspect separately
# -----------------------------------------------------------------------------

run_spearman <- function(data_session) {
  data_session |>
    group_by(task, cope, pipeline) |>
    summarise(
      rho_mean_z         = ifelse(sum(!is.na(mean_z) & !is.na(mean_fd)) >= 5,
                                  cor.test(mean_z, mean_fd, method = "spearman", exact = FALSE)$estimate, NA),
      p_mean_z           = ifelse(sum(!is.na(mean_z) & !is.na(mean_fd)) >= 5,
                                  cor.test(mean_z, mean_fd, method = "spearman", exact = FALSE)$p.value, NA),
      rho_mean_top1pct_z = ifelse(sum(!is.na(mean_top1pct_z) & !is.na(mean_fd)) >= 5,
                                  cor.test(mean_top1pct_z, mean_fd, method = "spearman", exact = FALSE)$estimate, NA),
      p_mean_top1pct_z   = ifelse(sum(!is.na(mean_top1pct_z) & !is.na(mean_fd)) >= 5,
                                  cor.test(mean_top1pct_z, mean_fd, method = "spearman", exact = FALSE)$p.value, NA),
      rho_n_suprathreshold = ifelse(sum(!is.na(n_suprathreshold) & !is.na(mean_fd)) >= 5,
                                    cor.test(n_suprathreshold, mean_fd, method = "spearman", exact = FALSE)$estimate, NA),
      p_n_suprathreshold   = ifelse(sum(!is.na(n_suprathreshold) & !is.na(mean_fd)) >= 5,
                                    cor.test(n_suprathreshold, mean_fd, method = "spearman", exact = FALSE)$p.value, NA),
      .groups = "drop"
    ) |>
    pivot_longer(
      cols      = starts_with("rho_") | starts_with("p_"),
      names_to  = c(".value", "outcome"),
      names_pattern = "(rho|p)_(.*)"
    )
}

# NOTE: if your tidyverse version doesn't support tibble-valued across(),
# replace the across() block with explicit calls, e.g.:
#   rho_mean_z = cor.test(mean_z, mean_fd, method="spearman", exact=FALSE)$estimate,
#   p_mean_z   = cor.test(mean_z, mean_fd, method="spearman", exact=FALSE)$p.value,
# etc. for each outcome.

spearman_s1 <- run_spearman(df_s1) |> mutate(session = 1)
spearman_s2 <- run_spearman(df_s2) |> mutate(session = 2)

spearman_all <- bind_rows(spearman_s1, spearman_s2) |>
  mutate(sig = p < .05)

spearman_all <- bind_rows(spearman_s1, spearman_s2) |>
  mutate(
    p_fdr = p.adjust(p, method = "fdr"),   # BH correction across all rows
    sig     = p     < .05,
    sig_fdr = p_fdr < .05
  )


print(spearman_all, n = 30)

spearman_all |> count(session)
# Should give roughly equal rows for session = 1 and session = 2

# -----------------------------------------------------------------------------
# 3. Visualise Spearman rho: heatmap (pipelines × outcomes), faceted by task/cope/session
# -----------------------------------------------------------------------------

p_heatmap <- spearman_all |>
  filter(outcome == "mean_z") |>
  filter(!(task == "oasis" & cope == "cope3"),
         !(task != "oasis" & cope != "cope3")) |>
  mutate(
    session = factor(session, labels = c("Session 1", "Session 2")),
    cope_label = case_when(
      task == "oasis" & cope == "cope4" ~ "OASIS (pos>neu)",
      task == "oasis" & cope == "cope5" ~ "OASIS (neg>neu)",
      task == "nback"  ~ "N-Back",
      task == "gonogo" ~ "Go-NoGo",
      task == "mid"    ~ "MID"
    ),
    cope_label = factor(cope_label, levels = c("N-Back", "Go-NoGo", "MID",
                                               "OASIS (pos>neu)", "OASIS (neg>neu)")),
    pipeline = gsub("_", "-", as.character(pipeline)),
    pipeline = factor(pipeline, levels = c(
      "6HMP",
      "24HMP",
      "6HMP-WM-CSF",
      "24HMP-WM-CSF",
      "6HMP-aCompCor",
      "24HMP-aCompCor",
      "6HMP-spike-reg",
      "24HMP-spike-reg",
      "6HMP-aCompCor-spike-reg",
      "24HMP-aCompCor-spike-reg",
      "6HMP-WM-CSF-spike-reg",
      "24HMP-WM-CSF-spike-reg"
    ))
  ) |> 
  ggplot(aes(x = cope_label, y = pipeline, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(sig_fdr, sprintf("%.2f*", rho), sprintf("%.2f", rho))),
            size = 2.5) +
  scale_fill_gradient2(low = "#3B6EB5", mid = "white", high = "#C84B31",
                       midpoint = 0, limits = c(-1, 1), name = "ρ") +
  scale_y_discrete(limits = rev) +
  facet_grid(session ~ .) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text  = element_text(size = 9)
  )
print(p_heatmap)

ggsave("spearman_heatmap.pdf", p_heatmap, width = 30, height = 10)


levels(df$task)
levels(df$cope)

# -----------------------------------------------------------------------------
# 4. Linear mixed models: outcome ~ FD * pipeline + (1|subject)
#    One model per session × task × cope × outcome
# -----------------------------------------------------------------------------

run_lmm <- function(data_session, outcome_var) {
  data_session |>
    group_by(task, cope) |>
    group_modify(~ {
      dat <- .x
      # Scale predictors for interpretability
      dat <- dat |>
        mutate(
          fd_z  = scale(mean_fd)[,1],
          outcome_val = .data[[outcome_var]]
        )
      # Fit model; catch convergence warnings
      m <- tryCatch(
        lmer(outcome_val ~ fd_z * pipeline + (1 | subject), data = dat,
             REML = FALSE,
             control = lmerControl(optimizer = "bobyqa")),
        error = function(e) NULL
      )
      if (is.null(m)) return(tibble(warning = "model failed"))
      
      # Fixed-effect table
      tidy(m, effects = "fixed") |>
        mutate(outcome = outcome_var)
    }) |>
    ungroup()
}

# Run for each outcome, each session
lmm_results <- map_dfr(
  outcomes,
  function(ov) {
    bind_rows(
      run_lmm(df_s1, ov) |> mutate(session = 1),
      run_lmm(df_s2, ov) |> mutate(session = 2)
    )
  }
)

# Tidy output with FDR correction within each task × cope × outcome × session block
lmm_results <- lmm_results |>
  group_by(task, cope, outcome, session) |>
  mutate(
    p_fdr   = p.adjust(p.value, method = "fdr"),
    sig_fdr = p_fdr < .05
  ) |>
  ungroup()

print(lmm_results |> filter(p_fdr < .05), n = 50)
View(lmm_results)
write_csv(lmm_results, "lmm_results.csv")


lmm_results |>
  filter(sig_fdr == TRUE, str_detect(term, "fd_z:pipeline")) |>
  select(session, task, cope, outcome, term, estimate, std.error, p_fdr) |>
  arrange(session, task, cope, outcome, estimate) |>
  View()

# -----------------------------------------------------------------------------
# 5. Extract and plot the FD × pipeline interaction terms specifically
# -----------------------------------------------------------------------------

interaction_terms <- lmm_results |>
  filter(str_detect(term, "fd_z:pipeline"))

p_interaction <- interaction_terms |>
  mutate(pipeline_name = str_remove(term, "fd_z:pipeline")) |>
  ggplot(aes(x = reorder(pipeline_name, estimate), y = estimate,
             colour = outcome, shape = factor(session))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(ymin = estimate - 1.96 * std.error,
                      ymax = estimate + 1.96 * std.error),
                  position = position_dodge(width = 0.5)) +
  facet_grid(outcome ~ task + cope, scales = "free_y", labeller = label_both) +
  coord_flip() +
  labs(
    title    = "FD × pipeline interaction coefficients",
    subtitle = "Negative = pipeline dampens motion-activation coupling relative to reference pipeline",
    x = "Pipeline (vs reference)", y = "Interaction estimate (β)",
    colour = "Outcome", shape = "Session"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 7))

print(p_interaction)
ggsave("lmm_interaction_plot.pdf", p_interaction, width = 16, height = 10)

# -----------------------------------------------------------------------------
# 6. Scatterplots: FD vs outcome, coloured by pipeline, per task/cope
#    (pick one outcome at a time to keep readable — here mean_z)
# -----------------------------------------------------------------------------

plot_scatter <- function(data, outcome_var, ses_label) {
  data |>
    ggplot(aes(x = mean_fd, y = .data[[outcome_var]], colour = pipeline)) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet_grid(task ~ cope, labeller = label_both) +
    labs(
      title    = glue::glue("{outcome_var} vs mean FD — session {ses_label}"),
      subtitle = "Flatter slopes indicate better motion mitigation by pipeline",
      x = "Mean framewise displacement (mm)", y = outcome_var,
      colour = "Pipeline"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
}

p_scatter_s1 <- plot_scatter(df_s1, "mean_z", 1)
p_scatter_s2 <- plot_scatter(df_s2, "mean_z", 2)

print(p_scatter_s1)
print
ggsave("scatter_mean_z_s1.pdf", p_scatter_s1, width = 12, height = 10)
ggsave("scatter_mean_z_s2.pdf", p_scatter_s2, width = 12, height = 10)

# -----------------------------------------------------------------------------
# 7. Session consistency check
#    Do the Spearman rhos replicate across sessions?
# -----------------------------------------------------------------------------

spearman_wide <- spearman_all |>
  select(task, cope, pipeline, outcome, session, rho) |>
  pivot_wider(names_from = session, values_from = rho,
              names_prefix = "session_")

p_session_consistency <- spearman_wide |>
  ggplot(aes(x = session_1, y = session_2, colour = pipeline)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.6) +
  facet_grid(outcome ~ task + cope, labeller = label_both) +
  labs(
    title    = "Session 1 vs session 2 Spearman ρ (FD–outcome)",
    subtitle = "Points on diagonal = perfect test-retest consistency",
    x = "Session 1 ρ", y = "Session 2 ρ", colour = "Pipeline"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(size = 7))

print(p_session_consistency)
ggsave("session_consistency.pdf", p_session_consistency, width = 14, height = 8)

# =============================================================================
# END
# =============================================================================
library(tidyverse)   # dplyr, purrr, stringr, readr
library(lme4)        # lmer
library(broom.mixed)
run_lmm <- function(data_session, outcome_var) {
  data_session |>
    group_by(task, cope) |>
    group_modify(~ {
      dat <- .x
      dat <- dat |>
        mutate(
          fd_z        = scale(mean_fd)[,1],
          outcome_val = .data[[outcome_var]]
        )
      m <- tryCatch(
        lmer(outcome_val ~ fd_z * pipeline + (1 | subject), data = dat,
             REML = FALSE,
             control = lmerControl(optimizer = "bobyqa")),
        error = function(e) NULL
      )
      if (is.null(m)) return(tibble(warning = "model failed"))
      
      tidy(m, effects = "fixed") |>
        mutate(outcome = outcome_var)
    }) |>
    ungroup()
}

# Run for mean_zstat only
lmm_results <- bind_rows(
  run_lmm(df_s1, "mean_z") |> mutate(session = 1),
  run_lmm(df_s2, "mean_z") |> mutate(session = 2)
)

# FDR correction
lmm_results <- lmm_results |>
  group_by(task, cope, outcome, session) |>
  mutate(
    p_fdr   = p.adjust(p.value, method = "fdr"),
    sig_fdr = p_fdr < .05
  ) |>
  ungroup()

print(lmm_results |> filter(p_fdr < .05), n = 50)
View (lmm_results |> filter(p_fdr < .05))
write_csv(lmm_results, "lmm_results.csv")

lmm_results |>
  filter(sig_fdr == TRUE, str_detect(term, "fd_z:pipeline")) |>
  select(session, task, cope, outcome, term, estimate, std.error, p_fdr) |>
  arrange(session, task, cope, outcome, estimate) |>
  View()

lmm_results |>
    filter(p_fdr < .05) |>
    count(task)
names(df_s1)
unique(df_s1$task)
unique(df_s2$task)
df_s1 |> count(task, cope)
lmm_results |> count(task)



fd_results <- lmm_results |>
  filter(str_detect(term, "fd_z"), sig_fdr == TRUE) |>
  select(session, task, cope, term, estimate, std.error, statistic, p.value, p_fdr) |>
  arrange(session, task, cope, term)

View(fd_results)
write_csv(fd_results, "fd_z_results.csv")
