library(tidyverse)
fd <- read_csv ("mean_fd.csv")

fd_plot <- fd %>%
  mutate(
    task    = recode(task,
                     "nback"  = "N-Back",
                     "mid"    = "MID",
                     "gonogo" = "Go-NoGo",
                     "oasis"  = "OASIS"),
    task    = factor(task, levels = c("N-Back", "MID", "Go-NoGo", "OASIS")),
    session = recode(as.character(session_id),
                     "1" = "Session 1",
                     "2" = "Session 2"),
    session = factor(session, levels = c("Session 1", "Session 2"))
  )

p_fd <- ggplot(fd_plot, aes(x = task, y = mean_fd, fill = session, colour = session)) +
  geom_boxplot(alpha = 0.4, width = 0.4, outlier.shape = NA, colour = "black",
               position = position_dodge(0.6)) +
  geom_point(shape = 16, size = 2, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)) +
  scale_fill_manual(values   = c("Session 1" = "#4C72B0", "Session 2" = "#C44E52"),
                    name = "Session") +
  scale_colour_manual(values = c("Session 1" = "#4C72B0", "Session 2" = "#C44E52"),
                      name = "Session") +
  labs(x = NULL, y = "Mean FD (mm)",
       title = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(colour = "black"),
    axis.ticks   = element_line(colour = "black"),
    axis.text.x  = element_text(size = 11),
    strip.text   = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_fd)
ggsave("mean_fd_plot.pdf", p_fd, width = 10, height = 5)
print(p_fd)
ggsave("mean_fd_plot.pdf", p_fd, width = 10, height = 5)



library(dplyr)

# t-tests per task, session 1 vs session 2
t_results <- fd_plot %>%
  group_by(task) %>%
  summarise(
    t_test = list(t.test(mean_fd[session == "Session 1"],
                         mean_fd[session == "Session 2"],
                         paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    t_stat  = sapply(t_test, function(x) round(x$statistic, 3)),
    df      = sapply(t_test, function(x) round(x$parameter, 1)),
    p_raw   = sapply(t_test, function(x) round(x$p.value, 4)),
    p_fdr   = round(p.adjust(p_raw, method = "fdr"), 4)
  ) %>%
  select(task, t_stat, df, p_raw, p_fdr)

print(t_results)


fd_plot %>%
  group_by(task, session) %>%
  summarise(mean = mean(mean_fd), sd = sd(mean_fd), .groups = "drop")

install.packages("effsize")
library(effsize)

effect_sizes <- fd_plot %>%
  group_by(task) %>%
  summarise(
    cohen_d = cohen.d(mean_fd[session == "Session 2"],
                      mean_fd[session == "Session 1"],
                      paired = TRUE)$estimate,
    .groups = "drop"
  )

t_results <- t_results %>%
  left_join(effect_sizes, by = "task")

print(t_results)



library(tidyverse)
library(readr)

fd <- read_csv("mean_fd.csv")
fd_plot <- fd %>%
  mutate(
    task    = recode(task,
                     "nback"  = "N-Back",
                     "mid"    = "MID",
                     "gonogo" = "Go-NoGo",
                     "oasis"  = "OASIS"),
    task    = factor(task, levels = c("N-Back", "MID", "Go-NoGo", "OASIS")),
    session = recode(as.character(session_id),
                     "1" = "Session 1",
                     "2" = "Session 2"),
    session = factor(session, levels = c("Session 1", "Session 2"))
  )

p_fd <- ggplot(fd_plot, aes(x = task, y = mean_fd, fill = session, colour = session)) +
  geom_boxplot(alpha = 0.4, width = 0.4, outlier.shape = NA, colour = "black",
               position = position_dodge(0.6)) +
  geom_point(shape = 16, size = 2.5, alpha = 0.7,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)) +
  scale_fill_manual(values = c("Session 1" = "#4C72B0", "Session 2" = "#C44E52"),
                    name = "Session") +
  scale_colour_manual(values = c("Session 1" = "#4C72B0", "Session 2" = "#C44E52"),
                      name = "Session") +
  scale_y_continuous(
    breaks = seq(0, 0.8, by = 0.1),
    minor_breaks = seq(0, 0.8, by = 0.05),
    limits = c(0, 0.8)
  ) +
  labs(x = NULL, y = "Mean FD (mm)", title = "") +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(colour = "black"),
    axis.ticks        = element_line(colour = "black"),
    axis.ticks.length = unit(4, "pt"),
    axis.text         = element_text(size = 13),
    axis.title.y      = element_text(size = 14),
    axis.text.x       = element_text(size = 13),
    strip.text        = element_text(face = "bold", size = 13),
    legend.position   = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 13)
  )

print(p_fd)
ggsave("mean_fd_plot.pdf", p_fd, width = 10, height = 5)
