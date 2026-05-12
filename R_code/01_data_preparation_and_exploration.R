################################################################################
# Probabilistic Index Models for Longitudinal Data: Small-Sample Inference
# Paper Analysis: Data Preparation and Exploratory Analysis
# 
# Dataset: NIMH Schizophrenia Collaborative Study
# - Two treatment groups (perphenazine vs risperidone)
# - 6 weeks follow-up with 5 measurement timepoints
# - BPRS (Brief Psychiatric Rating Scale) - ordinal outcome
# - N = ~40 subjects
################################################################################

# Clear workspace
rm(list = ls())

# Set seed for reproducibility
set.seed(2026)

################################################################################
# 1. LOAD REQUIRED PACKAGES
################################################################################

packages <- c("tidyverse", "lme4", "ordinal", "nlme", "ggplot2", "gridExtra", 
              "knitr", "kableExtra", "DescTools", "multcomp", "mlmRev")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

################################################################################
# 2. LOAD DATA: NIMH SCHIZOPHRENIA STUDY
################################################################################

# Load the Schizophrenia dataset from mlmRev package
data(Schizophrenia, package = "mlmRev")

# Create working dataset
dat_full <- Schizophrenia %>%
  as_tibble() %>%
  rename(
    Subject = subject,
    Week = week,
    Treatment = treatment,
    BPRS = outcome  # Brief Psychiatric Rating Scale (ordinal outcome)
  ) %>%
  mutate(
    Subject = as.factor(Subject),
    Treatment = factor(Treatment, levels = c("0", "1"), 
                       labels = c("Perphenazine", "Risperidone")),
    Week = as.numeric(as.character(Week))
  ) %>%
  arrange(Subject, Week)

# Quick check
head(dat_full, 12)
tail(dat_full, 12)

# Data dimensions
cat("Full Dataset Dimensions:\n")
cat("  Total observations:", nrow(dat_full), "\n")
cat("  Number of subjects:", n_distinct(dat_full$Subject), "\n")
cat("  Number of timepoints:", n_distinct(dat_full$Week), "\n")
cat("  Treatment groups:", n_distinct(dat_full$Treatment), "\n")

# Summary by treatment group
summary_by_group <- dat_full %>%
  group_by(Treatment) %>%
  summarise(
    N_subjects = n_distinct(Subject),
    N_obs = n(),
    BPRS_mean = mean(BPRS, na.rm = TRUE),
    BPRS_sd = sd(BPRS, na.rm = TRUE),
    BPRS_min = min(BPRS, na.rm = TRUE),
    BPRS_max = max(BPRS, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n\nSummary by Treatment Group:\n")
print(kable(summary_by_group, format = "rst"))

################################################################################
# 3. CREATE ORDINAL OUTCOME (CATEGORIZE CONTINUOUS BPRS)
################################################################################

# BPRS ranges from 18-126. Discretize into ordinal categories:
# Low (18-40), Moderate (41-75), High (76-100), Very High (101-126)
# This creates a 4-level ordinal outcome similar to typical rating scales

dat_full <- dat_full %>%
  mutate(
    BPRS_cat = cut(BPRS, 
                   breaks = c(0, 40, 75, 100, 126),
                   labels = c("Low", "Moderate", "High", "Very High"),
                   ordered = TRUE)
  )

# Check distribution of ordinal outcome
cat("\n\nDistribution of Ordinal BPRS Categories:\n")
print(table(dat_full$BPRS_cat, useNA = "ifany"))

cat("\nBy Treatment Group:\n")
print(table(dat_full$Treatment, dat_full$BPRS_cat, useNA = "ifany"))

################################################################################
# 4. CREATE REDUCED-SAMPLE SUBSET (20 per group)
################################################################################

# For demonstration of small-sample applicability
set.seed(2026)

dat_reduced <- dat_full %>%
  group_by(Treatment) %>%
  sample_n(size = 20, replace = FALSE) %>%
  ungroup() %>%
  arrange(Subject, Week)

cat("\n\nReduced Dataset Dimensions:\n")
cat("  Total observations:", nrow(dat_reduced), "\n")
cat("  Number of subjects:", n_distinct(dat_reduced$Subject), "\n")
cat("  Subjects per treatment:", 
    dat_reduced %>% group_by(Treatment) %>% 
    summarise(n = n_distinct(Subject)) %>% pull(n) %>% paste(collapse = " vs "), "\n")

################################################################################
# 5. EXPLORATORY DATA ANALYSIS - FULL SAMPLE
################################################################################

# Figure 1: Spaghetti plot (continuous BPRS by timepoint)
fig1 <- ggplot(dat_full, aes(x = Week, y = BPRS, group = Subject, color = Treatment)) +
  geom_line(alpha = 0.4, size = 0.8) +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~Treatment) +
  theme_minimal() +
  labs(title = "Individual Trajectories of BPRS Over Time (Full Sample)",
       x = "Week", y = "BPRS Score",
       color = "Treatment") +
  theme(legend.position = "bottom")

print(fig1)
ggsave("Figures/01_trajectories_full.png", fig1, width = 10, height = 6, dpi = 300)

# Figure 2: Mean trajectories with 95% CI
dat_summary <- dat_full %>%
  group_by(Week, Treatment) %>%
  summarise(
    Mean = mean(BPRS, na.rm = TRUE),
    SD = sd(BPRS, na.rm = TRUE),
    N = n(),
    SE = SD / sqrt(N),
    CI_low = Mean - 1.96 * SE,
    CI_high = Mean + 1.96 * SE,
    .groups = 'drop'
  )

fig2 <- ggplot(dat_summary, aes(x = Week, y = Mean, color = Treatment, fill = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Mean BPRS Trajectories with 95% CI (Full Sample)",
       x = "Week", y = "Mean BPRS Score",
       color = "Treatment", fill = "Treatment") +
  theme(legend.position = "bottom")

print(fig2)
ggsave("Figures/02_mean_trajectories_full.png", fig2, width = 10, height = 6, dpi = 300)

# Figure 3: Distribution of ordinal categories by timepoint
fig3_data <- dat_full %>%
  group_by(Week, Treatment, BPRS_cat) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  complete(Week, Treatment, BPRS_cat, fill = list(Count = 0)) %>%
  group_by(Week, Treatment) %>%
  mutate(Prop = Count / sum(Count))

fig3 <- ggplot(fig3_data, aes(x = factor(Week), fill = BPRS_cat)) +
  geom_bar(aes(y = Prop), stat = "identity") +
  facet_wrap(~Treatment) +
  theme_minimal() +
  labs(title = "Distribution of Ordinal BPRS Categories Over Time",
       x = "Week", y = "Proportion",
       fill = "BPRS Category") +
  scale_fill_brewer(palette = "RdYlGn", direction = -1) +
  theme(legend.position = "bottom")

print(fig3)
ggsave("Figures/03_ordinal_distribution.png", fig3, width = 10, height = 6, dpi = 300)

################################################################################
# 6. EXPLORATORY DATA ANALYSIS - REDUCED SAMPLE
################################################################################

# Figure 4: Trajectories for reduced sample
fig4 <- ggplot(dat_reduced, aes(x = Week, y = BPRS, group = Subject, color = Treatment)) +
  geom_line(alpha = 0.5, size = 0.8) +
  geom_point(alpha = 0.6, size = 2) +
  facet_wrap(~Treatment) +
  theme_minimal() +
  labs(title = "Individual Trajectories: Reduced Sample (n=20 per group)",
       x = "Week", y = "BPRS Score",
       color = "Treatment") +
  theme(legend.position = "bottom")

print(fig4)
ggsave("Figures/04_trajectories_reduced.png", fig4, width = 10, height = 6, dpi = 300)

# Figure 5: Box plots by timepoint and treatment
fig5 <- ggplot(dat_full, aes(x = factor(Week), y = BPRS, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.size = 2) +
  theme_minimal() +
  labs(title = "BPRS Distribution by Week and Treatment (Full Sample)",
       x = "Week", y = "BPRS Score",
       fill = "Treatment") +
  theme(legend.position = "bottom")

print(fig5)
ggsave("Figures/05_boxplots.png", fig5, width = 10, height = 6, dpi = 300)

################################################################################
# 7. DESCRIPTIVE STATISTICS TABLE
################################################################################

# Baseline characteristics (Week 0)
baseline <- dat_full %>%
  filter(Week == 0) %>%
  select(Subject, Treatment, BPRS) %>%
  group_by(Treatment) %>%
  summarise(
    N = n(),
    Mean_BPRS = mean(BPRS),
    SD_BPRS = sd(BPRS),
    Min_BPRS = min(BPRS),
    Max_BPRS = max(BPRS),
    .groups = 'drop'
  )

cat("\n\nBaseline Characteristics (Week 0):\n")
print(kable(baseline, digits = 2, format = "rst"))

# Change from baseline to last timepoint
change_score <- dat_full %>%
  select(Subject, Treatment, Week, BPRS) %>%
  pivot_wider(names_from = Week, values_from = BPRS, names_prefix = "Week_") %>%
  mutate(Change = Week_6 - Week_0) %>%
  group_by(Treatment) %>%
  summarise(
    N = n(),
    Mean_Change = mean(Change, na.rm = TRUE),
    SD_Change = sd(Change, na.rm = TRUE),
    Min_Change = min(Change, na.rm = TRUE),
    Max_Change = max(Change, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n\nChange from Baseline to Week 6:\n")
print(kable(change_score, digits = 2, format = "rst"))

################################################################################
# 8. SAVE PREPARED DATA
################################################################################

# Save full sample dataset
write.csv(dat_full, "data/NIMH_Schizophrenia_full_sample.csv", row.names = FALSE)

# Save reduced sample dataset
write.csv(dat_reduced, "data/NIMH_Schizophrenia_reduced_sample.csv", row.names = FALSE)

cat("\n\nData saved to:\n")
cat("  Full sample:    data/NIMH_Schizophrenia_full_sample.csv\n")
cat("  Reduced sample: data/NIMH_Schizophrenia_reduced_sample.csv\n")

################################################################################
# 9. SESSION INFO
################################################################################

cat("\n\nSession Information:\n")
print(sessionInfo())

################################################################################
# END OF SCRIPT
################################################################################
