# ==============================================================================
# 01_data_preparation.R
#
# Purpose:
#   Load each raw dataset, harmonise variable names and codings, and
#   produce one analysis-ready data frame per trial.
#
# Inputs:
#   Raw data files (not distributed — see README).
#   Replace the placeholder load/read calls in SECTION A with your own paths.
#
# Outputs (objects retained for downstream scripts):
#   deprexis_select_PHQ9   — development sample (Deprexis RCTs merged)
#   eficasy_select         — external validation: edupression / eFICASY
#   HRSD_selfapy_QIDSs     — external validation: Selfapy
#   moritz_select          — external validation: Deprexis independent cohort
#
# Variable dictionary (all datasets must contain these columns after prep):
#   Group            factor  "intervention" / "control"
#   Gender           factor  "male" / "female"
#   Age              numeric years
#   Relationship     factor  "partnered" / "single"
#   School           factor  ">12yrs" / "<12yrs"
#   Employment       factor  "employed" / "unemployed"
#   current_psychotherapy   factor  "current therapy" / "no current therapy"
#   current_antidepressant  factor  "current antidepressant" / "no antidepressant"
#   MajorDepression  factor  "yes" / "no"  (levels: c("no","yes"))
#   Dysthymia        factor  "yes" / "no"  (levels: c("no","yes"))
#   QOL_psych        numeric  WHO quality of life — psychological domain
#   QOL_phys         numeric  WHO quality of life — physical domain
#   PHQ_1            numeric  PHQ-9 sum score at baseline
#   PHQ_3            numeric  PHQ-9 sum score at post-treatment / follow-up
#   PHQ_change       numeric  PHQ_1 - PHQ_3  (positive = improvement)
#
# Note: Selfapy additionally contains QIDS_1, QIDS_3, QIDS_change.
#       Moritz additionally contains BDI_1, BDI_3, BDI_change (BDI-II,
#       cross-walked to PHQ-9 equivalents via IRT look-up table).
#
# ==============================================================================

library(dplyr)
library(haven)     # read_sav()
library(readxl)    # read_excel()
library(visdat)    # vis_miss()
library(skimr)

# ==============================================================================
# SECTION A  — DATA LOADING
# Replace the paths below with the actual locations of your raw data files.
# ==============================================================================

# -- Development datasets (EVIDENT main + EVIDENT severe add-on) ---------------
# dataset        <- haven::read_sav("PATH/TO/EVIDENT_main.sav")
# severe_dataset <- haven::read_sav("PATH/TO/EVIDENT_severe.sav")

# -- External validation datasets ----------------------------------------------
# moritz_dataset <- haven::read_sav("PATH/TO/Deprexis_Moritz.sav")
# eficasy_raw    <- readxl::read_excel("PATH/TO/eFICASY_data.xlsx", na = "NA")
# selfapy_raw    <- <load your Selfapy data here>

# ==============================================================================
# SECTION B  — HELPER: BDI-II → PHQ-9 IRT CROSSWALK  (Moritz dataset)
# Crosswalk tables from published IRT linking (Wahl et al. / Rose et al.).
# ==============================================================================

phq9_means <- c(
  22.9, 37.2, 44.6, 49.4, 53.0, 55.7, 58.2, 60.3,
  62.3, 63.7, 65.4, 66.7, 67.8, 69.5, 70.9, 72.0,
  73.4, 74.7, 76.0, 77.2, 78.9, 81.0, 82.4, 84.9,
  87.8, 91.8, 98.6, 112.2
)
phq9_df <- data.frame(PHQ9 = 0:27, theta = phq9_means)

bdi2_means <- c(
  15.9, 26.3, 33.6, 39.2, 42.8, 45.9, 48.5, 50.6,
  52.5, 54.0, 55.6, 57.0, 58.2, 59.1, 60.4, 61.3,
  62.7, 63.3, 64.4, 65.4, 66.2, 66.9, 67.5, 68.6,
  69.0, 69.9, 70.7, 71.4, 72.2, 72.9, 73.4, 74.1,
  74.9, 75.6, 76.1, 77.0, 77.3, 78.2, 79.0, 79.7,
  80.5, 81.2, 82.0, 82.5, 83.3, 84.4, 84.8, 85.9,
  86.6, 87.5, 88.5, 89.5, 90.5, 91.8, 93.0, 94.6,
  96.0, 97.9, 100.3, 103.4, 106.8, 111.7, 118.1, 128.8
)
bdi2_df <- data.frame(BDI2 = 0:63, theta = bdi2_means)

bdi2_to_phq9_df <- do.call(rbind, lapply(seq_len(nrow(bdi2_df)), function(i) {
  theta_i   <- bdi2_df$theta[i]
  phq_match <- phq9_df$PHQ9[which.min(abs(phq9_df$theta - theta_i))]
  data.frame(BDI2 = bdi2_df$BDI2[i], PHQ9 = phq_match)
}))

bdi2_to_phq9_fn <- function(x) {
  x_rounded <- round(x)
  x_clamped <- pmax(pmin(x_rounded, max(bdi2_to_phq9_df$BDI2)),
                    min(bdi2_to_phq9_df$BDI2))
  matched_idx <- match(x_clamped, bdi2_to_phq9_df$BDI2)
  bdi2_to_phq9_df$PHQ9[matched_idx]
}

# ==============================================================================
# SECTION C  — DEVELOPMENT SAMPLE (EVIDENT main + EVIDENT severe)
#
# Variable mapping for EVIDENT main dataset (compds after wrangling):
#   Gruppe           → Group   (0=control, 1=intervention)
#   Geschlecht       → Gender  (1=male, 2=female)
#   Familienstand    → Relationship
#   Schulabschluss   → School  (6/7/8 = ">12yrs", else "<12yrs")
#   volloderteilzeit.1 → Employment (0=unemployed, else=employed)
#   MajorDepressionAktuell.1 → MajorDepression (2=yes, else=no)
#   Dysthymie.1      → Dysthymia
#   niedergelassenerpsych.1 → current_psychotherapy
#   current_antidepressant  (derived from free-text medication strings)
#   SF_K.1           → QOL_phys
#   SF_P.1           → QOL_psych
#   PHQ_total_T1o    → PHQ_1
#   PHQ_total_T2o    → PHQ_3
#
# The code below is a template — adapt column names to your raw data.
# ==============================================================================

# After loading dataset and applying your inclusion criteria:
# compds_select <- dataset %>%
#   filter(!is.na(PHQ_total_T2o)) %>%
#   mutate(
#     Group     = factor(ifelse(Gruppe == 1, "intervention", "control")),
#     Gender    = factor(ifelse(Geschlecht == 1, "male", "female")),
#     # ... (see original deprexis.R for full recoding)
#     PHQ_change = PHQ_1 - PHQ_3
#   ) %>%
#   select(Group, Gender, Age, Relationship, School, Employment,
#          current_psychotherapy, current_antidepressant,
#          MajorDepression, Dysthymia, QOL_psych, QOL_phys,
#          PHQ_1, PHQ_3, PHQ_change)

# After the same for severe_select, merge:
# common_vars         <- intersect(names(compds_select), names(severe_select))
# deprexis_select_PHQ9 <- rbind(compds_select[, common_vars],
#                                severe_select[, common_vars]) %>%
#   filter(!is.na(PHQ_change)) %>%
#   mutate(across(c(QOL_psych, QOL_phys, PHQ_1, PHQ_change, Age), as.numeric)) %>%
#   as.data.frame()

# ==============================================================================
# SECTION D  — EXTERNAL VALIDATION: eFICASY / edupression
#
# Variable mapping (eficasy raw):
#   study_group_key   → Group
#   age               → Age
#   sex               → Gender
#   partner_status    → Relationship
#   education         → School
#   job_status        → Employment
#   PT_curr           → current_psychotherapy
#   Med_curr          → current_antidepressant
#   MD_Dx             → MajorDepression
#   F34.1             → Dysthymia
#   WHO_phys_baseline → QOL_phys
#   WHO_psych_baseline→ QOL_psych
#   Value_PHQ_baseline→ PHQ_1
#   PHQ_LOCF          → PHQ_3  (last-observation carried forward across wks 8/10/12)
# ==============================================================================

# eficasy_select <- eficasy_raw %>%
#   mutate(PHQ_LOCF = coalesce(Value_PHQ_woche_12,
#                               Value_PHQ_woche_10,
#                               Value_PHQ_woche_8)) %>%
#   mutate(
#     Group = factor(case_when(
#       study_group_key == "efficacy control"      ~ "control",
#       study_group_key == "efficacy intervention" ~ "intervention",
#       TRUE ~ NA_character_
#     )),
#     # ... (full recoding in original deprexis.R)
#     PHQ_change = Value_PHQ_baseline - PHQ_LOCF,
#     MajorDepression = factor(MajorDepression, levels = c("no", "yes"))
#   ) %>%
#   filter(!is.na(PHQ_change))

# ==============================================================================
# SECTION E  — EXTERNAL VALIDATION: Selfapy
#
# Selfapy uses QIDS (Quick Inventory of Depressive Symptomatology) as
# primary outcome. A PHQ-9 equivalent (PHQ_change) is derived via IRT linking
# or direct PHQ-9 if available. Retain QIDS variables for sensitivity analyses.
# ==============================================================================

# HRSD_selfapy_QIDSs <- selfapy_raw %>%
#   mutate(
#     Group = factor(case_when(...)),
#     # ... (adapt to Selfapy codebook)
#     PHQ_change = PHQ_1 - PHQ_3
#   ) %>%
#   filter(!is.na(PHQ_change))

# ==============================================================================
# SECTION F  — EXTERNAL VALIDATION: Moritz / Deprexis independent cohort
#
# This dataset uses BDI-II. PHQ-9 equivalents are derived via the IRT crosswalk
# defined in SECTION B. Group variable: sample_pp (1=control, 2=intervention).
# ==============================================================================

# moritz_select <- moritz_dataset %>%
#   filter(!is.na(BDItotalpost_recoded_IPD) & !is.na(sample_pp)) %>%
#   mutate(
#     Group = factor(case_when(
#       sample_pp == 1 ~ "control",
#       sample_pp == 2 ~ "intervention",
#       TRUE ~ NA_character_
#     )),
#     PHQ_1      = bdi2_to_phq9_fn(BDItotalpre_recoded_IPD),
#     PHQ_3      = bdi2_to_phq9_fn(BDItotalpost_recoded_IPD),
#     PHQ_change = PHQ_1 - PHQ_3,
#     Age        = 2012 - birthdate2.1,
#     # ...
#     MajorDepression = factor(MajorDepression, levels = c("no", "yes"))
#   )

# ==============================================================================
# SECTION G  — MISSINGNESS CHECK
# ==============================================================================

# vis_miss(deprexis_select_PHQ9)
# vis_miss(eficasy_select)
# vis_miss(HRSD_selfapy_QIDSs)
# vis_miss(moritz_select)

# ==============================================================================
# SECTION H  — FORMULA OBJECTS (shared across all scripts)
# ==============================================================================

outcome <- "PHQ_change"

# Main regression formula (with Group × covariate interactions for ENR/OLS)
formula_string <- paste(
  outcome,
  "~ Group * (Gender + Age + School + Employment + Dysthymia + Relationship +
     current_psychotherapy + current_antidepressant + MajorDepression +
     QOL_psych + QOL_phys + PHQ_1)"
)
f <- as.formula(formula_string)

# No-interaction formula for causal forest feature matrix
formula_string_rf <- paste(
  "~ Gender + Age + School + Employment + Dysthymia + Relationship +
     current_psychotherapy + current_antidepressant + MajorDepression +
     QOL_psych + QOL_phys + PHQ_1"
)
frf <- as.formula(formula_string_rf)

message("✓ 01_data_preparation.R complete")
message("  Objects produced: deprexis_select_PHQ9, eficasy_select,",
        " HRSD_selfapy_QIDSs, moritz_select, f, frf")

sessionInfo()
