# ==============================================================================
# 03_outcome_validation.R
#
# Purpose:
#   Assess outcome prediction performance (RMSE, MAE, R², calibration) of
#   ENR and OLS in each external validation sample (completer analysis, no MI).
#   Assembles performance_data_complete and writes the supplementary table.
#
# Requires (from earlier scripts):
#   full_model_deprexis, full_model_deprexis_OLS
#   optalpha, optlambda
#   sd_phq_dev
#   f, outcome
#   combined_metrics        (ENR CV metrics)
#   # OLS CV scalars: ols_cv_rmse_mean/sd, ols_cv_mae_mean/sd etc.
#   selfapy_imputed, eficasy_imputed, moritz_imputed
#
# Outputs:
#   validation_results              — list of per-dataset ENR/OLS metrics
#   performance_data_complete       — wide data frame for plotting
#   Supplementary_Table_Outcome_Performance_completer.csv
# ==============================================================================

library(dplyr)
library(caret)

# ==============================================================================
# SECTION 1  — VALIDATE ONE DATASET
# ==============================================================================

validate_one <- function(val_data, dataset_name) {

  # Ensure consistent factor levels
  if ("MajorDepression" %in% names(val_data))
    val_data$MajorDepression <- factor(val_data$MajorDepression,
                                        levels = c("no", "yes"))
  

  newx <- model.matrix(f, data = val_data)[, -1]

  pred_enr <- as.numeric(predict(full_model_deprexis,
                                  newx = newx, s = optlambda, type = "link"))
  pred_ols <- as.numeric(predict(full_model_deprexis_OLS,
                                  newx = newx, s = 0,         type = "link"))
  obs      <- val_data[[outcome]]

  perf_enr <- postResample(pred = pred_enr, obs = obs)
  perf_ols <- postResample(pred = pred_ols, obs = obs)
  cal_enr  <- lm(obs ~ pred_enr)
  cal_ols  <- lm(obs ~ pred_ols)

  list(
    ENR = list(
      RMSE          = perf_enr["RMSE"],
      MAE           = perf_enr["MAE"],
      R2            = perf_enr["Rsquared"],
      Cal_Intercept = coef(cal_enr)[1],
      Cal_Int_SE    = coef(summary(cal_enr))[1, 2],
      Cal_Slope     = coef(cal_enr)[2],
      Cal_Slope_SE  = coef(summary(cal_enr))[2, 2],
      SD_PHQ        = sd(obs, na.rm = TRUE)
    ),
    OLS = list(
      RMSE          = perf_ols["RMSE"],
      MAE           = perf_ols["MAE"],
      R2            = perf_ols["Rsquared"],
      Cal_Intercept = coef(cal_ols)[1],
      Cal_Int_SE    = coef(summary(cal_ols))[1, 2],
      Cal_Slope     = coef(cal_ols)[2],
      Cal_Slope_SE  = coef(summary(cal_ols))[2, 2],
      SD_PHQ        = sd(obs, na.rm = TRUE)
    )
  )
}

# ==============================================================================
# SECTION 2  — RUN EXTERNAL VALIDATION
# ==============================================================================

validation_datasets <- list(
  "Selfapy"     = selfapy_imputed,
  "Eficasy"     = eficasy_imputed,
  "Moritz"      = moritz_imputed
)

validation_results <- list()

for (dname in names(validation_datasets)) {
  message("\n========================================")
  message("EXTERNAL VALIDATION: ", dname)
  message("========================================")
  validation_results[[dname]] <- validate_one(validation_datasets[[dname]],
                                               dataset_name = dname)
  message("  ENR RMSE: ", round(validation_results[[dname]]$ENR$RMSE, 2),
          " | R²: ",      round(validation_results[[dname]]$ENR$R2,   3),
          " | Slope: ",   round(validation_results[[dname]]$ENR$Cal_Slope, 2))
}

# ==============================================================================
# SECTION 3  — ASSEMBLE performance_data_complete
# ==============================================================================

message("\nUsing development sample SD (", round(sd_phq_dev, 2),
        ") for standardisation across all samples")

performance_data_complete <- data.frame(
  Sample = rep(c("Deprexis (Dev)", "Selfapy", "edupression", "Deprexis (Val)"),
               each = 2),
  Model  = rep(c("ENR", "OLS"), 4),
  RMSE = c(
    enr_cv_rmse_mean,                       ols_cv_rmse_mean,
    validation_results$Selfapy$ENR$RMSE,    validation_results$Selfapy$OLS$RMSE,
    validation_results$Eficasy$ENR$RMSE,    validation_results$Eficasy$OLS$RMSE,
    validation_results$Moritz$ENR$RMSE,     validation_results$Moritz$OLS$RMSE
  ),
  RMSE_SD = c(enr_cv_rmse_sd, ols_cv_rmse_sd, NA, NA, NA, NA, NA, NA),
  MAE = c(
    enr_cv_mae_mean,                        ols_cv_mae_mean,
    validation_results$Selfapy$ENR$MAE,     validation_results$Selfapy$OLS$MAE,
    validation_results$Eficasy$ENR$MAE,     validation_results$Eficasy$OLS$MAE,
    validation_results$Moritz$ENR$MAE,      validation_results$Moritz$OLS$MAE
  ),
  MAE_SD = c(enr_cv_mae_sd, ols_cv_mae_sd, NA, NA, NA, NA, NA, NA),
  SD_PHQ = c(
    sd_phq_dev,                                  sd_phq_dev,
    validation_results$Selfapy$ENR$SD_PHQ,       validation_results$Selfapy$OLS$SD_PHQ,
    validation_results$Eficasy$ENR$SD_PHQ,       validation_results$Eficasy$OLS$SD_PHQ,
    validation_results$Moritz$ENR$SD_PHQ,        validation_results$Moritz$OLS$SD_PHQ
  ),
  SD_PHQ_dev = sd_phq_dev,
  R2 = c(
    enr_cv_r2_mean,                         ols_cv_r2_mean,
    validation_results$Selfapy$ENR$R2,      validation_results$Selfapy$OLS$R2,
    validation_results$Eficasy$ENR$R2,      validation_results$Eficasy$OLS$R2,
    validation_results$Moritz$ENR$R2,       validation_results$Moritz$OLS$R2
  ),
  R2_SD = c(enr_cv_r2_sd, ols_cv_r2_sd, NA, NA, NA, NA, NA, NA),
  Cal_Intercept = c(
    enr_cv_intercept_mean,                            ols_cv_intercept_mean,
    validation_results$Selfapy$ENR$Cal_Intercept,     validation_results$Selfapy$OLS$Cal_Intercept,
    validation_results$Eficasy$ENR$Cal_Intercept,     validation_results$Eficasy$OLS$Cal_Intercept,
    validation_results$Moritz$ENR$Cal_Intercept,      validation_results$Moritz$OLS$Cal_Intercept
  ),
  Cal_Int_SD = c(
    enr_cv_intercept_sd,                         ols_cv_intercept_sd,
    validation_results$Selfapy$ENR$Cal_Int_SE,   validation_results$Selfapy$OLS$Cal_Int_SE,
    validation_results$Eficasy$ENR$Cal_Int_SE,   validation_results$Eficasy$OLS$Cal_Int_SE,
    validation_results$Moritz$ENR$Cal_Int_SE,    validation_results$Moritz$OLS$Cal_Int_SE
  ),
  Cal_Slope = c(
    enr_cv_slope_mean,                           ols_cv_slope_mean,
    validation_results$Selfapy$ENR$Cal_Slope,    validation_results$Selfapy$OLS$Cal_Slope,
    validation_results$Eficasy$ENR$Cal_Slope,    validation_results$Eficasy$OLS$Cal_Slope,
    validation_results$Moritz$ENR$Cal_Slope,     validation_results$Moritz$OLS$Cal_Slope
  ),
  Cal_Slope_SD = c(
    enr_cv_slope_sd,                               ols_cv_slope_sd,
    validation_results$Selfapy$ENR$Cal_Slope_SE,   validation_results$Selfapy$OLS$Cal_Slope_SE,
    validation_results$Eficasy$ENR$Cal_Slope_SE,   validation_results$Eficasy$OLS$Cal_Slope_SE,
    validation_results$Moritz$ENR$Cal_Slope_SE,    validation_results$Moritz$OLS$Cal_Slope_SE
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Cal_Int_Lower   = Cal_Intercept - 1.96 * Cal_Int_SD,
    Cal_Int_Upper   = Cal_Intercept + 1.96 * Cal_Int_SD,
    Cal_Slope_Lower = Cal_Slope     - 1.96 * Cal_Slope_SD,
    Cal_Slope_Upper = Cal_Slope     + 1.96 * Cal_Slope_SD,
    RMSE_std        = RMSE    / SD_PHQ_dev,
    MAE_std         = MAE     / SD_PHQ_dev,
    RMSE_std_SD     = RMSE_SD / SD_PHQ_dev,
    MAE_std_SD      = MAE_SD  / SD_PHQ_dev,
    Sample = factor(Sample, levels = c("Deprexis (Dev)", "Selfapy",
                                        "edupression",   "Deprexis (Val)")),
    Model  = factor(Model,  levels = c("ENR", "OLS"))
  )

# ==============================================================================
# SECTION 4  — SUPPLEMENTARY TABLE
# ==============================================================================

supp_table <- performance_data_complete %>%
  mutate(
    `RMSE` = ifelse(!is.na(RMSE_SD),
                    sprintf("%.2f (%.2f)", RMSE, RMSE_SD),
                    sprintf("%.2f", RMSE)),
    `RMSE (SD units)*` = ifelse(!is.na(RMSE_SD),
                                sprintf("%.2f (%.2f)", RMSE_std, RMSE_std_SD),
                                sprintf("%.2f", RMSE_std)),
    `MAE` = ifelse(!is.na(MAE_SD),
                   sprintf("%.2f (%.2f)", MAE, MAE_SD),
                   sprintf("%.2f", MAE)),
    `R²` = ifelse(!is.na(R2_SD),
                  sprintf("%.2f (%.2f)", R2, R2_SD),
                  sprintf("%.2f", R2)),
    `Calibration Intercept (95% CI)` = sprintf("%.2f (%.2f to %.2f)",
                                                Cal_Intercept,
                                                Cal_Int_Lower, Cal_Int_Upper),
    `Calibration Slope (95% CI)`     = sprintf("%.2f (%.2f to %.2f)",
                                                Cal_Slope,
                                                Cal_Slope_Lower, Cal_Slope_Upper),
    `SD of PHQ-9 change` = sprintf("%.2f", SD_PHQ)
  ) %>%
  select(Sample, Model, RMSE, `RMSE (SD units)*`, MAE, `R²`,
         `Calibration Intercept (95% CI)`,
         `Calibration Slope (95% CI)`,
         `SD of PHQ-9 change`)

print(supp_table)
write.csv(supp_table,
          "Supplementary_Table_Outcome_Performance_completer.csv",
          row.names = FALSE)
message("✓ Supplementary outcome table saved.")

message("✓ 03_outcome_validation.R complete")
sessionInfo()
