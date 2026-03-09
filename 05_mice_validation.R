# ==============================================================================
# 05_mice_validation.R
#
# Purpose:
#   Sensitivity analysis repeating outcome validation and benefit calibration
#   using multiple imputation by chained equations (MICE, m = 20).
#   Applies Rubin's rules to pool estimates across imputations.
#
# Requires (from earlier scripts):
#   deprexis_select_PHQ9, eficasy_select, HRSD_selfapy_QIDSs, moritz_select
#   full_model_deprexis, full_model_deprexis_OLS, m3_final (or refitted below)
#   optalpha, optlambda, sd_phq_dev
#   f, frf, outcome
#   Combined CV metrics: enr_cv_rmse_mean/sd etc. (from 02)
#   OLS  CV metrics: ols_cv_rmse_mean/sd etc.
#
# Outputs:
#   imp_list_eficasy, imp_list_selfapy, imp_list_moritz — MI lists
#   validation_results (MI-pooled outcome performance)
#   benefit_prediction_results (MI-pooled benefit calibration)
#   pb_results_complete (MI-pooled population benefit)
#   Several CSV files and one TIFF/PNG figure
# ==============================================================================

library(dplyr)
library(mice)
library(glmnet)
library(grf)
library(caret)
library(ggplot2)

# Prevent dplyr / mice conflicts
summarise <- dplyr::summarise
summarize <- dplyr::summarize

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

rubin_scalar <- function(estimates, ses) {
  estimates <- as.numeric(estimates); ses <- as.numeric(ses)
  m     <- length(estimates)
  q_bar <- mean(estimates)
  U_bar <- mean(ses^2)
  B     <- ifelse(is.na(var(estimates)), 0, var(estimates))
  T_var <- U_bar + (1 + 1/m) * B
  list(estimate = q_bar, se = sqrt(T_var))
}

calculate_pooled_sd <- function(dat) {
  Y1 <- dat$Y[dat$t == 1]; Y0 <- dat$Y[dat$t == 0]
  n1 <- length(Y1);         n0 <- length(Y0)
  sqrt(((n1-1)*sd(Y1)^2 + (n0-1)*sd(Y0)^2) / (n1+n0-2))
}

label_sample <- function(ds) dplyr::case_when(
  ds == "Deprexis"    ~ "Deprexis (Dev)",
  ds == "edupression" ~ "edupression",
  ds == "Selfapy"     ~ "Selfapy",
  ds == "Moritz"      ~ "Deprexis (Val)"
)
sample_levels <- c("Deprexis (Dev)", "edupression", "Selfapy", "Deprexis (Val)")

# ==============================================================================
# SECTION 1  — MICE IMPUTATION  (m = 20 per dataset)
#
# Predictors used: Group, Gender, Age, School, Employment, Dysthymia,
#   Relationship, current_psychotherapy, current_antidepressant,
#   MajorDepression, QOL_psych, QOL_phys, PHQ_1
# PHQ_change is excluded from imputation (outcome).
# ==============================================================================

f_predictors <- c("Group", "Gender", "Age", "School", "Employment", "Dysthymia",
                   "Relationship", "current_psychotherapy",
                   "current_antidepressant", "MajorDepression",
                   "QOL_psych", "QOL_phys", "PHQ_1")

run_mice <- function(dat, missing_vars, methods_map, m = 20, seed = 42) {
  meth <- make.method(dat); meth[] <- ""
  pred <- make.predictorMatrix(dat); pred[] <- 0
  for (v in missing_vars) {
    meth[v] <- methods_map[[v]]
    avail    <- f_predictors[f_predictors %in% names(dat)]
    pred[v, avail] <- 1
    pred[v, v]     <- 0
    if ("PHQ_change" %in% colnames(pred)) pred[v, "PHQ_change"] <- 0
  }
  imp <- mice(dat, m = m, maxit = 20, method = meth,
               predictorMatrix = pred, seed = seed, printFlag = FALSE)
  lapply(seq_len(m), function(i) complete(imp, i))
}

# edupression / eFICASY: QOL_psych (14%), QOL_phys (14%) — continuous → pmm
imp_list_eficasy <- run_mice(
  eficasy_select,
  missing_vars  = c("QOL_psych", "QOL_phys"),
  methods_map   = list(QOL_psych = "pmm", QOL_phys = "pmm")
)

# Selfapy: Relationship (10%), School (20%), Employment (24%) — binary → logreg
imp_list_selfapy <- run_mice(
  HRSD_selfapy_QIDSs,
  missing_vars = c("Relationship", "School", "Employment"),
  methods_map  = list(Relationship = "logreg", School = "logreg",
                       Employment = "logreg")
)

# Moritz: current_psychotherapy (24%), current_antidepressant (10%) — logreg
imp_list_moritz  <- run_mice(
  moritz_select,
  missing_vars = c("current_psychotherapy", "current_antidepressant"),
  methods_map  = list(current_psychotherapy  = "logreg",
                       current_antidepressant = "logreg")
)

validation_mi_lists <- list(
  "edupression" = imp_list_eficasy,
  "Selfapy"     = imp_list_selfapy,
  "Moritz"      = imp_list_moritz
)

message("✓ MICE complete.  Imputed datasets: ",
        paste(sapply(validation_mi_lists, length), collapse = " / "))

# ==============================================================================
# SECTION 2  — MI-POOLED OUTCOME VALIDATION
# (mirrors 03_outcome_validation.R but pools across imputations)
# ==============================================================================

validate_one_mi <- function(val_data, dataset_name) {
  if ("MajorDepression" %in% names(val_data))
    val_data$MajorDepression <- factor(val_data$MajorDepression,
                                        levels = c("no", "yes"))

  newx     <- model.matrix(f, data = val_data)[, -1]
  pred_enr <- as.numeric(predict(full_model_deprexis,
                                  newx = newx, s = optlambda, type = "link"))
  pred_ols <- as.numeric(predict(full_model_deprexis_OLS,
                                  newx = newx, s = 0, type = "link"))
  obs      <- val_data[[outcome]]

  make_res <- function(pred, obs) {
    perf <- postResample(pred = pred, obs = obs)
    cal  <- lm(obs ~ pred)
    list(RMSE = perf["RMSE"], MAE = perf["MAE"], R2 = perf["Rsquared"],
         Cal_Intercept = coef(cal)[1],
         Cal_Int_SE    = coef(summary(cal))[1, 2],
         Cal_Slope     = coef(cal)[2],
         Cal_Slope_SE  = coef(summary(cal))[2, 2],
         SD_PHQ        = sd(obs, na.rm = TRUE))
  }
  list(ENR = make_res(pred_enr, obs), OLS = make_res(pred_ols, obs))
}

validation_results_mi <- list()

for (dname in names(validation_mi_lists)) {
  message("\nMI outcome validation: ", dname)
  imp_res <- lapply(validation_mi_lists[[dname]], validate_one_mi,
                    dataset_name = dname)
  pooled  <- list()
  for (mod in c("ENR", "OLS")) {
    get_vals <- function(fld) sapply(imp_res, function(r) r[[mod]][[fld]])
    pooled[[mod]] <- list(
      RMSE          = mean(get_vals("RMSE")),
      MAE           = mean(get_vals("MAE")),
      R2            = mean(get_vals("R2")),
      SD_PHQ        = mean(get_vals("SD_PHQ")),
      Cal_Intercept = rubin_scalar(get_vals("Cal_Intercept"),
                                   get_vals("Cal_Int_SE"))$estimate,
      Cal_Int_SE    = rubin_scalar(get_vals("Cal_Intercept"),
                                   get_vals("Cal_Int_SE"))$se,
      Cal_Slope     = rubin_scalar(get_vals("Cal_Slope"),
                                   get_vals("Cal_Slope_SE"))$estimate,
      Cal_Slope_SE  = rubin_scalar(get_vals("Cal_Slope"),
                                   get_vals("Cal_Slope_SE"))$se
    )
  }
  validation_results_mi[[dname]] <- pooled
}

# ==============================================================================
# SECTION 3  — MI-POOLED BENEFIT CALIBRATION
# (Rubin pooling across 20 imputations per dataset)
# ==============================================================================

# pool_performance_df: pools calibration regression metrics (Rubin's rules)
pool_performance_df <- function(perf_list) {
  combined <- bind_rows(perf_list, .id = "imp")

  rubin_ses <- combined %>%
    dplyr::group_by(Dataset, Method, N_bins) %>%
    dplyr::summarise(
      intercept_se_abs = {
        ests <- intercept_abs; ses <- intercept_se_abs; m_ <- dplyr::n()
        sqrt(mean(ses^2) + (1 + 1/m_) * var(ests))
      },
      slope_se_abs = {
        ests <- slope_abs; ses <- slope_se_abs; m_ <- dplyr::n()
        sqrt(mean(ses^2) + (1 + 1/m_) * var(ests))
      },
        .groups = "drop"
    )

  means <- combined %>%
    dplyr::group_by(Dataset, Method, N_bins) %>%
    dplyr::summarise(
      across(c(RMSE_abs, RMSE_abs_standardized,
                r2_abs, 
               intercept_abs, slope_abs),
             ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  means %>%
    dplyr::left_join(rubin_ses, by = c("Dataset", "Method", "N_bins")) %>%
    dplyr::mutate(
      intercept_lower_abs   = intercept_abs   - 1.96 * intercept_se_abs,
      intercept_upper_abs   = intercept_abs   + 1.96 * intercept_se_abs,
      slope_lower_abs       = slope_abs       - 1.96 * slope_se_abs,
      slope_upper_abs       = slope_abs       + 1.96 * slope_se_abs
    )
}

# process_one_imp: compute benefit predictions + calibration for one imputed dataset
process_one_imp <- function(val_data, dataset_name) {
  if (!is.numeric(val_data$Group))
    val_data <- val_data %>% mutate(Group = ifelse(Group == "intervention", 1, 0))
  if ("MajorDepression" %in% names(val_data))
    val_data$MajorDepression <- factor(val_data$MajorDepression,
                                        levels = c("no", "yes"))
 
  X_val    <- model.matrix(frf, data = val_data)[, -1]
  X_val_t1 <- model.matrix(f, data = transform(val_data, Group = 1))[, -1]
  X_val_t0 <- model.matrix(f, data = transform(val_data, Group = 0))[, -1]
  Y_val    <- val_data[[outcome]]
  T_val    <- val_data$Group

  benefit_enr <- as.numeric(predict(full_model_deprexis, newx = X_val_t1,
                                     s = optlambda) -
                               predict(full_model_deprexis, newx = X_val_t0,
                                       s = optlambda))
  benefit_ols <- as.numeric(predict(full_model_deprexis_OLS, newx = X_val_t1,
                                     s = 0) -
                               predict(full_model_deprexis_OLS, newx = X_val_t0,
                                       s = 0))
  benefit_cf  <- predict(m3_final, newdata = X_val)$predictions

  dat_enr <- data.frame(X_val, Y = Y_val, t = T_val, benefit = benefit_enr,
                         check.names = FALSE)
  dat_ols <- data.frame(X_val, Y = Y_val, t = T_val, benefit = benefit_ols,
                         check.names = FALSE)
  dat_cf  <- data.frame(X_val, Y = Y_val, t = T_val, benefit = benefit_cf,
                         check.names = FALSE)

  # Source calibrate_bins from 04_benefit_calibration.R (or define inline)
  # Here we call calibrate_bins which must be in scope
  Ngroups_val <- c( 3, 4, 5)
  cal_list    <- list(); perf_list <- list()
  for (mname in c("ENR", "OLS", "CF")) {
    dat_m <- list("ENR" = dat_enr, "OLS" = dat_ols, "CF" = dat_cf)[[mname]]
    for (k in Ngroups_val) {
      res <- calibrate_bins(dat_m, mname, dataset_name, k)
      key <- paste0(mname, "_bins", k)
      cal_list[[key]]  <- res$calibration
      perf_list[[key]] <- res$performance
    }
  }

  list(
    dats        = list(enr = dat_enr, ols = dat_ols, cf = dat_cf),
    calibration = bind_rows(cal_list),
    performance = bind_rows(perf_list)
  )
}

# compute_pb_one: population benefit for one imputation
compute_pb_one <- function(dat1, cutoff, dataset_name) {
  g11 <- dat1[dat1$benefit >  cutoff & dat1$t == 1, ]
  g12 <- dat1[dat1$benefit >  cutoff & dat1$t == 0, ]
  g13 <- dat1[dat1$benefit <= cutoff & dat1$t == 1, ]
  g14 <- dat1[dat1$benefit <= cutoff & dat1$t == 0, ]
  n11 <- nrow(g11); n12 <- nrow(g12)
  n13 <- nrow(g13); n14 <- nrow(g14)

  if (any(c(n11, n12, n13, n14) == 0)) {
    return(data.frame(Dataset = dataset_name, Method = "enr", Cutoff = cutoff,
                       n_G1 = n11, n_G2 = n12, n_G3 = n13, n_G4 = n14,
                       prop_optimal = (n11+n14)/nrow(dat1),
                       PB_adjusted = NA, PB_adj_se = NA,
                       PB1 = NA, PB1_se = NA))
  }

  dat1$agree <- as.numeric(sign(dat1$benefit) == sign(2 * dat1$t - 1))
  cov_cols   <- setdiff(names(dat1), c("Y", "t", "benefit", "agree"))
  dat_reg    <- data.frame(Y = dat1$Y, agree = dat1$agree,
                            dat1[, cov_cols, drop = FALSE])
  PB_adj <- SE_adj <- NA
  tryCatch({
    lm_adj  <- lm(Y ~ ., data = dat_reg)
    cag     <- summary(lm_adj)$coef["agree", ]
    PB_adj  <- cag[1]; SE_adj <- cag[2]
  }, error = function(e) NULL)

  PB1     <- sum(g14$Y)/(n11+n14) - sum(g13$Y)/(n11+n13) +
    sum(g11$Y) * (1/(n11+n14) - 1/(n11+n13))
  var_PB1 <- n14*var(g14$Y)/(n11+n14)^2 + n13*var(g13$Y)/(n11+n13)^2 +
    n11*var(g11$Y) * (1/(n11+n14) - 1/(n11+n13))^2

  data.frame(Dataset = dataset_name, Method = "enr", Cutoff = cutoff,
              n_G1 = n11, n_G2 = n12, n_G3 = n13, n_G4 = n14,
              prop_optimal = (n11+n14)/nrow(dat1),
              PB_adjusted = PB_adj, PB_adj_se = SE_adj,
              PB1 = PB1, PB1_se = sqrt(var_PB1))
}

# pool_pb_rubin: pool PB1 and PB_adjusted across imputations
pool_pb_rubin <- function(pb_list) {
  bind_rows(pb_list, .id = "imp") %>%
    dplyr::filter(!is.na(PB1_se)) %>%
    dplyr::group_by(Dataset, Method, Cutoff) %>%
    dplyr::summarise(
      n_G1         = mean(n_G1),
      n_G2         = mean(n_G2),
      n_G3         = mean(n_G3),
      n_G4         = mean(n_G4),
      prop_optimal = mean(prop_optimal, na.rm = TRUE),
      pb1_result   = list(rubin_scalar(as.numeric(PB1), as.numeric(PB1_se))),
      pbadj_result = list(
        if (any(!is.na(PB_adj_se)))
          rubin_scalar(as.numeric(PB_adjusted[!is.na(PB_adj_se)]),
                       as.numeric(PB_adj_se[!is.na(PB_adj_se)]))
        else list(estimate = NA_real_, se = NA_real_)
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      PB1         = sapply(pb1_result,   `[[`, "estimate"),
      PB1_se      = sapply(pb1_result,   `[[`, "se"),
      PB_adjusted = sapply(pbadj_result, `[[`, "estimate"),
      PB_adj_se   = sapply(pbadj_result, `[[`, "se"),
      PB1_lower    = PB1         - 1.96 * PB1_se,
      PB1_upper    = PB1         + 1.96 * PB1_se,
      PB_adj_lower = PB_adjusted - 1.96 * PB_adj_se,
      PB_adj_upper = PB_adjusted + 1.96 * PB_adj_se
    ) %>%
    dplyr::select(-pb1_result, -pbadj_result, -PB1_se, -PB_adj_se)
}

# ==============================================================================
# SECTION 4  — MAIN MI LOOP
# ==============================================================================

dat_storage_mi                <- list()
benefit_prediction_results_mi <- list()
pb_imp_d024                   <- list()
pb_imp_d050                   <- list()

for (dname in names(validation_mi_lists)) {
  message("\n=== MI benefit analysis: ", dname, " ===")
  imp_list    <- validation_mi_lists[[dname]]
  imp_results <- lapply(imp_list, process_one_imp, dataset_name = dname)
  m_imp       <- length(imp_list)

  # Average benefit predictions across imputations (Y and t are observed)
  dat_storage_mi[[dname]] <- list(
    enr = imp_results[[1]]$dats$enr %>%
      mutate(benefit = rowMeans(sapply(imp_results,
                                       function(r) r$dats$enr$benefit))),
    ols = imp_results[[1]]$dats$ols %>%
      mutate(benefit = rowMeans(sapply(imp_results,
                                       function(r) r$dats$ols$benefit))),
    cf  = imp_results[[1]]$dats$cf  %>%
      mutate(benefit = rowMeans(sapply(imp_results,
                                       function(r) r$dats$cf$benefit)))
  )

  # Pool calibration
  benefit_prediction_results_mi[[dname]]$calibration <-
    bind_rows(lapply(imp_results, `[[`, "calibration")) %>%
    dplyr::group_by(Dataset, Method, N_bins, benefit_group) %>%
    dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                     .groups = "drop")

  benefit_prediction_results_mi[[dname]]$performance <-
    pool_performance_df(lapply(imp_results, `[[`, "performance"))

  # PB at both thresholds
  psd_ds   <- calculate_pooled_sd(dat_storage_mi[[dname]]$enr)
  cutoff24 <- 0.24 * psd_ds
  cutoff50 <- 0.50 * psd_ds

  pb_imp_d024[[dname]] <- lapply(seq_len(m_imp), function(i)
    compute_pb_one(imp_results[[i]]$dats$enr, cutoff24, dname))
  pb_imp_d050[[dname]] <- lapply(seq_len(m_imp), function(i)
    compute_pb_one(imp_results[[i]]$dats$enr, cutoff50, dname))

  message("✓ Done: ", dname)
}

# Pool PB
pb_pooled_d024 <- bind_rows(lapply(names(pb_imp_d024),
                                    function(ds) pool_pb_rubin(pb_imp_d024[[ds]])))
pb_pooled_d050 <- bind_rows(lapply(names(pb_imp_d050),
                                    function(ds) pool_pb_rubin(pb_imp_d050[[ds]])))

# ==============================================================================
# SECTION 5  — SAVE OUTPUTS
# ==============================================================================

all_perf_mi <- bind_rows(
  performances_dev_final,   # reuse from 04 (development is completer)
  bind_rows(lapply(names(benefit_prediction_results_mi),
                   function(ds) benefit_prediction_results_mi[[ds]]$performance))
) %>%
  mutate(Sample = factor(label_sample(Dataset), levels = sample_levels),
         Method = factor(Method, levels = c("ENR", "OLS", "CF")))

supp_cal_mi <- all_perf_mi %>%
  filter((Dataset == "Deprexis" & N_bins %in% c(4, 11, 18)) |
           (Dataset != "Deprexis" & N_bins %in% c(4, 5))) %>%
  mutate(
    Intercept = sprintf("%.2f (%.2f to %.2f)",
                        intercept_abs, intercept_lower_abs, intercept_upper_abs),
    Slope     = sprintf("%.2f (%.2f to %.2f)",
                        slope_abs,     slope_lower_abs,     slope_upper_abs),
    RMSE      = round(RMSE_abs_standardized, 2),
    R2        = round(r2_abs, 2)
  ) %>%
  select(Sample, Method, `N groups` = N_bins, Intercept, Slope, RMSE, R2)

write.csv(pb_pooled_d024,  "population_benefit_d024_mice.csv",             row.names = FALSE)
write.csv(pb_pooled_d050,  "population_benefit_d050_mice.csv",             row.names = FALSE)
write.csv(supp_cal_mi,     "supplement_calibration_table_mice.csv",        row.names = FALSE)

message("✓ 05_mice_validation.R complete")
sessionInfo()
