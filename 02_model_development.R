# ==============================================================================
# 02_model_development.R
#
# Purpose:
#   (1)  Tune and cross-validate elastic net regression (ENR) and OLS on the
#        development sample (Deprexis) to obtain internal performance metrics.
#   (2)  Bootstrap coefficient stability (1 000 resamples).
#   (3)  Fit final full-sample ENR and OLS models.
#   (4)  Preprocess development data (mode/median imputation) and apply the
#        same preprocessing object to each validation dataset.
#
# Requires (from 01_data_preparation.R):
#   deprexis_select_PHQ9, eficasy_select, HRSD_selfapy_QIDSs, moritz_select
#   f, frf, outcome
#
# Outputs:
#   deprexis_d             — preprocessed development data frame (numeric Group)
#   full_model_deprexis    — final ENR glmnet object
#   full_model_deprexis_OLS— final OLS glmnet object (alpha=0, lambda=0)
#   optalpha, optlambda    — tuning parameters (medians across CV folds)
#   combined_metrics       — per-fold CV metrics (ENR)
#   pooled_results         — bootstrap coefficient stability table
#   selfapy_imputed, eficasy_imputed, moritz_imputed — preprocessed validation sets
#   sd_phq_dev             — mean training-fold SD (used for standardisation)
# ==============================================================================

library(dplyr)
library(caret)
library(glmnet)
library(doParallel)
library(foreach)

# -- assumes objects from 01_data_preparation.R are in environment:
# deprexis_select_PHQ9, eficasy_select, HRSD_selfapy_QIDSs, moritz_select, f, frf, outcome

# ==============================================================================
# SECTION 1  — HYPERPARAMETER GRID
# ==============================================================================

finegrid   <- 10^(seq(-1, -2.2, length.out = 50))
alpha_grid <- seq(0, 1, by = 0.1)
grids      <- expand.grid(alpha = alpha_grid, lambda = finegrid)

# ==============================================================================
# SECTION 2  — CROSS-VALIDATION  (10-fold × 10 repeats, parallelised)
#
# Produces combined_metrics: one row per repeat × fold with
#   TrainRMSE, TrainMAE, TrainRSquared,
#   TestRMSE,  TestMAE,  TestRSquared,
#   intercept, slope  (calibration regression: obs ~ pred),
#   alpha, lambda, sd_y_train, sd_y_test
# ==============================================================================

raw_data       <- as.data.frame(deprexis_select_PHQ9)
raw_data$row_id <- seq_len(nrow(raw_data))

split_ratio <- 0.9
n_splits    <- 10
n_repeats   <- 10   # increase to 100 for the published analysis

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

results <- foreach(
  repeat_split = seq_len(n_repeats),
  .combine      = rbind,
  .packages     = c("caret", "glmnet", "dplyr"),
  .errorhandling = "remove",
  .verbose       = FALSE
) %dopar% {

  set.seed(44 + repeat_split * 1044)
  seeds <- sample(1:5000, n_splits)

  inner_results <- lapply(seq_len(n_splits), function(split_index) {

    set.seed(seeds[split_index])

    trainIndex <- createDataPartition(raw_data[, outcome],
                                      p = split_ratio, list = FALSE)
    trainData  <- raw_data[ trainIndex, ]
    testData   <- raw_data[-trainIndex, ]

    # Median imputation fitted on training data only
    impute_model <- preProcess(trainData, method = "medianImpute")
    trainData    <- predict(impute_model, trainData)
    testData     <- predict(impute_model, testData)

    x     <- model.matrix(f, data = trainData)[, -1]
    y     <- trainData[, outcome]
    x_test <- model.matrix(f, data = testData)[, -1]
    y_test  <- testData[, outcome]

    ctrl <- trainControl(method = "repeatedcv", repeats = 10,
                         classProbs = FALSE)
    fit  <- train(x, y, method = "glmnet", trControl = ctrl,
                  family = "gaussian", tuneGrid = grids, metric = "MAE")

    test_preds <- predict(fit, newdata = x_test)
    lm_cal     <- lm(y_test ~ test_preds)

    data.frame(
      repeat_split  = repeat_split,
      split         = split_index,
      seed          = seeds[split_index],
      TrainRMSE     = getTrainPerf(fit)$TrainRMSE,
      TrainMAE      = getTrainPerf(fit)$TrainMAE,
      TrainRSquared = getTrainPerf(fit)$TrainRsquared,
      TestRMSE      = sqrt(mean((test_preds - y_test)^2)),
      TestMAE       = mean(abs(test_preds - y_test)),
      TestRSquared  = cor(test_preds, y_test)^2,
      intercept     = coef(lm_cal)[1],
      slope         = coef(lm_cal)[2],
      alpha         = fit$bestTune$alpha,
      lambda        = fit$bestTune$lambda,
      sd_y_train    = sd(y),
      sd_y_test     = sd(y_test)
    )
  })

  do.call(rbind, inner_results)
}

stopCluster(cl)
foreach::registerDoSEQ()

combined_metrics <- do.call(rbind, lapply(results[, 1], function(x) x))

# Summarise CV performance
enr_cv_rmse_mean      <- mean(combined_metrics$TestRMSE,     na.rm = TRUE)
enr_cv_rmse_sd        <- sd(combined_metrics$TestRMSE,       na.rm = TRUE)
enr_cv_mae_mean       <- mean(combined_metrics$TestMAE,      na.rm = TRUE)
enr_cv_mae_sd         <- sd(combined_metrics$TestMAE,        na.rm = TRUE)
enr_cv_r2_mean        <- mean(combined_metrics$TestRSquared, na.rm = TRUE)
enr_cv_r2_sd          <- sd(combined_metrics$TestRSquared,   na.rm = TRUE)
enr_cv_intercept_mean <- mean(combined_metrics$intercept,    na.rm = TRUE)
enr_cv_intercept_sd   <- sd(combined_metrics$intercept,      na.rm = TRUE)
enr_cv_slope_mean     <- mean(combined_metrics$slope,        na.rm = TRUE)
enr_cv_slope_sd       <- sd(combined_metrics$slope,          na.rm = TRUE)
sd_phq_dev            <- mean(combined_metrics$sd_y_train,   na.rm = TRUE)

optalpha  <- median(combined_metrics$alpha)
optlambda <- median(combined_metrics$lambda)

message("ENR — optimal alpha: ", round(optalpha, 3),
        "  lambda: ", round(optlambda, 4))

# ==============================================================================
# SECTION 3  — PREPROCESS FULL DEVELOPMENT SAMPLE
# ==============================================================================

# Mode imputation helpers
learn_modes <- function(data, vars) {
  setNames(
    lapply(vars, function(v) {
      if (v %in% names(data))
        names(sort(table(na.omit(data[[v]])), decreasing = TRUE))[1]
    }),
    vars
  )
}

apply_modes <- function(data, modes) {
  for (var in names(modes)) {
    if (!var %in% names(data) || is.null(modes[[var]])) next
    missing_idx <- is.na(data[[var]])
    if (!any(missing_idx)) next
    mode_val <- modes[[var]]
    if (is.factor(data[[var]])) {
      # Only impute if mode_val is already a valid level —
      # never add new levels to a validation factor
      if (mode_val %in% levels(data[[var]])) {
        data[[var]][missing_idx] <- mode_val
      } else {
        warning("Mode value '", mode_val, "' for variable '", var,
                "' is not a valid level in this dataset. Skipping imputation.")
      }
    } else {
      data[[var]][missing_idx] <- mode_val
    }
  }
  data
}

categorical_vars <- c(
  "Group", "Gender", "Relationship", "School", "Employment",
  "current_psychotherapy", "current_antidepressant", "Dysthymia"
)
continuous_vars <- c("Age", "QOL_phys", "QOL_psych", "PHQ_1")

interim        <- preProcess(deprexis_select_PHQ9[, continuous_vars, drop = FALSE],
                              method = "medianImpute")
deprexis_d     <- deprexis_select_PHQ9
deprexis_d[, continuous_vars] <- predict(interim,
                                          deprexis_select_PHQ9[, continuous_vars,
                                                                drop = FALSE])
mode_values    <- learn_modes(deprexis_d, categorical_vars)
deprexis_d     <- apply_modes(deprexis_d, mode_values)

# Ensure numeric Group for glmnet / causal forest
deprexis_d <- deprexis_d %>%
  mutate(Group = ifelse(Group == "intervention", 1, 0))

if ("MajorDepression" %in% names(deprexis_d))
  deprexis_d$MajorDepression <- factor(deprexis_d$MajorDepression,
                                        levels = c("no", "yes"))

message("✓ Development sample preprocessed. Rows: ", nrow(deprexis_d),
        "  Missing: ", sum(is.na(deprexis_d)))

# ==============================================================================
# SECTION 4  — FIT FINAL FULL-SAMPLE MODELS
# ==============================================================================

x <- model.matrix(f, data = deprexis_d)[, -1]
y <- deprexis_d[, outcome]

full_model_deprexis     <- glmnet(x, y, alpha = optalpha,  lambda = optlambda,
                                   family = "gaussian")
full_model_deprexis_OLS <- glmnet(x, y, alpha = 0,         lambda = 0,
                                   family = "gaussian")

message("✓ Final ENR model fitted.  Non-zero coefficients: ",
        sum(coef(full_model_deprexis) != 0) - 1)

# ==============================================================================
# SECTION 5  — BOOTSTRAP COEFFICIENT STABILITY  (1 000 resamples)
# ==============================================================================

n_bootstrap  <- 1000
selected_vars <- matrix(0, nrow = ncol(x) + 1, ncol = n_bootstrap)

set.seed(12345)
for (k in seq_len(n_bootstrap)) {

  idx        <- sample.int(nrow(x), replace = TRUE)
  data_boot  <- deprexis_select_PHQ9[idx, ]

  # Re-fit median imputation within bootstrap sample
  imp_boot   <- preProcess(data_boot[, continuous_vars, drop = FALSE],
                            method = "medianImpute")
  data_boot[, continuous_vars] <- predict(imp_boot,
                                          data_boot[, continuous_vars,
                                                     drop = FALSE])
  data_boot  <- apply_modes(data_boot, mode_values)
  data_boot  <- data_boot %>%
    mutate(Group = ifelse(Group == "intervention", 1, 0))
  if ("MajorDepression" %in% names(data_boot))
    data_boot$MajorDepression <- factor(data_boot$MajorDepression,
                                         levels = c("no", "yes"))

  x_boot <- model.matrix(f, data = data_boot)[, -1]
  y_boot <- as.numeric(data_boot[, outcome])

  boot_model          <- glmnet(x_boot, y_boot, alpha = optalpha,
                                 lambda = optlambda, family = "gaussian")
  selected_vars[, k]  <- as.numeric(coef(boot_model))
  rownames(selected_vars) <- rownames(coef(boot_model))
}

mean_coefs <- rowMeans(selected_vars)
freq       <- rowSums(selected_vars != 0) / n_bootstrap
consistency_neg <- rowSums(selected_vars < 0) /
  pmax(rowSums(selected_vars < 0) + rowSums(selected_vars > 0), 1) * 100
consistency_pos <- 100 - consistency_neg

pooled_results <- data.frame(
  Mean_Coefficients = round(mean_coefs, 4),
  lb_conf = round(mean_coefs - apply(selected_vars, 1,
                                      function(x) 1.96 * sd(x, na.rm = TRUE) /
                                        sqrt(sum(!is.na(x)))), 4),
  ub_conf = round(mean_coefs + apply(selected_vars, 1,
                                      function(x) 1.96 * sd(x, na.rm = TRUE) /
                                        sqrt(sum(!is.na(x)))), 4),
  frequency       = round(freq, 4),
  consistency_neg = round(consistency_neg, 4),
  consistency_pos = round(consistency_pos, 4)
)

write.csv(pooled_results, "bootstrap_ENR.csv", row.names = TRUE)
message("✓ Bootstrap complete. Results saved to bootstrap_ENR.csv")

# ==============================================================================
# SECTION 6  — PREPROCESS VALIDATION DATASETS
# (using imputation parameters learned from development data)
# ==============================================================================

preprocess_validation <- function(raw_val, median_imputer, mode_vals,
                                   cat_vars, cont_vars) {
  out <- raw_val
  out[, cont_vars] <- predict(median_imputer,
                               raw_val[, cont_vars, drop = FALSE])
  out <- apply_modes(out, mode_vals)
  # Drop any unused levels that may have accumulated during data wrangling
  out <- droplevels(out)
  if ("MajorDepression" %in% names(out))
    out$MajorDepression <- factor(out$MajorDepression, levels = c("no", "yes"))
  n_miss <- sum(is.na(out))
  if (n_miss > 0)
    warning("Remaining missing values after imputation: ", n_miss)
  out
}

selfapy_imputed  <- preprocess_validation(
  HRSD_selfapy_QIDSs, interim, mode_values,
  categorical_vars, continuous_vars
)
eficasy_imputed  <- preprocess_validation(
  eficasy_select, interim, mode_values,
  categorical_vars, continuous_vars
)
moritz_imputed   <- preprocess_validation(
  moritz_select, interim, mode_values,
  categorical_vars, continuous_vars
)

message("✓ Validation datasets preprocessed:")
message("  Selfapy N=",       nrow(selfapy_imputed),
        "  edupression N=",   nrow(eficasy_imputed),
        "  Deprexis-Val N=",  nrow(moritz_imputed))

# ==============================================================================
# SECTION 7  — OLS CV METRICS  (for supplementary table)
# Rerun cross-validation with alpha=0, lambda=0 to get OLS fold-level metrics.
# (Alternatively load pre-saved results: load("combined_metrics_OLS.RData"))
# ==============================================================================

# NOTE: The OLS CV is identical in structure to Section 2 above but with
# alpha = 0, lambda = 0 fixed (no tuning grid needed).  Produce
# ols_cv_rmse_mean/sd, ols_cv_mae_mean/sd, ols_cv_r2_mean/sd,
# ols_cv_intercept_mean/sd, ols_cv_slope_mean/sd.
#
# For brevity the loop is not repeated here; set up analogously or
# load from saved .RData.

message("✓ 02_model_development.R complete")
sessionInfo()
