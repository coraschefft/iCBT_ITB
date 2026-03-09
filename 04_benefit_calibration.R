library(devtools)
install_github("esm-ispm-unibe-ch/predieval")
install.packages("predieval")
library(predieval)
library(stringr)
library(caret)
library(glmnet)
library(grf)
library(patchwork)
library(Matching)
rm(list=ls()) # empty memory
set.seed(42) #the answer to life, the universe and everything 
# Calculate pooled SD from your data
calculate_pooled_sd <- function(dat) {
  Y1 <- dat$Y[dat$t == 1]
  Y0 <- dat$Y[dat$t == 0]
  n1 <- length(Y1)
  n0 <- length(Y0)
  s1 <- sd(Y1)
  s0 <- sd(Y0)
  
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
  return(pooled_sd)
}

#### simulate data ----

# was muss in dem data frame enthalten sein?
i<-7
outcome<-noquote(outcomes[i])
print(outcome)
formula_string <- paste(outcome, "~ Group*(Gender+Age+School+Employment+Relationship+
current_psychotherapy+current_antidepressant+MajorDepression+Dysthymia+
                        QOL_psych+QOL_phys+PHQ_1)")

f <- as.formula(formula_string)

formula_string_rf <- paste("~ Gender+Age+School+Employment+Dysthymia+Relationship+
current_psychotherapy+current_antidepressant+MajorDepression+
                        QOL_psych+QOL_phys+PHQ_1")

frf<-as.formula(formula_string_rf) #for random forests
#optimal_lambda<-median(combined_metrics$lambda)
#optimal_alpha<-median(combined_metrics$alpha)


is.numeric(deprexis_d$Group)
deprexis_d<-deprexis_d%>%mutate(Group=ifelse(Group=="intervention",1,0))
dat_efth<-deprexis_d

is.numeric(deprexis_d$Group) #must be numeric an 1 and 0 for this
##r### 10-fold CV repeated 100 times ------
k.folds<-10
repeats<-100
dat.CV<-list()
set.seed(1234)  # set global seed once
seeds <- sample.int(1e6, repeats)  # generate vector of seeds

for (k in 1:repeats) {
  set.seed(seeds[k])  # s
  flds <- createFolds(1:nrow(dat_efth), k = k.folds, list = TRUE, returnTrain = FALSE) #caret create folds
  dat.out.CV<-list() #liste hält outer fold predictions per repeat, das sollten dann 100 sein
  for (j in 1:k.folds){
    dat.in.CV<-dat_efth[-flds[[j]],] #inneres fold
    dat.out.CV[[j]]=dat_efth[flds[[j]],] #äußeres fold
    dat1<-dat.out.CV[[j]]; dat1$Group=1 #<-factor(1, levels=c("1","0")) äußeres fold, um Prädiktionen für t=1 zu erstellen
    dat0<-dat.out.CV[[j]]; dat0$Group=0 #<-factor(0, levels=c("1","0")) äußeres fold, um Prädiktionen für t=0 zu erstellen
    x_in<-model.matrix(f,data=dat.in.CV)[,-1]
    y_in<-dat.in.CV[,outcome]
    x_out1<-model.matrix(f,data=dat1)[,-1]
    x_out0<-model.matrix(f,data=dat0)[,-1]
    
    
   
   
    m1<-glmnet(x_in, y_in, alpha=optalpha,lambda = optlambda, family="gaussian")
    dat.out.CV[[j]]$m1.CV.treat1=predict(m1, newx = x_out1, s=optlambda, type="link") #hält die prediction für dat1
    dat.out.CV[[j]]$m1.CV.treat0=predict(m1, newx = x_out0, s=optlambda, type="link") 
    
    m2<-glmnet(x_in, y_in, alpha = 0, lambda = 0, family="gaussian")
    
    dat.out.CV[[j]]$m2.CV.treat1=predict(m2, newx = x_out1, s=0, type="link") #hält die prediction für dat1
    dat.out.CV[[j]]$m2.CV.treat0=predict(m2, newx = x_out0, s=0, type="link") 
    
    W <- as.integer(dat.in.CV$Group) 
    table(W)
    Y <- as.numeric(dat.in.CV[,outcome])
    X <- model.matrix(frf,data=dat.in.CV)[,-1]
    X_val<- model.matrix(frf,data=dat.out.CV[[j]])[,-1]
    
 #   Y.hat.mod<-regression_forest(X=X,
  #                               Y=Y,
   #                              num.trees=200,
   #                              ci.group.size = 1,
    #                             tune.parameters = "all")
    
  #  Y.hat.rf<-Y.hat.mod$predictions
  #  W.hat.mod<-regression_forest(X=X,
 #                                Y=W,
 #                                num.trees=200,
  #                               ci.group.size = 1,
  #                               tune.parameters = "all")
  #  W.hat.rf<-W.hat.mod$predictions
    
    m3 <- causal_forest(
      X = X,
      Y = Y,
      W = W,
      
    #  num.trees = 2000,
     # min.node.size=14,
    #  honesty.fraction=0.67,
    #  alpha=0.17,
    #  W.hat = W.hat.rf,
    #  Y.hat =Y.hat.rf,
       mtry = 4,
    #  sample.fraction = 0.11,
  #  imbalance.penalty = 0.333,
      seed = 1,
       #ci.group.size = 5,
    )
    dat.out.CV[[j]]$m3.benefit=predict(m3, newdata = X_val)$predictions
    }
 
  dat.CV[[k]]<-dat.out.CV[[1]] 
  for (j in 2:k.folds){dat.CV[[k]]<-rbind(dat.CV[[k]],dat.out.CV[[j]] )}#Combines all 10 folds for repeat k into a single data frame.
}

summary(unlist(lapply(dat.CV, function(df) sum(is.na(df$m3.benefit)))))
cate_all <- unlist(lapply(dat.CV, function(df) df$m3.benefit))
summary(cate_all)
hist(cate_all, breaks = 40, col = "lightblue", main = "Distribution of predicted CATEs", xlab = "CATE")
y1_preds <- unlist(lapply(dat.CV, function(df) df$m1.CV.treat1))
y0_preds <- unlist(lapply(dat.CV, function(df) df$m1.CV.treat0))
y1m2_preds <- unlist(lapply(dat.CV, function(df) df$m2.CV.treat1))
y0m2_preds <- unlist(lapply(dat.CV, function(df) df$m2.CV.treat0))
hist(y1_preds - y0_preds, breaks = 50, main = "Elastic Net Estimated Benefit (Treat1 - Treat0)", xlab = "Predicted benefit")
range(y1m2_preds - y0m2_preds)
range(y1_preds - y0_preds)
range(cate_all)

all_preds <- do.call(rbind, dat.CV)




for(i in 1:repeats){
  row.names(dat.CV[[i]])=as.numeric(str_remove(row.names(dat.CV[[i]]), "dat."))
  dat.CV[[i]]<-dat.CV[[i]][order(as.numeric(row.names(dat.CV[[i]]))),]}

dat.CV.all<-dat.CV[[1]]
for(k in 2:repeats){
  dat.CV.all[,c("m1.CV.treat1", "m1.CV.treat0",
                "m2.CV.treat1", "m2.CV.treat0","m3.benefit"
                )]<-
    dat.CV.all[,c("m1.CV.treat1", "m1.CV.treat0",
                  "m2.CV.treat1", "m2.CV.treat0","m3.benefit"
                  )] +
    dat.CV[[k]][,c("m1.CV.treat1", "m1.CV.treat0",
                   "m2.CV.treat1", "m2.CV.treat0","m3.benefit"
                   )] 
}

dat.CV.all[,c("m1.CV.treat1", "m1.CV.treat0",
              "m2.CV.treat1", "m2.CV.treat0","m3.benefit"
              )]<-
  dat.CV.all[,c("m1.CV.treat1", "m1.CV.treat0",
                "m2.CV.treat1", "m2.CV.treat0","m3.benefit"
                )]/repeats

all_preds$bin <- ntile(dat.CV.all$m3.benefit, 5)

calib_summary <- all_preds |>
  group_by(bin) |>
  dplyr::summarise(
    mean_pred = mean(m3.benefit),
    obs_effect = mean(PHQ_change[Group == 1 & bin == cur_group_id()]) -
      mean(PHQ_change[Group == 0 & bin == cur_group_id()])
  )



#bencalibr# Model M1
bencalibr(Ngroups=10, data=dat.CV.all, y.observed=PHQ_change,
          predicted.treat.1 = m1.CV.treat1,
          predicted.treat.0 = m1.CV.treat0, treat=Group, 
          , smoothing.function="lm")
bencalibr(Ngroups=18, data=dat.CV.all, y.observed=PHQ_change,
          predicted.treat.1 = m1.CV.treat1,
          predicted.treat.0 = m1.CV.treat0, treat=Group, 
          , smoothing.function="lm")

bencalibr(Ngroups=18, data=dat.CV.all, y.observed=PHQ_change,
          predicted.treat.1 = m2.CV.treat1,
          predicted.treat.0 = m2.CV.treat0, treat=Group, 
          , smoothing.function="lm")
# Model M2




#

X=as.data.frame(model.matrix(frf,data=dat.CV.all)[,-1])
Y=dat.CV.all[,outcome]
predicted.treat.1 = dat.CV.all$m1.CV.treat1
predicted.treat.0 = dat.CV.all$m1.CV.treat0
predicted.treat2.1 = dat.CV.all$m2.CV.treat1
predicted.treat2.0 = dat.CV.all$m2.CV.treat0
predicted.benefit=dat.CV.all$m3.benefit
treat=dat.CV.all$Group

  ################ continuous data --------------

   # dat2<-cbind(X,"Y"=c(Y),"t"=treat,"benefit"=predicted.treat.1-predicted.treat.0) 
    dat1_DS<-cbind(X,"Y"=c(Y),"t"=treat,"benefit"=predicted.treat.1-predicted.treat.0,"pred1"=predicted.treat.1,"pred0"=predicted.treat.0) 
    dat2_DS<-cbind(X,"Y"=c(Y),"t"=treat,"benefit"=predicted.treat2.1-predicted.treat2.0,"pred1"=predicted.treat2.1,"pred0"=predicted.treat2.0) 
  dat3_DS<-cbind(X,"Y"=c(Y),"t"=treat,"benefit"=predicted.benefit) 
    # calibration - mean bias
    mean.bias<-mean(dat3_DS$Y[dat3_DS$t==1])-mean(dat3_DS$Y[dat3_DS$t==0])-mean(dat3_DS$benefit) #The estimation of mean bias is straightforward for a continuous outcome, assuming randomization: we compare the observed mean benefit at the arm level (mean observed outcome in treatment minus control) to the mean-predicted benefit.


nrow(dat1_DS)

# ==============================================================================
###### EXTERNAL VALIDATION OF TREATMENT BENEFIT PREDICTION#######
# UPDATED: Consistent naming (edupression not Eficasy)
# ==============================================================================

# Set outcome


# ==============================================================================
# DEVELOPMENT SAMPLE — CALIBRATION + PREDICTION RANGE LOOP
# ==============================================================================

method_datasets <- list(
  enr = dat1_DS,
  ols = dat2_DS,
  cf  = dat3_DS
)

performances_dev_all   <- list()
calibration_dev_all    <- list()
prediction_summary_dev <- list()

for (method_name in names(method_datasets)) {
  
  cat("\n=== Processing method:", toupper(method_name), "===\n")
  
  datperf <- method_datasets[[method_name]]
  
  # ---- Prediction range summary (per method) --------------------------------
  benefit_vals <- datperf$benefit
  prediction_summary_dev[[method_name]] <- data.frame(
    Dataset = "Deprexis",
    Method  = toupper(method_name),
    Min     = min(benefit_vals,            na.rm = TRUE),
    Q25     = quantile(benefit_vals, 0.25, na.rm = TRUE),
    Median  = median(benefit_vals,         na.rm = TRUE),
    Mean    = mean(benefit_vals,           na.rm = TRUE),
    Q75     = quantile(benefit_vals, 0.75, na.rm = TRUE),
    Max     = max(benefit_vals,            na.rm = TRUE),
    SD      = sd(benefit_vals,             na.rm = TRUE)
  )
  
  # ---- Calibration bins loop ------------------------------------------------
  Ngroups <- c(3, 4, 5, 11, 18)
  performances_method    <- list()
  calibration_bins_method <- list()
  
  for (k in seq_along(Ngroups)) {
    
    n_bins <- Ngroups[k]
    cat("  Processing", n_bins, "bins...\n")
    
    dat_grouped <- datperf %>%
      mutate(benefit_group = ntile(benefit, n_bins))
    
    effect_sizes <- list()
    
    for (g in sort(unique(dat_grouped$benefit_group))) {
      
      df_sub <- dat_grouped %>% filter(benefit_group == g)
      
      m1 <- mean(df_sub$Y[df_sub$t == 1], na.rm = TRUE)
      m0 <- mean(df_sub$Y[df_sub$t == 0], na.rm = TRUE)
      s1 <- sd(df_sub$Y[df_sub$t == 1],   na.rm = TRUE)
      s0 <- sd(df_sub$Y[df_sub$t == 0],   na.rm = TRUE)
      n1 <- sum(df_sub$t == 1)
      n0 <- sum(df_sub$t == 0)
      
      se_diff   <- sqrt((s1^2 / n1) + (s0^2 / n0))
     # pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
     # d         <- (m1 - m0) / pooled_sd
     # d_pred    <- mean(df_sub$benefit) / pooled_sd
     # se_d      <- sqrt((n1 + n0) / (n1 * n0) + (d^2) / (2 * (n1 + n0)))
      
      effect_sizes[[as.character(g)]] <- data.frame(
        benefit_group          = g,
        n_total                = n1 + n0,
        n_treat                = n1,
        n_control              = n0,
        mean_predicted_benefit = mean(df_sub$benefit, na.rm = TRUE),
        sd_predicted_benefit   = sd(df_sub$benefit,   na.rm = TRUE),
     #   d_pred                 = d_pred,
        abs_diff_means         = (m1 - m0),
        lower_ci_means         = (m1 - m0) - se_diff,
        upper_ci_means         = (m1 - m0) + se_diff,
       pooled_sd              = pooled_sd#,
     #   cohen_d                = d,
     #   lower_ci_cohen         = d - se_d,
     #   upper_ci_cohen         = d + se_d
      )
    }
    
    calibration_summary <- bind_rows(effect_sizes)
    
    calibration_bins_method[[paste0("bins_", n_bins)]] <- calibration_summary %>%
      mutate(Dataset = "Deprexis", Method = toupper(method_name), N_bins = n_bins)
    
  #  lm_cohen <- lm(cohen_d      ~ d_pred,                 data = calibration_summary)
    lm_abs   <- lm(abs_diff_means ~ mean_predicted_benefit, data = calibration_summary)
    
  #  cohen_summary <- summary(lm_cohen)
    abs_summary   <- summary(lm_abs)
    
  #  intercept_se_cohen <- cohen_summary$coefficients["(Intercept)",          "Std. Error"]
   # slope_se_cohen     <- cohen_summary$coefficients["d_pred",               "Std. Error"]
    intercept_se_abs   <- abs_summary$coefficients["(Intercept)",            "Std. Error"]
    slope_se_abs       <- abs_summary$coefficients["mean_predicted_benefit", "Std. Error"]
    
    pooled_sd_dataset     <- mean(calibration_summary$pooled_sd, na.rm = TRUE)
   # RMSE_cohen            <- sqrt(mean((calibration_summary$d_pred - calibration_summary$cohen_d)^2))
    RMSE_abs              <- sqrt(mean((calibration_summary$mean_predicted_benefit -
                                          calibration_summary$abs_diff_means)^2))
    RMSE_abs_standardized <- RMSE_abs / pooled_sd_dataset
    
    performances_method[[paste0("bins_", n_bins)]] <- data.frame(
      Dataset  = "Deprexis",
      Method   = toupper(method_name),
      N_bins   = n_bins,
   #   RMSE_cohen            = RMSE_cohen,
    #  intercept_cohen       = coef(lm_cohen)[1],
  #    intercept_se_cohen    = intercept_se_cohen,
  #    intercept_lower_cohen = coef(lm_cohen)[1] - 1.96 * intercept_se_cohen,
  #    intercept_upper_cohen = coef(lm_cohen)[1] + 1.96 * intercept_se_cohen,
  #    slope_cohen           = coef(lm_cohen)[2],
  #    slope_se_cohen        = slope_se_cohen,
  #    slope_lower_cohen     = coef(lm_cohen)[2] - 1.96 * slope_se_cohen,
  #    slope_upper_cohen     = coef(lm_cohen)[2] + 1.96 * slope_se_cohen,
  #    r2_cohen              = cohen_summary$r.squared,
      RMSE_abs              = RMSE_abs,
      RMSE_abs_standardized = RMSE_abs_standardized,
      intercept_abs         = coef(lm_abs)[1],
      intercept_se_abs      = intercept_se_abs,
      intercept_lower_abs   = coef(lm_abs)[1] - 1.96 * intercept_se_abs,
      intercept_upper_abs   = coef(lm_abs)[1] + 1.96 * intercept_se_abs,
      slope_abs             = coef(lm_abs)[2],
      slope_se_abs          = slope_se_abs,
      slope_lower_abs       = coef(lm_abs)[2] - 1.96 * slope_se_abs,
      slope_upper_abs       = coef(lm_abs)[2] + 1.96 * slope_se_abs,
      r2_abs                = abs_summary$r.squared,
      pooled_sd             = pooled_sd_dataset
    )
  }
  
  performances_dev_all[[method_name]]    <- bind_rows(performances_method)
  calibration_dev_all[[method_name]]     <- bind_rows(calibration_bins_method)
}

performances_dev_final <- bind_rows(performances_dev_all)
calibration_dev_final  <- bind_rows(calibration_dev_all)

# Deposit dev prediction summaries into benefit_prediction_results
# (initialise if not yet created)
if (!exists("benefit_prediction_results")) benefit_prediction_results <- list()
benefit_prediction_results[["Deprexis"]]$prediction_summary <- bind_rows(prediction_summary_dev)

cat("\n=== Development Sample Results Summary ===\n")
print(performances_dev_final %>%
        filter(N_bins %in% c(4, 11)) %>%
        select(Method, N_bins, RMSE_abs_standardized, intercept_abs, slope_abs,  r2_abs))

# ==============================================================================
# FIT MODELS FOR EXTERNAL VALIDATION
# ==============================================================================

outcome <- "PHQ_change"

if (!is.numeric(deprexis_d$Group)) {
  deprexis_d <- deprexis_d %>%
    mutate(Group = ifelse(Group == "intervention", 1, 0))
}

W <- as.integer(deprexis_d$Group)
Y <- as.numeric(deprexis_d[, outcome])
X <- model.matrix(frf, data = deprexis_d)[, -1]

n_predictors <- ncol(X)
mtry         <- floor(sqrt(n_predictors))
m3           <- causal_forest(X, Y, W, num.trees = 2000, mtry = mtry, seed = 1)
message("✓ Causal forest fitted on development data")

# ==============================================================================
# EXTERNAL VALIDATION — LOOP OVER DATASETS
# ==============================================================================

validation_datasets <- list(
  "edupression" = list(data = eficasy_imputed,  formula = f, formula_rf = frf,
                       model_enr = full_model_deprexis, model_ols = full_model_deprexis_OLS,
                       model_cf  = m3),
  "Selfapy"     = list(data = selfapy_imputed,  formula = f, formula_rf = frf,
                       model_enr = full_model_deprexis, model_ols = full_model_deprexis_OLS,
                       model_cf  = m3),
  "Moritz"      = list(data = moritz_imputed,   formula = f, formula_rf = frf,
                       model_enr = full_model_deprexis, model_ols = full_model_deprexis_OLS,
                       model_cf  = m3)
)

dat_storage <- list()
benefit_prediction_results <- list()
for (dataset_name in names(validation_datasets)) {
  
  message("\n========================================")
  message("Processing: ", dataset_name)
  message("========================================")
  
  val_data      <- as.data.frame(validation_datasets[[dataset_name]]$data)
  formula_use   <- validation_datasets[[dataset_name]]$formula
  formula_rf_use <- validation_datasets[[dataset_name]]$formula_rf
  model_enr     <- validation_datasets[[dataset_name]]$model_enr
  model_ols     <- validation_datasets[[dataset_name]]$model_ols
  model_cf      <- validation_datasets[[dataset_name]]$model_cf
  
  if (!is.numeric(val_data$Group))
    val_data <- val_data %>% mutate(Group = ifelse(Group == "intervention", 1, 0))
  
  if ("MajorDepression" %in% names(val_data)) {
    val_data$MajorDepression <- factor(val_data$MajorDepression, levels = c("no", "yes"))
  }
  
  X_val    <- model.matrix(formula_rf_use, data = val_data)[, -1]
  X_val_t1 <- model.matrix(formula_use,   data = transform(val_data, Group = 1))[, -1]
  X_val_t0 <- model.matrix(formula_use,   data = transform(val_data, Group = 0))[, -1]
  
  Y_val     <- val_data[[outcome]]
  treat_val <- val_data$Group
  
  # Predict benefit
  benefit_enr <- as.numeric(predict(model_enr, newx = X_val_t1, s = optlambda) -
                              predict(model_enr, newx = X_val_t0, s = optlambda))
  benefit_ols <- as.numeric(predict(model_ols, newx = X_val_t1, s = 0) -
                              predict(model_ols, newx = X_val_t0, s = 0))
  benefit_cf  <- predict(model_cf, newdata = X_val)$predictions
  
  dat_enr <- data.frame(X_val, Y = Y_val, t = treat_val, benefit = benefit_enr, check.names = FALSE)
  dat_ols <- data.frame(X_val, Y = Y_val, t = treat_val, benefit = benefit_ols, check.names = FALSE)
  dat_cf  <- data.frame(X_val, Y = Y_val, t = treat_val, benefit = benefit_cf,  check.names = FALSE)
  
  dat_storage[[dataset_name]] <- list(enr = dat_enr, ols = dat_ols, cf = dat_cf)
  
  # ---- Prediction range summary --------------------------------------------
  prediction_summary <- data.frame(
    Dataset = dataset_name,
    Method  = c("ENR", "OLS", "CF"),
    Min     = c(min(benefit_enr),             min(benefit_ols),             min(benefit_cf)),
    Q25     = c(quantile(benefit_enr, 0.25),  quantile(benefit_ols, 0.25),  quantile(benefit_cf, 0.25)),
    Median  = c(median(benefit_enr),           median(benefit_ols),           median(benefit_cf)),
    Mean    = c(mean(benefit_enr),             mean(benefit_ols),             mean(benefit_cf)),
    Q75     = c(quantile(benefit_enr, 0.75),  quantile(benefit_ols, 0.75),  quantile(benefit_cf, 0.75)),
    Max     = c(max(benefit_enr),             max(benefit_ols),             max(benefit_cf)),
    SD      = c(sd(benefit_enr),              sd(benefit_ols),              sd(benefit_cf))
  )
  
  message("\nPrediction ranges:")
  print(prediction_summary %>% mutate(across(where(is.numeric), ~ round(., 2))))
  
  # ---- Calibration bins loop -----------------------------------------------
  methods <- list("ENR" = dat_enr, "OLS" = dat_ols, "CF" = dat_cf)
  
  calibration_list <- list()
  performance_list <- list()
  
  for (method_name in names(methods)) {
    
    datperf      <- methods[[method_name]]
    Ngroups_vec  <- c(2, 3, 4, 5)
    
    for (k in Ngroups_vec) {
      
      dat_grouped <- datperf %>% mutate(benefit_group = ntile(benefit, k))
      effect_sizes <- list()
      
      for (g in sort(unique(dat_grouped$benefit_group))) {
        
        df_sub <- dat_grouped %>% filter(benefit_group == g)
        
        m1 <- mean(df_sub$Y[df_sub$t == 1], na.rm = TRUE)
        m0 <- mean(df_sub$Y[df_sub$t == 0], na.rm = TRUE)
        s1 <- sd(df_sub$Y[df_sub$t == 1],   na.rm = TRUE)
        s0 <- sd(df_sub$Y[df_sub$t == 0],   na.rm = TRUE)
        n1 <- sum(df_sub$t == 1)
        n0 <- sum(df_sub$t == 0)
        
        se_diff   <- sqrt((s1^2 / n1) + (s0^2 / n0))
        pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
      #  d         <- (m1 - m0) / pooled_sd
      #  d_pred    <- mean(df_sub$benefit) / pooled_sd
      #  se_d      <- sqrt((n1 + n0) / (n1 * n0) + (d^2) / (2 * (n1 + n0)))
        
        effect_sizes[[as.character(g)]] <- data.frame(
          benefit_group          = g,
          n_total                = n1 + n0,
          n_treat                = n1,
          n_control              = n0,
          mean_predicted_benefit = mean(df_sub$benefit, na.rm = TRUE),
          sd_predicted_benefit   = sd(df_sub$benefit,   na.rm = TRUE),
      #    d_pred                 = d_pred,
          abs_diff_means         = (m1 - m0),
          lower_ci_means         = (m1 - m0) - se_diff,
          upper_ci_means         = (m1 - m0) + se_diff,
          pooled_sd              = pooled_sd#,
       #   cohen_d                = d,
       #   lower_ci_cohen         = d - se_d,
        #  upper_ci_cohen         = d + se_d
        )
      }
      
      calibration_summary <- bind_rows(effect_sizes)
      
   
      RMSE_abs   <- sqrt(mean((calibration_summary$mean_predicted_benefit -
                                 calibration_summary$abs_diff_means)^2))
      
   
      lm_abs   <- lm(abs_diff_means ~ mean_predicted_benefit, data = calibration_summary)
      

      abs_summary   <- summary(lm_abs)
      
     
      intercept_se_abs   <- abs_summary$coefficients["(Intercept)",            "Std. Error"]
      slope_se_abs       <- abs_summary$coefficients["mean_predicted_benefit", "Std. Error"]
      
      pooled_sd_dataset     <- mean(calibration_summary$pooled_sd, na.rm = TRUE)
      RMSE_abs_standardized <- RMSE_abs / pooled_sd_dataset
      
      key <- paste0(method_name, "_bins", k)
      
      calibration_list[[key]] <- calibration_summary %>%
        mutate(Dataset = dataset_name, Method = method_name, N_bins = k)
      
      performance_list[[key]] <- data.frame(
        Dataset  = dataset_name,
        Method   = method_name,
        N_bins   = k,
       
        RMSE_abs              = RMSE_abs,
        RMSE_abs_standardized = RMSE_abs_standardized,
        intercept_abs         = coef(lm_abs)[1],
        intercept_se_abs      = intercept_se_abs,
        intercept_lower_abs   = coef(lm_abs)[1] - 1.96 * intercept_se_abs,
        intercept_upper_abs   = coef(lm_abs)[1] + 1.96 * intercept_se_abs,
        slope_abs             = coef(lm_abs)[2],
        slope_se_abs          = slope_se_abs,
        slope_lower_abs       = coef(lm_abs)[2] - 1.96 * slope_se_abs,
        slope_upper_abs       = coef(lm_abs)[2] + 1.96 * slope_se_abs,
        r2_abs                = abs_summary$r.squared,
        pooled_sd             = pooled_sd_dataset
      )
    }
  }
  
  benefit_prediction_results[[dataset_name]] <- list(
    prediction_summary = prediction_summary,
    calibration        = bind_rows(calibration_list),
    performance        = bind_rows(performance_list)
  )
  
  message("\n✓ Completed analysis for ", dataset_name)
}

# ==============================================================================
# COMBINE ALL PREDICTION SUMMARIES (dev + validation)
# ==============================================================================

all_prediction_summaries <- bind_rows(
  lapply(names(benefit_prediction_results), function(ds) {
    benefit_prediction_results[[ds]]$prediction_summary %>%
      mutate(Dataset = ds)  # ensure Dataset col is set correctly
  })
) %>%
  # Consistent factor ordering for display
  mutate(
    Sample = case_when(
      Dataset == "Deprexis"    ~ "Deprexis (Dev)",
      Dataset == "edupression" ~ "edupression",
      Dataset == "Selfapy"     ~ "Selfapy",
      Dataset == "Moritz"      ~ "Deprexis (Val)"
    ),
    Sample = factor(Sample, levels = c("Deprexis (Dev)", "edupression",
                                       "Selfapy", "Deprexis (Val)")),
    Method = factor(Method, levels = c("ENR", "OLS", "CF"))
  ) %>%
  arrange(Sample, Method) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

# ==============================================================================
# COMBINE CALIBRATION + PERFORMANCE ACROSS ALL SAMPLES
# ==============================================================================

all_calibration_combined <- bind_rows(
  calibration_dev_final,
  bind_rows(lapply(names(benefit_prediction_results)[names(benefit_prediction_results) != "Deprexis"],
                   function(ds) benefit_prediction_results[[ds]]$calibration))
)

all_performance_combined <- bind_rows(
  performances_dev_final,
  bind_rows(lapply(names(benefit_prediction_results)[names(benefit_prediction_results) != "Deprexis"],
                   function(ds) benefit_prediction_results[[ds]]$performance))
) %>%
  mutate(
    Sample = case_when(
      Dataset == "Deprexis"    ~ "Deprexis (Dev)",
      Dataset == "edupression" ~ "edupression",
      Dataset == "Selfapy"     ~ "Selfapy",
      Dataset == "Moritz"      ~ "Deprexis (Val)"
    ),
    Sample = factor(Sample, levels = c("Deprexis (Dev)", "edupression",
                                       "Selfapy", "Deprexis (Val)")),
    Method = factor(Method, levels = c("ENR", "OLS", "CF"))
  )

all_calibration_combined <- all_calibration_combined %>%
  mutate(
    Sample = case_when(
      Dataset == "Deprexis"    ~ "Deprexis (Dev)",
      Dataset == "edupression" ~ "edupression",
      Dataset == "Selfapy"     ~ "Selfapy",
      Dataset == "Moritz"      ~ "Deprexis (Val)"
    ),
    Sample = factor(Sample, levels = c("Deprexis (Dev)", "edupression",
                                       "Selfapy", "Deprexis (Val)")),
    Method = factor(Method, levels = c("ENR", "OLS", "CF"))
  )

reduced_table_benefit <- all_performance_combined %>%
  select(Dataset, Method, N_bins,
         RMSE_abs, RMSE_abs_standardized,
         intercept_abs, intercept_se_abs, intercept_lower_abs, intercept_upper_abs,
         slope_abs, slope_se_abs, slope_lower_abs, slope_upper_abs,
         r2_abs, pooled_sd) %>%
  filter(N_bins > 2)


# ==============================================================================
# PRINT SUMMARY TABLES
# ==============================================================================

message("\n\n========================================")
message("PREDICTION RANGE SUMMARY — ALL SAMPLES")
message("========================================\n")
print(all_prediction_summaries %>%
        select(Sample, Method, Min, Q25, Median, Mean, Q75, Max, SD))

message("\n\n========================================")
message("CALIBRATION PERFORMANCE SUMMARY")
message("========================================\n")
print(all_performance_combined %>%
         filter((Dataset == "Deprexis" & N_bins %in% c(4, 11, 18)) |
                                  (Dataset != "Deprexis" & N_bins %in% c(3, 4, 5))) %>%
          select(Sample, Method, N_bins, RMSE_abs_standardized, r2_abs, slope_abs, intercept_abs))

message("\n\n========================================")
message("POPULATION BENEFIT")
message("========================================\n")

# ==============================================================================
# POPULATION BENEFIT  (cutoffs: d >= 0.24 and d >= 0.50 x pooled SD)
# ==============================================================================
# PB1 formula (Lesko et al.):
#   Groups: G1 = treated & benefit > cutoff,  G2 = control & benefit > cutoff
#           G3 = treated & benefit <= cutoff,  G4 = control & benefit <= cutoff
#   PB1 = mean(Y_G4) - mean(Y_G3) + mean(Y_G1) * (1/n(G1+G4) - 1/n(G1+G3))
# Uses ENR benefit predictions throughout (pre-specified primary model).

compute_pb1 <- function(dat, cutoff, dataset_name, method_name) {
  g1 <- dat[dat$benefit >  cutoff & dat$t == 1, ]
  g2 <- dat[dat$benefit >  cutoff & dat$t == 0, ]
  g3 <- dat[dat$benefit <= cutoff & dat$t == 1, ]
  g4 <- dat[dat$benefit <= cutoff & dat$t == 0, ]
  n1 <- nrow(g1); n2 <- nrow(g2); n3 <- nrow(g3); n4 <- nrow(g4)
  n_total <- n1 + n2 + n3 + n4

  if (any(c(n1, n2, n3, n4) == 0)) {
    message("  Skipping ", dataset_name, " [", method_name, "] cutoff=", round(cutoff,2),
            " — empty cell (n1=",n1," n2=",n2," n3=",n3," n4=",n4,")")
    return(data.frame(
      Dataset = dataset_name, Method = method_name, Cutoff = cutoff,
      n_total = n_total, n_G1 = n1, n_G2 = n2, n_G3 = n3, n_G4 = n4,
      prop_allocated_treatment = (n1 + n2) / n_total,
      PB1 = NA, PB1_lower = NA, PB1_upper = NA
    ))
  }

  PB1     <- sum(g4$Y) / (n1 + n4) - sum(g3$Y) / (n1 + n3) +
    sum(g1$Y) * (1 / (n1 + n4) - 1 / (n1 + n3))
  var_PB1 <- n4 * var(g4$Y) / (n1 + n4)^2 +
    n3 * var(g3$Y) / (n1 + n3)^2 +
    n1 * var(g1$Y) * (1 / (n1 + n4) - 1 / (n1 + n3))^2

  data.frame(
    Dataset  = dataset_name, Method = method_name, Cutoff = cutoff,
    n_total  = n_total,
    n_G1 = n1, n_G2 = n2, n_G3 = n3, n_G4 = n4,
    prop_allocated_treatment = (n1 + n2) / n_total,
    PB1       = PB1,
    PB1_lower = PB1 - 1.96 * sqrt(var_PB1),
    PB1_upper = PB1 + 1.96 * sqrt(var_PB1)
  )
}

# Build full dat_storage including development sample (uses CV benefit predictions)
dat_storage[["Deprexis"]] <- list(
  enr = data.frame(Y = dat.CV.all[, outcome], t = dat.CV.all$Group,
                   benefit = dat.CV.all$m1.CV.treat1 - dat.CV.all$m1.CV.treat0,
                   row.names = NULL),
  ols = data.frame(Y = dat.CV.all[, outcome], t = dat.CV.all$Group,
                   benefit = dat.CV.all$m2.CV.treat1 - dat.CV.all$m2.CV.treat0,
                   row.names = NULL),
  cf  = data.frame(Y = dat.CV.all[, outcome], t = dat.CV.all$Group,
                   benefit = dat.CV.all$m3.benefit,
                   row.names = NULL)
)

pb_results <- list()
counter    <- 1

for (dname in names(dat_storage)) {
  # Pooled SD from ENR frame (Y and t are identical across methods)
  enr_dat <- dat_storage[[dname]]$enr
  psd     <- calculate_pooled_sd(enr_dat)
  cutoff24 <- 0.24 * psd
  cutoff50 <- 0.50 * psd
  message("  ", dname, "  pooled SD = ", round(psd, 3),
          "  cutoff24 = ", round(cutoff24, 3),
          "  cutoff50 = ", round(cutoff50, 3))
  for (mname in c("enr", "ols", "cf")) {
    dat_m <- dat_storage[[dname]][[mname]]
    for (co in list(cutoff24, cutoff50)) {
      pb_results[[counter]] <- compute_pb1(dat_m, co, dname, toupper(mname))
      counter <- counter + 1
    }
  }
}

pb_results_df <- bind_rows(pb_results) %>%
  mutate(
    Sample = case_when(
      Dataset == "Deprexis"    ~ "Deprexis (Dev)",
      Dataset == "edupression" ~ "edupression",
      Dataset == "Selfapy"     ~ "Selfapy",
      Dataset == "Moritz"      ~ "Deprexis (Val)"
    ),
    Sample    = factor(Sample, levels = c("Deprexis (Dev)", "edupression",
                                           "Selfapy", "Deprexis (Val)")),
    Method    = factor(Method,  levels = c("ENR", "OLS", "CF")),
    Threshold = ifelse(round(Cutoff / ave(Cutoff, Dataset, FUN = max), 2) <= 0.70,
                       "d >= 0.24", "d >= 0.50"),
    across(c(PB1, PB1_lower, PB1_upper, prop_allocated_treatment),
           ~ round(., 3))
  ) %>%
  arrange(Sample, Method, Threshold)

print(pb_results_df %>%
        select(Sample, Method, Threshold, n_total, n_G1, n_G2, n_G3, n_G4,
               prop_allocated_treatment, PB1, PB1_lower, PB1_upper))


# ==============================================================================
# SAVE ALL OUTPUTS
# ==============================================================================

write.csv(all_prediction_summaries,  "benefit_prediction_summaries_all_completer.csv",  row.names = FALSE)
write.csv(all_calibration_combined,  "benefit_calibration_all_samples_completer.csv",   row.names = FALSE)
write.csv(all_performance_combined,  "benefit_performance_all_samples_completer.csv",   row.names = FALSE)
write.csv(reduced_table_benefit,     "benefit_performance_reduced_completer.csv",        row.names = FALSE)
write.csv(pb_results_df,             "population_benefit_completer.csv",                 row.names = FALSE)

message("\n✓ All results saved.")
message("  Naming: 'edupression' (not 'Eficasy'), methods 'ENR'/'OLS'/'CF' throughout")
# ==============================================================================
# SUPPLEMENT TABLE: Benefit Calibration Metrics Across Quantile Groupings
# ==============================================================================

# n is constant across bin groupings within a sample, so take it from N_bins == 4
# (present in all samples) to get one clean row per Sample x Method.
# We then join onto the full filtered table — n repeats across bin rows, which
# is correct since the same participants are just split differently each time.
n_by_group <- all_calibration_combined[all_calibration_combined$N_bins == 4, ]
n_by_group <- aggregate(
  cbind(n_treat, n_control) ~ Sample + Method,
  data = n_by_group,
  FUN  = sum
)
supplement_calibration_table <- all_performance_combined %>%
  filter((Dataset == "Deprexis" & N_bins %in% c(4, 11, 18)) |
           (Dataset != "Deprexis" & N_bins %in% c(3, 4, 5))) %>%
 
  mutate(
    Intercept = paste0(round(intercept_abs, 2), " [",
                       round(intercept_lower_abs, 2), ", ",
                       round(intercept_upper_abs, 2), "]"),
    Slope     = paste0(round(slope_abs, 2), " [",
                       round(slope_lower_abs, 2), ", ",
                       round(slope_upper_abs, 2), "]"),
    RMSE      = round(RMSE_abs_standardized, 2),
    R2        = round(r2_abs, 2)
  ) %>%
  select(Sample, Method, `N groups` = N_bins,
         
         Intercept, Slope, RMSE, R2) %>%
  arrange(Sample, Method, `N groups`)

message("\n\n========================================")
message("SUPPLEMENT TABLE: Benefit Calibration Metrics")
message("========================================\n")
print(supplement_calibration_table, row.names = FALSE)

write.csv(supplement_calibration_table,
          "supplement_calibration_table_completer.csv",
          row.names = FALSE)

message("\n✓ Saved to supplement_calibration_table.csv")
