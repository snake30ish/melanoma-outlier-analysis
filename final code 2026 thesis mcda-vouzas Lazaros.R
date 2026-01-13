# --- Melanoma Outlier & Influence Analysis Script ---
# Author: Lazaros vouzas
# Date: Jan 2026

# Φόρτωση απαραίτητων εργαλείων
libs <- c("MASS", "FNN", "dplyr", "survival", "flexsurv", "survminer")
invisible(lapply(libs, library, character.only = TRUE))

set.seed(123)

# Εισαγωγή δεδομένων
raw_data <- read.csv("clinical_melanoma_clinical_plus_scores outlies + non outliers.csv", check.names = FALSE)

# --- 1. Ανάλυση Outliers στα Βιολογικά Μονοπάτια (Pathways) ---
path_features <- c("MAPK_Score", "CellCycle_Score", "TumorSuppressor_Score",
                   "Differentiation_Score", "AntigenPresentation_Score",
                   "ImmuneCheckpoint_Score", "RTK_PI3K_Score", "Metabolic_Score",
                   "WNT_Score", "Migration_Score")

# Καθαρισμός και κανονικοποίηση
path_subset <- na.omit(raw_data[, path_features])
X_scaled <- scale(path_subset)
n_obs <- nrow(X_scaled)
orig_idx <- as.numeric(rownames(X_scaled))

# Υπολογισμός αποστάσεων (Mahalanobis & kNN)
m_dist <- mahalanobis(X_scaled, colMeans(X_scaled), cov(X_scaled))
k_neighbors <- floor(sqrt(n_obs))
knn_dist <- rowMeans(get.knn(X_scaled, k = k_neighbors)$nn.dist)

# Bootstrap εκτίμηση ορίων (Thresholding)
B_reps <- 10000; conf_level <- 0.95
m_boot_limits <- numeric(B_reps); knn_boot_limits <- numeric(B_reps)

for(b in 1:B_reps){
  boot_sample <- X_scaled[sample(1:n_obs, n_obs, replace = TRUE), ]
  m_boot_limits[b] <- quantile(mahalanobis(boot_sample, colMeans(boot_sample), cov(boot_sample)), conf_level)
  knn_boot_limits[b] <- quantile(rowMeans(get.knn(boot_sample, k = k_neighbors)$nn.dist), conf_level)
}

m_limit <- mean(m_boot_limits)
knn_limit <- mean(knn_boot_limits)

# Φιλτράρισμα Pathway Outliers
raw_data$flag_pathway <- FALSE
raw_data$flag_pathway[orig_idx] <- (m_dist > m_limit) & (knn_dist > knn_limit)

# Προετοιμασία βοηθητικού πίνακα για το reporting
pathway_results <- data.frame(
  cases.submitter_id = raw_data$cases.submitter_id[orig_idx],
  Mahalanobis_Dist = m_dist,
  kNN_Dist = knn_dist,
  is_outlier = raw_data$flag_pathway[orig_idx]
) %>% filter(is_outlier == TRUE)

# --- 2. Μοντελοποίηση Επιβίωσης (Cox Proportional Hazards) ---
# Διόρθωση οριακών τιμών χρόνου
raw_data$time <- ifelse(raw_data$time == 0, 0.1, raw_data$time)
treatment_vars <- c("had_pharmaceutical", "had_radiation", "had_chemotherapy", "had_immunotherapy", "had_surgery")
raw_data[treatment_vars] <- lapply(raw_data[treatment_vars], factor)

# Fit Cox model & Stepwise Selection
cox_input <- raw_data[, c("time", "event", path_features, treatment_vars, "diagnoses.age_at_diagnosis")]
fit_cox_full <- coxph(Surv(time, event) ~ ., data = cox_input)
fit_cox_final <- step(fit_cox_full, direction = "backward", trace = 0)

# Martingale Residuals & Outliers
raw_data$martingale_res <- residuals(fit_cox_final, type = "martingale")
m_cutoff <- quantile(abs(raw_data$martingale_res), 0.95, na.rm = TRUE)
raw_data$flag_cox <- abs(raw_data$martingale_res) > m_cutoff

# DFBETA (Influential Cases)
dfb_matrix <- residuals(fit_cox_final, type = "dfbeta")
dfb_limit <- 2 / sqrt(nrow(dfb_matrix))
raw_data$DFBETA_Outlier <- apply(abs(dfb_matrix), 1, function(val) any(val > dfb_limit))

# --- 3. Παραμετρική Ανάλυση (FHT / Log-Logistic) ---
fit_llogis <- flexsurvreg(Surv(time, event) ~ MAPK_Score + CellCycle_Score + 
                            ImmuneCheckpoint_Score + RTK_PI3K_Score + Metabolic_Score + 
                            had_pharmaceutical + had_immunotherapy, 
                          data = raw_data, dist = "llogis")

# PIT Residuals
s_estimates <- summary(fit_llogis, type = "survival", t = raw_data$time, tidy = TRUE)$est
raw_data$pit_val <- 1 - s_estimates
raw_data$flag_survival_param <- (raw_data$pit_val < 0.05 | raw_data$pit_val > 0.95)

# Hazard Visualization
grid_t <- seq(1, quantile(raw_data$time, 0.95, na.rm = TRUE), length.out = 300)
haz_data <- summary(fit_llogis, type = "hazard", t = grid_t, tidy = TRUE)
plot(grid_t, haz_data$est, type = "l", col = "red", lwd = 2,
     xlab = "Days", ylab = "Hazard Rate", main = "Log-Logistic Hazard Profile")

# --- 4. Leave-One-Out (LOO) Influence Metrics ---
n_total <- nrow(raw_data)
s_formula <- formula(fit_llogis)
inf_aft <- numeric(n_total)
inf_fht <- numeric(n_total)

for(i in 1:n_total){
  fit_loo <- try(flexsurvreg(s_formula, data = raw_data[-i, ], dist = "llogis"), silent = TRUE)
  if(!inherits(fit_loo, "try-error")){
    # AFT Parametric shift
    inf_aft[i] <- mean(abs(exp(fit_loo$res[-(1:2), 1]) - exp(fit_llogis$res[-(1:2), 1])))
    # FHT/PIT stability
    pit_loo <- 1 - summary(fit_loo, type = "survival", t = raw_data$time[-i], tidy = TRUE)$est
    inf_fht[i] <- mean(abs(pit_loo - raw_data$pit_val[-i]))
  }
}

# --- 5. Εξαγωγή Αναφορών (CSV Generation) ---

clinical_cols <- c("cases.submitter_id", "hospital", "demographic.gender", 
                   "diagnoses.age_at_diagnosis", "time", "event", 
                   "demographic.race", "demographic.ethnicity",
                   treatment_vars, "had_radiation")

# Pathway Outliers Report
report_pathway <- pathway_results %>%
  left_join(raw_data %>% select(all_of(clinical_cols), flag_pathway), by = "cases.submitter_id") %>%
  arrange(desc(Mahalanobis_Dist)) %>%
  select(hospital, diagnoses.age_at_diagnosis, demographic.gender, everything())

write.csv(report_pathway, "pathway_outliers_18_clinical_profile.csv", row.names = FALSE)

# Cox Outliers Report
report_cox <- raw_data %>%
  filter(flag_cox | DFBETA_Outlier) %>%
  select(all_of(clinical_cols), martingale_res, DFBETA_Outlier) %>%
  arrange(desc(abs(martingale_res)))

write.csv(report_cox, "cox_outliers_profile.csv", row.names = FALSE)

# FHT/AFT Outliers Report
report_fht <- raw_data %>%
  filter(flag_survival_param) %>%
  select(all_of(clinical_cols), pit_val) %>%
  arrange(desc(abs(pit_val - 0.5)))

write.csv(report_fht, "fht_outliers_profile.csv", row.names = FALSE)


cat("Ανάλυση ολοκληρώθηκε. Τα αρχεία αποθηκεύτηκαν επιτυχώς.\n")
