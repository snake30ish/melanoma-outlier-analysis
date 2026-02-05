# --- Melanoma Outlier & Influence Analysis Script ---
# Author: Lazaros Vouzas
# Date: Jan 2026
# Description: Multi-level detection of biological and survival outliers.
# Περιγραφή: Πολυεπίπεδος εντοπισμός βιολογικών και κλινικών εκτόπων τιμών.

# --- 0. Load Libraries | Φόρτωση Βιβλιοθηκών ---
libs <- c("MASS", "FNN", "dplyr", "survival", "flexsurv", "survminer")
invisible(lapply(libs, library, character.only = TRUE))

set.seed(123)

# --- 1. Data Import | Εισαγωγή Δεδομένων ---
raw_data <- read.csv("clinical_melanoma_clinical_plus_scores outlies + non outliers.csv", check.names = FALSE)

# --- 2. Pathway Outliers (Mahalanobis + kNN) | Βιολογικοί Outliers ---
path_features <- c("MAPK_Score", "CellCycle_Score", "TumorSuppressor_Score",
                   "Differentiation_Score", "AntigenPresentation_Score",
                   "ImmuneCheckpoint_Score", "RTK_PI3K_Score", "Metabolic_Score",
                   "WNT_Score", "Migration_Score")

# Clean and Scale data | Καθαρισμός και Κανονικοποίηση
path_subset <- na.omit(raw_data[, path_features])
X_scaled <- scale(path_subset)
n_obs <- nrow(X_scaled)
orig_idx <- as.numeric(rownames(X_scaled))

# Compute distances | Υπολογισμός αποστάσεων
m_dist <- mahalanobis(X_scaled, colMeans(X_scaled), cov(X_scaled))
k_neighbors <- floor(sqrt(n_obs))
knn_dist <- rowMeans(get.knn(X_scaled, k = k_neighbors)$nn.dist)

# Bootstrap estimation for thresholds | Bootstrap εκτίμηση ορίων
B_reps <- 10000; conf_level <- 0.95
m_boot_limits <- numeric(B_reps); knn_boot_limits <- numeric(B_reps)

for(b in 1:B_reps){
  boot_sample <- X_scaled[sample(1:n_obs, n_obs, replace = TRUE), ]
  m_boot_limits[b] <- quantile(mahalanobis(boot_sample, colMeans(boot_sample), cov(boot_sample)), conf_level)
  knn_boot_limits[b] <- quantile(rowMeans(get.knn(boot_sample, k = k_neighbors)$nn.dist), conf_level)
}

m_limit <- mean(m_boot_limits)
knn_limit <- mean(knn_boot_limits)

# Identify biological outliers | Ταυτοποίηση βιολογικών εκτόπων τιμών
raw_data$flag_pathway <- FALSE
raw_data$flag_pathway[orig_idx] <- (m_dist > m_limit) & (knn_dist > knn_limit)

# Prep intermediate results | Προετοιμασία ενδιάμεσων αποτελεσμάτων
pathway_results <- data.frame(
  cases.submitter_id = raw_data$cases.submitter_id[orig_idx],
  Mahalanobis_Dist = m_dist,
  kNN_Dist = knn_dist,
  is_outlier = raw_data$flag_pathway[orig_idx]
) %>% filter(is_outlier == TRUE)

# --- 3. Cox Model Diagnostics | Μοντελοποίηση Cox και Διαγνωστικά ---
# Adjust time zeros | Διόρθωση μηδενικών τιμών χρόνου
raw_data$time <- ifelse(raw_data$time == 0, 0.1, raw_data$time)
treatment_vars <- c("had_pharmaceutical", "had_radiation", "had_chemotherapy", "had_immunotherapy", "had_surgery")
raw_data[treatment_vars] <- lapply(raw_data[treatment_vars], factor)

# Fit Cox model & Selection | Επιλογή και Προσαρμογή Μοντέλου Cox
cox_input <- raw_data[, c("time", "event", path_features, treatment_vars, "diagnoses.age_at_diagnosis")]
fit_cox_full <- coxph(Surv(time, event) ~ ., data = cox_input)
fit_cox_final <- step(fit_cox_full, direction = "backward", trace = 0)

# Martingale Residuals | Υπόλοιπα Martingale (Σφάλμα Πρόβλεψης)
raw_data$martingale_res <- residuals(fit_cox_final, type = "martingale")
m_cutoff <- quantile(abs(raw_data$martingale_res), 0.95, na.rm = TRUE)
raw_data$flag_cox <- abs(raw_data$martingale_res) > m_cutoff

# DFBETA Influence | Δείκτες DFBETA (Στατιστική Επιρροή)
dfb_matrix <- residuals(fit_cox_final, type = "dfbeta")
dfb_limit <- 2 / sqrt(nrow(dfb_matrix))
raw_data$DFBETA_Outlier <- apply(abs(dfb_matrix), 1, function(val) any(val > dfb_limit))

# --- 4.1. Προσαρμογή του μοντέλου με τις 7 μεταβλητές ---
fit_llogis_7 <- flexsurvreg(
  Surv(time, event) ~ MAPK_Score + CellCycle_Score + 
    ImmuneCheckpoint_Score + RTK_PI3K_Score + Metabolic_Score + 
    had_pharmaceutical + had_immunotherapy, 
  data = raw_data, 
  dist = "llogis"
)

# --- 4.2. Εμφάνιση του κλασικού summary της flexsurv ---
# Αυτό σου δίνει τα αποτελέσματα σε log-scale (τα beta)
print(fit_llogis_7)

# --- 4.3. Δημιουργία καθαρού πίνακα αποτελεσμάτων (R Data Frame) ---
# Περιλαμβάνει τα Beta, τα Time Ratios (TR) και τα p-values
results_table <- as.data.frame(fit_llogis_7$res)
results_table$Variable <- rownames(results_table)

# Υπολογισμός Time Ratios και P-values για τις συμμεταβλητές
results_table <- results_table %>%
  mutate(
    Time_Ratio = ifelse(!Variable %in% c("shape", "scale"), exp(est), NA),
    z_value = est / se,
    p_value = 2 * (1 - pnorm(abs(z_value)))
  ) %>%
  select(Variable, est, Time_Ratio, se, p_value)

# --- 5. Influence Analysis (LOO) | Ανάλυση Επιρροής Leave-One-Out ---
n_total <- nrow(raw_data)
s_formula <- formula(fit_llogis)
inf_aft <- numeric(n_total)
inf_fht <- numeric(n_total)

for(i in 1:n_total){
  fit_loo <- try(flexsurvreg(s_formula, data = raw_data[-i, ], dist = "llogis"), silent = TRUE)
  if(!inherits(fit_loo, "try-error")){
    # Parametric shift in AFT | Μεταβολή παραμέτρων AFT
    inf_aft[i] <- mean(abs(exp(fit_loo$res[-(1:2), 1]) - exp(fit_llogis$res[-(1:2), 1])))
    # PIT stability for FHT | Σταθερότητα PIT για FHT
    pit_loo <- 1 - summary(fit_loo, type = "survival", t = raw_data$time[-i], tidy = TRUE)$est
    inf_fht[i] <- mean(abs(pit_loo - raw_data$pit_val[-i]))
  }
}

# --- 6. Export Reports | Εξαγωγή Αναφορών σε CSV ---
clinical_cols <- c("cases.submitter_id", "hospital", "demographic.gender", 
                   "diagnoses.age_at_diagnosis", "time", "event", 
                   "demographic.race", "demographic.ethnicity",
                   treatment_vars, "had_radiation")

# 6.1 Biological Outliers | Βιολογικοί Outliers
report_pathway <- pathway_results %>%
  left_join(raw_data %>% select(all_of(clinical_cols), flag_pathway), by = "cases.submitter_id") %>%
  arrange(desc(Mahalanobis_Dist)) %>%
  select(hospital, diagnoses.age_at_diagnosis, demographic.gender, everything())

write.csv(report_pathway, "pathway_outliers_18_clinical_profile.csv", row.names = FALSE)

# 6.2 Cox Outliers | Outliers του Μοντέλου Cox
report_cox <- raw_data %>%
  filter(flag_cox | DFBETA_Outlier) %>%
  select(all_of(clinical_cols), martingale_res, DFBETA_Outlier) %>%
  arrange(desc(abs(martingale_res)))

write.csv(report_cox, "cox_outliers_profile.csv", row.names = FALSE)

# 6.3 FHT/AFT Outliers | Στοχαστικοί Outliers
report_fht <- raw_data %>%
  filter(flag_survival_param) %>%
  select(all_of(clinical_cols), pit_val) %>%
  arrange(desc(abs(pit_val - 0.5)))

write.csv(report_fht, "fht_outliers_profile.csv", row.names = FALSE)

cat("Analysis complete. Reports exported successfully. | Η ανάλυση ολοκληρώθηκε επιτυχώς.\n")

