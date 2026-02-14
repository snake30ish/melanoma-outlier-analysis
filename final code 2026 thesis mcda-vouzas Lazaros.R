############################################################
# MELANOMA ANALYSIS: AFT LOG-LOGISTIC & OUTLIERS
# Περιγραφή: Εντοπισμός εκτόπων τιμών μέσω AFT μοντέλου
############################################################

# --- 1. Βιβλιοθήκες & Ρυθμίσεις ---
libs <- c("MASS", "FNN", "dplyr", "survival", "flexsurv", "survminer")
invisible(lapply(libs, library, character.only = TRUE))

set.seed(123) 

# --- 2. Εισαγωγή & Προετοιμασία Δεδομένων ---
raw_data <- read.csv("clinical_melanoma_clinical_plus_scores outlies + non outliers.csv", check.names = FALSE)

# Διόρθωση χρόνου (ο λογάριθμος της Log-logistic απαιτεί t > 0)
raw_data$time <- ifelse(raw_data$time == 0, 0.1, raw_data$time)
treatment_vars <- c("had_pharmaceutical", "had_radiation", "had_chemotherapy", "had_immunotherapy", "had_surgery")
raw_data[treatment_vars] <- lapply(raw_data[treatment_vars], factor)

# --- 3. Εντοπισμός Βιολογικών Outliers (Pathways) ---
path_features <- c("MAPK_Score", "CellCycle_Score", "TumorSuppressor_Score",
                   "Differentiation_Score", "AntigenPresentation_Score",
                   "ImmuneCheckpoint_Score", "RTK_PI3K_Score", "Metabolic_Score",
                   "WNT_Score", "Migration_Score")

path_subset <- na.omit(raw_data[, path_features])
X_scaled <- scale(path_subset)
n_obs <- nrow(X_scaled)
orig_idx <- as.numeric(rownames(X_scaled))

# Mahalanobis & kNN Αποστάσεις
m_dist <- mahalanobis(X_scaled, colMeans(X_scaled), cov(X_scaled))
k_neighbors <- floor(sqrt(n_obs))
knn_dist <- rowMeans(get.knn(X_scaled, k = k_neighbors)$nn.dist)

# Bootstrap για καθορισμό ορίων (Thresholds)
B_reps <- 10000; conf_level <- 0.95 
m_boot_limits <- replicate(B_reps, quantile(mahalanobis(X_scaled[sample(1:n_obs, replace=T),], colMeans(X_scaled), cov(X_scaled)), conf_level))
knn_boot_limits <- replicate(B_reps, quantile(rowMeans(get.knn(X_scaled[sample(1:n_obs, replace=T),], k=k_neighbors)$nn.dist), conf_level))

m_limit <- mean(m_boot_limits)
knn_limit <- mean(knn_boot_limits)

raw_data$flag_pathway <- FALSE
raw_data$flag_pathway[orig_idx] <- (m_dist > m_limit) & (knn_dist > knn_limit)

# --- 4. Μοντελοποίηση Cox & Κλινικοί Outliers ---
cox_input <- raw_data[, c("time", "event", path_features, treatment_vars, "diagnoses.age_at_diagnosis")]
fit_cox_full <- coxph(Surv(time, event) ~ ., data = cox_input)
fit_cox_final <- step(fit_cox_full, direction = "backward", trace = 0)

# Martingale Outliers
raw_data$martingale_res <- residuals(fit_cox_final, type = "martingale")
m_cutoff <- quantile(abs(raw_data$martingale_res), 0.95, na.rm = TRUE)
raw_data$flag_cox <- abs(raw_data$martingale_res) > m_cutoff

# DFBETAs Outliers (Επιρροή στις παραμέτρους)
dfb_matrix <- residuals(fit_cox_final, type = "dfbeta")
dfb_limit <- 2 / sqrt(nrow(dfb_matrix))
raw_data$DFBETA_Outlier <- apply(abs(dfb_matrix), 1, function(val) any(val > dfb_limit))

# --- 5. Παραμετρική Ανάλυση AFT (Log-Logistic) & Outliers ---
# Προσαρμογή του μοντέλου AFT
fit_llogis <- flexsurvreg(Surv(time, event) ~ MAPK_Score + CellCycle_Score + 
                            ImmuneCheckpoint_Score + RTK_PI3K_Score + Metabolic_Score + 
                            had_pharmaceutical + had_immunotherapy, 
                          data = raw_data, dist = "llogis")

# Υπολογισμός PIT Residuals (Probability Integral Transform)
# Τα PIT residuals (1 - Survival Probability) δείχνουν πόσο "αναμενόμενο" ήταν το event
s_estimates <- summary(fit_llogis, type = "survival", t = raw_data$time, tidy = TRUE)$est
raw_data$aft_pit_val <- 1 - s_estimates

# Εντοπισμός AFT Outliers (Ακραίο 5% σε κάθε πλευρά)
raw_data$flag_aft_outlier <- (raw_data$aft_pit_val < 0.05 | raw_data$aft_pit_val > 0.95)

# --- 6. Εξαγωγή Αποτελεσμάτων ---
# Φιλτράρουμε τους ασθενείς που είναι outliers σε τουλάχιστον μία μέθοδο
final_outliers_report <- raw_data %>% 
  filter(flag_pathway | flag_cox | flag_aft_outlier | DFBETA_Outlier)

# Καταμέτρηση
cat("Σύνοψη Outliers:\n")
cat("- Pathway Outliers:", sum(raw_data$flag_pathway), "\n")
cat("- Cox Martingale Outliers:", sum(raw_data$flag_cox), "\n")
cat("- AFT Log-Logistic Outliers:", sum(raw_data$flag_aft_outlier), "\n")
cat("- DFBETA Outliers:", sum(raw_data$DFBETA_Outlier), "\n")

write.csv(final_outliers_report, "melanoma_aft_outliers_final.csv", row.names = FALSE)
cat("\nΗ ανάλυση ολοκληρώθηκε. Το αρχείο 'melanoma_aft_outliers_final.csv' δημιουργήθηκε.")
