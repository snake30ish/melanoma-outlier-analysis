# --- Melanoma Outlier & Influence Analysis Script ---
# Author: Lazaros Vouzas | Date: Jan 2026

# --- 0. Load Libraries ---
libs <- c("MASS", "FNN", "dplyr", "survival", "flexsurv", "survminer", "truncnorm", "tibble")
invisible(lapply(libs, library, character.only = TRUE))

set.seed(123)

# --- 1. Data Import & Preprocessing ---
raw_data <- read.csv("clinical_melanoma_clinical_plus_scores outlies + non outliers.csv", check.names = FALSE)

# Διόρθωση χρόνων και μετατροπή σε numeric για τους υπολογισμούς του simulation
raw_data$time <- ifelse(raw_data$time == 0, 0.1, raw_data$time)
treatment_vars <- c("had_pharmaceutical", "had_radiation", "had_chemotherapy", "had_immunotherapy", "had_surgery")
raw_data[treatment_vars] <- lapply(raw_data[treatment_vars], function(x) as.numeric(as.character(x)))

path_features <- c("MAPK_Score", "CellCycle_Score", "TumorSuppressor_Score",
                   "Differentiation_Score", "AntigenPresentation_Score",
                   "ImmuneCheckpoint_Score", "RTK_PI3K_Score", "Metabolic_Score",
                   "WNT_Score", "Migration_Score")

# --- 2. Pathway Outliers (Mahalanobis + kNN) ---
path_subset <- na.omit(raw_data[, path_features])
X_scaled <- scale(path_subset)
n_obs <- nrow(X_scaled)
orig_idx <- as.numeric(rownames(X_scaled))

m_dist <- mahalanobis(X_scaled, colMeans(X_scaled), cov(X_scaled))
k_neighbors <- floor(sqrt(n_obs))
knn_dist <- rowMeans(get.knn(X_scaled, k = k_neighbors)$nn.dist)

# Thresholds
m_limit <- mean(replicate(500, quantile(mahalanobis(X_scaled[sample(1:n_obs, replace=T),], colMeans(X_scaled), cov(X_scaled)), 0.95)))
knn_limit <- mean(replicate(500, quantile(rowMeans(get.knn(X_scaled[sample(1:n_obs, replace=T),], k=k_neighbors)$nn.dist), 0.95)))

raw_data$flag_pathway <- FALSE
raw_data$flag_pathway[orig_idx] <- (m_dist > m_limit) & (knn_dist > knn_limit)

# --- 3. Cox Model Outliers ---
cox_input <- raw_data[, c("time", "event", path_features, treatment_vars, "diagnoses.age_at_diagnosis")]
fit_cox_full <- coxph(Surv(time, event) ~ ., data = cox_input)
fit_cox_final <- step(fit_cox_full, direction = "backward", trace = 0)

# Martingale & DFBETA
raw_data$martingale_res <- residuals(fit_cox_final, type = "martingale")
raw_data$flag_cox <- abs(raw_data$martingale_res) > quantile(abs(raw_data$martingale_res), 0.95, na.rm=T)
dfb_matrix <- residuals(fit_cox_final, type = "dfbeta")
raw_data$DFBETA_Outlier <- apply(abs(dfb_matrix), 1, function(val) any(val > (2/sqrt(nrow(dfb_matrix)))))

# --- 4. FHT / AFT Outliers ---
fit_llogis <- flexsurvreg(Surv(time, event) ~ MAPK_Score + CellCycle_Score + 
                            ImmuneCheckpoint_Score + RTK_PI3K_Score + Metabolic_Score + 
                            had_pharmaceutical + had_immunotherapy, 
                          data = raw_data, dist = "llogis")

# PIT Residuals
s_estimates <- summary(fit_llogis, type = "survival", t = raw_data$time, tidy = TRUE)$est
raw_data$pit_val <- 1 - s_estimates
raw_data$flag_survival_param <- (raw_data$pit_val < 0.05 | raw_data$pit_val > 0.95)

# --- 5. Εξαγωγή Outliers σε CSV ---

# 1. Cox Outliers
report_cox <- raw_data %>%
  filter(flag_cox | DFBETA_Outlier) %>%
  arrange(desc(abs(martingale_res)))
write.csv(report_cox, "cox_outliers_report.csv", row.names = FALSE)

# 2. FHT/AFT Outliers
report_fht <- raw_data %>%
  filter(flag_survival_param) %>%
  arrange(desc(abs(pit_val - 0.5)))
write.csv(report_fht, "fht_aft_outliers_report.csv", row.names = FALSE)

# --- 6. FHT Simulation (Wiener Process) ---
intercept <- 7.679087
shape_alpha <- 1.4488916
betas <- c(MAPK_Score = 0.293875, CellCycle_Score = -0.202272,
           ImmuneCheckpoint_Score = 0.529326, RTK_PI3K_Score = -0.229661,
           Metabolic_Score = -0.230366, had_pharmaceutical = 0.297245,
           had_immunotherapy = 0.911868)

matrix_vars <- names(betas)
# Εξασφάλιση ότι οι στήλες υπάρχουν πριν τον πολλαπλασιασμό
matrix_data <- as.matrix(raw_data[, matrix_vars])
raw_data$mu0_val <- exp(intercept + matrix_data %*% betas)
raw_data$sigma0_val <- raw_data$mu0_val / shape_alpha

# Επιλογή 3 ασθενών για το γράφημα
target_patients <- raw_data %>% filter(event == 1) %>% sample_n(3)
target_patients$y0_sampled <- rtruncnorm(3, a = 0, b = Inf, mean = target_patients$mu0_val, sd = target_patients$sigma0_val)

# Plotting
km_median <- summary(survfit(Surv(time, event) ~ 1, data = raw_data))$table["median"]
plot(NULL, xlim = c(0, 5000), ylim = c(0, max(target_patients$y0_sampled) + 500),
     xlab = "Days", ylab = "Health State", main = "FHT Simulation")
abline(h = 0, col = "red", lty = 2)

for(i in 1:3) {
  y0 <- target_patients$y0_sampled[i]
  drift <- -(y0 / km_median)
  times <- 0:5000
  path <- y0 + drift * times + cumsum(rnorm(length(times), 0, 10))
  end_idx <- which(path <= 0)[1]
  if(is.na(end_idx)) end_idx <- 5001
  lines(times[1:end_idx], path[1:end_idx], col = i+1, lwd = 2)
}

cat("Τα αρχεία 'cox_outliers_report.csv' και 'fht_aft_outliers_report.csv' δημιουργήθηκαν επιτυχώς.\n")
