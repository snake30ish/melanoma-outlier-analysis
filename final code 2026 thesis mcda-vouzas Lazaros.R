############################################################
# MELANOMA MULTI-LEVEL ANALYSIS: OUTLIERS & FHT SIMULATION
# Συγγραφέας: Lazaros Vouzas & Gemini
# Περιγραφή: Πλήρης ροή ανάλυσης από τον εντοπισμό εκτόπων τιμών 
# έως τη στοχαστική προσομοίωση φθοράς υγείας (Wiener Process).
############################################################

# --- 1. Βιβλιοθήκες & Ρυθμίσεις ---
libs <- c("MASS", "FNN", "dplyr", "survival", "flexsurv", "survminer", "truncnorm")
invisible(lapply(libs, library, character.only = TRUE))

set.seed(123) # Για αναπαραγωγιμότητα

# --- 2. Εισαγωγή & Προετοιμασία Δεδομένων ---
raw_data <- read.csv("clinical_melanoma_clinical_plus_scores outlies + non outliers.csv", check.names = FALSE)

# Διόρθωση χρόνου και μετατροπή μεταβλητών θεραπείας σε factors
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

# Διαγνωστικά Cox (Martingale & DFBETAs)
raw_data$martingale_res <- residuals(fit_cox_final, type = "martingale")
m_cutoff <- quantile(abs(raw_data$martingale_res), 0.95, na.rm = TRUE)
raw_data$flag_cox <- abs(raw_data$martingale_res) > m_cutoff

dfb_matrix <- residuals(fit_cox_final, type = "dfbeta")
dfb_limit <- 2 / sqrt(nrow(dfb_matrix))
raw_data$DFBETA_Outlier <- apply(abs(dfb_matrix), 1, function(val) any(val > dfb_limit))

# --- 5. Παραμετρική Ανάλυση (Log-Logistic AFT/FHT) ---
fit_llogis <- flexsurvreg(Surv(time, event) ~ MAPK_Score + CellCycle_Score + 
                            ImmuneCheckpoint_Score + RTK_PI3K_Score + Metabolic_Score + 
                            had_pharmaceutical + had_immunotherapy, 
                          data = raw_data, dist = "llogis")

# Στοχαστικοί Outliers (PIT Residuals)
s_estimates <- summary(fit_llogis, type = "survival", t = raw_data$time, tidy = TRUE)$est
raw_data$pit_val <- 1 - s_estimates
raw_data$flag_survival_param <- (raw_data$pit_val < 0.05 | raw_data$pit_val > 0.95)

# --- 6. Προσομοίωση FHT (Wiener Process) ---
# Εξαγωγή παραμέτρων για το Drift & Diffusion
params <- fit_llogis$res[, "est"]
shape_alpha <- as.numeric(params["shape"])
scale_lambda <- as.numeric(params["scale"]) 

# Επιλογή 3 τυχαίων ασθενών με event
target_patients <- raw_data %>% filter(event == 1) %>% sample_n(3)

# Υπολογισμός αρχικής "Υγείας" y0 (Latent Health)
target_patients$y0_sampled <- rtruncnorm(3, a = 0, b = Inf, 
                                         mean = scale_lambda, 
                                         sd = scale_lambda / shape_alpha)

# Plotting Προσομοίωσης
colors <- c("#E41A1C", "#377EB8", "#4DAF4A")
t_max <- 6000; sigma_w <- 10
km_median <- summary(survfit(Surv(time, event) ~ 1, data = raw_data))$table["median"]

plot(NULL, xlim = c(0, t_max), ylim = c(0, max(target_patients$y0_sampled) + 500),
     xlab = "Ημέρες (Time)", ylab = "Latent Health (y_t)",
     main = "FHT Model: Στοχαστικές Τροχιές 3 Ασθενών")
abline(h = 0, col = "red", lwd = 2, lty = 2)

for(i in 1:nrow(target_patients)) {
  y0 <- target_patients$y0_sampled[i]
  drift <- -(y0 / km_median)
  times <- seq(0, t_max, by = 1)
  noise_path <- cumsum(rnorm(length(times), mean = 0, sd = sigma_w))
  health_path <- y0 + drift * times + noise_path
  
  death_idx <- border <- which(health_path <= 0)[1]
  valid_time <- if(!is.na(death_idx)) times[1:death_idx] else times
  valid_path <- if(!is.na(death_idx)) health_path[1:death_idx] else health_path
  
  lines(valid_time, valid_path, col = colors[i], lwd = 2.5)
  if(!is.na(death_idx)) points(max(valid_time), 0, col = colors[i], pch = 19)
}

# --- 7. Εξαγωγή Αποτελεσμάτων ---
write.csv(raw_data %>% filter(flag_pathway | flag_cox | flag_survival_param), 
          "melanoma_combined_outliers.csv", row.names = FALSE)

cat("Η ανάλυση ολοκληρώθηκε. Τα αρχεία εξήχθησαν.")
