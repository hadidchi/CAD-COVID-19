# Install and load required package
install.packages("survival")
library("survival")

# Load data
attach(CAD)  # Replace with your actual dataset if different
data <- CAD

# Define outcomes of interest
columns_of_interest <- c("DEATH", "CMP_MACE")

# Define outcome labels
column_meanings <- c(
  "DEATH" = "All-Cause Mortality",
  "CMP_MACE" = "Major Adverse Cardiovascular Events"
)

# Define whether to exclude pre-existing conditions (1 = exclude, 0 = don't exclude)
exclude_pre <- c(
  "DEATH" = 0,
  "CMP_MACE" = 0
)

# Covariates to adjust for
covariates <- c("COVID")  # <- Adjust this list as needed

# Loop through outcomes
for (column in columns_of_interest) {
  
  outcome_name <- column_meanings[column]
  subset_condition <- if (exclude_pre[column] == 1) {
    data[[paste0("PRE_", column)]] == 0
  } else {
    rep(TRUE, nrow(data))
  }
  
  # Initialize results lists
  univariate_results <- list()
  multivariate_results <- list()
  
  # ---- UNIVARIATE MODELS ----
  for (cov in covariates) {
    formula_uni <- as.formula(paste0("Surv(", column, "_TIME, ", column, ") ~ ", cov))
    cox_uni <- coxph(formula_uni, data = data[subset_condition, ])
    
    sum_uni <- summary(cox_uni)
    hr <- sprintf("%.2f", sum_uni$coef[1, "exp(coef)"])
    ci_lower <- sprintf("%.2f", sum_uni$conf.int[,"lower .95"][1])
    ci_upper <- sprintf("%.2f", sum_uni$conf.int[,"upper .95"][1])
    p <- sum_uni$coef[1, "Pr(>|z|)"]
    p_fmt <- ifelse(p < 0.005, "<0.005", sprintf("%.3f", p))
    hr_ci <- paste0(hr, " [", ci_lower, ", ", ci_upper, "]")
    
    univariate_results[[length(univariate_results) + 1]] <- data.frame(
      Outcome = outcome_name,
      Model = "Univariate",
      Covariate = cov,
      HR_95CI = hr_ci,
      P_value = p_fmt
    )
  }
  
  # ---- MULTIVARIATE MODEL ----
  formula_multi <- as.formula(paste0("Surv(", column, "_TIME, ", column, ") ~ ", paste(covariates, collapse = " + ")))
  cox_multi <- coxph(formula_multi, data = data[subset_condition, ])
  sum_multi <- summary(cox_multi)
  
  for (i in seq_along(covariates)) {
    cov <- covariates[i]
    hr <- sprintf("%.2f", sum_multi$coef[i, "exp(coef)"])
    ci_lower <- sprintf("%.2f", sum_multi$conf.int[i, "lower .95"])
    ci_upper <- sprintf("%.2f", sum_multi$conf.int[i, "upper .95"])
    p <- sum_multi$coef[i, "Pr(>|z|)"]
    p_fmt <- ifelse(p < 0.005, "<0.005", sprintf("%.3f", p))
    hr_ci <- paste0(hr, " [", ci_lower, ", ", ci_upper, "]")
    
    multivariate_results[[length(multivariate_results) + 1]] <- data.frame(
      Outcome = outcome_name,
      Model = "Multivariate",
      Covariate = cov,
      HR_95CI = hr_ci,
      P_value = p_fmt
    )
  }
  
  # Combine and write to CSV
  all_results <- do.call(rbind, c(univariate_results, multivariate_results))
  write.csv(all_results, paste0("coxph_", column, "_results.csv"), row.names = FALSE)
}