# ---------- Load Packages ----------
install.packages(c("cmprsk", "nnet", "survival"))
library(cmprsk)
library(nnet)
library(survival)

# ---------- Attach Data ----------
attach(CAD)
data <- CAD
data$ID <- 1:nrow(data)

# ---------- Setup ----------
columns_of_interest <- c("NEW_MI", "NEW_CHF", "NEW_STROKE", "CMP_MACE", "DEATH")
column_meanings <- c(
  "NEW_MI" = "Myocardial Infarction",
  "NEW_CHF" = "Congestive Heart Failure",
  "NEW_STROKE" = "Ischemic or Hemorrhagic Stroke",
  "CMP_MACE" = "Major Adverse Cardiovascular Events",
  "DEATH" = "All-Cause Mortality"
)

covariates <- c("AGE", "MALE", "BLACK", "ASIAN", "OTHER_RACE", "HISPANIC",
                "HTN", "DIABETES_II", "COPD", "ASTHMA", "CKD", "LIVER", "SMOKING",
                "Medicaid", "Medicare", "Self.Pay", "Bottom.25q", "X25.50q", "X50.75q")

# ---------- Toggle Pre-Exclusion for Each Outcome ----------
apply_pre_exclusion <- c(
  "NEW_MI" = FALSE,
  "NEW_CHF" = FALSE,
  "NEW_STROKE" = FALSE,
  "CMP_MACE" = FALSE,
  "DEATH" = FALSE
)

# ---------- Choose Plot Type for Each Outcome ----------
plot_type_per_outcome <- c(
  "NEW_MI" = "CIF",
  "NEW_CHF" = "CIF",
  "NEW_STROKE" = "CIF",
  "CMP_MACE" = "KM",
  "DEATH" = "KM"
)

# ---------- Propensity Score Weights ----------
data$covid <- as.factor(data$covid)
ps_model <- multinom(covid ~ ., data = data[, c("covid", covariates)], trace = FALSE)
ps_probs <- predict(ps_model, type = "probs")

actual_treatment <- as.integer(as.character(data$covid))
treatment_dist <- prop.table(table(actual_treatment))[as.character(actual_treatment)]
data$ipw <- mapply(function(i, t) {
  pi_t <- ps_probs[i, as.character(t)]
  pi_pop <- treatment_dist[as.character(t)]
  pi_pop / pi_t
}, i = 1:nrow(data), t = actual_treatment)

# ---------- Loop Through Outcomes ----------
for (outcome in columns_of_interest) {
  
  label <- column_meanings[outcome]
  time_col <- paste0(outcome, "_TIME")
  status_col <- outcome
  group_var <- "covid"
  plot_type <- plot_type_per_outcome[[outcome]]
  
  subset_condition <- if (apply_pre_exclusion[[outcome]]) {
    data[[paste0("PRE_", outcome)]] == 0
  } else {
    rep(TRUE, nrow(data))
  }
  
  data_sub <- data[subset_condition, ]
  ipw_rounded <- round(data_sub$ipw)
  data_weighted <- data_sub[rep(1:nrow(data_sub), ipw_rounded), ]
  
  if (plot_type == "CIF") {
    # ---------- CIF ----------
    cif_out <- cuminc(
      ftime = data_weighted[[time_col]],
      fstatus = data_weighted[[status_col]],
      group = data_weighted[[group_var]],
      cencode = 0
    )
    cif_primary <- cif_out[grep("1$", names(cif_out))]
    
    plot(NULL, xlim = c(0, 45), ylim = c(0, 0.15),
         xlab = "Time Since Index Date (Months)",
         ylab = paste("Cumulative Incidence of", label),
         main = paste("IPW-Adjusted Cumulative Incidence Function for", label))
    
    for (i in seq_along(cif_primary)) {
      lwd_val <- if (i == 1) 1 else if (i == 2) 2 else 4
      lines(cif_primary[[i]]$time, cif_primary[[i]]$est, col = "black", lwd = lwd_val)
    }
    
    legend("topleft",
           legend = c("COVID+ Hospitalized", "COVID+ Non-Hospitalized", "COVID–"),
           lwd = c(4, 2, 1),
           col = "black",
           bty = "n")
    
  } else if (plot_type == "KM") {
    # ---------- KM ----------
    covid_levels <- sort(unique(as.integer(as.character(data_weighted$covid))))
    plotted <- FALSE
    ymin <- 1
    ymax <- 0
    
    for (i in covid_levels) {
      df_group <- data_weighted[data_weighted$covid == i, ]
      if (sum(df_group[[status_col]] == 1, na.rm = TRUE) == 0) next
      
      surv_obj <- Surv(df_group[[time_col]], df_group[[status_col]] == 1)
      km_fit <- survfit(surv_obj ~ 1)
      
      ymin <- min(ymin, min(km_fit$surv, na.rm = TRUE))
      ymax <- max(ymax, max(km_fit$surv, na.rm = TRUE))
      
      lwd_val <- if (i == 0) 1 else if (i == 1) 2 else 4
      
      if (!plotted) {
        plot(km_fit$time, km_fit$surv, type = "s", col = "black", lwd = lwd_val,
             xlab = "Time Since Index Date (Months)",
             ylab = paste("Survival Probability for", label),
             main = paste("IPW-Adjusted Kaplan-Meier Curve for", label),
             xlim = c(0, 45), ylim = c(0.90, 1))
        plotted <- TRUE
      } else {
        lines(km_fit$time, km_fit$surv, col = "black", lwd = lwd_val)
      }
    }
    
    if (plotted) {
      legend("bottomleft",
             legend = c("COVID+ Hospitalized", "COVID+ Non-Hospitalized", "COVID–"),
             lwd = c(4, 2, 1),
             col = "black",
             bty = "n")
    }
  }
}