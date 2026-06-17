# Modified: Return per-simulation, per-ID raw correlations
extract_all_within_person_correlations <- function(simulated_data_list, target_var, vars) {
  
  sim_correlations_list <- list()
  
  for (sim_idx in seq_along(simulated_data_list)) {
    sim_data <- simulated_data_list[[sim_idx]] %>% mutate(ID = as.character(ID))
    
    sim_correlations <- sim_data %>%
      group_by(ID) %>%
      summarize(across(
        all_of(vars),
        ~ cor(.x, !!sym(target_var), use = "complete.obs"),
        .names = "{.col}"
      ), .groups = "drop") %>%
      mutate(simulation = sim_idx)
    
    sim_correlations_list[[sim_idx]] <- sim_correlations
  }
  
  # Return list of data frames
  return(sim_correlations_list)
}

target_variable <- "depressedmood_state"
vars_to_use <- c("loneliness_state_pmc","socintsatisfaction_state_pmc","responsiveness_state_pmc",
                 "selfdisclosure_state_pmc","otherdisclosure_state_pmc")

raw_correlation_list <- extract_all_within_person_correlations(
  simulated_data_list = sim_data_list_VAR_trendlarge,
  target_var = target_variable,
  vars = vars_to_use
)

library(dplyr)

# Arrange each data frame in the list by ID (numeric sort)
raw_correlation_list <- lapply(raw_correlation_list, function(df) {
  df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID)
})

compare_raw_vs_sigma <- function(sigma_correlations, raw_correlations, correlation_threshold = 0.2) {
  
  # Define important values based on the threshold
  sigma_important <- ifelse(abs(sigma_correlations) > correlation_threshold, 1, 0)  # True important values
  raw_important <- ifelse(abs(raw_correlations) > correlation_threshold, 1, 0)      # Detected important values
  
  # Calculate True Positives, False Negatives, False Positives, and True Negatives
  TP <- as.numeric(sum(sigma_important == 1 & raw_important == 1))  # True Positives
  FN <- as.numeric(sum(sigma_important == 1 & raw_important == 0))  # False Negatives
  FP <- as.numeric(sum(sigma_important == 0 & raw_important == 1))  # False Positives
  TN <- as.numeric(sum(sigma_important == 0 & raw_important == 0))  # True Negatives
  
  # Calculate Sensitivity and Specificity
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  
  # Calculate F1 Score
  f1_score <- ifelse(
    !is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0,
    2 * (precision * sensitivity) / (precision + sensitivity),
    NA_real_
  )
  
  # Calculate MCC
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  mcc <- ifelse(
    is.finite(mcc_denom) & mcc_denom > 0,
    ((TP * TN) - (FP * FN)) / mcc_denom,
    NA_real_
  )
  
  # Return results in a data frame
  results <- data.frame(
    Metric = c("Sensitivity", "Specificity", "Precision", "F1 Score", "MCC"),
    Value = c(sensitivity, specificity, precision, f1_score, mcc)
  )
  
  return(results)
}

compare_raw_vs_sigma_all_sims <- function(sigma_correlations, raw_cor_list, corr_threshold = 0.2) {
  results_list <- list()
  
  variable_names <- setdiff(names(sigma_correlations), "ID")
  
  for (sim in seq_along(raw_cor_list)) {
    raw_df <- raw_cor_list[[sim]] %>%
      select(all_of(c("ID", variable_names)))
    
    # Align and flatten
    raw_vec <- as.numeric(as.matrix(raw_df %>% select(-ID)))
    sigma_vec <- as.numeric(as.matrix(sigma_correlations %>% select(-ID)))
    
    perf <- compare_raw_vs_sigma(
      sigma_correlations = sigma_vec,
      raw_correlations = raw_vec,
      correlation_threshold = corr_threshold
    )
    
    perf_wide <- tibble(simulation = sim) %>%
      bind_cols(t(perf$Value) %>% as.data.frame())
    colnames(perf_wide)[-1] <- perf$Metric
    results_list[[sim]] <- perf_wide
    
  }
  
  # Combine and return
  return(bind_rows(results_list))
}

# New column names (adjust to match your variables)
new_colnames <- c("ID", "loneliness","socintsatisfaction","responsiveness","selfdisclosure","otherdisclosure")

# Apply renaming across the list
raw_correlation_list <- lapply(raw_correlation_list, function(df) {
  colnames(df) <- new_colnames
  return(df)
})

raw_vs_sigma_performance <- compare_raw_vs_sigma_all_sims(
  sigma_correlations = sigmacorrelations,
  raw_cor_list = raw_correlation_list,
  corr_threshold = 0.2
)

average_results <- raw_vs_sigma_performance %>%
  select(-simulation) %>%        # Remove the simulation column
  summarise(across(everything(), mean, na.rm = TRUE))  # Compute average for each metric


print(average_results)

library(dplyr)
library(tidyr)

# 1. Remove the 'simulation' column to focus on metric columns
metrics_only <- raw_vs_sigma_performance %>% select(-simulation)

# 2. Calculate mean and 95% confidence intervals
average_metrics <- metrics_only %>%
  summarise(across(everything(), list(
    mean = ~mean(.),
    lower_CI = ~mean(.) - 1.96 * sd(.) / sqrt(n()),
    upper_CI = ~mean(.) + 1.96 * sd(.) / sqrt(n())
  )))

# 3. Optional: reshape to clearer format
# Combine results into a long-format summary
average_metrics_long <- average_metrics %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Metric", ".value"),
    names_sep = "_"
  )

# 4. View the result
print(average_metrics_long)

saveRDS(imputed_realdata, file = "imputed_realdata.rds")
