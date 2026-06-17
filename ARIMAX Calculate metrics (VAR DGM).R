library(dplyr)
library(tidyr)

df <- df %>%
  mutate(ID = row_number())

library(stringr)

iarimax_VAR_detrended <- lapply(iarimax_VAR_detrended, function(df) {
  names(df) <- names(df) %>%
    str_replace("loneliness_state_pmc", "IV_1") %>%
    str_replace("socintsatisfaction_state_pmc", "IV_2") %>%
    str_replace("responsiveness_state_pmc", "IV_3") %>%
    str_replace("selfdisclosure_state_pmc", "IV_4") %>%
    str_replace("otherdisclosure_state_pmc", "IV_5")
  df
})

# Function to calculate t-statistics for one variable in one simulation
calculate_t_statistics <- function(sim_results, variable_name) {
  tryCatch({
    xreg_coefficients <- sim_results$results_df$xreg
    standard_errors <- sim_results$results_df$stderr_xreg
    
    t_statistics <- xreg_coefficients / standard_errors
    
    return(data.frame(
      ID = sim_results$results_df$ID,
      variable = variable_name,
      t_stat = t_statistics
    ))
  }, error = function(e) {
    message(paste("Error calculating t-statistics for variable:", variable_name))
    return(NULL)
  })
}

calculate_t_statistics_wide <- function(iarimax_results_list, variable_names) {
  all_results <- list()
  
  for (sim in seq_along(iarimax_results_list)) {
    sim_results <- iarimax_results_list[[sim]]
    
    sim_t_stats <- lapply(seq_along(variable_names), function(j) {
      variable_result <- sim_results[[paste0("Sim_", sim, "_IV_", j)]]
      if (!is.null(variable_result)) {
        t_stat_data <- calculate_t_statistics(variable_result, variable_names[j])
        return(t_stat_data)
      } else {
        return(NULL)
      }
    })
    
    sim_t_stats_combined <- do.call(rbind, sim_t_stats)
    
    if (!is.null(sim_t_stats_combined)) {
      sim_t_stats_wide <- sim_t_stats_combined %>%
        pivot_wider(names_from = variable, values_from = t_stat) %>%
        mutate(simulation = sim)
      
      all_results[[sim]] <- sim_t_stats_wide
    }
  }
  
  final_wide_t_stats <- bind_rows(all_results)
  
  # ✅ Preserve original ID, just sort (don't relabel it!)
  final_wide_t_stats <- final_wide_t_stats %>%
    arrange(simulation, ID)
  
  return(final_wide_t_stats)
}

# Define variable names (full original names)
variable_names <- c("loneliness", "socintsatisfaction", "responsiveness", "selfdisclosure", "otherdisclosure")

# Run function
t_statistics_wide <- calculate_t_statistics_wide(iarimax_VAR_detrended, variable_names)

# Rename from IV_* to full names
t_statistics_wide <- t_statistics_wide %>%
  rename(
    loneliness_state = IV_1,
    social_interaction_satisfaction_state = IV_2,
    perceived_responsiveness_state = IV_3,
    self_disclosure_state = IV_4,
    other_disclosure_state = IV_5
  )

t_statistics_wide <- t_statistics_wide %>%
  mutate(ID = as.numeric(ID)) %>%           # Ensure ID is numeric
  arrange(simulation, ID)                   # Sort by simulation and ID

calculate_simulation_performance <- function(expected_correlations, t_statistics,
                                             correlation_threshold = 0.2,
                                             tstat_threshold = 1.96) {
  
  expected_important <- ifelse(abs(expected_correlations) > correlation_threshold, 1, 0)
  detected_important <- ifelse(abs(t_statistics) > tstat_threshold, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  f1_score    <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * (precision * sensitivity) / (precision + sensitivity)
  } else {
    NA_real_
  }
  
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  if (is.finite(mcc_denom) && mcc_denom > 0) {
    mcc <- ((TP * TN) - (FP * FN)) / mcc_denom
  } else {
    mcc <- NA_real_
  }
  
  c(Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1Score = f1_score,
    MCC = mcc)
}

calculate_performance_all_simulations <- function(t_stats_wide, expected_correlation_df, t_threshold = 1.96, corr_threshold = 0.2) {
  results_list <- list()
  
  # Identify variable names only (exclude ID/simulation)
  variable_names <- setdiff(colnames(expected_correlation_df), "ID")
  
  for (sim in unique(t_stats_wide$simulation)) {
    # Extract t-stats for the current simulation
    sim_data <- t_stats_wide %>%
      filter(simulation == sim)
    
    # Extract only variable columns (drop ID and simulation)
    sim_t_stats <- sim_data %>% select(all_of(variable_names))
    
    # Extract expected correlations (drop ID only)
    expected_corrs <- expected_correlation_df %>% select(all_of(variable_names))
    
    # Flatten to vectors
    t_vec <- as.numeric(as.matrix(sim_t_stats))
    cor_vec <- as.numeric(as.matrix(expected_corrs))
    
    # Calculate performance
    perf <- calculate_simulation_performance(
      expected_correlations = cor_vec,
      t_statistics = t_vec,
      correlation_threshold = corr_threshold,
      tstat_threshold = t_threshold
    )
    
    results_list[[sim]] <- as.data.frame(t(perf)) %>% mutate(simulation = sim)
  }
  
  performance_df <- bind_rows(results_list)
  return(performance_df)
}


performance_per_sim <- calculate_performance_all_simulations(
  t_stats_wide = t_statistics_wide,
  expected_correlation_df = sigmacorrelations,
  corr_threshold = 0.2
)

average_results <- performance_per_sim %>%
  select(-simulation) %>%        # Remove the simulation column
  summarise(across(everything(), mean, na.rm = TRUE))  # Compute average for each metric


print(average_results)

library(dplyr)

performance_per_sim %>%
  summarise(across(
    .cols = -simulation,
    .fns = list(
      mean = ~mean(.x, na.rm = TRUE),
      lower_CI = ~mean(.x, na.rm = TRUE) - 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
      upper_CI = ~mean(.x, na.rm = TRUE) + 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
    ),
    .names = "{.col}_{.fn}"
  ))


head(performance_per_sim)
print(average_results)

sim54_tstats <- t_statistics_wide %>%
  filter(simulation == 11) %>%
  arrange(ID) %>%
  select(-simulation)  # remove simulation column to match expected correlations

t_vec <- as.numeric(as.matrix(sim54_tstats %>% select(-ID)))
cor_vec <- as.numeric(as.matrix(sim54_expected %>% select(-ID)))

manual_perf_54 <- calculate_simulation_performance(
  expected_correlations = cor_vec,
  t_statistics = t_vec,
  correlation_threshold = 0.2,
  tstat_threshold = 1.96
)
print(manual_perf_54)
performance_per_sim %>% filter(simulation == 11)

# 1. Remove the 'simulation' column to focus on metric columns
metrics_only <- performance_per_sim %>% select(-simulation)

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
