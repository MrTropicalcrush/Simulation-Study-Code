library(forecast)
extract_arimax_results <- function(model_list) {
  extracted_results <- list()
  
  for (id in names(model_list)) {
    model <- model_list[[id]]
    model_summary <- summary(model)
    
    # Try extracting the coefficients table
    coef_matrix <- model_summary$coef
    
    # If the coefficients table doesn't include "Std. Error", compute them from var.coef.
    if (is.null(dim(coef_matrix)) || !("Std. Error" %in% colnames(coef_matrix))) {
      # If coef_matrix is a vector, use it directly for estimates.
      estimates <- if (is.null(dim(coef_matrix))) coef_matrix else coef_matrix[, "Estimate"]
      std_errors <- sqrt(diag(model$var.coef))
    } else {
      estimates <- coef_matrix[, "Estimate"]
      std_errors <- coef_matrix[, "Std. Error"]
    }
    
    # Calculate t-statistics for each coefficient
    t_stats <- estimates / std_errors
    
    # Store the results in the list
    extracted_results[[id]] <- list(
      coefficients = estimates,
      std_errors   = std_errors,
      t_stats      = t_stats
    )
  }
  
  return(extracted_results)
}

# Function to process 100 simulations
extract_arimax_results_from_simulations <- function(simulations_list) {
  all_results <- list()
  
  # Loop over each simulation
  for (sim in seq_along(simulations_list)) {
    simulation_data <- simulations_list[[sim]]
    
    # Extract results for this simulation
    extracted_results <- extract_arimax_results(simulation_data)
    
    # Store the results for this simulation
    all_results[[as.character(sim)]] <- extracted_results
  }
  
  return(all_results)
}

# Assuming your 100 simulations are stored in a list called simulations_list
# Each element of simulations_list is a list of 102 individual ARIMAX results
all_simulations_results <- extract_arimax_results_from_simulations(MVarimax_autocormedium)

# Define the exogenous variables of interest
exog_vars <- c("loneliness_state_pmc_PSD", "socintsatisfaction_state_pmc_PSD", 
                 "responsiveness_state_pmc_PSD", "selfdisclosure_state_pmc_PSD", "otherdisclosure_state_pmc_PSD")

extract_all_tstats_long_df <- function(all_simulations_results, exog_vars) {
  all_rows <- list()
  
  for (sim_index in seq_along(all_simulations_results)) {
    sim_data <- all_simulations_results[[sim_index]]
    
    for (id_index in seq_along(sim_data)) {
      individual_result <- sim_data[[id_index]]
      
      # Check if t_stats exists
      if (!is.null(individual_result$t_stats)) {
        t_stats <- individual_result$t_stats
        
        # Clean coefficient names (remove "xreg." prefix)
        names(t_stats) <- sub("^xreg\\.", "", names(t_stats))
        
        # Ensure all exogenous variables are present
        if (all(exog_vars %in% names(t_stats))) {
          filtered_stats <- t_stats[exog_vars]
          
          # Store as a row in the list
          all_rows[[length(all_rows) + 1]] <- data.frame(
            Simulation = sim_index,
            ID = id_index,
            as.list(filtered_stats)
          )
        }
      }
    }
  }
  
  # Combine all rows into a single dataframe
  if (length(all_rows) > 0) {
    final_df <- do.call(rbind, all_rows)
    return(final_df[, c("Simulation", "ID", exog_vars)])
  } else {
    warning("No valid t-statistics found.")
    return(data.frame())
  }
}

# Run the function
mv_arimax_tstats_long_df <- extract_all_tstats_long_df(all_simulations_results, exog_vars)

mv_arimax_tstats_long_df <- mv_arimax_tstats_long_df %>%
  rename(
    loneliness = loneliness_state_pmc_PSD,
    socintsatisfaction = socintsatisfaction_state_pmc_PSD,
    responsiveness = responsiveness_state_pmc_PSD,
    selfdisclosure = selfdisclosure_state_pmc_PSD,
    otherdisclosure = otherdisclosure_state_pmc_PSD,
    simulation = Simulation
    # Add more as needed
  )

# Function to calculate performance metrics for one simulation
calculate_simulation_performance <- function(expected_correlations, t_statistics, correlation_threshold = 0.2, tstat_threshold = 1.96) {
  expected_important <- ifelse(abs(expected_correlations) > correlation_threshold, 1, 0)
  detected_important <- ifelse(abs(t_statistics) > tstat_threshold, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  f1_score    <- ifelse((!is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0),
                        2 * (precision * sensitivity) / (precision + sensitivity), NA_real_)
  
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  mcc <- ifelse(
    is.finite(mcc_denom) & mcc_denom > 0,
    ((TP * TN) - (FP * FN)) / mcc_denom,
    NA_real_
  )
  
  return(c(Sensitivity = sensitivity,
           Specificity = specificity,
           Precision = precision,
           F1Score = f1_score,
           MCC = mcc))
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
  t_stats_wide = mv_arimax_tstats_long_df,
  expected_correlation_df = sigmacorrelations,
  corr_threshold = 0.2
)

library(dplyr)

average_results <- performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(average_results)



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



























