generate_confirmed_results <- function(BORUTA_results_all) {
  # Initialize a list to store data frames for each simulation
  confirmed_results_simulations <- vector("list", length(BORUTA_results_all))
  
  # Loop through each simulation
  for (sim in seq_along(BORUTA_results_all)) {
    
    # Get the number of individuals in the current simulation
    num_individuals <- length(BORUTA_results_all[[sim]])
    
    # Initialize a data frame to store the results for all individuals in this simulation
    confirmed_results_individuals <- data.frame(
      ID = numeric(num_individuals)
    )
    
    # Loop through each individual in the simulation
    for (i in seq_along(BORUTA_results_all[[sim]])) {
      # Add individual ID (assuming it is present in the data)
      confirmed_results_individuals$ID[i] <- i  # Replace this with actual IDs if available
      
      # Check and record the final decision for each predictor
      confirmed_results_individuals$loneliness_state_pmc[i] <- ifelse(
        BORUTA_results_all[[sim]][[i]]$finalDecision["loneliness_state_pmc"] == "Confirmed", 1, 0)
      
      confirmed_results_individuals$socintsatisfaction_state_pmc[i] <- ifelse(
        BORUTA_results_all[[sim]][[i]]$finalDecision["socintsatisfaction_state_pmc"] == "Confirmed", 1, 0)
      
      confirmed_results_individuals$responsiveness_state_pmc[i] <- ifelse(
        BORUTA_results_all[[sim]][[i]]$finalDecision["responsiveness_state_pmc"] == "Confirmed", 1, 0)
      
      confirmed_results_individuals$selfdisclosure_state_pmc[i] <- ifelse(
        BORUTA_results_all[[sim]][[i]]$finalDecision["selfdisclosure_state_pmc"] == "Confirmed", 1, 0)
      
      confirmed_results_individuals$otherdisclosure_state_pmc[i] <- ifelse(
        BORUTA_results_all[[sim]][[i]]$finalDecision["otherdisclosure_state_pmc"] == "Confirmed", 1, 0)
    
    }
    
    # Store the results for this simulation
    confirmed_results_simulations[[sim]] <- confirmed_results_individuals
  }
  
  # Return the list of data frames
  return(confirmed_results_simulations)
}

# Example usage with BORUTA_results_all_50
BORUTA_Confirmed_condition <- generate_confirmed_results(BORUTA_results_all_250)

# New column names (adjust to match your variables)BORUTA_results_residuals_4var# New column names (adjust to match your variables)
new_colnames <- c("ID", "loneliness","socintsatisfaction","responsiveness","selfdisclosure","otherdisclosure")

# Apply renaming across the list
BORUTA_Confirmed_condition <- lapply(BORUTA_Confirmed_condition, function(df) {
  colnames(df) <- new_colnames
  return(df)
})

calculate_boruta_performance <- function(expected_correlations, boruta_decisions,
                                         correlation_threshold = 0.2) {
  
  if (length(expected_correlations) != length(boruta_decisions)) {
    stop("Expected correlations and BORUTA decisions must have the same length.")
  }
  
  # IMPORTANT LOGIC unchanged
  expected_important <- ifelse(abs(expected_correlations) > correlation_threshold, 1, 0)
  detected_important <- ifelse(boruta_decisions == 1, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  
  f1_score <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * (precision * sensitivity) / (precision + sensitivity)
  } else {
    NA_real_
  }
  
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  mcc <- if (is.finite(mcc_denom) && mcc_denom > 0) {
    ((TP * TN) - (FP * FN)) / mcc_denom
  } else {
    NA_real_
  }
  
  c(
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1Score = f1_score,
    MCC = mcc
  )
}

calculate_boruta_performance_all_sims <- function(boruta_confirmed_list, expected_correlation_df,
                                                  corr_threshold = 0.2) {
  results_list <- list()
  
  variable_names <- setdiff(colnames(expected_correlation_df), "ID")
  
  for (sim in seq_along(boruta_confirmed_list)) {
    sim_data <- boruta_confirmed_list[[sim]] %>%
      select(all_of(c("ID", variable_names)))
    
    expected_corrs <- expected_correlation_df %>%
      select(all_of(c("ID", variable_names)))
    
    # KEEP ORIGINAL ORDER
    boruta_vec <- as.numeric(as.matrix(sim_data %>% select(-ID)))
    cor_vec <- as.numeric(as.matrix(expected_corrs %>% select(-ID)))
    
    perf <- calculate_boruta_performance(
      expected_correlations = cor_vec,
      boruta_decisions = boruta_vec,
      correlation_threshold = corr_threshold
    )
    
    # only fix output formatting
    perf_transposed <- as.data.frame(t(perf))
    perf_transposed$simulation <- sim
    
    results_list[[sim]] <- perf_transposed
  }
  
  performance_df <- dplyr::bind_rows(results_list)
  return(performance_df)
}

library(dplyr)

# Example usage
boruta_performance_per_sim <- calculate_boruta_performance_all_sims(
  boruta_confirmed_list = BORUTA_Confirmed_condition,
  expected_correlation_df = sigma_correlations,
  corr_threshold = 0.2
)

average_results <- boruta_performance_per_sim %>%
  select(-simulation) %>%        # Remove the simulation column
  summarise(across(everything(), mean, na.rm = TRUE))  # Compute average for each metric


print(average_results)



# 1. Extract Boruta confirBORUTA_Confirmed_trendlow1# 1. Extract Boruta confirmed results for simulation 1
boruta_sim1 <- BORUTA_Confirmed_all[[15]]  # Replace with your actual list name if different

# 2. Extract expected correlations (assumed same for all simulations)
expected_corrs <- effect_size_wide  # Replace with your actual expected correlation dataframe

# 3. Make sure column order matches and only includes variables (exclude ID)
# Assume both have same variable names: c("X", "Z") for example
variable_names <- setdiff(colnames(expected_corrs), "ID")

# 4. Convert to numeric vectors for comparison
boruta_vec <- as.numeric(as.matrix(boruta_sim1 %>% select(all_of(variable_names))))
cor_vec <- as.numeric(as.matrix(expected_corrs %>% select(all_of(variable_names))))

# 5. Run performance function
manual_check_sim1 <- calculate_boruta_performance(
  expected_correlations = cor_vec,
  boruta_decisions = boruta_vec,
  correlation_threshold = 0.1  # Or your threshold
)

# 6. Print results
print(manual_check_sim1)

library(dplyr)
library(tidyr)
# 1. Remove the 'simulation' column to focus on metric columns
metrics_only <- boruta_performance_per_sim %>% select(-simulation)

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

