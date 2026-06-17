#########################################################
#### Calculating performance metrics lagged (VAR DGM)####
#########################################################

# i-ARIMAX #
library(dplyr)
library(tidyr)

iarimax_results_VAR_lagged<- lapply(iarimax_results_VAR_lagged, function(df) {
  names(df) <- names(df) %>%
    str_replace("loneliness_lag", "IV_1") %>%
    str_replace("socintsatisfaction_lag", "IV_2") %>%
    str_replace("responsiveness_lag", "IV_3") %>%
    str_replace("selfdisclosure_lag", "IV_4") %>%
    str_replace("otherdisclosure_lag", "IV_5")
  df
})

# =========================================================
# 1. Function to calculate t-statistics for one variable
# =========================================================

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
# =========================================================
# 2. Extract lagged t-statistics into wide format
# =========================================================

calculate_t_statistics_wide_lag <- function(iarimax_results_list, variable_names) {
  all_results <- list()
  
  for (sim in seq_along(iarimax_results_list)) {
    sim_results <- iarimax_results_list[[sim]]
    
    sim_t_stats <- lapply(seq_along(variable_names), function(j) {
      
      # Change this if your lagged objects are named differently
      object_name <- paste0("Sim_", sim, "_IV_", j)
      
      variable_result <- sim_results[[object_name]]
      
      if (!is.null(variable_result)) {
        t_stat_data <- calculate_t_statistics(variable_result, variable_names[j])
        return(t_stat_data)
      } else {
        message(paste("Missing:", object_name))
        return(NULL)
      }
    })
    
    sim_t_stats_combined <- do.call(rbind, sim_t_stats)
    
    if (!is.null(sim_t_stats_combined) && nrow(sim_t_stats_combined) > 0) {
      sim_t_stats_wide <- sim_t_stats_combined %>%
        pivot_wider(names_from = variable, values_from = t_stat) %>%
        mutate(simulation = sim)
      
      all_results[[sim]] <- sim_t_stats_wide
    }
  }
  
  final_wide_t_stats <- bind_rows(all_results) %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(simulation, ID)
  
  return(final_wide_t_stats)
}

variable_names_lag <- c(
  "loneliness_state_pmc",
  "socintsatisfaction_state_pmc",
  "responsiveness_state_pmc",
  "selfdisclosure_state_pmc",
  "otherdisclosure_state_pmc"
)

t_statistics_lag_wide <- calculate_t_statistics_wide_lag(
  iarimax_results_list = iarimax_results_VAR_lagged,
  variable_names = variable_names_lag
)

calculate_simulation_performance_lag <- function(expected_phi,
                                                 t_statistics,
                                                 phi_threshold = 0.05,
                                                 tstat_threshold = 1.96) {
  
  valid_idx <- !is.na(expected_phi) & !is.na(t_statistics)
  
  expected_phi <- expected_phi[valid_idx]
  t_statistics <- t_statistics[valid_idx]
  
  expected_important <- ifelse(abs(expected_phi) > phi_threshold, 1, 0)
  detected_important <- ifelse(abs(t_statistics) > tstat_threshold, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  
  f1_score <- if (!is.na(precision) && !is.na(sensitivity) && 
                  (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
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

calculate_performance_all_simulations_lag <- function(t_stats_wide,
                                                      expected_phi_df,
                                                      t_threshold = 1.96,
                                                      phi_threshold = 0.05) {
  results_list <- list()
  
  variable_names <- setdiff(colnames(expected_phi_df), "ID")
  
  expected_phi_df <- expected_phi_df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID)
  
  for (sim in unique(t_stats_wide$simulation)) {
    
    sim_data <- t_stats_wide %>%
      filter(simulation == sim) %>%
      mutate(ID = as.numeric(ID)) %>%
      arrange(ID)
    
    sim_t_stats <- sim_data %>%
      select(all_of(variable_names))
    
    expected_phi <- expected_phi_df %>%
      select(all_of(variable_names))
    
    t_vec <- as.numeric(as.matrix(sim_t_stats))
    phi_vec <- as.numeric(as.matrix(expected_phi))
    
    perf <- calculate_simulation_performance_lag(
      expected_phi = phi_vec,
      t_statistics = t_vec,
      phi_threshold = phi_threshold,
      tstat_threshold = t_threshold
    )
    
    results_list[[as.character(sim)]] <- as.data.frame(t(perf)) %>%
      mutate(simulation = sim)
  }
  
  bind_rows(results_list)
}

performance_per_sim_lag <- calculate_performance_all_simulations_lag(
  t_stats_wide = t_statistics_lag_wide,
  expected_phi_df = phi_coefficients,
  t_threshold = 1.96,
  phi_threshold = 0.12
)

average_results_lag <- performance_per_sim_lag %>%
  select(-simulation) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

print(average_results_lag)

# 1. Remove the 'simulation' column to focus on metric columns
metrics_only <- performance_per_sim_lag %>% select(-simulation)

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

#############################################################################

# MV i-ARIMAX #
library(forecast)
library(dplyr)
library(tidyr)

# =========================================================
# 1. Extract coefficients, SEs, and t-stats
# =========================================================

extract_arimax_results <- function(model_list) {
  extracted_results <- list()
  
  for (id in names(model_list)) {
    model <- model_list[[id]]
    model_summary <- summary(model)
    
    coef_matrix <- model_summary$coef
    
    if (is.null(dim(coef_matrix)) || !("Std. Error" %in% colnames(coef_matrix))) {
      estimates <- if (is.null(dim(coef_matrix))) coef_matrix else coef_matrix[, "Estimate"]
      std_errors <- sqrt(diag(model$var.coef))
    } else {
      estimates <- coef_matrix[, "Estimate"]
      std_errors <- coef_matrix[, "Std. Error"]
    }
    
    t_stats <- estimates / std_errors
    
    extracted_results[[id]] <- list(
      coefficients = estimates,
      std_errors   = std_errors,
      t_stats      = t_stats
    )
  }
  
  return(extracted_results)
}

extract_arimax_results_from_simulations <- function(simulations_list) {
  all_results <- list()
  
  for (sim in seq_along(simulations_list)) {
    simulation_data <- simulations_list[[sim]]
    extracted_results <- extract_arimax_results(simulation_data)
    all_results[[as.character(sim)]] <- extracted_results
  }
  
  return(all_results)
}

# Replace this with your actual MV i-ARIMAX lagged result object
all_simulations_results_lag <- extract_arimax_results_from_simulations(
  MV_ARIMAX_lagged_VAR
)

# =========================================================
# 2. Define lagged exogenous variables
# =========================================================

exog_vars_lag <- c(
  "loneliness_lag_PSD",
  "socintsatisfaction_lag_PSD",
  "responsiveness_lag_PSD",
  "selfdisclosure_lag_PSD",
  "otherdisclosure_lag_PSD"
)

# =========================================================
# 3. Extract lagged t-stats into wide dataframe
# =========================================================

extract_all_tstats_long_df_lag <- function(all_simulations_results, exog_vars) {
  all_rows <- list()
  
  for (sim_index in seq_along(all_simulations_results)) {
    sim_data <- all_simulations_results[[sim_index]]
    
    for (id_index in seq_along(sim_data)) {
      individual_result <- sim_data[[id_index]]
      
      if (!is.null(individual_result$t_stats)) {
        t_stats <- individual_result$t_stats
        
        # Clean coefficient names
        names(t_stats) <- sub("^xreg\\.", "", names(t_stats))
        
        available_vars <- intersect(exog_vars, names(t_stats))
        
        if (length(available_vars) > 0) {
          
          row_values <- setNames(
            as.list(rep(NA_real_, length(exog_vars))),
            exog_vars
          )
          
          row_values[available_vars] <- as.list(as.numeric(t_stats[available_vars]))
          
          all_rows[[length(all_rows) + 1]] <- data.frame(
            Simulation = sim_index,
            ID = id_index,
            row_values,
            check.names = FALSE
          )
        }
      }
    }
  }
  
  if (length(all_rows) > 0) {
    final_df <- bind_rows(all_rows)
    return(final_df[, c("Simulation", "ID", exog_vars)])
  } else {
    warning("No valid lagged t-statistics found.")
    return(data.frame())
  }
}

mv_arimax_tstats_lag_df <- extract_all_tstats_long_df_lag(
  all_simulations_results_lag,
  exog_vars_lag
)

# =========================================================
# 4. Rename lagged t-stat columns to match phi_coefficients
# =========================================================

mv_arimax_tstats_lag_df <- mv_arimax_tstats_lag_df %>%
  rename(
    loneliness_state_pmc = loneliness_lag_PSD,
    socintsatisfaction_state_pmc = socintsatisfaction_lag_PSD,
    responsiveness_state_pmc = responsiveness_lag_PSD,
    selfdisclosure_state_pmc = selfdisclosure_lag_PSD,
    otherdisclosure_state_pmc = otherdisclosure_lag_PSD,
    simulation = Simulation
  ) %>%
  mutate(ID = as.numeric(ID)) %>%
  arrange(simulation, ID)

# =========================================================
# 5. Performance function for lagged Phi truth
# =========================================================

calculate_simulation_performance_lag <- function(expected_phi,
                                                 t_statistics,
                                                 phi_threshold = 0.12,
                                                 tstat_threshold = 1.96) {
  
  valid_idx <- !is.na(expected_phi) & !is.na(t_statistics)
  
  expected_phi <- expected_phi[valid_idx]
  t_statistics <- t_statistics[valid_idx]
  
  expected_important <- ifelse(abs(expected_phi) > phi_threshold, 1, 0)
  detected_important <- ifelse(abs(t_statistics) > tstat_threshold, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  
  f1_score <- ifelse(
    (!is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0),
    2 * (precision * sensitivity) / (precision + sensitivity),
    NA_real_
  )
  
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  mcc <- ifelse(
    is.finite(mcc_denom) & mcc_denom > 0,
    ((TP * TN) - (FP * FN)) / mcc_denom,
    NA_real_
  )
  
  return(c(
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1Score = f1_score,
    MCC = mcc
  ))
}

# =========================================================
# 6. Calculate lagged performance across simulations
# =========================================================

calculate_performance_all_simulations_lag <- function(t_stats_wide,
                                                      expected_phi_df,
                                                      t_threshold = 1.96,
                                                      phi_threshold = 0.12) {
  results_list <- list()
  
  variable_names <- setdiff(colnames(expected_phi_df), "ID")
  
  expected_phi_df <- expected_phi_df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID)
  
  for (sim in unique(t_stats_wide$simulation)) {
    
    sim_data <- t_stats_wide %>%
      filter(simulation == sim) %>%
      mutate(ID = as.numeric(ID)) %>%
      arrange(ID)
    
    sim_t_stats <- sim_data %>%
      select(all_of(variable_names))
    
    expected_phi <- expected_phi_df %>%
      select(all_of(variable_names))
    
    t_vec <- as.numeric(as.matrix(sim_t_stats))
    phi_vec <- as.numeric(as.matrix(expected_phi))
    
    perf <- calculate_simulation_performance_lag(
      expected_phi = phi_vec,
      t_statistics = t_vec,
      phi_threshold = phi_threshold,
      tstat_threshold = t_threshold
    )
    
    results_list[[as.character(sim)]] <- as.data.frame(t(perf)) %>%
      mutate(simulation = sim)
  }
  
  performance_df <- bind_rows(results_list)
  return(performance_df)
}

# =========================================================
# 7. Run performance calculation using phi_coefficients
# =========================================================

performance_per_sim_lag <- calculate_performance_all_simulations_lag(
  t_stats_wide = mv_arimax_tstats_lag_df,
  expected_phi_df = phi_coefficients,
  t_threshold = 1.96,
  phi_threshold = 0.12
)

# =========================================================
# 8. Average performance
# =========================================================

average_results_lag <- performance_per_sim_lag %>%
  select(-simulation) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(average_results_lag)

# =========================================================
# 9. Average performance + 95% CIs
# =========================================================

metrics_only_lag <- performance_per_sim_lag %>%
  select(-simulation)

average_metrics_lag <- metrics_only_lag %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    lower_CI = ~ mean(.x, na.rm = TRUE) - 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
    upper_CI = ~ mean(.x, na.rm = TRUE) + 1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
  )))

average_metrics_lag_long <- average_metrics_lag %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Metric", ".value"),
    names_sep = "_"
  )

print(average_metrics_lag_long)

###############################################################################

# Boruta #
library(dplyr)

generate_confirmed_results_lag <- function(BORUTA_results_all) {
  
  confirmed_results_simulations <- vector("list", length(BORUTA_results_all))
  
  for (sim in seq_along(BORUTA_results_all)) {
    
    num_individuals <- length(BORUTA_results_all[[sim]])
    
    confirmed_results_individuals <- data.frame(
      ID = numeric(num_individuals)
    )
    
    for (i in seq_along(BORUTA_results_all[[sim]])) {
      
      confirmed_results_individuals$ID[i] <- i
      
      fd <- BORUTA_results_all[[sim]][[i]]$finalDecision
      
      confirmed_results_individuals$loneliness_state_pmc[i] <- ifelse(
        fd["loneliness_lag"] == "Confirmed", 1, 0
      )
      
      confirmed_results_individuals$socintsatisfaction_state_pmc[i] <- ifelse(
        fd["socintsatisfaction_lag"] == "Confirmed", 1, 0
      )
      
      confirmed_results_individuals$responsiveness_state_pmc[i] <- ifelse(
        fd["responsiveness_lag"] == "Confirmed", 1, 0
      )
      
      confirmed_results_individuals$selfdisclosure_state_pmc[i] <- ifelse(
        fd["selfdisclosure_lag"] == "Confirmed", 1, 0
      )
      
      confirmed_results_individuals$otherdisclosure_state_pmc[i] <- ifelse(
        fd["otherdisclosure_lag"] == "Confirmed", 1, 0
      )
    }
    
    confirmed_results_simulations[[sim]] <- confirmed_results_individuals
  }
  
  return(confirmed_results_simulations)
}

BORUTA_Confirmed_lag_condition <- generate_confirmed_results_lag(
  tsBoruta_results_VAR_lagged
)

calculate_boruta_performance_lag <- function(expected_phi,
                                             boruta_decisions,
                                             phi_threshold = 0.12) {
  
  if (length(expected_phi) != length(boruta_decisions)) {
    stop("Expected Phi values and Boruta decisions must have the same length.")
  }
  
  expected_important <- ifelse(abs(expected_phi) > phi_threshold, 1, 0)
  detected_important <- ifelse(boruta_decisions == 1, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  
  f1_score <- if (!is.na(precision) && !is.na(sensitivity) && 
                  (precision + sensitivity) > 0) {
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

calculate_boruta_performance_all_sims_lag <- function(boruta_confirmed_list,
                                                      expected_phi_df,
                                                      phi_threshold = 0.12) {
  results_list <- list()
  
  variable_names <- setdiff(colnames(expected_phi_df), "ID")
  
  expected_phi_df <- expected_phi_df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID)
  
  for (sim in seq_along(boruta_confirmed_list)) {
    
    sim_data <- boruta_confirmed_list[[sim]] %>%
      mutate(ID = as.numeric(ID)) %>%
      arrange(ID) %>%
      select(all_of(c("ID", variable_names)))
    
    expected_phi <- expected_phi_df %>%
      select(all_of(c("ID", variable_names)))
    
    boruta_vec <- as.numeric(as.matrix(sim_data %>% select(-ID)))
    phi_vec    <- as.numeric(as.matrix(expected_phi %>% select(-ID)))
    
    perf <- calculate_boruta_performance_lag(
      expected_phi = phi_vec,
      boruta_decisions = boruta_vec,
      phi_threshold = phi_threshold
    )
    
    perf_transposed <- as.data.frame(t(perf))
    perf_transposed$simulation <- sim
    
    results_list[[sim]] <- perf_transposed
  }
  
  performance_df <- dplyr::bind_rows(results_list)
  return(performance_df)
}

boruta_lag_performance_per_sim <- calculate_boruta_performance_all_sims_lag(
  boruta_confirmed_list = BORUTA_Confirmed_lag_condition,
  expected_phi_df = phi_coefficients,
  phi_threshold = 0.12
)

average_results_boruta_lag <- boruta_lag_performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

print(average_results_boruta_lag)

average_metrics_boruta_lag <- boruta_lag_performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    lower_CI = ~ mean(.x, na.rm = TRUE) -
      1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
    upper_CI = ~ mean(.x, na.rm = TRUE) +
      1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
  )))

average_metrics_boruta_lag_long <- average_metrics_boruta_lag %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("Metric", ".value"),
    names_sep = "_"
  )

print(average_metrics_boruta_lag_long)

################################################################################

# GIMME #
library(dplyr)
library(tidyr)

# Define the parent directory containing the GIMME output folders
parent_directory <- "C:/Users/WillLi/Documents/mlvar simulation/output_directory_70timeGIMME"

# Function to process all simulation folders for lagged paths
process_simulation_folders_lag <- function(directory) {
  
  subfolders <- list.dirs(directory, recursive = FALSE)
  all_simulation_results <- list()
  
  for (subfolder in subfolders) {
    
    csv_file <- file.path(subfolder, "indivPathEstimates.csv")
    
    if (file.exists(csv_file)) {
      
      data <- read.csv(csv_file)
      
      relevant_paths <- data %>%
        filter(lhs == "depressedmood_state") %>%
        filter(rhs %in% c(
          "loneliness_state_pmclag",
          "socintsatisfaction_state_pmclag",
          "responsiveness_state_pmclag",
          "selfdisclosure_state_pmclag",
          "otherdisclosure_state_pmclag"
        )) %>%
        select(ID = file, rhs, beta) %>%
        mutate(simulation = basename(subfolder))
      
      all_simulation_results[[subfolder]] <- relevant_paths
    }
  }
  
  do.call(rbind, all_simulation_results)
}

# Process all folders
all_paths_lag <- process_simulation_folders_lag(parent_directory)

# Lagged variable names as they appear in GIMME output
all_lag_variables <- c(
  "loneliness_state_pmclag",
  "socintsatisfaction_state_pmclag",
  "responsiveness_state_pmclag",
  "selfdisclosure_state_pmclag",
  "otherdisclosure_state_pmclag"
)

# Extract unique IDs and simulations
all_IDs <- unique(all_paths_lag$ID)
all_sims <- unique(all_paths_lag$simulation)

# Create full grid
full_grid_lag <- expand.grid(
  simulation = all_sims,
  ID = all_IDs,
  rhs = all_lag_variables,
  stringsAsFactors = FALSE
)

# Join and create binary variable
full_paths_binary_lag <- full_grid_lag %>%
  left_join(all_paths_lag, by = c("simulation", "ID", "rhs")) %>%
  mutate(path_present = ifelse(!is.na(beta), 1, 0)) %>%
  select(simulation, ID, rhs, path_present)

# Wide binary GIMME lagged path dataframe
gimme_paths_wide_binary_lag <- full_paths_binary_lag %>%
  pivot_wider(
    names_from = rhs,
    values_from = path_present
  ) %>%
  arrange(simulation, ID)

# Rename lagged GIMME variables back to match phi_coefficients
gimme_paths_wide_binary_lag <- gimme_paths_wide_binary_lag %>%
  rename(
    loneliness_state_pmc = loneliness_state_pmclag,
    socintsatisfaction_state_pmc = socintsatisfaction_state_pmclag,
    responsiveness_state_pmc = responsiveness_state_pmclag,
    selfdisclosure_state_pmc = selfdisclosure_state_pmclag,
    otherdisclosure_state_pmc = otherdisclosure_state_pmclag
  )

# Convert simulation number
gimme_paths_wide_binary_lag$simulation <- as.numeric(
  sub("^Sim", "", gimme_paths_wide_binary_lag$simulation)
)

# Make sure ID is numeric
gimme_paths_wide_binary_lag <- gimme_paths_wide_binary_lag %>%
  mutate(ID = as.numeric(ID)) %>%
  arrange(simulation, ID)

calculate_simulation_performance_binary_lag <- function(expected_phi,
                                                        binary_decisions,
                                                        phi_threshold = 0.12) {
  
  expected_important <- ifelse(abs(expected_phi) > phi_threshold, 1, 0)
  detected_important <- ifelse(binary_decisions == 1, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  
  f1_score <- ifelse(
    (!is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0),
    2 * (precision * sensitivity) / (precision + sensitivity),
    NA_real_
  )
  
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  mcc <- ifelse(
    is.finite(mcc_denom) & mcc_denom > 0,
    ((TP * TN) - (FP * FN)) / mcc_denom,
    NA_real_
  )
  
  return(c(
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1Score = f1_score,
    MCC = mcc
  ))
}

calculate_gimme_performance_all_simulations_lag <- function(binary_wide_df,
                                                            expected_phi_df,
                                                            phi_threshold = 0.12) {
  results_list <- list()
  
  variable_names <- setdiff(colnames(expected_phi_df), "ID")
  
  expected_phi_df <- expected_phi_df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID)
  
  for (sim in unique(binary_wide_df$simulation)) {
    
    sim_data <- binary_wide_df %>%
      filter(simulation == sim) %>%
      mutate(ID = as.numeric(ID)) %>%
      arrange(ID)
    
    binary_decisions <- sim_data %>%
      select(all_of(variable_names))
    
    expected_phi <- expected_phi_df %>%
      select(all_of(variable_names))
    
    binary_vec <- as.numeric(as.matrix(binary_decisions))
    phi_vec <- as.numeric(as.matrix(expected_phi))
    
    perf <- calculate_simulation_performance_binary_lag(
      expected_phi = phi_vec,
      binary_decisions = binary_vec,
      phi_threshold = phi_threshold
    )
    
    results_list[[as.character(sim)]] <- as.data.frame(t(perf)) %>%
      mutate(simulation = sim)
  }
  
  performance_df <- bind_rows(results_list)
  return(performance_df)
}

gimme_lag_performance_per_sim <- calculate_gimme_performance_all_simulations_lag(
  binary_wide_df = gimme_paths_wide_binary_lag,
  expected_phi_df = phi_coefficients,
  phi_threshold = 0.12
)

gimme_lag_average_results <- gimme_lag_performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

print(gimme_lag_average_results)

metrics_only_lag <- gimme_lag_performance_per_sim %>%
  select(-simulation)

average_metrics_lag <- metrics_only_lag %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    lower_CI = ~ mean(.x, na.rm = TRUE) -
      1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
    upper_CI = ~ mean(.x, na.rm = TRUE) +
      1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
  )))

average_metrics_lag_long <- average_metrics_lag %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Metric", ".value"),
    names_sep = "_"
  )

print(average_metrics_lag_long)

###################################################################################

# indSEM #
library(dplyr)
library(tidyr)

# =========================================================
# STEP 1 - Process indSEM folders for LAGGED paths
# =========================================================

process_indSEM_simulations_lag <- function(directory) {
  
  sim_folders <- list.dirs(directory, recursive = FALSE)
  all_simulation_results <- list()
  column_names <- NULL
  
  for (sim_folder in sim_folders) {
    
    individual_folder <- file.path(sim_folder, "individual")
    
    if (dir.exists(individual_folder)) {
      message(paste("Processing folder:", individual_folder))
      
      beta_files <- list.files(
        individual_folder,
        pattern = "\\d+BetasStd\\.csv$",
        full.names = TRUE
      )
      
      if (length(beta_files) == 0) {
        message(paste("No beta files found in:", individual_folder))
      }
      
      for (beta_file in beta_files) {
        message(paste("Processing file:", beta_file))
        
        individual_id <- sub("BetasStd\\.csv$", "", basename(beta_file))
        
        beta_matrix <- tryCatch({
          read.csv(beta_file, row.names = 1, check.names = FALSE)
        }, error = function(e) {
          message(paste("Error reading file:", beta_file, ":", e$message))
          NULL
        })
        
        if (is.null(beta_matrix)) next
        
        if (is.null(column_names)) {
          column_names <- colnames(beta_matrix)
          message("Column names detected:", paste(column_names, collapse = ", "))
        }
        
        if ("depressedmood_state" %in% rownames(beta_matrix)) {
          
          depressedmood_state_row <- beta_matrix["depressedmood_state", ]
          
          beta_data <- data.frame(
            ID = as.numeric(individual_id),
            rhs = names(depressedmood_state_row),
            beta = as.numeric(depressedmood_state_row),
            simulation = basename(sim_folder),
            stringsAsFactors = FALSE
          )
          
          all_simulation_results[[paste0(sim_folder, "_", individual_id)]] <- beta_data
          
        } else {
          message(paste("Row 'depressedmood_state' not found in file:", beta_file))
        }
      }
      
    } else {
      message(paste("Folder does not exist:", individual_folder))
    }
  }
  
  if (length(all_simulation_results) > 0) {
    
    combined_data <- do.call(rbind, all_simulation_results)
    
    if (!is.null(column_names)) {
      combined_data <- combined_data %>%
        mutate(rhs = factor(rhs, levels = column_names))
    }
    
    return(combined_data)
    
  } else {
    message("No valid data found in the simulations.")
    return(NULL)
  }
}

# =========================================================
# STEP 2 - Run extraction
# =========================================================

parent_directory <- "C:/Users/WillLi/Documents/mlvar simulation/output_directory_70timeindSEM"

all_paths_indSEM_lag <- process_indSEM_simulations_lag(parent_directory)

# =========================================================
# STEP 3 - Define lagged variables as named in indSEM output
# =========================================================

all_lag_variables <- c(
  "loneliness_state_pmclag",
  "socintsatisfaction_state_pmclag",
  "responsiveness_state_pmclag",
  "selfdisclosure_state_pmclag",
  "otherdisclosure_state_pmclag"
)

# =========================================================
# STEP 4 - Filter lagged paths, convert to binary, pivot wide
# =========================================================

indSEM_paths_wide_binary_lag <- all_paths_indSEM_lag %>%
  filter(rhs %in% all_lag_variables) %>%
  mutate(path_present = ifelse(beta != 0, 1L, 0L)) %>%
  select(simulation, ID, rhs, path_present) %>%
  pivot_wider(
    names_from = rhs,
    values_from = path_present,
    values_fill = 0
  ) %>%
  arrange(simulation, ID)

# Convert simulation labels like Sim1 to numeric
indSEM_paths_wide_binary_lag$simulation <- as.numeric(
  sub("^Sim", "", indSEM_paths_wide_binary_lag$simulation)
)

# Rename lagged variables back to match phi_coefficients column names
indSEM_paths_wide_binary_lag <- indSEM_paths_wide_binary_lag %>%
  rename(
    loneliness_state_pmc = loneliness_state_pmclag,
    socintsatisfaction_state_pmc = socintsatisfaction_state_pmclag,
    responsiveness_state_pmc = responsiveness_state_pmclag,
    selfdisclosure_state_pmc = selfdisclosure_state_pmclag,
    otherdisclosure_state_pmc = otherdisclosure_state_pmclag
  ) %>%
  mutate(ID = as.numeric(ID)) %>%
  arrange(simulation, ID)

# =========================================================
# STEP 5 - Performance function using Phi threshold
# =========================================================

calculate_simulation_performance_binary_lag <- function(expected_phi,
                                                        binary_decisions,
                                                        phi_threshold = 0.12) {
  
  expected_important <- ifelse(abs(expected_phi) > phi_threshold, 1, 0)
  detected_important <- ifelse(binary_decisions == 1, 1, 0)
  
  TP <- as.numeric(sum(expected_important == 1 & detected_important == 1))
  FN <- as.numeric(sum(expected_important == 1 & detected_important == 0))
  FP <- as.numeric(sum(expected_important == 0 & detected_important == 1))
  TN <- as.numeric(sum(expected_important == 0 & detected_important == 0))
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  precision   <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  
  f1_score <- ifelse(
    (!is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0),
    2 * (precision * sensitivity) / (precision + sensitivity),
    NA_real_
  )
  
  mcc_denom <- as.numeric(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  mcc <- ifelse(
    is.finite(mcc_denom) & mcc_denom > 0,
    ((TP * TN) - (FP * FN)) / mcc_denom,
    NA_real_
  )
  
  return(c(
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1Score = f1_score,
    MCC = mcc
  ))
}

# =========================================================
# STEP 6 - Calculate performance across simulations
# =========================================================

calculate_indSEM_performance_all_simulations_lag <- function(binary_wide_df,
                                                             expected_phi_df,
                                                             phi_threshold = 0.12) {
  results_list <- list()
  
  variable_names <- setdiff(colnames(expected_phi_df), "ID")
  
  expected_phi_df <- expected_phi_df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID)
  
  for (sim in unique(binary_wide_df$simulation)) {
    
    sim_data <- binary_wide_df %>%
      filter(simulation == sim) %>%
      mutate(ID = as.numeric(ID)) %>%
      arrange(ID)
    
    binary_decisions <- sim_data %>%
      select(all_of(variable_names))
    
    expected_phi <- expected_phi_df %>%
      select(all_of(variable_names))
    
    binary_vec <- as.numeric(as.matrix(binary_decisions))
    phi_vec <- as.numeric(as.matrix(expected_phi))
    
    perf <- calculate_simulation_performance_binary_lag(
      expected_phi = phi_vec,
      binary_decisions = binary_vec,
      phi_threshold = phi_threshold
    )
    
    results_list[[as.character(sim)]] <- as.data.frame(t(perf)) %>%
      mutate(simulation = sim)
  }
  
  performance_df <- bind_rows(results_list)
  return(performance_df)
}

# =========================================================
# STEP 7 - Run lagged indSEM performance using phi_coefficients
# =========================================================

indSEM_lag_performance_per_sim <- calculate_indSEM_performance_all_simulations_lag(
  binary_wide_df = indSEM_paths_wide_binary_lag,
  expected_phi_df = phi_coefficients,
  phi_threshold = 0.12
)

# =========================================================
# STEP 8 - Average results
# =========================================================

indSEM_lag_average_results <- indSEM_lag_performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))

print(indSEM_lag_average_results)

# =========================================================
# STEP 9 - Average results + 95% CIs
# =========================================================

metrics_only_lag <- indSEM_lag_performance_per_sim %>%
  select(-simulation)

average_metrics_lag <- metrics_only_lag %>%
  summarise(across(everything(), list(
    mean = ~ mean(.x, na.rm = TRUE),
    lower_CI = ~ mean(.x, na.rm = TRUE) -
      1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
    upper_CI = ~ mean(.x, na.rm = TRUE) +
      1.96 * sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
  )))

average_metrics_lag_long <- average_metrics_lag %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Metric", ".value"),
    names_sep = "_"
  )

print(average_metrics_lag_long)

