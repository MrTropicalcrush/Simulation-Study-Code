library(dplyr)
library(tidyr)
# STEP 1 - Process Sim folders #
# Define the parent directory containing the batch folder
parent_directory <- "C:/Users/WillLi/Documents/mlvar simulation/output_directory_GIMMEdepressedhigh"  # Replace with your actual path

# Function to process all simulation folders
process_simulation_folders <- function(directory) {
  # List all subfolders in the batch directory
  subfolders <- list.dirs(directory, recursive = FALSE)
  
  # Initialize a list to store results for all simulations
  all_simulation_results <- list()
  
  for (subfolder in subfolders) {
    # Construct the full path to the CSV file
    csv_file <- file.path(subfolder, "indivPathEstimates.csv")
    
    # Check if the file exists (skip if missing)
    if (file.exists(csv_file)) {
      # Read the CSV file
      data <- read.csv(csv_file)
      
      # Extract relevant paths where lhs is "depressedmood_state"
      relevant_paths <- data %>%
        filter(lhs == "depressedmood_state") %>%  # Filter for the desired dependent variable
        select(ID = file, rhs, beta)  # Select relevant columns (adjust if needed)
      
      # Add simulation identifier
      relevant_paths <- relevant_paths %>%
        mutate(simulation = basename(subfolder))  # Use subfolder name as simulation ID
      
      # Append to results list
      all_simulation_results[[subfolder]] <- relevant_paths
    }
  }
  
  # Combine results from all simulations into one dataframe
  do.call(rbind, all_simulation_results)
}

# Process all folders in the batch directory
all_paths <- process_simulation_folders(parent_directory)

all_variables <- c("loneliness_state_pmc", "socintsatisfaction_state_pmc", 
                   "responsiveness_state_pmc", "selfdisclosure_state_pmc", "otherdisclosure_state_pmc")

# Extract unique IDs and simulations
all_IDs <- unique(all_paths$ID)
all_sims <- unique(all_paths$simulation)

# Create a full combination grid
full_grid <- expand.grid(
  simulation = all_sims,
  ID = all_IDs,
  rhs = all_variables,
  stringsAsFactors = FALSE
)

# Join and create binary variable
full_paths_binary <- full_grid %>%
  left_join(all_paths, by = c("simulation", "ID", "rhs")) %>%
  mutate(path_present = ifelse(!is.na(beta), 1, 0)) %>%  # <- THIS is the key change
  select(simulation, ID, rhs, path_present)

gimme_paths_wide_binary <- full_paths_binary %>%
  pivot_wider(
    names_from = rhs,
    values_from = path_present
  ) %>%
  arrange(simulation, ID)

gimme_paths_wide_binary <- gimme_paths_wide_binary %>%
  rename(
    loneliness = loneliness_state_pmc,
    socintsatisfaction = socintsatisfaction_state_pmc,
    responsiveness = responsiveness_state_pmc,
    selfdisclosure = selfdisclosure_state_pmc,
    otherdisclosure = otherdisclosure_state_pmc
    # Add more as needed
  )

gimme_paths_wide_binary$simulation <- as.numeric(sub("^Sim", "", gimme_paths_wide_binary$simulation))

# Modified version of the simulation performance function
calculate_simulation_performance_binary <- function(expected_correlations, binary_decisions, correlation_threshold = 0.2) {
  expected_important <- ifelse(abs(expected_correlations) > correlation_threshold, 1, 0)
  detected_important <- ifelse(binary_decisions == 1, 1, 0)  # already binary
  
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

calculate_gimme_performance_all_simulations <- function(binary_wide_df, expected_correlation_df, corr_threshold = 0.2) {
  results_list <- list()
  
  # Identify variable names (excluding ID)
  variable_names <- setdiff(colnames(expected_correlation_df), "ID")
  
  for (sim in unique(binary_wide_df$simulation)) {
    # Filter data for current simulation
    sim_data <- binary_wide_df %>%
      filter(simulation == sim)
    
    binary_decisions <- sim_data %>%
      select(all_of(variable_names))
    
    expected_corrs <- expected_correlation_df %>%
      select(all_of(variable_names))
    
    binary_vec <- as.numeric(as.matrix(binary_decisions))
    cor_vec <- as.numeric(as.matrix(expected_corrs))
    
    perf <- calculate_simulation_performance_binary(
      expected_correlations = cor_vec,
      binary_decisions = binary_vec,
      correlation_threshold = corr_threshold
    )
    
    results_list[[sim]] <- as.data.frame(t(perf)) %>%
      mutate(simulation = sim)
  }
  
  performance_df <- bind_rows(results_list)
  return(performance_df)
}


gimme_performance_per_sim <- calculate_gimme_performance_all_simulations(
  binary_wide_df = gimme_paths_wide_binary,
  expected_correlation_df = sigmacorrelations,
  corr_threshold = 0.2
)


# Average across all simulations
gimme_average_results <- gimme_performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(gimme_average_results)

# 1. Remove the 'simulation' column to focus on metric columns
metrics_only <- gimme_performance_per_sim %>% select(-simulation)

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
###########################################################################################################
# STEP 1 - Process indSEM Folders #
process_indSEM_simulations <- function(directory) {
  # List all simulation folders
  sim_folders <- list.dirs(directory, recursive = FALSE)
  
  # Initialize a list to store results for all simulations
  all_simulation_results <- list()
  column_names <- NULL  # Initialize column names as NULL
  
  for (sim_folder in sim_folders) {
    # Construct the path to the individual folder
    individual_folder <- file.path(sim_folder, "individual")
    
    # Check if the individual folder exists
    if (dir.exists(individual_folder)) {
      message(paste("Processing folder:", individual_folder))
      
      # List all individual beta files (e.g., 1BetasStd.csv, 2BetasStd.csv, etc.)
      beta_files <- list.files(individual_folder, pattern = "\\d+BetasStd\\.csv$", full.names = TRUE)
      
      if (length(beta_files) == 0) {
        message(paste("No beta files found in:", individual_folder))
      }
      
      # Process each beta file
      for (beta_file in beta_files) {
        message(paste("Processing file:", beta_file))
        
        # Extract the individual ID from the filename
        individual_id <- sub("BetasStd\\.csv$", "", basename(beta_file))
        
        # Attempt to read the beta matrix
        beta_matrix <- tryCatch({
          read.csv(beta_file, row.names = 1)
        }, error = function(e) {
          message(paste("Error reading file:", beta_file, ":", e$message))
          NULL
        })
        
        # Skip processing if beta_matrix is NULL
        if (is.null(beta_matrix)) next
        
        # Extract column names from the first valid beta matrix
        if (is.null(column_names)) {
          column_names <- colnames(beta_matrix)
          message("Column names detected:", paste(column_names, collapse = ", "))
        }
        
        # Extract the row for depressedmood_state
        if ("depressedmood_state" %in% rownames(beta_matrix)) {
          depressedmood_state_row <- beta_matrix["depressedmood_state", ]
          
          # Convert the row to a long format
          beta_data <- data.frame(
            ID = as.numeric(individual_id),
            rhs = names(depressedmood_state_row),   # Names of the columns (variables)
            beta = as.numeric(depressedmood_state_row),  # Beta coefficients
            simulation = basename(sim_folder)  # Add simulation identifier
          )
          
          # Append to results
          all_simulation_results[[paste0(sim_folder, "_", individual_id)]] <- beta_data
        } else {
          message(paste("Row 'DV' not found in file:", beta_file))
        }
      }
    } else {
      message(paste("Folder does not exist:", individual_folder))
    }
  }
  
  # Combine results from all simulations into one dataframe
  if (length(all_simulation_results) > 0) {
    combined_data <- do.call(rbind, all_simulation_results)
    
    # Ensure column names match the BetasStd.csv variables
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

# Process all folders in the batch directory
parent_directory <- "C:/Users/WillLi/Documents/mlvar simulation/output_directory_indSEMdepressedhigh"
all_paths_indSEM <- process_indSEM_simulations(parent_directory)

library(dplyr)
library(tidyr)

# Define the variables you care about
all_variables <- c("loneliness_state_pmc", "socintsatisfaction_state_pmc", 
                   "responsiveness_state_pmc", "selfdisclosure_state_pmc", "otherdisclosure_state_pmc")


# Filter, convert beta to binary, then pivot wide
indSEM_paths_wide_binary <- all_paths_indSEM %>%
  filter(rhs %in% all_variables) %>%
  mutate(path_present = ifelse(beta != 0, 1, 0)) %>%
  select(simulation, ID, rhs, path_present) %>%
  tidyr::pivot_wider(
    names_from = rhs,
    values_from = path_present
  ) %>%
  arrange(simulation, ID)

indSEM_paths_wide_binary$simulation <- as.numeric(sub("^Sim", "", indSEM_paths_wide_binary$simulation))

indSEM_paths_wide_binary <- indSEM_paths_wide_binary %>%
  rename(
    loneliness = loneliness_state_pmc,
    socintsatisfaction = socintsatisfaction_state_pmc,
    responsiveness = responsiveness_state_pmc,
    selfdisclosure = selfdisclosure_state_pmc,
    otherdisclosure = otherdisclosure_state_pmc
    # Add more as needed
  )

# Modified version of the simulation performance function
calculate_simulation_performance_binary <- function(expected_correlations, binary_decisions, correlation_threshold = 0.2) {
  expected_important <- ifelse(abs(expected_correlations) > correlation_threshold, 1, 0)
  detected_important <- ifelse(binary_decisions == 1, 1, 0)  # already binary
  
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

calculate_gimme_performance_all_simulations <- function(binary_wide_df, expected_correlation_df, corr_threshold = 0.2) {
  results_list <- list()
  
  # Identify variable names (excluding ID)
  variable_names <- setdiff(colnames(expected_correlation_df), "ID")
  
  for (sim in unique(binary_wide_df$simulation)) {
    # Filter data for current simulation
    sim_data <- binary_wide_df %>%
      filter(simulation == sim)
    
    binary_decisions <- sim_data %>%
      select(all_of(variable_names))
    
    expected_corrs <- expected_correlation_df %>%
      select(all_of(variable_names))
    
    binary_vec <- as.numeric(as.matrix(binary_decisions))
    cor_vec <- as.numeric(as.matrix(expected_corrs))
    
    perf <- calculate_simulation_performance_binary(
      expected_correlations = cor_vec,
      binary_decisions = binary_vec,
      correlation_threshold = corr_threshold
    )
    
    results_list[[sim]] <- as.data.frame(t(perf)) %>%
      mutate(simulation = sim)
  }
  
  performance_df <- bind_rows(results_list)
  return(performance_df)
}


gimme_performance_per_sim <- calculate_gimme_performance_all_simulations(
  binary_wide_df = indSEM_paths_wide_binary,
  expected_correlation_df = sigmacorrelations,
  corr_threshold = 0.2
)


# Average across all simulations
gimme_average_results <- gimme_performance_per_sim %>%
  select(-simulation) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

print(gimme_average_results)

# 1. Remove the 'simulation' column to focus on metric columns
metrics_only <- gimme_performance_per_sim %>% select(-simulation)

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
