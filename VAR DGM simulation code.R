#####################################################
#### Simulating linear data (based on real data) ####
#####################################################
#### Cleaning real data ####
rawsimdata$pid <- match(rawsimdata$pid, unique(rawsimdata$pid))
realdata <- Projectpingemadata

realdata$pid <- match(realdata$pid, unique(realdata$pid))

# Create data frame with only variables that we need
rawsimdata <- rawsimdata[, c("pid", "pingTotal", "depressedmood_state", "loneliness_state_pmc", "socintsatisfaction_state_pmc", "responsiveness_state_pmc", "selfdisclosure_state_pmc", "otherdisclosure_state_pmc")]
realdata <- realdata[, c("pid", "pingTotal", "depressedmood_state", "loneliness_state_pmc", "socintsatisfaction_state_pmc", "responsiveness_state_pmc", "selfdisclosure_state_pmc", "otherdisclosure_state_pmc")]

library(dplyr)

# Convert columns to numeric and ensure NA values are present where needed
cols_of_interest <- c("depressedmood_state", "loneliness_state_pmc", "socintsatisfaction_state_pmc", "responsiveness_state_pmc", "selfdisclosure_state_pmc", "otherdisclosure_state_pmc")
realdata <- realdata %>%
  mutate(across(all_of(cols_of_interest), ~ as.numeric(as.character(.))))

# Fill missing values using idionomics package
install.packages("devtools")
devtools::install_github("cristobalehc/idionomics")
library(idionomics)

imputed <- imputatron_2000(data = realdata, id_col = "pid", time_col = "pingTotal", cols_of_interest = c("depressedmood_state", "loneliness_state_pmc", "socintsatisfaction_state_pmc", "responsiveness_state_pmc", "selfdisclosure_state_pmc", "otherdisclosure_state_pmc"))

imputed_realdata <- imputed$imputed_df

# Fit a VAR model on real data
# Load necessary package
library(MTS)

# Split by individuals into a list
realdata_list <- split(imputed_realdata, imputed_realdata$pid)

# Reorder the list numerically
realdata_list <- realdata_list[order(as.numeric(names(realdata_list)))]


#### Fitting VAR model to data and them simulating new data ####
# Function to fit a VAR model (p=1) and extract Phi and Sigma
fit_VAR_individual <- function(realdata_list) {
  # Remove ID and time columns
  individual_matrix <- realdata_list[, -c(1, 2)]
  
  # Convert to matrix
  individual_matrix <- as.matrix(individual_matrix)
  
  # Fit a simple VAR model (p=1)
  fit <- MTS::VAR(individual_matrix, p = 1)
  
  # Extract AR coefficient matrix (Phi) and error covariance matrix (Sigma)
  phi <- fit$Phi
  sigma <- fit$Sigma
  
  # Return Phi and Sigma for this individual
  return(list(phi = phi, sigma = sigma))
}

# Apply the function to each individual's data and store results
VAR_fit_results <- lapply(realdata_list, fit_VAR_individual)

# Simulate data using VARMA model only 
simulate_individual <- function(var_fit, nTime = 250, id, seed = NULL) {
  # Set the seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Step 1: Simulate temporal effects using VARMA model
  var_sim_data <- tryCatch({
    VARMAsim(n = nTime, arlags = c(1), phi = var_fit$phi, sigma = var_fit$sigma)$series
  }, error = function(e) {
    stop(paste("Error simulating VAR model for temporal effects:", e$message))
  })
  
  # Convert to data frame and add ID and time columns
  simulated_data_with_id_time <- as.data.frame(var_sim_data)
  simulated_data_with_id_time$ID <- id  # Add ID column
  simulated_data_with_id_time$time <- 1:nTime  # Add time column
  
  return(simulated_data_with_id_time)
}

# Function to simulate new data using pre-fitted Phi and Sigma
simulate_individual <- function(phi, sigma, n = 70) {
  # Simulate new data using the provided Phi and Sigma matrices
  sim_data <- VARMAsim(n = n, arlags = c(1), malags = NULL, 
                       phi = phi, sigma = sigma, skip = 200)
  
  simulated_series <- sim_data$series
  
  # Return the simulated series
  return(simulated_series)
}

# Simulate multiple datasets (e.g., 100 datasets) for each individual
simulate_multiple_datasets <- function(VAR_fit_results, nTime = 70, nSim = 10) {
  # Create an empty list to store all simulations
  all_simulated_data <- list()
  
  # Loop through each simulation iteration
  for (sim_num in 1:nSim) {
    # Set a unique seed for the entire dataset for each simulation
    set.seed(12345 + sim_num)
    
    # Simulate data for each individual in the current iteration
    simulated_data_list <- lapply(1:102, function(i) {
      simulate_individual(phi = VAR_fit_results[[i]]$phi, 
                          sigma = VAR_fit_results[[i]]$sigma, 
                          n = nTime)
    })
    
    # Combine into one large data frame for further analysis
    combined_simulated_data <- do.call(rbind, lapply(1:length(simulated_data_list), function(i) {
      individual_df <- as.data.frame(simulated_data_list[[i]])
      
      # Ensure appropriate column names
      colnames(individual_df) <- paste0("V", 1:ncol(individual_df))
      
      # Add ID and time columns
      individual_df$ID <- i
      individual_df$time <- 1:nrow(individual_df)
      
      return(individual_df)
    }))
    
    # Store each simulation result in the list with the simulation number
    all_simulated_data[[sim_num]] <- combined_simulated_data
  }
  
  return(all_simulated_data)
}

# Example usage:
# Assuming 'combined_matrices' contains fitted Phi and Sigma for 102 individuals
nSimulations <- 100  # Number of simulations to run
simulated_data_list <- simulate_multiple_datasets(VAR_fit_results, nTime = 70, nSim = nSimulations)

# Function to rename columns for all datasets in the list
rename_columns_in_simulated_data <- function(simulated_data_list, new_names) {
  # Loop through each dataset in the list
  renamed_simulated_data_list <- lapply(simulated_data_list, function(dataset) {
    # Rename columns of the current dataset
    dataset <- dataset %>% 
      rename_at(vars(names(dataset)), ~ new_names)
    
    return(dataset)
  })
  
  return(renamed_simulated_data_list)
}

# Example usage:
# Assuming you want to rename columns V1, V2, V3, etc. to more meaningful names
new_column_names <- c("depressedmood_state", "loneliness_state_pmc", "socintsatisfaction_state_pmc", 
                      "responsiveness_state_pmc", "selfdisclosure_state_pmc","otherdisclosure_state_pmc","ID", "time")

# Apply the renaming function to your simulated_data_list100
simulated_data_list <- rename_columns_in_simulated_data(simulated_data_list, new_column_names)

#####################################################################
# Calculate simulated contemporaneous correlations using sigma matrix 
library(dplyr)

# Extract the Sigma matrices of each individual 
sigma_matrices <- lapply(VAR_fit_results, function(individual_data) {
  individual_data$sigma  # Adjust if the element name for sigma is different
})

# Function to calculate correlation between predictors and outcome for each individual based on the sigma matrix used during simulation
validate_correlations_range <- function(sigma_matrices) {
  # Initialize an empty list to store data for each individual
  correlation_data <- list()
  
  for (id in names(sigma_matrices)) {
    sigma_matrix <- sigma_matrices[[id]]
    cor_matrix <- cov2cor(sigma_matrix)
    cor_values <- cor_matrix[1, 2:6]  # Correlations of variables 2-6 with variable 1
    
    # Store the correlations along with the ID in a list
    correlation_data[[id]] <- c(ID = id, cor_values)
  }
  
  # Combine all individual results into a data frame
  correlation_df <- do.call(rbind, correlation_data)
  correlation_df <- as.data.frame(correlation_df, stringsAsFactors = FALSE)
  
  # Convert the correlation columns to numeric
  correlation_df[-1] <- lapply(correlation_df[-1], as.numeric)
  
  # Rename the columns for clarity
  colnames(correlation_df) <- c("ID", paste0("Correlation_with_V", 2:6))
  
  return(correlation_df)
}

# Example usage
sigmacorrelations <- validate_correlations_range(sigma_matrices)

# Rename columns
sigmacorrelations <- sigmacorrelations %>%
  rename(loneliness = Correlation_with_V2,
         socintsatisfaction = Correlation_with_V3,
         responsiveness = Correlation_with_V4,
         selfdisclosure = Correlation_with_V5,
         otherdisclosure = Correlation_with_V6)

# Sigma correlations are the known contemporaneous relationships between predictors and outcome
