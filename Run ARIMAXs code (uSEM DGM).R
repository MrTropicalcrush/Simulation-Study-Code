#########################################
#### MV-iARIMAX for uSEM DGM Parallel ####
#########################################
run_simulated_arimax_parallel_uSEM <- function(sim_data_list, impute_vars, id_col = "id") {
  
  library(parallel)
  library(forecast)
  library(idionomics)
  
  num_cores <- detectCores() - 4
  cl <- makeCluster(num_cores)
  
  clusterExport(cl, varlist = c("i_standarbot_300", "impute_vars", "id_col"), envir = environment())
  
  clusterEvalQ(cl, {
    library(forecast)
    library(idionomics)
  })
  
  simulation_results <- parLapply(cl, sim_data_list, function(sim_data) {
    
    # Split by individual
    individual_list <- split(sim_data, sim_data[[id_col]])
    individual_list <- individual_list[order(as.numeric(names(individual_list)))]
    
    # Run MV-ARIMAX for each individual
    arimax_results <- lapply(individual_list, function(ind_data) {
      
      # Standardize using i_standarbot_300
      standardized_data <- i_standarbot_300(ind_data, impute_vars, id_col, explanation = TRUE)
      
      # Convert dependent variable to time series
      standardized_data$V1_PSD <- ts(standardized_data$V1_PSD, frequency = 1)
      
      # Predictor matrix
      xreg <- standardized_data[, c("V2_PSD", "V3_PSD", "V4_PSD", "V5_PSD", "V6_PSD")]
      xreg <- as.matrix(xreg)
      
      # Run ARIMAX
      model <- auto.arima(standardized_data$V1_PSD, xreg = xreg)
      
      return(model)
    })
    
    return(arimax_results)
  })
  
  stopCluster(cl)
  
  return(simulation_results)
}

# Variables to standardize
imputeThese_uSEM <- c("V1", "V2", "V3", "V4", "V5", "V6")

# Run function
MV_ARIMAX_uSEM_DGM <- run_simulated_arimax_parallel_uSEM(
  sim_data_list,
  imputeThese_uSEM,
  id_col = "id"
)

############################################
############ i-ARIMAX ######################
############################################
library(dplyr)
library(idionomics)
library(MTS)
library(parallel)

# Function to run i-ARIMAX on uSEM DGM simulated data in parallel
run_iarimax_on_uSEM_parallel <- function(sim_data_list,
                                         IV = c("V2", "V3", "V4", "V5", "V6"),
                                         id_col = "id",
                                         time_col = "time",
                                         outcome = "V1") {
  
  # Set up number of cores
  num_cores <- max(1, detectCores() - 4)
  cl <- makeCluster(num_cores)
  
  # Load required libraries on workers
  clusterEvalQ(cl, {
    library(dplyr)
    library(idionomics)
    library(MTS)
    NULL
  })
  
  # Export objects/functions to workers
  clusterExport(
    cl,
    c("sim_data_list", "IV", "id_col", "time_col", "outcome",
      "i_standarbot_300", "IARIMAXoid_Pro"),
    envir = environment()
  )
  
  # Run simulations in parallel
  all_iarimax_results <- parLapply(cl, seq_along(sim_data_list), function(i) {
    
    current_data <- as.data.frame(sim_data_list[[i]])
    
    # Basic checks
    needed_cols <- c(id_col, time_col, outcome, IV)
    missing_cols <- setdiff(needed_cols, names(current_data))
    
    if (length(missing_cols) > 0) {
      return(list(
        status = "error",
        reason = paste("Missing columns:", paste(missing_cols, collapse = ", "))
      ))
    }
    
    # Keep only relevant columns
    current_data <- current_data[, needed_cols, drop = FALSE]
    
    # Rename to match IARIMAXoid_Pro workflow
    names(current_data)[names(current_data) == id_col] <- "ID"
    names(current_data)[names(current_data) == time_col] <- "Time"
    names(current_data)[names(current_data) == outcome] <- "Y"
    
    # Standardization
    imputeThese <- c("Y", IV)
    zDatasim <- i_standarbot_300(current_data, imputeThese, "ID", explanation = TRUE)
    
    # Store i-ARIMAX results for each predictor
    iarimax_results <- list()
    
    for (j in seq_along(IV)) {
      current_iv <- IV[[j]]
      model_name <- paste0("Sim_", i, "_", current_iv)
      
      # Create temporary standardized X variable
      zData_temp <- zDatasim
      zData_temp$X_PSD <- zData_temp[[paste0(current_iv, "_PSD")]]
      
      # Run one-predictor i-ARIMAX using standardized variables
      iarimax_results[[model_name]] <- tryCatch(
        {
          IARIMAXoid_Pro(
            zData_temp,
            x_series = "X_PSD",
            y_series = "Y_PSD",
            id_var = "ID",
            hlm_compare = FALSE,
            timevar = "Time",
            metaanalysis = TRUE
          )
        },
        error = function(e) {
          list(
            status = "error",
            reason = e$message,
            predictor = current_iv
          )
        }
      )
    }
    
    return(iarimax_results)
  })
  
  stopCluster(cl)
  
  names(all_iarimax_results) <- names(sim_data_list)
  
  return(all_iarimax_results)
}

IV_uSEM <- c("V2", "V3", "V4", "V5", "V6")

iarimax_results_uSEM_DGM <- run_iarimax_on_uSEM_parallel(
  sim_data_list = sim_data_list,
  IV = IV_uSEM,
  id_col = "id",
  time_col = "time",
  outcome = "V1"
)