#########################################
#### MV-iARIMAX for AR=1 Parallel ####
#########################################
run_simulated_arimax_parallel_AR1 <- function(simulated_data_list_AR1, impute_vars, id_col = "ID") {
  
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
  
  simulation_results <- parLapply(cl, simulated_data_list_AR1, function(sim_data) {
    
    # Split by individual
    individual_list <- split(sim_data, sim_data[[id_col]])
    individual_list <- individual_list[order(as.numeric(names(individual_list)))]
    
    # Run MV-ARIMAX for each individual
    arimax_results <- lapply(individual_list, function(ind_data) {
      
      # Standardize using i_standarbot_300
      standardized_data <- i_standarbot_300(ind_data, impute_vars, id_col, explanation = TRUE)
      
      # Convert dependent variable to time series
      standardized_data$depressedmood_state_PSD <- ts(
        standardized_data$depressedmood_state_PSD,
        frequency = 1
      )
      
      # Predictor matrix
      xreg <- standardized_data[, c(
        "loneliness_state_pmc_PSD",
        "socintsatisfaction_state_pmc_PSD",
        "responsiveness_state_pmc_PSD",
        "selfdisclosure_state_pmc_PSD",
        "otherdisclosure_state_pmc_PSD"
      )]
      xreg <- as.matrix(xreg)
      
      # Run ARIMAX
      model <- auto.arima(standardized_data$depressedmood_state_PSD, xreg = xreg)
      
      return(model)
    })
    
    return(arimax_results)
  })
  
  stopCluster(cl)
  
  return(simulation_results)
}

# Variables to standardize
imputeThese_AR1 <- c(
  "loneliness_state_pmc",
  "socintsatisfaction_state_pmc",
  "responsiveness_state_pmc",
  "depressedmood_state",
  "selfdisclosure_state_pmc",
  "otherdisclosure_state_pmc"
)

# Run function
MV_ARIMAX_AR1 <- run_simulated_arimax_parallel_AR1(
  simulated_data_list_AR1,
  imputeThese_AR1,
  id_col = "ID"
)

############################################
############ i-ARIMAX for AR=1 #############
############################################
library(dplyr)
library(idionomics)
library(MTS)
library(parallel)

run_iarimax_on_AR1_parallel <- function(simulated_data_list_AR1,
                                        IV = c("loneliness_state_pmc",
                                               "socintsatisfaction_state_pmc",
                                               "responsiveness_state_pmc",
                                               "selfdisclosure_state_pmc",
                                               "otherdisclosure_state_pmc"),
                                        id_col = "ID",
                                        time_col = "time",
                                        outcome = "depressedmood_state") {
  
  num_cores <- max(1, detectCores() - 4)
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, {
    library(dplyr)
    library(idionomics)
    library(MTS)
    NULL
  })
  
  clusterExport(
    cl,
    c("simulated_data_list_AR1", "IV", "id_col", "time_col", "outcome",
      "i_standarbot_300", "IARIMAXoid_Pro"),
    envir = environment()
  )
  
  all_iarimax_results <- parLapply(cl, seq_along(simulated_data_list_AR1), function(i) {
    
    current_data <- as.data.frame(simulated_data_list_AR1[[i]])
    
    needed_cols <- c(id_col, time_col, outcome, IV)
    missing_cols <- setdiff(needed_cols, names(current_data))
    
    if (length(missing_cols) > 0) {
      return(list(
        status = "error",
        reason = paste("Missing columns:", paste(missing_cols, collapse = ", "))
      ))
    }
    
    current_data <- current_data[, needed_cols, drop = FALSE]
    
    # Rename to match IARIMAXoid_Pro expectations
    names(current_data)[names(current_data) == id_col] <- "ID"
    names(current_data)[names(current_data) == time_col] <- "Time"
    names(current_data)[names(current_data) == outcome] <- "Y"
    
    imputeThese <- c("Y", IV)
    zDatasim <- i_standarbot_300(current_data, imputeThese, "ID", explanation = TRUE)
    
    iarimax_results <- list()
    
    for (j in seq_along(IV)) {
      current_iv <- IV[[j]]
      model_name <- paste0("Sim_", i, "_", current_iv)
      
      zData_temp <- zDatasim
      zData_temp$X_PSD <- zData_temp[[paste0(current_iv, "_PSD")]]
      
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
  
  names(all_iarimax_results) <- names(simulated_data_list_AR1)
  
  return(all_iarimax_results)
}

# Define IVs
IV_AR1 <- c(
  "loneliness_state_pmc",
  "socintsatisfaction_state_pmc",
  "responsiveness_state_pmc",
  "selfdisclosure_state_pmc",
  "otherdisclosure_state_pmc"
)

# Run i-ARIMAX on AR=1 data
iarimax_results_AR1 <- run_iarimax_on_AR1_parallel(
  simulated_data_list_AR1 = simulated_data_list_AR1,
  IV = IV_AR1,
  id_col = "ID",
  time_col = "time",
  outcome = "depressedmood_state"
)