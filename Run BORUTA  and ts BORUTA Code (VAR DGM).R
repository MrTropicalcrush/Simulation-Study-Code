# =========================================================
# Convert each rep dataframe into a list of individual dataframes by ID
# =========================================================

convert_to_individual_list <- function(simulated_data_list, id_col = "ID") {
  
  if (!is.list(simulated_data_list)) {
    stop("simulated_data_list must be a list.")
  }
  
  individual_simulation_list <- lapply(simulated_data_list, function(df) {
    
    if (!is.data.frame(df)) {
      stop("Each element of simulated_data_list must be a data.frame.")
    }
    
    if (!(id_col %in% names(df))) {
      stop(paste0("Each dataframe must contain an '", id_col, "' column."))
    }
    
    split_list <- split(df, df[[id_col]])
    split_list <- split_list[order(as.numeric(names(split_list)))]
    
    return(split_list)
  })
  
  return(individual_simulation_list)
}

# Apply to AR=1 simulated data
BORUTA_AR1 <- convert_to_individual_list(simulated_data_list_AR1, id_col = "ID")

# Check first rep
names(BORUTA_AR1[[1]])
length(BORUTA_AR1[[1]])

#######################################
#### BORUTA cluster Parallel Cores ####
#######################################

library(Boruta)
library(parallel)

apply_boruta_parallel <- function(individual_simulation_list) {
  
  num_cores <- max(1, detectCores() - 6)
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, {
    library(Boruta)
    NULL
  })
  
  clusterExport(cl, varlist = c("individual_simulation_list"), envir = environment())
  
  BORUTA_results_simulations <- parLapply(cl, seq_along(individual_simulation_list), function(sim_index) {
    
    current_simulation <- individual_simulation_list[[sim_index]]
    
    BORUTAsim <- vector("list", length(current_simulation))
    names(BORUTAsim) <- names(current_simulation)
    
    for (i in seq_along(current_simulation)) {
      
      individual_data <- current_simulation[[i]]
      
      BORUTAsim[[i]] <- tryCatch(
        {
          Boruta(
            depressedmood_state ~ loneliness_state_pmc +
              socintsatisfaction_state_pmc +
              responsiveness_state_pmc +
              selfdisclosure_state_pmc +
              otherdisclosure_state_pmc,
            data = individual_data,
            maxRuns = 500,
            doTrace = 0
          )
        },
        error = function(e) {
          list(error = e$message)
        }
      )
    }
    
    BORUTAsim
  })
  
  names(BORUTA_results_simulations) <- names(individual_simulation_list)
  
  stopCluster(cl)
  
  return(BORUTA_results_simulations)
}

# Run BORUTA on AR=1 data
BORUTA_results_AR1 <- apply_boruta_parallel(BORUTA_AR1)

#############################################
#### Parallel residualisation for AR=1 ######
#############################################

library(parallel)
library(forecast)

num_cores <- max(1, detectCores() - 4)
cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(forecast)
  NULL
})

clusterExport(cl, varlist = c("simulated_data_list_AR1"), envir = environment())

residual_AR1 <- parLapply(cl, seq_along(simulated_data_list_AR1), function(dataset_index) {
  
  dataset <- simulated_data_list_AR1[[dataset_index]]
  residuals_list <- list()
  
  id_col <- if ("ID" %in% names(dataset)) {
    "ID"
  } else if ("id" %in% names(dataset)) {
    "id"
  } else {
    stop(paste("Rep", dataset_index, "has no ID/id column"))
  }
  
  time_col <- if ("time" %in% names(dataset)) {
    "time"
  } else if ("Time" %in% names(dataset)) {
    "Time"
  } else {
    NULL
  }
  
  vars_to_model <- c(
    "depressedmood_state",
    "loneliness_state_pmc",
    "socintsatisfaction_state_pmc",
    "responsiveness_state_pmc",
    "selfdisclosure_state_pmc",
    "otherdisclosure_state_pmc"
  )
  
  vars_to_model <- vars_to_model[vars_to_model %in% names(dataset)]
  
  for (id in unique(dataset[[id_col]])) {
    
    individual_data <- dataset[dataset[[id_col]] == id, ]
    
    if (nrow(individual_data) == 0) next
    
    individual_residuals <- data.frame(ID = rep(id, nrow(individual_data)))
    
    if (!is.null(time_col)) {
      individual_residuals[[time_col]] <- individual_data[[time_col]]
    }
    
    for (var in vars_to_model) {
      
      series <- individual_data[[var]]
      
      if (length(unique(na.omit(series))) <= 1) {
        individual_residuals[[var]] <- rep(NA_real_, nrow(individual_data))
        next
      }
      
      fit <- tryCatch(
        auto.arima(series),
        error = function(e) NULL
      )
      
      if (!is.null(fit)) {
        resids <- as.numeric(residuals(fit))
        
        if (length(resids) < nrow(individual_data)) {
          resids <- c(resids, rep(NA_real_, nrow(individual_data) - length(resids)))
        }
        if (length(resids) > nrow(individual_data)) {
          resids <- resids[1:nrow(individual_data)]
        }
        
        individual_residuals[[var]] <- resids
      } else {
        individual_residuals[[var]] <- rep(NA_real_, nrow(individual_data))
      }
    }
    
    residuals_list[[as.character(id)]] <- individual_residuals
  }
  
  residuals_list <- residuals_list[order(as.numeric(names(residuals_list)))]
  
  return(residuals_list)
})

names(residual_AR1) <- names(simulated_data_list_AR1)

stopCluster(cl)

#######################################
#### Running tsBoruta on AR=1 data ####
#######################################

library(Boruta)
library(parallel)

apply_tsboruta_parallel <- function(residual_data_list) {
  
  num_cores <- max(1, detectCores() - 4)
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, {
    library(Boruta)
    NULL
  })
  
  clusterExport(cl, varlist = c("residual_data_list"), envir = environment())
  
  tsBORUTA_results_simulations <- parLapply(cl, seq_along(residual_data_list), function(sim_index) {
    
    current_simulation <- residual_data_list[[sim_index]]
    BORUTAsim <- vector("list", length(current_simulation))
    names(BORUTAsim) <- names(current_simulation)
    
    for (i in seq_along(current_simulation)) {
      
      individual_data <- current_simulation[[i]]
      
      required_vars <- c(
        "depressedmood_state",
        "loneliness_state_pmc",
        "socintsatisfaction_state_pmc",
        "responsiveness_state_pmc",
        "selfdisclosure_state_pmc",
        "otherdisclosure_state_pmc"
      )
      
      if (is.data.frame(individual_data) && all(required_vars %in% colnames(individual_data))) {
        
        individual_data <- na.omit(individual_data)
        
        if (nrow(individual_data) > 10) {
          BORUTAsim[[i]] <- tryCatch(
            {
              Boruta(
                depressedmood_state ~ loneliness_state_pmc +
                  socintsatisfaction_state_pmc +
                  responsiveness_state_pmc +
                  selfdisclosure_state_pmc +
                  otherdisclosure_state_pmc,
                data = individual_data,
                maxRuns = 500,
                doTrace = 0
              )
            },
            error = function(e) {
              list(error = e$message)
            }
          )
        } else {
          BORUTAsim[[i]] <- list(error = "Too few rows after na.omit")
        }
        
      } else {
        BORUTAsim[[i]] <- list(error = "Missing required variables or invalid data frame")
      }
    }
    
    return(BORUTAsim)
  })
  
  names(tsBORUTA_results_simulations) <- names(residual_data_list)
  
  stopCluster(cl)
  
  return(tsBORUTA_results_simulations)
}

# Run tsBoruta on residualised AR=1 data
tsBoruta_results_AR1 <- apply_tsboruta_parallel(residual_AR1)

