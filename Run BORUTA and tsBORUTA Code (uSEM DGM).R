# BORUTA
# =========================================================
# Convert each rep dataframe into a list of individual dataframes by ID
# =========================================================

convert_to_individual_list <- function(simulated_data_list) {
  
  # Check input is a list
  if (!is.list(simulated_data_list)) {
    stop("simulated_data_list must be a list.")
  }
  
  # Loop through each rep
  individual_simulation_list <- lapply(simulated_data_list, function(df) {
    
    # Check each rep is a dataframe
    if (!is.data.frame(df)) {
      stop("Each element of simulated_data_list must be a data.frame.")
    }
    
    # Check ID column exists
    if (!("id" %in% names(df))) {
      stop("Each dataframe must contain an 'ID' column.")
    }
    
    # Split by ID
    split_list <- split(df, df$id)
    
    # Reorder individuals numerically by ID
    split_list <- split_list[order(as.numeric(names(split_list)))]
    
    return(split_list)
  })
  
  return(individual_simulation_list)
}

# Apply to your new uSEM DGM simulated data
BORUTA_uSEM_DGM <- convert_to_individual_list(sim_data_list)

# Check first rep
names(BORUTA_uSEM_DGM[[1]])
length(BORUTA_uSEM_DGM[[1]])

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
            V1 ~ V2 + V3 + V4 + V5 + V6,   # 👈 UPDATED FOR uSEM DGM
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

# Run it
BORUTA_results_uSEM_DGM <- apply_boruta_parallel(BORUTA_uSEM_DGM)

#################################################################################
########## tsBoruta #############################################################
#################################################################################

# Parallel residualisation for uSEM DGM (robust ID handling)

library(parallel)
library(forecast)

num_cores <- max(1, detectCores() - 4)
cl <- makeCluster(num_cores)

clusterEvalQ(cl, {
  library(forecast)
  NULL
})

clusterExport(cl, varlist = c("sim_data_list"), envir = environment())

residual_uSEM_DGM <- parLapply(cl, seq_along(sim_data_list), function(dataset_index) {
  
  dataset <- sim_data_list[[dataset_index]]
  residuals_list <- list()
  
  # 🔥 Detect ID column automatically
  id_col <- if ("id" %in% names(dataset)) {
    "id"
  } else if ("ID" %in% names(dataset)) {
    "ID"
  } else {
    stop(paste("Rep", dataset_index, "has no id/ID column"))
  }
  
  # Detect time column (optional)
  time_col <- if ("time" %in% names(dataset)) {
    "time"
  } else if ("Time" %in% names(dataset)) {
    "Time"
  } else {
    NULL
  }
  
  # Only use V variables
  vars_to_model <- names(dataset)[grepl("^V[0-9]+$", names(dataset))]
  
  # Loop through individuals
  for (id in unique(dataset[[id_col]])) {
    
    individual_data <- dataset[dataset[[id_col]] == id, ]
    
    if (nrow(individual_data) == 0) next
    
    # Create residual dataframe
    individual_residuals <- data.frame(id = rep(id, nrow(individual_data)))
    
    if (!is.null(time_col)) {
      individual_residuals[[time_col]] <- individual_data[[time_col]]
    }
    
    # Residualise each V variable
    for (var in vars_to_model) {
      
      series <- individual_data[[var]]
      
      # Handle constant series
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
        
        # Pad or trim
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
  
  # Order individuals numerically
  residuals_list <- residuals_list[order(as.numeric(names(residuals_list)))]
  
  return(residuals_list)
})

names(residual_uSEM_DGM) <- names(sim_data_list)

stopCluster(cl)

#### Running tsBoruta on residual data ####

library(Boruta)
library(parallel)

apply_tsboruta_parallel <- function(residual_uSEM_DGM) {
  
  # Set up number of cores
  num_cores <- max(1, detectCores() - 4)
  cl <- makeCluster(num_cores)
  
  # Load Boruta on workers
  clusterEvalQ(cl, {
    library(Boruta)
    NULL
  })
  
  # Export residual data to workers
  clusterExport(cl, varlist = c("residual_uSEM_DGM"), envir = environment())
  
  # Run each rep in parallel
  tsBORUTA_results_simulations <- parLapply(cl, seq_along(residual_uSEM_DGM), function(sim_index) {
    
    current_simulation <- residual_uSEM_DGM[[sim_index]]
    num_individuals <- length(current_simulation)
    
    # Store results for this rep
    BORUTAsim <- vector("list", num_individuals)
    names(BORUTAsim) <- names(current_simulation)
    
    # Loop through each individual
    for (i in seq_along(current_simulation)) {
      
      individual_data <- current_simulation[[i]]
      
      # Convert ts object if needed
      if (inherits(individual_data, "ts")) {
        individual_data <- as.data.frame(individual_data)
        colnames(individual_data) <- paste0("V", seq_len(ncol(individual_data)))
      }
      
      # Check required uSEM columns exist
      required_vars <- c("V1", "V2", "V3", "V4", "V5", "V6")
      
      if (is.data.frame(individual_data) && all(required_vars %in% colnames(individual_data))) {
        
        # Remove rows with missing data
        individual_data <- na.omit(individual_data)
        
        # Only run if enough rows remain
        if (nrow(individual_data) > 10) {
          BORUTAsim[[i]] <- tryCatch(
            {
              Boruta(
                V1 ~ V2 + V3 + V4 + V5 + V6,
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
        BORUTAsim[[i]] <- list(error = "Missing required V columns or invalid data frame")
      }
    }
    
    return(BORUTAsim)
  })
  
  # Keep rep names
  names(tsBORUTA_results_simulations) <- names(residual_uSEM_DGM)
  
  stopCluster(cl)
  
  return(tsBORUTA_results_simulations)
}

# Run tsBoruta on residualised uSEM DGM data
tsBoruta_results_uSEM_DGM <- apply_tsboruta_parallel(residual_uSEM_DGM)

