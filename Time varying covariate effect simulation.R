library(dplyr)

create_VAR_time_varying_loneliness <- function(data,
                                               id_var = "ID",
                                               time_var = "time",
                                               outcome_var = "depressedmood_state",
                                               predictor_var = "loneliness_state_pmc",
                                               early_effect = 0.00,
                                               late_effect = 0.40) {
  
  data_tve <- data %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[time_var]], .by_group = TRUE) %>%
    mutate(
      row_num = row_number(),
      n_obs = n(),
      split_point = floor(n_obs / 2),
      
      loneliness_timevarying_beta = ifelse(
        row_num <= split_point,
        early_effect,
        late_effect
      ),
      
      "{outcome_var}" := .data[[outcome_var]] +
        loneliness_timevarying_beta * .data[[predictor_var]]
    ) %>%
    ungroup() %>%
    select(-row_num, -n_obs, -split_point)
  
  return(data_tve)
}

VAR_timevarying_loneliness_list <- lapply(seq_along(sim_data_list_VAR_baseline), function(i) {
  create_VAR_time_varying_loneliness(
    data = sim_data_list_VAR_baseline[[i]],
    id_var = "ID",
    time_var = "time",
    outcome_var = "depressedmood_state",
    predictor_var = "loneliness_state_pmc",
    early_effect = 0.00,
    late_effect = 0.40
  )
})

VAR_timevarying_loneliness_list[[1]] %>%
  group_by(ID) %>%
  summarise(
    n_obs = n(),
    early_beta = first(loneliness_timevarying_beta),
    late_beta = last(loneliness_timevarying_beta),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(ID))

#########################################################
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
BORUTA_timevarying <- convert_to_individual_list(VAR_timevarying_loneliness_list, id_col = "ID")

# Check first rep
names(BORUTA_AR1[[1]])
length(BORUTA_AR1[[1]])

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

clusterExport(cl, varlist = c("VAR_timevarying_loneliness_list"), envir = environment())

residual_AR1 <- parLapply(cl, seq_along(VAR_timevarying_loneliness_list), function(dataset_index) {
  
  dataset <- VAR_timevarying_loneliness_list[[dataset_index]]
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

names(residual_AR1) <- names(VAR_timevarying_loneliness_list)

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
tsBoruta_results_timevarying <- apply_tsboruta_parallel(residual_timevarying)


