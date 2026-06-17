#########################################################
#### Calculting lagged peformance metrics (uSEM DGM)####
########################################################

# i-ARIMAX #
library(dplyr)
library(tidyr)
library(stringr)

# Function to calculate t-statistics for one variable in one simulation
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

calculate_t_statistics_wide <- function(iarimax_results_list, variable_names) {
  all_results <- list()
  
  for (sim in seq_along(iarimax_results_list)) {
    sim_results <- iarimax_results_list[[sim]]
    
    sim_t_stats <- lapply(seq_along(variable_names), function(j) {
      object_name <- paste0("Sim_", sim, "_", variable_names[j])
      variable_result <- sim_results[[object_name]]
      
      if (!is.null(variable_result)) {
        t_stat_data <- calculate_t_statistics(variable_result, variable_names[j])
        return(t_stat_data)
      } else {
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
    arrange(simulation, ID)
  
  return(final_wide_t_stats)
}

# Lagged predictor names
variable_names_lag <- c("V2_lag", "V3_lag", "V4_lag", "V5_lag", "V6_lag")

# Extract lagged t-stats from results
iarimax_tstats_uSEM_lag_wide <- calculate_t_statistics_wide(
  iarimax_results_list = iarimax_results_uSEM_lagged,
  variable_names = variable_names_lag
)

iarimax_tstats_uSEM_lag_wide <- iarimax_tstats_uSEM_lag_wide %>%
  mutate(ID = as.numeric(ID)) %>%
  arrange(simulation, ID)

# Convert lagged t-stats to binary estimates
generate_iarimax_est_lag_df <- function(t_stats_wide_df, t_threshold = 1.96) {
  est_df <- t_stats_wide_df %>%
    pivot_longer(
      cols = c("V2_lag", "V3_lag", "V4_lag", "V5_lag", "V6_lag"),
      names_to = "predictor",
      values_to = "t_stat"
    ) %>%
    mutate(
      predictor = str_remove(predictor, "_lag")
    ) %>%
    transmute(
      rep = simulation,
      id = ID,
      predictor = predictor,
      est_lag = ifelse(!is.na(t_stat) & abs(t_stat) > t_threshold, 1L, 0L)
    ) %>%
    arrange(rep, id, predictor)
  
  return(est_df)
}

iarimax_est_uSEMDGM_lag <- generate_iarimax_est_lag_df(
  t_stats_wide_df = iarimax_tstats_uSEM_lag_wide,
  t_threshold = 1.96
)

# Score lagged metrics
score_all_reps_lag <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
  rep_scores <- lapply(rep_ids, function(r) {
    truth_rep <- subset(truth_df, rep == r)
    est_rep   <- subset(est_df, rep == r)
    
    merged <- merge(
      truth_rep,
      est_rep,
      by = c("rep", "id", "predictor"),
      all.x = TRUE,
      sort = FALSE
    )
    
    merged$est_lag[is.na(merged$est_lag)] <- 0L
    
    metrics <- calc_binary_metrics(merged$true_lag, merged$est_lag)
    data.frame(rep = r, metrics, row.names = NULL)
  })
  
  rep_scores <- do.call(rbind, rep_scores)
  rownames(rep_scores) <- NULL
  
  avg_scores <- aggregate(
    cbind(TP, TN, FP, FN, Sensitivity, Specificity, Precision, F1, MCC) ~ 1,
    data = rep_scores,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  avg_scores <- avg_scores[, -1, drop = FALSE]
  
  list(
    rep_scores = rep_scores,
    avg_scores = avg_scores
  )
}

calc_binary_metrics <- function(truth, estimate) {
  stopifnot(length(truth) == length(estimate))
  
  tp <- as.numeric(sum(truth == 1 & estimate == 1, na.rm = TRUE))
  tn <- as.numeric(sum(truth == 0 & estimate == 0, na.rm = TRUE))
  fp <- as.numeric(sum(truth == 0 & estimate == 1, na.rm = TRUE))
  fn <- as.numeric(sum(truth == 1 & estimate == 0, na.rm = TRUE))
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  
  f1 <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else {
    NA_real_
  }
  
  denom <- as.numeric(sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
  mcc <- if (is.finite(denom) && denom > 0) {
    ((tp * tn) - (fp * fn)) / denom
  } else {
    NA_real_
  }
  
  data.frame(
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    MCC = mcc
  )
}

iarimax_scores_uSEM_lag <- score_all_reps_lag(
  truth_df = truth_df,
  est_df = iarimax_est_uSEMDGM_lag
)

iarimax_scores_uSEM_lag$avg_scores

calculate_metric_cis <- function(rep_scores_df,
                                 metric_cols = c("Sensitivity", "Specificity", "Precision", "F1", "MCC"),
                                 conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  
  ci_results <- rep_scores_df %>%
    select(all_of(metric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      n = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(n),
      t_crit = qt(1 - alpha/2, df = n - 1),
      Lower_CI = Mean - t_crit * SE,
      Upper_CI = Mean + t_crit * SE,
      .groups = "drop"
    )
  
  return(ci_results)
}

arimax_ci_results <- calculate_metric_cis(iarimax_scores_uSEM_lag$rep_scores)

print(arimax_ci_results)

###########################################################################################

# MV i-ARIMAX #
library(forecast)
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

mv_iarimax_uSEM_extracted <- extract_arimax_results_from_simulations(MV_ARIMAX_uSEM_lagged)

generate_mv_iarimax_est_lag_df <- function(extracted_results_list,
                                           predictor_map = c(
                                             "V2_lag_PSD" = "V2",
                                             "V3_lag_PSD" = "V3",
                                             "V4_lag_PSD" = "V4",
                                             "V5_lag_PSD" = "V5",
                                             "V6_lag_PSD" = "V6"
                                           ),
                                           t_threshold = 1.96) {
  est_list <- list()
  
  raw_predictors <- names(predictor_map)
  
  for (rep in seq_along(extracted_results_list)) {
    rep_results <- extracted_results_list[[rep]]
    person_rows <- list()
    
    for (person_id in names(rep_results)) {
      person_tstats <- rep_results[[person_id]]$t_stats
      
      available_preds <- intersect(raw_predictors, names(person_tstats))
      
      person_df <- data.frame(
        rep = rep,
        id = as.numeric(person_id),
        predictor = unname(predictor_map[available_preds]),
        est_lag = ifelse(abs(person_tstats[available_preds]) > t_threshold, 1L, 0L),
        stringsAsFactors = FALSE
      )
      
      person_rows[[person_id]] <- person_df
    }
    
    est_list[[rep]] <- do.call(rbind, person_rows)
  }
  
  est_df <- do.call(rbind, est_list)
  rownames(est_df) <- NULL
  est_df
}

MViarimax_est_uSEMDGM_lag <- generate_mv_iarimax_est_lag_df(
  extracted_results_list = mv_iarimax_uSEM_extracted,
  predictor_map = c(
    "V2_lag_PSD" = "V2",
    "V3_lag_PSD" = "V3",
    "V4_lag_PSD" = "V4",
    "V5_lag_PSD" = "V5",
    "V6_lag_PSD" = "V6"
  ),
  t_threshold = 1.96
)

score_one_rep_mv_iarimax_lag <- function(truth_df, est_df, rep_id) {
  truth_rep <- subset(truth_df, rep == rep_id)
  est_rep   <- subset(est_df, rep == rep_id)
  
  merged <- merge(
    truth_rep,
    est_rep,
    by = c("rep", "id", "predictor"),
    all.x = TRUE,
    sort = FALSE
  )
  
  if (!"est_lag" %in% names(merged)) {
    merged$est_lag <- 0L
  }
  
  merged$est_lag[is.na(merged$est_lag)] <- 0L
  
  list(
    lag = calc_binary_metrics(merged$true_lag, merged$est_lag),
    merged = merged
  )
}

score_all_reps_mv_iarimax_lag <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
  rep_scores <- lapply(rep_ids, function(r) {
    res <- score_one_rep_mv_iarimax_lag(truth_df, est_df, rep_id = r)
    
    data.frame(
      rep = r,
      target = "lag",
      res$lag,
      row.names = NULL
    )
  })
  
  rep_scores <- do.call(rbind, rep_scores)
  rownames(rep_scores) <- NULL
  
  avg_scores <- aggregate(
    cbind(TP, TN, FP, FN, Sensitivity, Specificity, Precision, F1, MCC) ~ target,
    data = rep_scores,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  list(
    rep_scores = rep_scores,
    avg_scores = avg_scores
  )
}

MViarimax_scores_uSEM_lag <- score_all_reps_mv_iarimax_lag(
  truth_df = truth_df,
  est_df = MViarimax_est_uSEMDGM_lag
)

MViarimax_scores_uSEM_lag$avg_scores

calculate_metric_cis <- function(rep_scores_df,
                                 metric_cols = c("Sensitivity", "Specificity", "Precision", "F1", "MCC"),
                                 conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  
  ci_results <- rep_scores_df %>%
    select(all_of(metric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      n = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(n),
      t_crit = qt(1 - alpha/2, df = n - 1),
      Lower_CI = Mean - t_crit * SE,
      Upper_CI = Mean + t_crit * SE,
      .groups = "drop"
    )
  
  return(ci_results)
}

MVarimax_ci_results <- calculate_metric_cis(MViarimax_scores_uSEM_lag$rep_scores)

print(MVarimax_ci_results)

######################################################################################

# iBoruta/tsBoruta #
generate_boruta_est_lag_df <- function(BORUTA_results_all) {
  est_list <- vector("list", length(BORUTA_results_all))
  
  lag_predictors <- c("V2_lag", "V3_lag", "V4_lag", "V5_lag", "V6_lag")
  truth_predictors <- c("V2", "V3", "V4", "V5", "V6")
  
  for (sim in seq_along(BORUTA_results_all)) {
    sim_rows <- list()
    
    for (i in seq_along(BORUTA_results_all[[sim]])) {
      fd <- BORUTA_results_all[[sim]][[i]]$finalDecision
      
      person_df <- data.frame(
        rep = sim,
        id = i,
        predictor = truth_predictors,
        est_lag = c(
          ifelse(fd["V2_lag"] == "Confirmed", 1L, 0L),
          ifelse(fd["V3_lag"] == "Confirmed", 1L, 0L),
          ifelse(fd["V4_lag"] == "Confirmed", 1L, 0L),
          ifelse(fd["V5_lag"] == "Confirmed", 1L, 0L),
          ifelse(fd["V6_lag"] == "Confirmed", 1L, 0L)
        ),
        stringsAsFactors = FALSE
      )
      
      sim_rows[[i]] <- person_df
    }
    
    est_list[[sim]] <- do.call(rbind, sim_rows)
  }
  
  est_df <- do.call(rbind, est_list)
  rownames(est_df) <- NULL
  est_df
}

score_all_reps_boruta_lag <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
  rep_scores <- lapply(rep_ids, function(r) {
    truth_rep <- subset(truth_df, rep == r)
    est_rep   <- subset(est_df, rep == r)
    
    merged <- merge(
      truth_rep,
      est_rep,
      by = c("rep", "id", "predictor"),
      all.x = TRUE,
      sort = FALSE
    )
    
    merged$est_lag[is.na(merged$est_lag)] <- 0L
    
    metrics <- calc_binary_metrics(merged$true_lag, merged$est_lag)
    data.frame(rep = r, metrics, row.names = NULL)
  })
  
  rep_scores <- do.call(rbind, rep_scores)
  rownames(rep_scores) <- NULL
  
  avg_scores <- aggregate(
    cbind(TP, TN, FP, FN, Sensitivity, Specificity, Precision, F1, MCC) ~ 1,
    data = rep_scores,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  avg_scores <- avg_scores[, -1, drop = FALSE]
  
  list(
    rep_scores = rep_scores,
    avg_scores = avg_scores
  )
}

TSboruta_est_uSEMDGM_lag <- generate_boruta_est_lag_df(
  BORUTA_results_uSEM_lagged
)

TSboruta_uSEM_lag <- score_all_reps_boruta_lag(
  truth_df = truth_df,
  est_df = TSboruta_est_uSEMDGM_lag
)

TSboruta_uSEM_lag$avg_scores

calculate_metric_cis <- function(rep_scores_df,
                                 metric_cols = c("Sensitivity", "Specificity", "Precision", "F1", "MCC"),
                                 conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  
  ci_results <- rep_scores_df %>%
    select(all_of(metric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      n = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(n),
      t_crit = qt(1 - alpha/2, df = n - 1),
      Lower_CI = Mean - t_crit * SE,
      Upper_CI = Mean + t_crit * SE,
      .groups = "drop"
    )
  
  return(ci_results)
}

boruta_ci_results <- calculate_metric_cis(TSboruta_uSEM_lag$rep_scores)

print(boruta_ci_results)

######################################################################################

# GIMME #
library(dplyr)

process_gimme_simulation_folders_lag <- function(directory,
                                                 outcome_var = "V1",
                                                 lag_predictor_vars = c("V2lag", "V3lag", "V4lag", "V5lag", "V6lag")) {
  subfolders <- list.dirs(directory, recursive = FALSE)
  all_simulation_results <- list()
  
  for (subfolder in subfolders) {
    csv_file <- file.path(subfolder, "indivPathEstimates.csv")
    
    if (file.exists(csv_file)) {
      data <- read.csv(csv_file)
      
      relevant_paths <- data %>%
        filter(lhs == outcome_var, rhs %in% lag_predictor_vars) %>%
        select(ID = file, rhs, beta) %>%
        mutate(
          rep = as.numeric(sub("^Sim", "", basename(subfolder))),
          predictor = gsub("lag", "", rhs)
        )
      
      all_simulation_results[[subfolder]] <- relevant_paths
    }
  }
  
  bind_rows(all_simulation_results)
}

all_paths_uSEM_lag <- process_gimme_simulation_folders_lag(
  directory = "C:/Users/WillLi/Documents/mlvar simulation/output_directory_GIMME_uSEM_DGM",
  outcome_var = "V1",
  lag_predictor_vars = c("V2lag", "V3lag", "V4lag", "V5lag", "V6lag")
)

predictor_vars <- c("V2", "V3", "V4", "V5", "V6")

all_IDs <- sort(unique(all_paths_uSEM_lag$ID))
all_reps <- sort(unique(all_paths_uSEM_lag$rep))

full_grid <- expand.grid(
  rep = all_reps,
  id = all_IDs,
  predictor = predictor_vars,
  stringsAsFactors = FALSE
)

GIMME_est_uSEMDGM_lag <- full_grid %>%
  left_join(
    all_paths_uSEM_lag %>%
      transmute(
        rep = rep,
        id = ID,
        predictor = predictor,
        beta = beta
      ),
    by = c("rep", "id", "predictor")
  ) %>%
  mutate(
    est_lag = ifelse(!is.na(beta), 1L, 0L)
  ) %>%
  select(rep, id, predictor, est_lag) %>%
  arrange(rep, id, predictor)

GIMME_est_uSEMDGM_lag <- GIMME_est_uSEMDGM_lag %>%
  mutate(id = as.numeric(id))

truth_df <- truth_df %>%
  mutate(id = as.numeric(id))

GIMME_scores_uSEM_lag <- score_all_reps_lag(
  truth_df = truth_df,
  est_df = GIMME_est_uSEMDGM_lag
)

GIMME_scores_uSEM_lag$avg_scores

calculate_metric_cis <- function(rep_scores_df,
                                 metric_cols = c("Sensitivity", "Specificity", "Precision", "F1", "MCC"),
                                 conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  
  ci_results <- rep_scores_df %>%
    select(all_of(metric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      n = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(n),
      t_crit = qt(1 - alpha/2, df = n - 1),
      Lower_CI = Mean - t_crit * SE,
      Upper_CI = Mean + t_crit * SE,
      .groups = "drop"
    )
  
  return(ci_results)
}

GIMME_ci_results <- calculate_metric_cis(GIMME_scores_uSEM_lag$rep_scores)

print(GIMME_ci_results)

##################################################################################

# indSEM #
library(dplyr)
library(tidyr)

process_indSEM_simulations_lag <- function(directory,
                                           outcome_row = "V1",
                                           lag_predictor_vars = c("V2lag", "V3lag", "V4lag", "V5lag", "V6lag")) {
  
  sim_folders <- list.dirs(directory, recursive = FALSE)
  all_simulation_results <- list()
  
  for (sim_folder in sim_folders) {
    individual_folder <- file.path(sim_folder, "individual")
    
    if (!dir.exists(individual_folder)) {
      message(paste("Folder does not exist:", individual_folder))
      next
    }
    
    message(paste("Processing folder:", individual_folder))
    
    beta_files <- list.files(
      individual_folder,
      pattern = "\\d+BetasStd\\.csv$",
      full.names = TRUE
    )
    
    if (length(beta_files) == 0) {
      message(paste("No beta files found in:", individual_folder))
      next
    }
    
    sim_results <- list()
    
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
      
      if (!(outcome_row %in% rownames(beta_matrix))) {
        message(paste("Row", outcome_row, "not found in file:", beta_file))
        message(paste("Available row names:", paste(rownames(beta_matrix), collapse = ", ")))
        next
      }
      
      available_predictors <- intersect(lag_predictor_vars, colnames(beta_matrix))
      
      if (length(available_predictors) == 0) {
        message(paste("No matching lagged predictors found in file:", beta_file))
        message(paste("Available columns:", paste(colnames(beta_matrix), collapse = ", ")))
        next
      }
      
      outcome_row_values <- beta_matrix[outcome_row, available_predictors, drop = TRUE]
      
      beta_data <- data.frame(
        rep = as.numeric(sub("^Sim", "", basename(sim_folder))),
        id = as.numeric(individual_id),
        predictor_lag = available_predictors,
        predictor = gsub("lag", "", available_predictors),
        beta = as.numeric(outcome_row_values),
        stringsAsFactors = FALSE
      )
      
      sim_results[[basename(beta_file)]] <- beta_data
    }
    
    if (length(sim_results) > 0) {
      all_simulation_results[[basename(sim_folder)]] <- bind_rows(sim_results)
    }
  }
  
  if (length(all_simulation_results) == 0) {
    message("No valid lagged data found in the simulations.")
    return(NULL)
  }
  
  bind_rows(all_simulation_results) %>%
    arrange(rep, id, predictor)
}

parent_directory <- "C:/Users/WillLi/Documents/mlvar simulation/output_directory_GIMME_uSEM_DGM"

all_paths_indSEM_lag <- process_indSEM_simulations_lag(
  directory = parent_directory,
  outcome_row = "V1",
  lag_predictor_vars = c("V2lag", "V3lag", "V4lag", "V5lag", "V6lag")
)

generate_indSEM_est_lag_df <- function(all_paths_indSEM_lag,
                                       predictor_vars = c("V2", "V3", "V4", "V5", "V6")) {
  
  all_IDs  <- sort(unique(all_paths_indSEM_lag$id))
  all_reps <- sort(unique(all_paths_indSEM_lag$rep))
  
  full_grid <- expand.grid(
    rep = all_reps,
    id = all_IDs,
    predictor = predictor_vars,
    stringsAsFactors = FALSE
  )
  
  indSEM_est_df <- full_grid %>%
    left_join(
      all_paths_indSEM_lag %>%
        select(rep, id, predictor, beta),
      by = c("rep", "id", "predictor")
    ) %>%
    mutate(
      est_lag = ifelse(!is.na(beta) & beta != 0, 1L, 0L)
    ) %>%
    select(rep, id, predictor, est_lag) %>%
    arrange(rep, id, predictor)
  
  return(indSEM_est_df)
}

indSEM_est_uSEMDGM_lag <- generate_indSEM_est_lag_df(
  all_paths_indSEM_lag,
  predictor_vars = c("V2", "V3", "V4", "V5", "V6")
)

score_all_reps_lag <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
  rep_scores <- lapply(rep_ids, function(r) {
    truth_rep <- subset(truth_df, rep == r)
    est_rep   <- subset(est_df, rep == r)
    
    merged <- merge(
      truth_rep,
      est_rep,
      by = c("rep", "id", "predictor"),
      all.x = TRUE,
      sort = FALSE
    )
    
    merged$est_lag[is.na(merged$est_lag)] <- 0L
    
    metrics <- calc_binary_metrics(merged$true_lag, merged$est_lag)
    data.frame(rep = r, metrics, row.names = NULL)
  })
  
  rep_scores <- do.call(rbind, rep_scores)
  rownames(rep_scores) <- NULL
  
  avg_scores <- aggregate(
    cbind(TP, TN, FP, FN, Sensitivity, Specificity, Precision, F1, MCC) ~ 1,
    data = rep_scores,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  avg_scores <- avg_scores[, -1, drop = FALSE]
  
  list(
    rep_scores = rep_scores,
    avg_scores = avg_scores
  )
}

indSEM_scores_uSEM_lag <- score_all_reps_lag(
  truth_df = truth_df,
  est_df = indSEM_est_uSEMDGM_lag
)

indSEM_scores_uSEM_lag$avg_scores
head(indSEM_scores_uSEM_lag$rep_scores)

calculate_metric_cis <- function(rep_scores_df,
                                 metric_cols = c("Sensitivity", "Specificity", "Precision", "F1", "MCC"),
                                 conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  
  ci_results <- rep_scores_df %>%
    select(all_of(metric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    summarise(
      n = sum(!is.na(Value)),
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(n),
      t_crit = qt(1 - alpha/2, df = n - 1),
      Lower_CI = Mean - t_crit * SE,
      Upper_CI = Mean + t_crit * SE,
      .groups = "drop"
    )
  
  return(ci_results)
}

indSEM_ci_results <- calculate_metric_cis(indSEM_scores_uSEM_lag$rep_scores)

print(indSEM_ci_results)

