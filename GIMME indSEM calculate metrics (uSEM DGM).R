library(dplyr)

process_gimme_simulation_folders <- function(directory,
                                             outcome_var = "V1",
                                             predictor_vars = c("V2", "V3", "V4", "V5", "V6")) {
  subfolders <- list.dirs(directory, recursive = FALSE)
  all_simulation_results <- list()
  
  for (subfolder in subfolders) {
    csv_file <- file.path(subfolder, "indivPathEstimates.csv")
    
    if (file.exists(csv_file)) {
      data <- read.csv(csv_file)
      
      relevant_paths <- data %>%
        filter(lhs == outcome_var, rhs %in% predictor_vars) %>%
        select(ID = file, rhs, beta) %>%
        mutate(rep = as.numeric(sub("^Sim", "", basename(subfolder))))
      
      all_simulation_results[[subfolder]] <- relevant_paths
    }
  }
  
  bind_rows(all_simulation_results)
}

all_paths_uSEM <- process_gimme_simulation_folders(
  directory = "C:/Users/WillLi/Documents/mlvar simulation/output_directory_GIMME_uSEM_DGM",   # replace with your actual uSEM GIMME folder
  outcome_var = "V1",
  predictor_vars = c("V2", "V3", "V4", "V5", "V6")
)

predictor_vars <- c("V2", "V3", "V4", "V5", "V6")

all_IDs <- sort(unique(all_paths_uSEM$ID))
all_reps <- sort(unique(all_paths_uSEM$rep))

full_grid <- expand.grid(
  rep = all_reps,
  id = all_IDs,
  predictor = predictor_vars,
  stringsAsFactors = FALSE
)

GIMME_est_uSEMDGM <- full_grid %>%
  left_join(
    all_paths_uSEM %>%
      transmute(
        rep = rep,
        id = ID,
        predictor = rhs,
        beta = beta
      ),
    by = c("rep", "id", "predictor")
  ) %>%
  mutate(
    est_contemp = ifelse(!is.na(beta), 1L, 0L)
  ) %>%
  select(rep, id, predictor, est_contemp) %>%
  arrange(rep, id, predictor)

GIMME_est_uSEMDGM <- GIMME_est_uSEMDGM %>%
  mutate(id = as.numeric(id))

truth_df <- truth_df %>%
  mutate(id = as.numeric(id))


score_all_reps_contemp <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
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
    
    merged$est_contemp[is.na(merged$est_contemp)] <- 0L
    
    metrics <- calc_binary_metrics(merged$true_contemp, merged$est_contemp)
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

GIMME_scores_uSEM <- score_all_reps_contemp(
  truth_df = truth_df,
  est_df = GIMME_est_uSEMDGM
)

GIMME_scores_uSEM$avg_scores
head(GIMME_scores_uSEM$rep_scores)

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

GIMME_ci_results <- calculate_metric_cis(GIMME_scores_uSEM$rep_scores)

print(GIMME_ci_results)

################################################################################
##### indSEM #####
###############################################################################
library(dplyr)
library(tidyr)

process_indSEM_simulations <- function(directory,
                                       outcome_row = "V1",
                                       predictor_vars = c("V2", "V3", "V4", "V5", "V6")) {
  
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
      
      available_predictors <- intersect(predictor_vars, colnames(beta_matrix))
      
      if (length(available_predictors) == 0) {
        message(paste("No matching predictors found in file:", beta_file))
        message(paste("Available columns:", paste(colnames(beta_matrix), collapse = ", ")))
        next
      }
      
      outcome_row_values <- beta_matrix[outcome_row, available_predictors, drop = TRUE]
      
      beta_data <- data.frame(
        rep = as.numeric(sub("^Sim", "", basename(sim_folder))),
        id = as.numeric(individual_id),
        predictor = available_predictors,
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
    message("No valid data found in the simulations.")
    return(NULL)
  }
  
  bind_rows(all_simulation_results) %>%
    arrange(rep, id, predictor)
}

# Run extract
all_paths_indSEM <- process_indSEM_simulations(
  directory = parent_directory,
  outcome_row = "V1",
  predictor_vars = c("V2", "V3", "V4", "V5", "V6")
)

# Covert to binary
generate_indSEM_est_contemp_df <- function(all_paths_indSEM,
                                           predictor_vars = c("V2", "V3", "V4", "V5", "V6")) {
  
  all_IDs  <- sort(unique(all_paths_indSEM$id))
  all_reps <- sort(unique(all_paths_indSEM$rep))
  
  full_grid <- expand.grid(
    rep = all_reps,
    id = all_IDs,
    predictor = predictor_vars,
    stringsAsFactors = FALSE
  )
  
  indSEM_est_df <- full_grid %>%
    left_join(all_paths_indSEM, by = c("rep", "id", "predictor")) %>%
    mutate(
      est_contemp = ifelse(!is.na(beta) & beta != 0, 1L, 0L)
    ) %>%
    select(rep, id, predictor, est_contemp) %>%
    arrange(rep, id, predictor)
  
  return(indSEM_est_df)
  
  
}

parent_directory <- "C:/Users/WillLi/Documents/mlvar simulation/output_directory_GIMME_uSEM_DGM"

# Run binary
indSEM_est_uSEMDGM <- generate_indSEM_est_contemp_df(
  all_paths_indSEM,
  predictor_vars = c("V2", "V3", "V4", "V5", "V6")
)


score_all_reps_contemp <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
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
    
    merged$est_contemp[is.na(merged$est_contemp)] <- 0L
    
    metrics <- calc_binary_metrics(merged$true_contemp, merged$est_contemp)
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

indSEM_scores_uSEM <- score_all_reps_contemp(
  truth_df = truth_df,
  est_df = indSEM_est_uSEMDGM
)

indSEM_scores_uSEM$avg_scores

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

indSEM_ci_results <- calculate_metric_cis(indSEM_scores_uSEM$rep_scores)

print(indSEM_ci_results)
