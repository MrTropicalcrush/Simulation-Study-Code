library(dplyr)
library(tidyr)

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

variable_names <- c("V2", "V3", "V4", "V5", "V6")

# Extract tstats from results
iarimax_tstats_uSEM_wide <- calculate_t_statistics_wide(
  iarimax_results_list = iarimax_results_uSEM_DGM,
  variable_names = variable_names
)

iarimax_tstats_uSEM_wide<- iarimax_tstats_uSEM_wide %>%
  mutate(ID = as.numeric(ID)) %>%           # Ensure ID is numeric
  arrange(simulation, ID)                   # Sort by simulation and ID

generate_iarimax_est_contemp_df <- function(t_stats_wide_df, t_threshold = 1.96) {
  est_df <- t_stats_wide_df %>%
    pivot_longer(
      cols = c("V2", "V3", "V4", "V5", "V6"),
      names_to = "predictor",
      values_to = "t_stat"
    ) %>%
    transmute(
      rep = simulation,
      id = ID,
      predictor = predictor,
      est_contemp = ifelse(!is.na(t_stat) & abs(t_stat) > t_threshold, 1L, 0L)
    ) %>%
    arrange(rep, id, predictor)
  
  return(est_df)
}

# Convert tstats to binary 
iarimax_est_uSEMDGM <- generate_iarimax_est_contemp_df(
  t_stats_wide_df = iarimax_tstats_uSEM_wide,
  t_threshold = 1.96
)

# Score metrics and average across simulations
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

iarimax_scores_uSEM <- score_all_reps_contemp(
  truth_df = truth_df,
  est_df = iarimax_est_uSEMDGM
)

iarimax_scores_uSEM$avg_scores

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

arimax_ci_results <- calculate_metric_cis(iarimax_scores_uSEM$rep_scores)

print(arimax_ci_results)
