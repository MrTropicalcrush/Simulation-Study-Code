library(dplyr)

extract_all_within_person_correlations_uSEM <- function(simulated_data_list,
                                                        target_var = "V1",
                                                        vars = c("V2", "V3", "V4", "V5", "V6")) {
  
  sim_correlations_list <- list()
  
  for (sim_idx in seq_along(simulated_data_list)) {
    sim_data <- simulated_data_list[[sim_idx]] %>%
      mutate(id = as.character(id))
    
    sim_correlations <- sim_data %>%
      group_by(id) %>%
      summarize(across(
        all_of(vars),
        ~ cor(.x, .data[[target_var]], use = "complete.obs"),
        .names = "{.col}"
      ), .groups = "drop") %>%
      mutate(simulation = sim_idx)
    
    sim_correlations_list[[sim_idx]] <- sim_correlations
  }
  
  return(sim_correlations_list)
}

raw_correlation_list_uSEM <- extract_all_within_person_correlations_uSEM(
  simulated_data_list = sim_data_list,   # replace with your actual uSEM dataset object
  target_var = "V1",
  vars = c("V2", "V3", "V4", "V5", "V6")
)

raw_correlation_list_uSEM <- lapply(raw_correlation_list_uSEM, function(df) {
  df %>%
    mutate(id = as.numeric(id)) %>%
    arrange(id)
})

library(tidyr)

generate_rawcor_est_contemp_df <- function(raw_correlation_list,
                                           corr_threshold = 0.2) {
  
  est_list <- lapply(seq_along(raw_correlation_list), function(sim) {
    raw_correlation_list[[sim]] %>%
      pivot_longer(
        cols = c("V2", "V3", "V4", "V5", "V6"),
        names_to = "predictor",
        values_to = "raw_r"
      ) %>%
      transmute(
        rep = simulation,
        id = as.numeric(id),
        predictor = predictor,
        est_contemp = ifelse(!is.na(raw_r) & abs(raw_r) > corr_threshold, 1L, 0L)
      )
  })
  
  est_df <- bind_rows(est_list) %>%
    arrange(rep, id, predictor)
  
  return(est_df)
}

rawcor_est_uSEMDGM <- generate_rawcor_est_contemp_df(
  raw_correlation_list = raw_correlation_list_uSEM,
  corr_threshold = 0.2
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

RawCor_scores_uSEM <- score_all_reps_contemp(
  truth_df = truth_df,
  est_df = rawcor_est_uSEMDGM
)

RawCor_scores_uSEM$avg_scores
head(RawCor_scores_uSEM$rep_scores)

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

rawcor_ci_results <- calculate_metric_cis(RawCor_scores_uSEM$rep_scores)

print(rawcor_ci_results)
