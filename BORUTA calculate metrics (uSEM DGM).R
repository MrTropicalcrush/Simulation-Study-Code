generate_boruta_est_contemp_df <- function(BORUTA_results_all) {
  est_list <- vector("list", length(BORUTA_results_all))
  
  predictors <- c("V2", "V3", "V4", "V5", "V6")
  
  for (sim in seq_along(BORUTA_results_all)) {
    sim_rows <- list()
    
    for (i in seq_along(BORUTA_results_all[[sim]])) {
      fd <- BORUTA_results_all[[sim]][[i]]$finalDecision
      
      person_df <- data.frame(
        rep = sim,
        id = i,
        predictor = predictors,
        est_contemp = c(
          ifelse(fd["V2"] == "Confirmed", 1L, 0L),
          ifelse(fd["V3"] == "Confirmed", 1L, 0L),
          ifelse(fd["V4"] == "Confirmed", 1L, 0L),
          ifelse(fd["V5"] == "Confirmed", 1L, 0L),
          ifelse(fd["V6"] == "Confirmed", 1L, 0L)
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


score_all_reps_boruta_contemp <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
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

TSboruta_est_uSEMDGM <- generate_boruta_est_contemp_df(tsBoruta_results_uSEM_DGM)

TSboruta_uSEM_baseline <- score_all_reps_boruta_contemp(
  truth_df = truth_df,
  est_df = TSboruta_est_uSEMDGM
)

TSboruta_uSEM_baseline$avg_scores


library(dplyr)
library(tidyr)

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

boruta_ci_results <- calculate_metric_cis(TSboruta_uSEM_baseline$rep_scores)

print(boruta_ci_results)

