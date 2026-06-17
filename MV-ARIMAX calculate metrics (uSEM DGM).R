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

mv_iarimax_uSEM_extracted <- extract_arimax_results_from_simulations(MV_ARIMAX_uSEM_DGM)

generate_mv_iarimax_est_contemp_df <- function(extracted_results_list,
                                               predictor_map = c(
                                                 "V2_PSD" = "V2",
                                                 "V3_PSD" = "V3",
                                                 "V4_PSD" = "V4",
                                                 "V5_PSD" = "V5",
                                                 "V6_PSD" = "V6"
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
        est_contemp = ifelse(abs(person_tstats[available_preds]) > t_threshold, 1L, 0L),
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


MViarimax_est_uSEMDGM <- generate_mv_iarimax_est_contemp_df(
  extracted_results_list = mv_iarimax_uSEM_extracted,
  predictor_map = c(
    "V2_PSD" = "V2",
    "V3_PSD" = "V3",
    "V4_PSD" = "V4",
    "V5_PSD" = "V5",
    "V6_PSD" = "V6"
  ),
  t_threshold = 1.96
)

score_one_rep_mv_iarimax_contemp <- function(truth_df, est_df, rep_id) {
  truth_rep <- subset(truth_df, rep == rep_id)
  est_rep   <- subset(est_df, rep == rep_id)
  
  merged <- merge(
    truth_rep,
    est_rep,
    by = c("rep", "id", "predictor"),
    all.x = TRUE,
    sort = FALSE
  )
  
  if (!"est_contemp" %in% names(merged)) {
    merged$est_contemp <- 0L
  }
  merged$est_contemp[is.na(merged$est_contemp)] <- 0L
  
  list(
    contemp = calc_binary_metrics(merged$true_contemp, merged$est_contemp),
    merged = merged
  )
}

score_all_reps_mv_iarimax_contemp <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
  rep_scores <- lapply(rep_ids, function(r) {
    res <- score_one_rep_mv_iarimax_contemp(truth_df, est_df, rep_id = r)
    
    data.frame(
      rep = r,
      target = "contemp",
      res$contemp,
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

MViarimax_scores_uSEM <- score_all_reps_mv_iarimax_contemp(
  truth_df = truth_df,
  est_df = MViarimax_est_uSEMDGM
)

MViarimax_scores_uSEM$avg_scores

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

MVarimax_ci_results <- calculate_metric_cis(MViarimax_scores_uSEM$rep_scores)

print(MVarimax_ci_results)

