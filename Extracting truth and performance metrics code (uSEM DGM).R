# Extract truth for one individual 
extract_truth_individual <- function(true_mat,
                                     rep_id,
                                     person_id,
                                     outcome_var = 1,
                                     predictor_vars = 2:6) {
  v <- nrow(true_mat)
  
  Phi <- true_mat[, 1:v, drop = FALSE]
  A   <- true_mat[, (v + 1):(2 * v), drop = FALSE]
  
  out <- data.frame(
    rep       = rep_id,
    id        = person_id,
    predictor = paste0("V", predictor_vars),
    true_lag  = as.integer(Phi[outcome_var, predictor_vars] != 0),
    true_contemp = as.integer(A[outcome_var, predictor_vars] != 0),
    stringsAsFactors = FALSE
  )
  
  out$true_any <- as.integer(out$true_lag == 1 | out$true_contemp == 1)
  out
}

# Extract truth for one replication
extract_truth_rep <- function(true_rep_list,
                              rep_id,
                              outcome_var = 1,
                              predictor_vars = 2:6) {
  rep_truth <- lapply(seq_along(true_rep_list), function(person_id) {
    extract_truth_individual(
      true_mat = true_rep_list[[person_id]],
      rep_id = rep_id,
      person_id = person_id,
      outcome_var = outcome_var,
      predictor_vars = predictor_vars
    )
  })
  
  do.call(rbind, rep_truth)
}

# Extract truth for all replications
extract_truth_all <- function(true_data_list,
                              outcome_var = 1,
                              predictor_vars = 2:6) {
  all_truth <- lapply(seq_along(true_data_list), function(rep_id) {
    extract_truth_rep(
      true_rep_list = true_data_list[[rep_id]],
      rep_id = rep_id,
      outcome_var = outcome_var,
      predictor_vars = predictor_vars
    )
  })
  
  truth_df <- do.call(rbind, all_truth)
  rownames(truth_df) <- NULL
  truth_df
}

# Run it and save the master truth table
truth_df <- extract_truth_all(true_data_list)

saveRDS(truth_df, file.path(base_dir, "truth_df.rds"))
write.csv(truth_df, file.path(base_dir, "truth_df.csv"), row.names = FALSE)

#####################################################################################
# Generic metric function
calc_binary_metrics <- function(truth, estimate) {
  stopifnot(length(truth) == length(estimate))
  
  tp <- sum(truth == 1 & estimate == 1, na.rm = TRUE)
  tn <- sum(truth == 0 & estimate == 0, na.rm = TRUE)
  fp <- sum(truth == 0 & estimate == 1, na.rm = TRUE)
  fn <- sum(truth == 1 & estimate == 0, na.rm = TRUE)
  
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  
  f1 <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else {
    NA_real_
  }
  
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- if (denom > 0) ((tp * tn) - (fp * fn)) / denom else NA_real_
  
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

# Scoring one person against truth
score_one_person <- function(truth_df, est_df, rep_id, person_id) {
  truth_sub <- subset(truth_df, rep == rep_id & id == person_id)
  est_sub   <- subset(est_df,   rep == rep_id & id == person_id)
  
  merged <- merge(
    truth_sub,
    est_sub,
    by = c("rep", "id", "predictor"),
    all.x = TRUE,
    sort = FALSE
  )
  
  # Fill missing estimates with 0
  for (nm in c("est_lag", "est_contemp", "est_any")) {
    if (!nm %in% names(merged)) merged[[nm]] <- 0L
    merged[[nm]][is.na(merged[[nm]])] <- 0L
  }
  
  list(
    lag = calc_binary_metrics(merged$true_lag, merged$est_lag),
    contemp = calc_binary_metrics(merged$true_contemp, merged$est_contemp),
    any = calc_binary_metrics(merged$true_any, merged$est_any),
    merged = merged
  )
}

score_one_person(truth_df, est_df, rep_id = 1, person_id = 1)

# Scoring an entire replication
# If your estimated results for a whole replication are in one dataframe with columns:
score_one_rep <- function(truth_df, est_rep_df, rep_id) {
  truth_rep <- subset(truth_df, rep == rep_id)
  est_rep   <- subset(est_rep_df, rep == rep_id)
  
  merged <- merge(
    truth_rep,
    est_rep,
    by = c("rep", "id", "predictor"),
    all.x = TRUE,
    sort = FALSE
  )
  
  for (nm in c("est_lag", "est_contemp", "est_any")) {
    if (!nm %in% names(merged)) merged[[nm]] <- 0L
    merged[[nm]][is.na(merged[[nm]])] <- 0L
  }
  
  list(
    lag = calc_binary_metrics(merged$true_lag, merged$est_lag),
    contemp = calc_binary_metrics(merged$true_contemp, merged$est_contemp),
    any = calc_binary_metrics(merged$true_any, merged$est_any),
    merged = merged
  )
}

# Scoring all replications and averaging
score_all_reps <- function(truth_df, est_df, rep_ids = sort(unique(truth_df$rep))) {
  rep_scores <- lapply(rep_ids, function(r) {
    res <- score_one_rep(truth_df, est_df, rep_id = r)
    
    data.frame(
      rep = r,
      target = c("lag", "contemp", "any"),
      rbind(res$lag, res$contemp, res$any),
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