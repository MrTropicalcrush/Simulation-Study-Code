# Detrending baseline conditions
# fitting a separate linear regression within each person and each variable
detrend_by_id <- function(df, cols, id_var = "id", timevar = "time",
                          min_n_subject = 20, minvar = 0.01,
                          append = TRUE) {
  
  if (!all(c(cols, id_var, timevar) %in% names(df))) {
    stop("Some requested columns are missing from the dataframe.")
  }
  
  if (!is.numeric(df[[timevar]])) {
    stop("timevar must be numeric.")
  }
  
  if (any(is.na(df[[timevar]]))) {
    stop("timevar cannot contain NA values.")
  }
  
  out <- df
  
  for (col in cols) {
    
    new_col <- paste0(col, "_dt")
    out[[new_col]] <- NA_real_
    
    for (id in unique(out[[id_var]])) {
      
      rows <- which(out[[id_var]] == id)
      y <- out[[col]][rows]
      time <- out[[timevar]][rows]
      
      n_obs <- sum(!is.na(y))
      pre_var <- var(y, na.rm = TRUE)
      
      if (all(is.na(y)) ||
          n_obs < min_n_subject ||
          is.na(pre_var) ||
          pre_var < minvar) {
        
        out[[new_col]][rows] <- NA_real_
        
      } else {
        
        fit <- lm(y ~ time, na.action = na.exclude)
        resid_vals <- residuals(fit)
        post_var <- var(resid_vals, na.rm = TRUE)
        
        if (is.na(post_var) || post_var < minvar) {
          out[[new_col]][rows] <- NA_real_
        } else {
          out[[new_col]][rows] <- resid_vals
        }
      }
    }
  }
  
  if (append) {
    return(out)
  } else {
    return(out[, c(id_var, timevar, paste0(cols, "_dt"))])
  }
}

# Run detrending function
sim_data_list_dt <- lapply(sim_data_list, function(df) {
  detrend_by_id(
    df = df,
    cols = c("V1", "V2", "V3", "V4", "V5", "V6"),
    id_var = "id",
    timevar = "time",
    append = TRUE
  )
})

# Replace variables with detrended residuals
replace_with_detrended <- function(df) {
  
  df_dt <- detrend_by_id(
    df = df,
    cols = c("V1", "V2", "V3", "V4", "V5", "V6"),
    id_var = "id",
    timevar = "time",
    append = TRUE
  )
  
  df_dt$V1 <- df_dt$V1_dt
  df_dt$V2 <- df_dt$V2_dt
  df_dt$V3 <- df_dt$V3_dt
  df_dt$V4 <- df_dt$V4_dt
  df_dt$V5 <- df_dt$V5_dt
  df_dt$V6 <- df_dt$V6_dt
  
  df_dt <- df_dt[, c("id", "time", "V1", "V2", "V3", "V4", "V5", "V6")]
  
  return(df_dt)
}

sim_data_list_detrended_uSEM <- lapply(sim_data_list, replace_with_detrended)

saveRDS(sim_data_list_detrended_uSEM, file = "sim_data_list_detrended_uSEM.RDS")


###########################################################################################
var_cols <- c(
  "depressedmood_state",
  "loneliness_state_pmc",
  "socintsatisfaction_state_pmc",
  "responsiveness_state_pmc",
  "selfdisclosure_state_pmc",
  "otherdisclosure_state_pmc"
)

detrend_by_id <- function(df, cols, id_var = "ID", timevar = "time",
                          min_n_subject = 20, minvar = 0.01,
                          append = TRUE) {
  
  if (!all(c(cols, id_var, timevar) %in% names(df))) {
    missing <- c(cols, id_var, timevar)[!c(cols, id_var, timevar) %in% names(df)]
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  out <- df
  
  for (col in cols) {
    
    new_col <- paste0(col, "_dt")
    out[[new_col]] <- NA_real_
    
    for (id in unique(out[[id_var]])) {
      
      rows <- which(out[[id_var]] == id)
      y <- out[[col]][rows]
      time <- out[[timevar]][rows]
      
      n_obs <- sum(!is.na(y))
      pre_var <- var(y, na.rm = TRUE)
      
      if (all(is.na(y)) ||
          n_obs < min_n_subject ||
          is.na(pre_var) ||
          pre_var < minvar) {
        
        out[[new_col]][rows] <- NA_real_
        
      } else {
        
        fit <- lm(y ~ time, na.action = na.exclude)
        resid_vals <- residuals(fit)
        post_var <- var(resid_vals, na.rm = TRUE)
        
        if (is.na(post_var) || post_var < minvar) {
          out[[new_col]][rows] <- NA_real_
        } else {
          out[[new_col]][rows] <- resid_vals
        }
      }
    }
  }
  
  if (append) {
    return(out)
  } else {
    return(out[, c(id_var, timevar, paste0(cols, "_dt"))])
  }
}

replace_with_detrended_VAR <- function(df) {
  
  df_dt <- detrend_by_id(
    df = df,
    cols = var_cols,
    id_var = "ID",
    timevar = "time",
    append = TRUE
  )
  
  for (v in var_cols) {
    df_dt[[v]] <- df_dt[[paste0(v, "_dt")]]
  }
  
  df_dt <- df_dt[, c("ID", "time", var_cols)]
  
  return(df_dt)
}

VAR_baseline_detrended <- lapply(sim_data_list_VAR_baseline, function(df) {
  
  df_dt <- detrend_by_id(
    df = df,
    cols = var_cols,
    id_var = "ID",
    timevar = "time",
    append = TRUE
  )
  
  # Replace original variables
  for (v in var_cols) {
    df_dt[[v]] <- df_dt[[paste0(v, "_dt")]]
  }
  
  df_dt <- df_dt[, c("ID", "time", var_cols)]
  
  return(df_dt)
})

saveRDS(VAR_baseline_detrended, "VAR_baseline_detrended.RDS")

