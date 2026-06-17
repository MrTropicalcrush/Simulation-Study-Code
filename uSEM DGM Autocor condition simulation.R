# Simulating uSEM DGM Autocorrelation conditions
# Note - baseline condition is already small AR (0.2) condition

mat.generate.custom <- function(nvar = 6,
                                AR_v1 = 0.2,
                                AR_other = 0.2,
                                con_vals = c(0.35, -0.35),
                                lag_vals = c(-0.30, 0.30)) {
  
  # nvar should be 6: V1 = outcome, V2-V6 = predictors
  A   <- matrix(0, nrow = nvar, ncol = nvar)   # contemporaneous
  Phi <- matrix(0, nrow = nvar, ncol = nvar)   # lagged
  
  # Baseline AR/self-lag for all variables
  diag(Phi) <- AR_other
  
  # Overwrite outcome AR only
  Phi[1, 1] <- AR_v1
  
  # -----------------------------
  # Shared contemporaneous edges
  # -----------------------------
  A[1, 2] <- con_vals[1]   # V2 -> V1
  A[1, 3] <- con_vals[2]   # V3 -> V1
  
  # ----------------------
  # Shared lagged edges
  # ----------------------
  Phi[1, 4] <- lag_vals[1] # V4(t-1) -> V1(t)
  Phi[1, 5] <- lag_vals[2] # V5(t-1) -> V1(t)
  
  # -----------------------------
  # Shared predictor network edges
  # -----------------------------
  Phi[2, 3] <- 0.25   # V3(t-1) -> V2(t)
  Phi[3, 4] <- -0.25  # V4(t-1) -> V3(t)
  
  A[4, 5] <- 0.30     # V5 -> V4
  A[5, 6] <- -0.25    # V6 -> V5
  
  # Check contemporaneous matrix validity
  A.test <- 1 * (A != 0)
  A.test <- A.test + t(A.test)
  
  if (max(A.test) == 2) {
    stop("Bidirectional contemporaneous paths detected in A.")
  }
  
  if (max(abs(eigen(A, only.values = TRUE)$values)) >= 1) {
    stop("A matrix is unstable: max eigenvalue >= 1.")
  }
  
  all <- cbind(Phi, A)
  
  # levels matrix: mark all nonzero shared paths as 'grp'
  all.lvl <- matrix(NA, nrow = nrow(all), ncol = ncol(all))
  all.lvl[all != 0] <- "grp"
  diag(all.lvl[, 1:nvar]) <- "grp"
  
  res <- list(
    sub1 = all,
    lvl1 = all.lvl
  )
  
  return(res)
}

ts.generate.custom <- function(mat, lvl, t,
                               extra_edge_probs = c(0.20, 0.45, 0.35),
                               con_noise_sd = 0.05,
                               lag_noise_sd = 0.05,
                               ar_noise_sd  = 0.05) {
  
  repeat {
    
    v <- ncol(mat) / 2
    Phi <- mat[, 1:v]
    A   <- mat[, (v + 1):(v * 2)]
    
    extra_pool <- data.frame(
      type = c("A","Phi","Phi","A",
               "Phi","Phi","A","A","Phi","A"),
      row  = c(1,1,3,4,
               2,3,4,5,6,5),
      col  = c(6,6,2,5,
               3,4,5,6,2,3),
      val  = c(0.30,-0.30,0.25,0.25,
               0.20,-0.20,0.25,-0.25,0.20,0.20)
    )
    
    n_extra <- sample(0:2, size = 1, prob = extra_edge_probs)
    
    if (n_extra > 0) {
      chosen <- extra_pool[sample(1:nrow(extra_pool), n_extra, replace = FALSE), ]
      
      for (j in 1:nrow(chosen)) {
        if (chosen$type[j] == "A") {
          A[chosen$row[j], chosen$col[j]] <- chosen$val[j]
        } else if (chosen$type[j] == "Phi") {
          Phi[chosen$row[j], chosen$col[j]] <- chosen$val[j]
        }
      }
    }
    
    # Noise for nonzero contemporaneous paths
    noise.inds.A <- which(A != 0, arr.ind = TRUE)
    if (nrow(noise.inds.A) > 0) {
      A[noise.inds.A] <- A[noise.inds.A] + rnorm(nrow(noise.inds.A), 0, con_noise_sd)
    }
    
    # Noise for nonzero lagged paths excluding AR diagonal
    noise.inds.Phi <- which(Phi != 0, arr.ind = TRUE)
    noise.inds.Phi.offdiag <- noise.inds.Phi[noise.inds.Phi[,1] != noise.inds.Phi[,2], , drop = FALSE]
    if (nrow(noise.inds.Phi.offdiag) > 0) {
      Phi[noise.inds.Phi.offdiag] <- Phi[noise.inds.Phi.offdiag] +
        rnorm(nrow(noise.inds.Phi.offdiag), 0, lag_noise_sd)
    }
    
    # Noise for AR/self-lag terms
    noise.inds.Phi.diag <- noise.inds.Phi[noise.inds.Phi[,1] == noise.inds.Phi[,2], , drop = FALSE]
    if (nrow(noise.inds.Phi.diag) > 0) {
      Phi[noise.inds.Phi.diag] <- Phi[noise.inds.Phi.diag] +
        rnorm(nrow(noise.inds.Phi.diag), 0, ar_noise_sd)
    }
    
    # Check contemporaneous matrix validity
    A.test <- 1 * (A != 0)
    A.test <- A.test + t(A.test)
    
    if ((max(A.test) != 2) &&
        (max(abs(eigen(A, only.values = TRUE)$values)) < 1)) {
      break
    }
  }
  
  repeat {
    st <- t + 50
    noise <- matrix(rnorm(v * st, 0, 1), nrow = v)
    I <- diag(v)
    time <- matrix(0, nrow = v, ncol = st + 1)
    time1 <- matrix(0, nrow = v, ncol = st)
    
    for (i in 1:st) {
      time1[, i]  <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
      time[, i+1] <- time1[, i]
    }
    
    time1  <- time1[, 51:(50 + t)]
    series <- t(time1)
    paths  <- cbind(Phi, A)
    
    if (abs(max(series, na.rm = TRUE)) < 20 &&
        abs(min(series, na.rm = TRUE)) > .01 &&
        abs(min(series, na.rm = TRUE)) < 20) {
      break
    }
  }
  
  lvl_new <- lvl
  lvl_new[is.na(lvl_new) & paths != 0] <- "ind"
  
  out <- list(
    series = series,
    paths  = paths,
    levels = lvl_new
  )
  
  return(out)
}

############################################################
#### 2. Simulate Data: Medium + Large #########
############################################################

v   <- 6
n   <- 102
t   <- 70
rep <- 1:100

# Baseline AR for V2-V6
AR_other <- 0.2

# V1 AR conditions only
ar_conditions <- data.frame(
  condition = c("medium", "large"),
  AR_v1 = c(0.5, 0.8)
)

all <- expand.grid(
  rep = rep,
  condition_row = 1:nrow(ar_conditions)
)

all$condition <- ar_conditions$condition[all$condition_row]
all$AR_v1     <- ar_conditions$AR_v1[all$condition_row]
all$AR_other  <- AR_other
all$t         <- t
all$n         <- n
all$v         <- v

all$folder <- paste0("AR_", all$condition, "_rep_", all$rep)

base_dir <- file.path(getwd(), "Sim_Study", "gimme_style_condition_V1_AR_extra")

data.path  <- file.path(base_dir, "data")
true.path  <- file.path(base_dir, "true")
level.path <- file.path(base_dir, "levels")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data.path, showWarnings = FALSE)
dir.create(true.path, showWarnings = FALSE)
dir.create(level.path, showWarnings = FALSE)

for (i in 1:nrow(all)) {
  dir.create(file.path(data.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(true.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(level.path, all$folder[i]), recursive = TRUE, showWarnings = FALSE)
}

sim_data_list   <- vector("list", length = nrow(all))
true_data_list  <- vector("list", length = nrow(all))
level_data_list <- vector("list", length = nrow(all))

names(sim_data_list)   <- all$folder
names(true_data_list)  <- all$folder
names(level_data_list) <- all$folder

for (i in 1:nrow(all)) {
  
  res <- mat.generate.custom(
    nvar     = all$v[i],
    AR_v1    = all$AR_v1[i],
    AR_other = all$AR_other[i],
    con_vals = c(0.35, -0.35),
    lag_vals = c(-0.30, 0.30)
  )
  
  rep_data_list  <- vector("list", length = all$n[i])
  rep_true_list  <- vector("list", length = all$n[i])
  rep_level_list <- vector("list", length = all$n[i])
  
  for (a in 1:all$n[i]) {
    
    out <- ts.generate.custom(
      mat = res$sub1,
      lvl = res$lvl1,
      t   = all$t[i],
      extra_edge_probs = c(0.35, 0.40, 0.25),
      con_noise_sd = 0.05,
      lag_noise_sd = 0.05,
      ar_noise_sd  = 0.05
    )
    
    out$series <- round(out$series, digits = 5)
    
    write.csv(
      out$series,
      file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    colnames(out$paths) <- c(
      paste0("lag_V", 1:6),
      paste0("contemp_V", 1:6)
    )
    rownames(out$paths) <- paste0("V", 1:6)
    
    write.csv(
      out$paths,
      file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    colnames(out$levels) <- c(
      paste0("lag_V", 1:6),
      paste0("contemp_V", 1:6)
    )
    rownames(out$levels) <- paste0("V", 1:6)
    
    write.csv(
      out$levels,
      file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    df <- as.data.frame(out$series)
    df$id <- a
    df$time <- 1:nrow(df)
    df <- df[, c("id", "time", setdiff(names(df), c("id", "time")))]
    
    rep_data_list[[a]]  <- df
    rep_true_list[[a]]  <- out$paths
    rep_level_list[[a]] <- out$levels
  }
  
  sim_data_list[[i]]   <- do.call(rbind, rep_data_list)
  true_data_list[[i]]  <- rep_true_list
  level_data_list[[i]] <- rep_level_list
}

ms.path <- file.path(base_dir, "outputMS")
ar.path <- file.path(base_dir, "outputAR")

dir.create(ms.path, recursive = TRUE, showWarnings = FALSE)
dir.create(ar.path, recursive = TRUE, showWarnings = FALSE)

folders <- all$folder

for (i in 1:nrow(all)) {
  dir.create(file.path(ms.path, folders[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ar.path, folders[i]), recursive = TRUE, showWarnings = FALSE)
}

sim_data_ARmedium <- sim_data_list[grepl("^AR_medium_rep_", names(sim_data_list))]
sim_data_ARlarge  <- sim_data_list[grepl("^AR_large_rep_",  names(sim_data_list))]

saveRDS(sim_data_ARmedium, file = "sim_data_ARmedium.RDS")
saveRDS(sim_data_ARlarge, file = "sim_data_ARlarge.RDS")

true_data_ARmedium <- true_data_list[grepl("^AR_medium_rep_", names(true_data_list))]
true_data_ARlarge  <- true_data_list[grepl("^AR_large_rep_",  names(true_data_list))]



##################################################################
#### AR = 1 Condition ####
#################################################################
mat.generate.custom <- function(nvar = 6,
                                AR_v1 = 1.0,
                                AR_other = 0.2,
                                con_vals = c(0.35, -0.35),
                                lag_vals = c(-0.30, 0.30)) {
  
  # nvar should be 6: V1 = outcome, V2-V6 = predictors
  A   <- matrix(0, nrow = nvar, ncol = nvar)   # contemporaneous
  Phi <- matrix(0, nrow = nvar, ncol = nvar)   # lagged
  
  # Baseline AR/self-lag for all variables
  diag(Phi) <- AR_other
  
  # Overwrite outcome AR only
  Phi[1, 1] <- AR_v1
  
  # -----------------------------
  # Shared contemporaneous edges
  # -----------------------------
  A[1, 2] <- con_vals[1]   # V2 -> V1
  A[1, 3] <- con_vals[2]   # V3 -> V1
  
  # ----------------------
  # Shared lagged edges
  # ----------------------
  Phi[1, 4] <- lag_vals[1] # V4(t-1) -> V1(t)
  Phi[1, 5] <- lag_vals[2] # V5(t-1) -> V1(t)
  
  # -----------------------------
  # Shared predictor network edges
  # -----------------------------
  Phi[2, 3] <- 0.25   # V3(t-1) -> V2(t)
  Phi[3, 4] <- -0.25  # V4(t-1) -> V3(t)
  
  A[4, 5] <- 0.30     # V5 -> V4
  A[5, 6] <- -0.25    # V6 -> V5
  
  # Check contemporaneous matrix validity only
  A.test <- 1 * (A != 0)
  A.test <- A.test + t(A.test)
  
  if (max(A.test) == 2) {
    stop("Bidirectional contemporaneous paths detected in A.")
  }
  
  if (max(abs(eigen(A, only.values = TRUE)$values)) >= 1) {
    stop("A matrix is unstable: max eigenvalue >= 1.")
  }
  
  all <- cbind(Phi, A)
  
  all.lvl <- matrix(NA, nrow = nrow(all), ncol = ncol(all))
  all.lvl[all != 0] <- "grp"
  diag(all.lvl[, 1:nvar]) <- "grp"
  
  list(
    sub1 = all,
    lvl1 = all.lvl
  )
}

ts.generate.AR1 <- function(mat, lvl, t,
                            extra_edge_probs = c(0.20, 0.45, 0.35),
                            con_noise_sd = 0.05,
                            lag_noise_sd = 0.05,
                            ar_noise_sd  = 0.05) {
  
  # -----------------------------------
  # Build person-specific model once
  # -----------------------------------
  repeat {
    
    v <- ncol(mat) / 2
    Phi <- mat[, 1:v]
    A   <- mat[, (v + 1):(v * 2)]
    
    extra_pool <- data.frame(
      type = c("A","Phi","Phi","A",
               "Phi","Phi","A","A","Phi","A"),
      row  = c(1,1,3,4,
               2,3,4,5,6,5),
      col  = c(6,6,2,5,
               3,4,5,6,2,3),
      val  = c(0.30,-0.30,0.25,0.25,
               0.20,-0.20,0.25,-0.25,0.20,0.20)
    )
    
    n_extra <- sample(0:2, size = 1, prob = extra_edge_probs)
    
    if (n_extra > 0) {
      chosen <- extra_pool[sample(1:nrow(extra_pool), n_extra, replace = FALSE), ]
      
      for (j in 1:nrow(chosen)) {
        if (chosen$type[j] == "A") {
          A[chosen$row[j], chosen$col[j]] <- chosen$val[j]
        } else if (chosen$type[j] == "Phi") {
          Phi[chosen$row[j], chosen$col[j]] <- chosen$val[j]
        }
      }
    }
    
    # Add noise to nonzero contemporaneous paths
    noise.inds.A <- which(A != 0, arr.ind = TRUE)
    if (nrow(noise.inds.A) > 0) {
      A[noise.inds.A] <- A[noise.inds.A] + rnorm(nrow(noise.inds.A), 0, con_noise_sd)
    }
    
    # Add noise to nonzero lagged off-diagonal paths
    noise.inds.Phi <- which(Phi != 0, arr.ind = TRUE)
    noise.inds.Phi.offdiag <- noise.inds.Phi[noise.inds.Phi[,1] != noise.inds.Phi[,2], , drop = FALSE]
    if (nrow(noise.inds.Phi.offdiag) > 0) {
      Phi[noise.inds.Phi.offdiag] <- Phi[noise.inds.Phi.offdiag] +
        rnorm(nrow(noise.inds.Phi.offdiag), 0, lag_noise_sd)
    }
    
    # Add noise to AR/self-lag terms
    noise.inds.Phi.diag <- noise.inds.Phi[noise.inds.Phi[,1] == noise.inds.Phi[,2], , drop = FALSE]
    if (nrow(noise.inds.Phi.diag) > 0) {
      Phi[noise.inds.Phi.diag] <- Phi[noise.inds.Phi.diag] +
        rnorm(nrow(noise.inds.Phi.diag), 0, ar_noise_sd)
    }
    
    # Only check contemporaneous validity
    A.test <- 1 * (A != 0)
    A.test <- A.test + t(A.test)
    
    if ((max(A.test) != 2) &&
        (max(abs(eigen(A, only.values = TRUE)$values)) < 1)) {
      break
    }
  }
  
  # -----------------------------------
  # Simulate ONCE only (no rejection loop)
  # -----------------------------------
  st <- t + 50
  noise <- matrix(rnorm(v * st, 0, 1), nrow = v)
  I <- diag(v)
  time <- matrix(0, nrow = v, ncol = st + 1)
  time1 <- matrix(0, nrow = v, ncol = st)
  
  for (i in 1:st) {
    time1[, i]  <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
    time[, i+1] <- time1[, i]
  }
  
  time1  <- time1[, 51:(50 + t)]
  series <- t(time1)
  paths  <- cbind(Phi, A)
  
  lvl_new <- lvl
  lvl_new[is.na(lvl_new) & paths != 0] <- "ind"
  
  list(
    series = series,
    paths  = paths,
    levels = lvl_new
  )
}

#### Simulate Data: AR1 condition only ###############
v   <- 6
n   <- 102
t   <- 70
rep <- 1:100

AR_other <- 0.2
AR_v1    <- 1.0

all <- data.frame(
  t = t,
  n = n,
  v = v,
  rep = rep
)

all$folder <- paste0("AR_AR1_rep_", all$rep)

base_dir <- file.path(getwd(), "Sim_Study", "gimme_style_condition_AR1_only")

data.path  <- file.path(base_dir, "data")
true.path  <- file.path(base_dir, "true")
level.path <- file.path(base_dir, "levels")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data.path, showWarnings = FALSE)
dir.create(true.path, showWarnings = FALSE)
dir.create(level.path, showWarnings = FALSE)

for (i in 1:nrow(all)) {
  dir.create(file.path(data.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(true.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(level.path, all$folder[i]), recursive = TRUE, showWarnings = FALSE)
}

sim_data_AR1   <- vector("list", length = nrow(all))
true_data_AR1  <- vector("list", length = nrow(all))
level_data_AR1 <- vector("list", length = nrow(all))

names(sim_data_AR1)   <- all$folder
names(true_data_AR1)  <- all$folder
names(level_data_AR1) <- all$folder

for (i in 1:nrow(all)) {
  
  message("Running replication: ", i, " / ", nrow(all))
  
  res <- mat.generate.custom(
    nvar     = all$v[i],
    AR_v1    = AR_v1,
    AR_other = AR_other,
    con_vals = c(0.35, -0.35),
    lag_vals = c(-0.30, 0.30)
  )
  
  rep_data_list  <- vector("list", length = all$n[i])
  rep_true_list  <- vector("list", length = all$n[i])
  rep_level_list <- vector("list", length = all$n[i])
  
  for (a in 1:all$n[i]) {
    
    out <- ts.generate.AR1(
      mat = res$sub1,
      lvl = res$lvl1,
      t   = all$t[i],
      extra_edge_probs = c(0.35, 0.40, 0.25),
      con_noise_sd = 0.05,
      lag_noise_sd = 0.05,
      ar_noise_sd  = 0.05
    )
    
    out$series <- round(out$series, digits = 5)
    
    write.csv(
      out$series,
      file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    colnames(out$paths) <- c(
      paste0("lag_V", 1:6),
      paste0("contemp_V", 1:6)
    )
    rownames(out$paths) <- paste0("V", 1:6)
    
    write.csv(
      out$paths,
      file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    colnames(out$levels) <- c(
      paste0("lag_V", 1:6),
      paste0("contemp_V", 1:6)
    )
    rownames(out$levels) <- paste0("V", 1:6)
    
    write.csv(
      out$levels,
      file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    df <- as.data.frame(out$series)
    df$id <- a
    df$time <- 1:nrow(df)
    df <- df[, c("id", "time", setdiff(names(df), c("id", "time")))]
    
    rep_data_list[[a]]  <- df
    rep_true_list[[a]]  <- out$paths
    rep_level_list[[a]] <- out$levels
  }
  
  sim_data_AR1[[i]]   <- do.call(rbind, rep_data_list)
  true_data_AR1[[i]]  <- rep_true_list
  level_data_AR1[[i]] <- rep_level_list
}

ms.path <- file.path(base_dir, "outputMS")
ar.path <- file.path(base_dir, "outputAR")

dir.create(ms.path, recursive = TRUE, showWarnings = FALSE)
dir.create(ar.path, recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(all)) {
  dir.create(file.path(ms.path, all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ar.path, all$folder[i]), recursive = TRUE, showWarnings = FALSE)
}

saveRDS(sim_data_AR1,   file.path(base_dir, "sim_data_AR1.rds"))


##################################################################
#### AR = 1 Null Predictor-to-Outcome Condition ####
##################################################################

mat.generate.custom.AR1.null <- function(nvar = 6,
                                         AR_v1 = 1.0,
                                         AR_other = 0.2) {
  
  # nvar should be 6: V1 = outcome, V2-V6 = predictors
  A   <- matrix(0, nrow = nvar, ncol = nvar)   # contemporaneous
  Phi <- matrix(0, nrow = nvar, ncol = nvar)   # lagged
  
  # Baseline AR/self-lag for all variables
  diag(Phi) <- AR_other
  
  # Outcome AR = 1
  Phi[1, 1] <- AR_v1
  
  # ----------------------------------------------------
  # IMPORTANT:
  # No predictor -> outcome effects
  # ----------------------------------------------------
  # A[1, 2:6]   = 0 by default
  # Phi[1, 2:6] = 0 by default
  
  # -----------------------------
  # Keep predictor network edges
  # -----------------------------
  # These do NOT point to outcome V1, so they can stay
  Phi[2, 3] <- 0.25   # V3(t-1) -> V2(t)
  Phi[3, 4] <- -0.25  # V4(t-1) -> V3(t)
  
  A[4, 5] <- 0.30     # V5 -> V4
  A[5, 6] <- -0.25    # V6 -> V5
  
  # Check contemporaneous matrix validity only
  A.test <- 1 * (A != 0)
  A.test <- A.test + t(A.test)
  
  if (max(A.test) == 2) {
    stop("Bidirectional contemporaneous paths detected in A.")
  }
  
  if (max(abs(eigen(A, only.values = TRUE)$values)) >= 1) {
    stop("A matrix is unstable: max eigenvalue >= 1.")
  }
  
  all <- cbind(Phi, A)
  
  all.lvl <- matrix(NA, nrow = nrow(all), ncol = ncol(all))
  all.lvl[all != 0] <- "grp"
  diag(all.lvl[, 1:nvar]) <- "grp"
  
  list(
    sub1 = all,
    lvl1 = all.lvl
  )
}

ts.generate.AR1.null <- function(mat, lvl, t,
                                 extra_edge_probs = c(0.20, 0.45, 0.35),
                                 con_noise_sd = 0.05,
                                 lag_noise_sd = 0.05,
                                 ar_noise_sd  = 0.05) {
  
  repeat {
    
    v <- ncol(mat) / 2
    Phi <- mat[, 1:v]
    A   <- mat[, (v + 1):(v * 2)]
    
    # ----------------------------------------------------
    # Extra edges allowed ONLY among predictors
    # No row 1 edges from predictors to outcome
    # ----------------------------------------------------
    extra_pool <- data.frame(
      type = c("Phi", "Phi", "A", "A", "Phi", "A"),
      row  = c(2, 3, 4, 5, 6, 5),
      col  = c(3, 4, 5, 6, 2, 3),
      val  = c(0.20, -0.20, 0.25, -0.25, 0.20, 0.20)
    )
    
    n_extra <- sample(0:2, size = 1, prob = extra_edge_probs)
    
    if (n_extra > 0) {
      chosen <- extra_pool[sample(1:nrow(extra_pool), n_extra, replace = FALSE), ]
      
      for (j in 1:nrow(chosen)) {
        if (chosen$type[j] == "A") {
          A[chosen$row[j], chosen$col[j]] <- chosen$val[j]
        } else if (chosen$type[j] == "Phi") {
          Phi[chosen$row[j], chosen$col[j]] <- chosen$val[j]
        }
      }
    }
    
    # Add noise to nonzero contemporaneous paths
    noise.inds.A <- which(A != 0, arr.ind = TRUE)
    if (nrow(noise.inds.A) > 0) {
      A[noise.inds.A] <- A[noise.inds.A] + rnorm(nrow(noise.inds.A), 0, con_noise_sd)
    }
    
    # Add noise to nonzero lagged off-diagonal paths
    noise.inds.Phi <- which(Phi != 0, arr.ind = TRUE)
    noise.inds.Phi.offdiag <- noise.inds.Phi[
      noise.inds.Phi[, 1] != noise.inds.Phi[, 2], ,
      drop = FALSE
    ]
    
    if (nrow(noise.inds.Phi.offdiag) > 0) {
      Phi[noise.inds.Phi.offdiag] <- Phi[noise.inds.Phi.offdiag] +
        rnorm(nrow(noise.inds.Phi.offdiag), 0, lag_noise_sd)
    }
    
    # Add noise to AR/self-lag terms
    noise.inds.Phi.diag <- noise.inds.Phi[
      noise.inds.Phi[, 1] == noise.inds.Phi[, 2], ,
      drop = FALSE
    ]
    
    if (nrow(noise.inds.Phi.diag) > 0) {
      Phi[noise.inds.Phi.diag] <- Phi[noise.inds.Phi.diag] +
        rnorm(nrow(noise.inds.Phi.diag), 0, ar_noise_sd)
    }
    
    # Force outcome AR exactly 1 after noise
    Phi[1, 1] <- 1
    
    # Force no predictor -> outcome paths after all changes/noise
    Phi[1, 2:v] <- 0
    A[1, 2:v]   <- 0
    
    # Check contemporaneous validity
    A.test <- 1 * (A != 0)
    A.test <- A.test + t(A.test)
    
    if ((max(A.test) != 2) &&
        (max(abs(eigen(A, only.values = TRUE)$values)) < 1)) {
      break
    }
  }
  
  # -----------------------------------
  # Simulate ONCE only
  # -----------------------------------
  st <- t + 50
  noise <- matrix(rnorm(v * st, 0, 1), nrow = v)
  I <- diag(v)
  time <- matrix(0, nrow = v, ncol = st + 1)
  time1 <- matrix(0, nrow = v, ncol = st)
  
  for (i in 1:st) {
    time1[, i]  <- solve(I - A) %*% (Phi %*% time[, i] + noise[, i])
    time[, i+1] <- time1[, i]
  }
  
  time1  <- time1[, 51:(50 + t)]
  series <- t(time1)
  paths  <- cbind(Phi, A)
  
  lvl_new <- lvl
  lvl_new[is.na(lvl_new) & paths != 0] <- "ind"
  
  list(
    series = series,
    paths  = paths,
    levels = lvl_new
  )
}

##################################################################
#### Simulate Data: AR = 1 Null Predictor-to-Outcome Condition ####
##################################################################

v   <- 6
n   <- 102
t   <- 70
rep <- 1:100

AR_other <- 0.2
AR_v1    <- 1.0

all <- data.frame(
  t = t,
  n = n,
  v = v,
  rep = rep
)

all$folder <- paste0("AR1_null_predictors_rep_", all$rep)

base_dir <- file.path(
  getwd(),
  "Sim_Study",
  "gimme_style_condition_AR1_null_predictors"
)

data.path  <- file.path(base_dir, "data")
true.path  <- file.path(base_dir, "true")
level.path <- file.path(base_dir, "levels")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data.path, recursive = TRUE, showWarnings = FALSE)
dir.create(true.path, recursive = TRUE, showWarnings = FALSE)
dir.create(level.path, recursive = TRUE, showWarnings = FALSE)

for (i in 1:nrow(all)) {
  dir.create(file.path(data.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(true.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(level.path, all$folder[i]), recursive = TRUE, showWarnings = FALSE)
}

sim_data_AR1_nullpredictors   <- vector("list", length = nrow(all))
true_data_AR1_nullpredictors  <- vector("list", length = nrow(all))
level_data_AR1_nullpredictors <- vector("list", length = nrow(all))

names(sim_data_AR1_nullpredictors)   <- all$folder
names(true_data_AR1_nullpredictors)  <- all$folder
names(level_data_AR1_nullpredictors) <- all$folder

for (i in 1:nrow(all)) {
  
  message("Running replication: ", i, " / ", nrow(all))
  
  res <- mat.generate.custom.AR1.null(
    nvar     = all$v[i],
    AR_v1    = AR_v1,
    AR_other = AR_other
  )
  
  rep_data_list  <- vector("list", length = all$n[i])
  rep_true_list  <- vector("list", length = all$n[i])
  rep_level_list <- vector("list", length = all$n[i])
  
  for (a in 1:all$n[i]) {
    
    out <- ts.generate.AR1.null(
      mat = res$sub1,
      lvl = res$lvl1,
      t   = all$t[i],
      extra_edge_probs = c(0.35, 0.40, 0.25),
      con_noise_sd = 0.05,
      lag_noise_sd = 0.05,
      ar_noise_sd  = 0.05
    )
    
    out$series <- round(out$series, digits = 5)
    
    write.csv(
      out$series,
      file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    colnames(out$paths) <- c(
      paste0("lag_V", 1:6),
      paste0("contemp_V", 1:6)
    )
    
    rownames(out$paths) <- paste0("V", 1:6)
    
    write.csv(
      out$paths,
      file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    colnames(out$levels) <- c(
      paste0("lag_V", 1:6),
      paste0("contemp_V", 1:6)
    )
    
    rownames(out$levels) <- paste0("V", 1:6)
    
    write.csv(
      out$levels,
      file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
      row.names = FALSE
    )
    
    df <- as.data.frame(out$series)
    colnames(df) <- paste0("V", 1:v)
    
    df$id <- a
    df$time <- 1:nrow(df)
    
    df <- df[, c("id", "time", setdiff(names(df), c("id", "time")))]
    
    rep_data_list[[a]]  <- df
    rep_true_list[[a]]  <- out$paths
    rep_level_list[[a]] <- out$levels
  }
  
  sim_data_AR1_nullpredictors[[i]]   <- do.call(rbind, rep_data_list)
  true_data_AR1_nullpredictors[[i]]  <- rep_true_list
  level_data_AR1_nullpredictors[[i]] <- rep_level_list
}

##################################################################
#### Save objects ####
##################################################################

saveRDS(
  sim_data_AR1_nullpredictors,
  file.path(base_dir, "sim_data_AR1_nullpredictors_uSEM.rds")
)

saveRDS(
  true_data_AR1_nullpredictors,
  file.path(base_dir, "true_data_AR1_nullpredictors_uSEM.rds")
)

saveRDS(
  level_data_AR1_nullpredictors,
  file.path(base_dir, "level_data_AR1_nullpredictors.rds")
)