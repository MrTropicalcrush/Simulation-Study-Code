mat.generate.custom <- function(nvar = 6, AR = 0.3,
                                con_vals = c(0.35, -0.35),
                                lag_vals = c(-0.30, 0.30)) {
  
  # nvar should be 6: V1 = outcome, V2-V6 = predictors
  A   <- matrix(0, nrow = nvar, ncol = nvar)   # contemporaneous
  Phi <- matrix(0, nrow = nvar, ncol = nvar)   # lagged
  
  # AR/self-lag for all variables
  diag(Phi) <- AR
  
  # -----------------------------
  # Shared contemporaneous edges
  # -----------------------------
  # V2 -> V1
  A[1, 2] <- con_vals[1]
  
  # V3 -> V1
  A[1, 3] <- con_vals[2]
  
  # ----------------------
  # Shared lagged edges
  # ----------------------
  # V4(t-1) -> V1(t)
  Phi[1, 4] <- lag_vals[1]
  
  # V5(t-1) -> V1(t)
  Phi[1, 5] <- lag_vals[2]
  
  # -----------------------------
  # Shared predictor network edges
  # -----------------------------
  
  # Lagged network edges
  Phi[2, 3] <- 0.25   # V3(t-1) -> V2(t)
  Phi[3, 4] <- -0.25  # V4(t-1) -> V3(t)
  
  # Contemporaneous network edges
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
  
  # combine
  all <- cbind(Phi, A)
  
  # levels matrix: mark all nonzero shared paths as 'grp'
  all.lvl <- matrix(NA, nrow = nrow(all), ncol = ncol(all))
  all.lvl[all != 0] <- "grp"
  diag(all.lvl[, 1:nvar]) <- "grp"   # AR terms also shared
  
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
    
    # ------------------------------------------------
    # Define pool of optional extra person-specific edges
    # type = "A" for contemporaneous, "Phi" for lagged
    # row = target variable at time t
    # col = predictor variable (same time for A, t-1 for Phi)
    # ------------------------------------------------
    extra_pool <- data.frame(
      type = c("A","Phi","Phi","A",   # outcome edges
               "Phi","Phi","A","A","Phi","A"),  # network edges
      row  = c(1,1,3,4,
               2,3,4,5,6,5),
      col  = c(6,6,2,5,
               3,4,5,6,2,3),
      val  = c(0.30,-0.30,0.25,0.25,
               0.20,-0.20,0.25,-0.25,0.20,0.20)
    )
    # Decide how many extra edges this person gets: 0, 1, or 2
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
    
    # -----------------------------------
    # Add small person-level coefficient noise
    # -----------------------------------
    
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
    
    # -----------------------------------
    # Check contemporaneous matrix validity
    # -----------------------------------
    A.test <- 1 * (A != 0)
    A.test <- A.test + t(A.test)
    
    if ((max(A.test) != 2) &&
        (max(abs(eigen(A, only.values = TRUE)$values)) < 1)) {
      break
    }
  }
  
  # -----------------------------------
  # Simulate time series
  # -----------------------------------
  repeat {
    st <- t + 50   # burn-in
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
  
  # Update levels matrix
  lvl_new <- lvl
  lvl_new[is.na(lvl_new) & paths != 0] <- "ind"
  
  out <- list(
    series = series,
    paths  = paths,
    levels = lvl_new
  )
  
  return(out)
}

#######################################################
#### 2. Simulate Data: Single GIMME-style Condition ###
#######################################################

# -----------------------------
# Single condition parameters
# -----------------------------
v   <- 6         # 1 outcome + 5 predictors
n   <- 102       # number of individuals
t   <- 70       # number of time points
ar  <- 0.2       # small AR/self-lag
rep <- 1:100     # 100 replications

all <- data.frame(
  t = t,
  n = n,
  v = v,
  ar = ar,
  rep = rep
)

# Folder name for each replication
all$folder <- paste0("rep_", all$rep)

# -----------------------------
# DIRECTORY SETUP
# -----------------------------

base_dir <- file.path(getwd(), "Sim_Study", "gimme_style_condition")

data.path  <- file.path(base_dir, "data")
true.path  <- file.path(base_dir, "true")
level.path <- file.path(base_dir, "levels")

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data.path, showWarnings = FALSE)
dir.create(true.path, showWarnings = FALSE)
dir.create(level.path, showWarnings = FALSE)

# Create folders for each replication
for (i in 1:nrow(all)) {
  dir.create(file.path(data.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(true.path,  all$folder[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(level.path, all$folder[i]), recursive = TRUE, showWarnings = FALSE)
}

sim_data_list  <- vector("list", length = nrow(all))
true_data_list <- vector("list", length = nrow(all))
level_data_list <- vector("list", length = nrow(all))

names(sim_data_list)  <- all$folder
names(true_data_list) <- all$folder
names(level_data_list) <- all$folder

# -----------------------------
# Run the simulations
# -----------------------------
for (i in 1:nrow(all)) {
  
  # Generate shared structure for this replication
  res <- mat.generate.custom(
    nvar = all$v[i],
    AR   = all$ar[i],
    con_vals = c(0.35, -0.35),
    lag_vals = c(-0.30, 0.30)
  )
  
  # Storage for this replication
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
    
    # Save csv files
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
    
    # Save in-memory objects
    df <- as.data.frame(out$series)
    df$id <- a
    df$time <- 1:nrow(df)
    df <- df[, c("id", "time", setdiff(names(df), c("id", "time")))]
    
    rep_data_list[[a]]  <- df
    rep_true_list[[a]]  <- out$paths
    rep_level_list[[a]] <- out$levels
  }
  
  # Combine all individuals' time series for this replication
  sim_data_list[[i]]  <- do.call(rbind, rep_data_list)
  true_data_list[[i]] <- rep_true_list
  level_data_list[[i]] <- rep_level_list
}

# -----------------------------
# Output folders for model fits
# -----------------------------
# Use the SAME base_dir as before
ms.path <- file.path(base_dir, "outputMS")
ar.path <- file.path(base_dir, "outputAR")

dir.create(ms.path, recursive = TRUE, showWarnings = FALSE)
dir.create(ar.path, recursive = TRUE, showWarnings = FALSE)

# Create replication subfolders
folders <- all$folder

for (i in 1:nrow(all)) {
  dir.create(file.path(ms.path, folders[i]), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ar.path, folders[i]), recursive = TRUE, showWarnings = FALSE)
}

saveRDS(sim_data_list, file.path(base_dir, "sim_data_list.rds"))
saveRDS(true_data_list, file.path(base_dir, "true_data_list.rds"))
saveRDS(level_data_list, file.path(base_dir, "level_data_list.rds"))


