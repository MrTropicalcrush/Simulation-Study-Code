########################################
#### Format uSEM DGM data for GIMME ####
########################################

# Function to restructure each rep into a list of individuals
transform_to_individual_lists <- function(simulations_list, id_col = "id") {
  
  transformed_list <- lapply(simulations_list, function(simulation_data) {
    
    if (!is.data.frame(simulation_data)) {
      stop("Each rep must be a data.frame.")
    }
    
    if (!(id_col %in% names(simulation_data))) {
      stop(paste0(
        "Column '", id_col, "' not found. Available columns are: ",
        paste(names(simulation_data), collapse = ", ")
      ))
    }
    
    split_list <- split(simulation_data, simulation_data[[id_col]])
    
    # Reorder numerically by individual ID
    split_list <- split_list[order(as.numeric(names(split_list)))]
    
    split_list
  })
  
  return(transformed_list)
}

# Split raw uSEM DGM data by individual
transformed_data_uSEM <- transform_to_individual_lists(sim_data_list, id_col = "id")


# Function to remove columns for every individual in every rep
remove_columns_from_simulations <- function(simulations_list, columns_to_remove = c("id", "time")) {
  
  lapply(simulations_list, function(simulation) {
    lapply(simulation, function(individual_data) {
      individual_data[, !colnames(individual_data) %in% columns_to_remove, drop = FALSE]
    })
  })
}

# Remove id and time so GIMME only gets V variables
GIMME_uSEM_DGM <- remove_columns_from_simulations(
  transformed_data_uSEM,
  columns_to_remove = c("id", "time")
)

####################################################
#### RUNNING GIMME SIMULTANEOUSLY USING SCRIPTS ####
####################################################

# Save the function to a file
writeLines(c(
  "run_gimme_sequential <- function(simulations_list, base_output_dir) {",
  "  if (!dir.exists(base_output_dir)) {",
  "    dir.create(base_output_dir, recursive = TRUE)",
  "  }",
  "  for (i in seq_along(simulations_list)) {",
  "    current_simulation <- simulations_list[[i]]",
  "    sim_output_dir <- file.path(base_output_dir, paste0('Sim', i))",
  "    if (!dir.exists(sim_output_dir)) {",
  "      dir.create(sim_output_dir, recursive = TRUE)",
  "    }",
  "    tryCatch({",
  "      gimme(",
  "        data = current_simulation,",
  "        out = sim_output_dir,",
  "        subgroup = FALSE,",
  "        groupcutoff = 0.75,",
  "        outcome = 'V1'",
  "      )",
  "      message(sprintf('Simulation %d completed successfully', i))",
  "    }, error = function(e) {",
  "      message(sprintf('Simulation %d failed with error: %s', i, e$message))",
  "    })",
  "  }",
  "}"
), "run_gimme_sequential.R")

library(parallel)

# Split simulations into batches
split_into_batches <- function(simulations_list, batch_size) {
  split(simulations_list, ceiling(seq_along(simulations_list) / batch_size))
}

# Save each batch as an R script
save_gimme_batch_script <- function(batch, batch_index, output_dir, start_number) {
  
  simulation_names <- sprintf("Sim%d", seq(start_number, start_number + length(batch) - 1))
  named_batch <- setNames(batch, simulation_names)
  
  batch_file <- file.path(output_dir, paste0("batch_", batch_index, ".rds"))
  saveRDS(named_batch, batch_file)
  
  script_content <- paste0(
    "library(gimme)\n",
    "source('run_gimme_sequential.R')\n",
    "start_time <- Sys.time()\n",
    "batch <- readRDS('", batch_file, "')\n",
    "run_gimme_sequential(batch, '", file.path(output_dir, paste0("batch_", batch_index)), "')\n",
    "end_time <- Sys.time()\n",
    "cat(sprintf('Batch ", batch_index, " completed in %s seconds\\n', as.numeric(difftime(end_time, start_time, units = 'secs'))))\n"
  )
  
  writeLines(script_content, paste0("gimme_batch_", batch_index, ".R"))
  
  return(named_batch)
}

# Run batch scripts in parallel
run_gimme_batch_scripts <- function(num_cores, num_batches) {
  cl <- makeCluster(num_cores)
  results <- parLapply(cl, seq_len(num_batches), function(i) {
    start_time <- Sys.time()
    system(paste0("Rscript gimme_batch_", i, ".R > gimme_batch_", i, ".log 2>&1"))
    end_time <- Sys.time()
    elapsed_time <- difftime(end_time, start_time, units = "secs")
    cat(sprintf("Batch %d completed in %s seconds\n", i, as.numeric(elapsed_time)))
    return(as.numeric(elapsed_time))
  })
  stopCluster(cl)
  total_time <- sum(unlist(results))
  cat(sprintf("Total time for all GIMME batches: %s seconds\n", total_time))
}

# =========================
# Example usage for uSEM DGM
# =========================

batch_size <- 10
batches <- split_into_batches(GIMME_uSEM_DGM, batch_size)

output_dir <- "output_directory_GIMME_uSEM_DGM"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

start_number <- 1
batch_results <- lapply(seq_along(batches), function(i) {
  batch <- save_gimme_batch_script(
    batches[[i]],
    batch_index = i,
    output_dir = output_dir,
    start_number = start_number
  )
  start_number <<- start_number + length(batch)
  return(batch)
})

num_cores <- 48
start_time <- Sys.time()
run_gimme_batch_scripts(num_cores, length(batches))
end_time <- Sys.time()

total_execution_time <- difftime(end_time, start_time, units = "secs")
cat(sprintf("Total execution time for all GIMME scripts: %s seconds\n", total_execution_time))

