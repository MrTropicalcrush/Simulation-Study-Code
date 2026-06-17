###############################
#### Adding linear trend by ID ####
###############################
library(dplyr)

add_trend_to_variable_by_id <- function(data, 
                                        variable_name, 
                                        trend_size_range,
                                        id_var = "ID",
                                        time_var = "time") {
  
  if (!variable_name %in% names(data)) {
    stop(paste("Variable", variable_name, "not found in dataset."))
  }
  
  if (!id_var %in% names(data)) {
    stop(paste("ID variable", id_var, "not found in dataset."))
  }
  
  if (!time_var %in% names(data)) {
    stop(paste("Time variable", time_var, "not found in dataset."))
  }
  
  data %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[time_var]], .by_group = TRUE) %>%
    mutate(
      linear_trend = seq(
        from = 0,
        to = runif(1, trend_size_range[1], trend_size_range[2]),
        length.out = n()
      ),
      "{variable_name}" := .data[[variable_name]] + linear_trend
    ) %>%
    select(-linear_trend) %>%
    ungroup()
}


add_trend_to_simulations_by_id <- function(simulation_list, 
                                           variable_name, 
                                           trend_size_range,
                                           id_var = "ID",
                                           time_var = "time") {
  
  lapply(simulation_list, function(sim_data) {
    add_trend_to_variable_by_id(
      data = sim_data,
      variable_name = variable_name,
      trend_size_range = trend_size_range,
      id_var = id_var,
      time_var = time_var
    )
  })
}

#### Small trend condition ####

sim_data_list_VAR_trendsmall <- add_trend_to_simulations_by_id(
  simulation_list = sim_data_list_VAR_baseline,
  variable_name = "depressedmood_state",
  trend_size_range = c(0.1, 0.3),
  id_var = "ID",
  time_var = "time"
)



#### Medium trend condition ####

sim_data_list_VAR_trendmedium <- add_trend_to_simulations_by_id(
  simulation_list = sim_data_list_VAR_baseline,
  variable_name = "depressedmood_state",
  trend_size_range = c(0.4, 0.6),
  id_var = "ID",
  time_var = "time"
)


#### Large trend condition ####

sim_data_list_VAR_trendlarge <- add_trend_to_simulations_by_id(
  simulation_list = sim_data_list_VAR_baseline,
  variable_name = "depressedmood_state",
  trend_size_range = c(0.7, 0.9),
  id_var = "ID",
  time_var = "time"
)

library(dplyr)

sort_simulation_by_id_time <- function(simulation_list,
                                       id_var = "ID",
                                       time_var = "time") {
  
  lapply(simulation_list, function(df) {
    
    df %>%
      mutate(
        ID_numeric_for_sorting = as.numeric(.data[[id_var]])
      ) %>%
      arrange(ID_numeric_for_sorting, .data[[time_var]]) %>%
      select(-ID_numeric_for_sorting)
  })
}

sim_data_list_VAR_trendsmall <- sort_simulation_by_id_time(
  sim_data_list_VAR_trendsmall,
  id_var = "ID",
  time_var = "time"
)

sim_data_list_VAR_trendmedium <- sort_simulation_by_id_time(
  sim_data_list_VAR_trendmedium,
  id_var = "ID",
  time_var = "time"
)

sim_data_list_VAR_trendlarge <- sort_simulation_by_id_time(
  sim_data_list_VAR_trendlarge,
  id_var = "ID",
  time_var = "time"
)

saveRDS(sim_data_list_VAR_trendsmall,  file = "sim_data_list_VAR_trendsmall.rds")
saveRDS(sim_data_list_VAR_trendmedium, file = "sim_data_list_VAR_trendmedium.rds")
saveRDS(sim_data_list_VAR_trendlarge,  file = "sim_data_list_VAR_trendlarge.rds")



df <- sim_data_list_VAR_trendlarge[[1]]

trend_check <- df %>%
  group_by(ID) %>%
  summarise(
    n_time = n(),
    slope = coef(lm(depressedmood_state ~ time))[2],
    first_value = depressedmood_state[which.min(time)],
    last_value = depressedmood_state[which.max(time)],
    change = last_value - first_value,
    mean_dep = mean(depressedmood_state, na.rm = TRUE),
    .groups = "drop"
  )

head(trend_check)
summary(trend_check$slope)
summary(trend_check$change)

################################
#### Adding nonlinear trend ####
################################
############################################
#### Add nonlinear trend to VAR DGM data ####
############################################
library(dplyr)

add_nonlinear_trend_within_person_VAR <- function(data,
                                                  id_var = "ID",
                                                  time_var = "time",
                                                  variable_name = "depressedmood_state",
                                                  trend_strength_range = c(0.7, 0.9)) {
  
  if (!variable_name %in% colnames(data)) {
    stop(paste("Variable", variable_name, "not found in dataset."))
  }
  
  data %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[time_var]], .by_group = TRUE) %>%
    mutate(
      trend_strength = runif(1, trend_strength_range[1], trend_strength_range[2]),
      time_seq = seq(0, 1, length.out = n()),
      nonlinear_trend = (time_seq^2) * trend_strength * (n() - 1),
      "{variable_name}" := .data[[variable_name]] + nonlinear_trend
    ) %>%
    ungroup() %>%
    select(-trend_strength, -time_seq, -nonlinear_trend)
}

VAR_nonlineartrend_large <- lapply(
  sim_data_list_VAR,
  function(dat) {
    add_nonlinear_trend_within_person_VAR(
      data = dat,
      id_var = "ID",
      time_var = "time",
      variable_name = "depressedmood_state",
      trend_strength_range = c(0.7, 0.9)
    )
  }
)

before <- sim_data_list_VAR[[1]]
after  <- VAR_nonlineartrend_large[[1]]

before_sorted <- before %>%
  mutate(ID = as.numeric(ID)) %>%
  arrange(ID, time)

after_sorted <- after %>%
  mutate(ID = as.numeric(ID)) %>%
  arrange(ID, time)

check_trend_VAR <- after_sorted
check_trend_VAR$diff <- after_sorted$depressedmood_state - before_sorted$depressedmood_state

check_trend_VAR %>%
  group_by(ID) %>%
  summarise(
    min_diff = min(diff),
    max_diff = max(diff),
    mean_diff = mean(diff),
    .groups = "drop"
  ) %>%
  arrange(ID)

test <- VAR_nonlineartrend_large[[1]]

VAR_nonlineartrend_large <- lapply(VAR_nonlineartrend_large, function(df) {
  df %>%
    mutate(ID = as.numeric(ID)) %>%
    arrange(ID, time)
})

saveRDS(VAR_nonlineartrend_large, file = "VAR_nonlineartrend_large.RDS")
##################################
#### Nonlinear trend uSEM DGM ####
##################################
library(dplyr)

add_nonlinear_trend_within_person_uSEM <- function(data,
                                                   id_var = "id",
                                                   time_var = "time",
                                                   variable_name = "V1",
                                                   trend_strength_range = c(0.7, 0.9)) {
  
  if (!variable_name %in% colnames(data)) {
    stop(paste("Variable", variable_name, "not found in dataset."))
  }
  
  data %>%
    group_by(.data[[id_var]]) %>%
    arrange(.data[[time_var]], .by_group = TRUE) %>%
    mutate(
      trend_strength = runif(1, trend_strength_range[1], trend_strength_range[2]),
      time_seq = seq(0, 1, length.out = n()),
      nonlinear_trend = (time_seq^2) * trend_strength * (n() - 1),
      "{variable_name}" := .data[[variable_name]] + nonlinear_trend
    ) %>%
    ungroup() %>%
    select(-trend_strength, -time_seq, -nonlinear_trend)
}

uSEM_nonlineartrend_large <- lapply(
  sim_data_list_uSEM,
  function(dat) {
    add_nonlinear_trend_within_person_uSEM(
      data = dat,
      id_var = "id",
      time_var = "time",
      variable_name = "V1",
      trend_strength_range = c(0.7, 0.9)
    )
  }
)

before <- sim_data_list_uSEM[[1]]
after  <- uSEM_nonlineartrend_large[[1]]

check_trend_uSEM <- after
check_trend_uSEM$diff <- after$V1 - before$V1

check_trend_uSEM %>%
  group_by(id) %>%
  summarise(
    min_diff = min(diff),
    max_diff = max(diff),
    mean_diff = mean(diff),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(id))

saveRDS(uSEM_nonlineartrend_large, file = "uSEM_nonlineartrend_large.RDS")
