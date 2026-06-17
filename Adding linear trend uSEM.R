library(dplyr)

#################################################
#### Add linear trend to uSEM baseline by ID ####
#################################################

add_trend_to_variable_by_id <- function(data, 
                                        variable_name, 
                                        trend_size_range,
                                        id_var = "id",
                                        time_var = "time") {
  
  # Check variable exists
  if (!variable_name %in% names(data)) {
    stop(paste("Variable", variable_name, "not found in dataset."))
  }
  
  # Check ID variable exists
  if (!id_var %in% names(data)) {
    stop(paste("ID variable", id_var, "not found in dataset."))
  }
  
  # Check time variable exists
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
                                           id_var = "id",
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

###############################
#### Small trend condition ####
###############################

sim_data_list_uSEM_trendsmall <- add_trend_to_simulations_by_id(
  simulation_list = sim_data_list_uSEM,
  variable_name = "V1",
  trend_size_range = c(0.1, 0.3),
  id_var = "id",
  time_var = "time"
)


###############################
#### Medium trend condition ####
###############################

sim_data_list_uSEM_trendmedium <- add_trend_to_simulations_by_id(
  simulation_list = sim_data_list_uSEM,
  variable_name = "V1",
  trend_size_range = c(0.4, 0.6),
  id_var = "id",
  time_var = "time"
)


###############################
#### Large trend condition ####
###############################

sim_data_list_uSEM_trendlarge <- add_trend_to_simulations_by_id(
  simulation_list = sim_data_list_uSEM,
  variable_name = "V1",
  trend_size_range = c(0.7, 0.9),
  id_var = "id",
  time_var = "time"
)


library(ggplot2)

baseline_df <- sim_data_list_uSEM[[1]]
trended_df  <- sim_data_list_uSEM_trendlarge[[1]]

baseline_df <- baseline_df %>%
  mutate(id_num = as.numeric(id)) %>%
  arrange(id_num, time) %>%
  select(-id_num)

trended_df <- trended_df %>%
  mutate(id_num = as.numeric(id)) %>%
  arrange(id_num, time) %>%
  select(-id_num)

trend_check_df <- baseline_df %>%
  select(id, time, baseline_V1 = V1) %>%
  left_join(
    trended_df %>%
      select(id, time, trended_V1 = V1),
    by = c("id", "time")
  ) %>%
  mutate(
    added_trend = trended_V1 - baseline_V1
  )

trend_check_df %>%
  filter(id %in% c(1, 2, 3, 4, 5)) %>%
  ggplot(aes(x = time, y = added_trend)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ id) +
  theme_minimal() +
  labs(
    title = "Added linear trend in V1: uSEM DGM",
    x = "Time",
    y = "Trended V1 - baseline V1"
  )


trend_summary <- trend_check_df %>%
  group_by(id) %>%
  summarise(
    first_added = added_trend[which.min(time)],
    last_added  = added_trend[which.max(time)],
    total_added = last_added - first_added,
    min_added   = min(added_trend, na.rm = TRUE),
    max_added   = max(added_trend, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(id_num = as.numeric(id)) %>%
  arrange(id_num)

head(trend_summary, 10)
summary(trend_summary$total_added)

saveRDS(sim_data_list_uSEM_trendsmall,  file = "sim_data_list_uSEM_trendsmall.rds")
saveRDS(sim_data_list_uSEM_trendmedium, file = "sim_data_list_uSEM_trendmedium.rds")
saveRDS(sim_data_list_uSEM_trendlarge,  file = "sim_data_list_uSEM_trendlarge.rds")
