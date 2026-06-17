# Creating lagged variables for VAR DGM
library(dplyr)

add_lagged_vars_VAR <- function(data,
                                id_var = "ID",
                                time_var = "time") {
  
  data %>%
    mutate(
      ID_num_for_sort = as.numeric(.data[[id_var]]),
      time_num_for_sort = as.numeric(.data[[time_var]])
    ) %>%
    arrange(ID_num_for_sort, time_num_for_sort) %>%
    group_by(ID_num_for_sort) %>%
    mutate(
      loneliness_lag = dplyr::lag(loneliness_state_pmc, 1),
      socintsatisfaction_lag = dplyr::lag(socintsatisfaction_state_pmc, 1),
      responsiveness_lag = dplyr::lag(responsiveness_state_pmc, 1),
      selfdisclosure_lag = dplyr::lag(selfdisclosure_state_pmc, 1),
      otherdisclosure_lag = dplyr::lag(otherdisclosure_state_pmc, 1)
    ) %>%
    filter(!is.na(loneliness_lag)) %>%   # removes time 1 for each ID
    ungroup() %>%
    select(-ID_num_for_sort, -time_num_for_sort)
}

sim_data_list_VAR_lagged <- lapply(
  sim_data_list_VAR,
  add_lagged_vars_VAR,
  id_var = "ID",
  time_var = "time"
)

original_df <- sim_data_list_VAR[[1]] %>%
  mutate(
    ID_num = as.numeric(ID),
    time_num = as.numeric(time)
  ) %>%
  arrange(ID_num, time_num)

lagged_df <- sim_data_list_VAR_lagged[[1]] %>%
  mutate(
    ID_num = as.numeric(ID),
    time_num = as.numeric(time)
  ) %>%
  arrange(ID_num, time_num)

original_df %>%
  filter(ID_num == 1) %>%
  select(ID, time, loneliness_state_pmc) %>%
  head(5)

lagged_df %>%
  filter(ID_num == 1) %>%
  select(ID, time, loneliness_state_pmc, loneliness_lag) %>%
  head(5)

# First 5 simulations
sim_data_list_lagged_5 <- head(sim_data_list_VAR_lagged, 5)

#################################################################
# Creating lagged for uSEM DGM
library(dplyr)

add_lagged_vars_uSEM <- function(data,
                                 id_var = "id",
                                 time_var = "time") {
  
  data %>%
    mutate(
      id_num_for_sort = as.numeric(.data[[id_var]]),
      time_num_for_sort = as.numeric(.data[[time_var]])
    ) %>%
    arrange(id_num_for_sort, time_num_for_sort) %>%
    group_by(id_num_for_sort) %>%
    mutate(
      V2_lag = dplyr::lag(V2, 1),
      V3_lag = dplyr::lag(V3, 1),
      V4_lag = dplyr::lag(V4, 1),
      V5_lag = dplyr::lag(V5, 1),
      V6_lag = dplyr::lag(V6, 1)
    ) %>%
    filter(!is.na(V2_lag)) %>%   # removes time 1 for each id
    ungroup() %>%
    select(-id_num_for_sort, -time_num_for_sort)
}

sim_data_list_uSEM_lagged <- lapply(
  sim_data_list_uSEM,
  add_lagged_vars_uSEM,
  id_var = "id",
  time_var = "time"
)

original_df_uSEM <- sim_data_list_uSEM[[1]] %>%
  mutate(
    id_num = as.numeric(id),
    time_num = as.numeric(time)
  ) %>%
  arrange(id_num, time_num)

lagged_df_uSEM <- sim_data_list_uSEM_lagged[[1]] %>%
  mutate(
    id_num = as.numeric(id),
    time_num = as.numeric(time)
  ) %>%
  arrange(id_num, time_num)

original_df_uSEM %>%
  filter(id_num == 1) %>%
  select(id, time, V2) %>%
  head(5)

lagged_df_uSEM %>%
  filter(id_num == 1) %>%
  select(id, time, V2, V2_lag) %>%
  head(5)
