# fit a Bayesian model-based estimate of R_effective over time, quantifying the
# impacts of both quarantine and physical distancing measures.

source("R/lib.R")

set.seed(2020-04-29)
source("R/functions.R")

# sync up the case data
sync_nndss()

# prepare data for Reff modelling
data <- reff_model_data()

data$dates$linelist

# save the key dates for Freya and David to read in, and tabulated local cases
# data for the Robs
 write_reff_key_dates(data)
 write_local_cases(data)

# format and write out any new linelists to the past_cases folder for Rob H
#update_past_cases()

# define the model (and greta arrays) for Reff, and sample until convergence
fitted_model <- fit_reff_model(data)

# save the fitted model object
saveRDS(fitted_model, "outputs/vl_full_ext/fitted_reff_model.RDS")
# # fitted_model <- readRDS("outputs/fitted_reff_model.RDS")
# 
# # output Reff trajectory draws for Rob M
 write_reff_sims(fitted_model, dir = "outputs/vl_full_ext/projection")
# 
# # visual checks of model fit
# #plot_reff_checks(fitted_model)
# 
timeseries <- readRDS("outputs/vaccine_timeseries.RDS")
# 
# # do plots for main period
reff_plotting(fitted_model, dir = "outputs/vl_full_ext")
# 
# fitted_model_ext <- fitted_model
# fitted_model_ext$data$dates$latest_mobility <- mobility_data$date %>% max
# 
# fitted_model_ext$data$dates$mobility <- seq(fitted_model$data$dates$mobility %>% min, mobility_data$date %>% max, by = 1)

# and for projected part
reff_plotting(fitted_model,
              dir = "outputs/vl_full_ext/projection",
              max_date = fitted_model$data$dates$latest_project,
              mobility_extrapolation_rectangle = FALSE,
              #projection_date = fitted_model$data$dates$latest_mobility
              projection_date = as.Date("2021-06-07"))


# produce simulations where proportion VOC is zero throughout period for reporting
# in common operating picture table
#simulate_wild_type()
