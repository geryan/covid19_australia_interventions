# mock up new, simmpler internals for greta neut model that enable
# immunity/variant interactions and other constraints

# define a regualrised horseshoe prior on a vector of parameters:
# https://forum.greta-stats.org/t/is-the-regularised-horseshoe-prior-example-correct/264/4
regularized_horseshoe <- function (tau = 1,  c = 1, dim = NULL) {
  stopifnot(c > 0)
  stopifnot(tau > 0)
  lambda <- cauchy(0, 1, truncation = c(0, Inf), dim = dim)
  sd <- tau * c * lambda / sqrt(c^2 + tau^2 * lambda^2)
  normal(0, sd, dim = dim)
}

# given a character vector (or matrix or dataframe with multiple character
# columns) of values, and a matching vaector (or matrix or dataframe) of
# 'options', find the numeric index of each of the values to the options. This
# is essentially just match(), but can match calues on mulltiple columns
get_index <- function(
  values,
  options
) {
  
  if (NCOL(values) > 1) {
    
    if (NCOL(options) != NCOL(values)) {
      stop("values and options have different numbers of columns")
    }
    
    # paste together the columns into a unique single code
    values <- apply(values, 1, paste, collapse = "_X_")
    options <- apply(options, 1, paste, collapse = "_X_")
    
  }
  
  match(values, options)
  
}

# construct a matrix encoding the possible state inheritance structure of the
# transition matrix, from a tibble of the inter-state transitions that are
# possible (ignoring self-transitions, which are added automatically)
manual_inheritance_matrix <- function(possible_transitions,
                                      self_inherit = FALSE,
                                      states = NULL) {
  if (is.null(states)) {
    states <- unique(c(possible_transitions$from, possible_transitions$to))
  }
  
  possible_transitions_index <- as.matrix(possible_transitions[, c("to", "from")])
  inheritance <- empty_matrix(states)
  inheritance[possible_transitions_index] <- 1
  if(self_inherit) {
    diag(inheritance) <- 1
  }
  inheritance
}



# define the names of the different types of immunity (need not all be observed)
immunity_type_name <- c("wt_inf",
                        "vaccine_1",
                        "vaccine_2",
                        "delta_inf",
                        "delta_inf_+_vaccine_1")
# enable structure in these, so that combinations must not decrease immunity?

# define the names of the different variants against which people may be immune
# (need not all be observed)
variant_name <- c("wt", "delta", "omicron")

# compute all possible combinations of immunity type and variant
immunity_variant_interaction_name <- expand.grid(
  immunity_type = immunity_type_name,
  variant = variant_name
)

# mock up some observed (relative, log10) neut values
n_neut_obs <- 100
neut_obs <- data.frame(
  immunity_type = sample(immunity_type_name, n_neut_obs, replace = TRUE),
  variant = sample(variant_name, n_neut_obs, replace = TRUE),
  neuts = rnorm(n_neut_obs)
)

# create greta arrays for the variables corresponding to each of these
n_immunity_types <- length(immunity_type_name)
n_variants <- length(variant_name)
n_immunity_variant_interactions <- n_immunity_types * n_variants

library(greta)

# note that we don't want to use priors for these any more, as we will provide
# neut study info as data instead

# minimally-informative priors on these
immunity_type_values <- normal(0, 10, dim = n_immunity_types)
# enable structure in these, so that combinations must not decrease immunity?

variant_values <- normal(0, 10, dim = n_variants)

# shrinkage prior (horseshoe) on the interaction coefficients
immunity_variant_interaction_values <- regularized_horseshoe(
  dim = n_immunity_variant_interactions
)

# observation error on neut values (constant on log10 scale implies it scales
# linearly with the absolute neut values)
sd_obs <- normal(0, 1, truncation = c(0, Inf))

# pull out expected neut value for observed neuts
variant_idx <- get_index(neut_obs$variant, variant_name)
immunity_type_idx <- get_index(neut_obs$immunity_type, immunity_type_name)
immunity_variant_interaction_idx <- get_index(neut_obs[, c("immunity_type", "variant")], immunity_variant_interaction_name)

expected_neut_obs <- variant_values[variant_idx] +
  immunity_type_values[immunity_type_idx] +
  immunity_variant_interaction_values[immunity_variant_interaction_idx]

# define likelihood over observed neut data:
distribution(neut_obs$neuts) <- normal(expected_neut_obs, sd_obs)

m <- model(immunity_type_values, variant_values, immunity_variant_interaction_values)
draws <- mcmc(m)
plot(draws)





# to do:

# 1. enable strictly increasing immunity from hybrid immune conditions
# (vaccine+infection cannot give lower neuts than either vaccine or infection; 2
# doses of vaccine cannot give lower neuts than 1 dose, etc.)

# 2. constrain parameters to be relative to wt infection (convalescent serum),
# with log10 neut parameter for wt fixed at 0


# 1.1
# how to constrain combinations of immunity to be increasing?
# can define positive offsets (added to highest of components)
a <- normal(0, 1)
b <- normal(0, 1)
a_plus_b_extra <- normal(0, 1, truncation = c(0, Inf))
a_plus_b <- max(a, b) + a_plus_b_extra

# check this will work in a model
obs <- rnorm(10)
distribution(obs) <- normal(a_plus_b, 1)
m <- model(a, b, a_plus_b)
draws <- mcmc(m)
plot(draws)

# check it's always higher
do.call(cbind, calculate(a, b, a_plus_b, nsim = 10))

# 1. define the inheritance structure for immunity types


# given a dataframe with each column givin a different types of baseline
# immunity, and the values giving the number of those immune events, return a
# vector with number of elements equal to the number of rows, describing that
# level of immunity
get_hybrid_immunity_names <- function(immunity_combinations) {
  
  make_name <- function(counts, names) {
    exists <- counts > 0
    paste(paste0(names, counts)[exists], collapse = "_")
  }
  
  apply(immunity_combinations,
        1,
        make_name,
        colnames(immunity_combinations))
}

# need to define inheritance structure from this (inheriting from another
# combination if that has exactly one fewer of one type of vaccine)
find_inheritance <- function(immunity_combinations) {
  n <- nrow(immunity_combinations)
  result <- matrix(0, n, n)
  rownames(result) <- colnames(result) <- get_hybrid_immunity_names(immunity_combinations) 
  for (j in 1:n) {
    for (i in j:n) {
      result[i, j] <- inherits_from(immunity_combinations[i, ],
                                    immunity_combinations[j, ])
    }
  }
  result
  
}

inherits_from <- function(parent, child) {
  # is the child the same as the parent, but with exactly 1 less in one column?
  # ie. can you get from the child to the parent by adding one vaccine dose?
  diff <- parent - child
  max(diff) == 1 & min(diff) == 0 & sum(diff) == 1
}


# for a vector of immunity type names, get a list of the same length, each
# element giving a character vector of the immunity types they inherit from
find_parent_names <- function(names, inheritance) {
  vectors <- apply(inheritance[names, ] == 1, 1, which)
  lapply(vectors, names)
}

# given a character vector of names of immunity types, a vector of values for
# those, and the corrspondng vector of names, return the maximum of al the
# immunity values in that set
get_max_neuts_from_set <- function(names, values, value_names) {
  max(values[get_index(names, value_names)])
}

# given a vector of names of second dose immunity types, an inheritance
# structure, and the values and names of the first dose immunity types, return a
# vector of the minimum neuts (the max of the parents)
get_min_neuts <- function(second_dose_names, inheritance, first_dose_values, first_dose_names) {
  
  # for each of these second dose immunities, find the parents (first dose immunities)
  parent_names <- find_parent_names(second_dose_names,
                                    inheritance)
  
  # for each of these get the minimum value of the hybrid immunity (max of the
  # parents)
  minima_list <- lapply(parent_names,
                        get_max_neuts_from_set,
                        first_dose_values,
                        first_dose_names)
  
  # combine these into a vector
  do.call(c, minima_list)
  
}

# start with vaccines only

library(tidyverse)

# get all possible vaccination combinations, with AZ only available up to dose
# 2, and a max of 4 doses across all types
vaccine_combinations <- expand_grid(
  AZ = 0:2,
  Pf = 0:4,
  Mod = 0:4
) %>%
  filter(
    rowSums(.) > 0,
    rowSums(.) <= 4
  )

# number of vaccinations for each of these combos
vaccine_counts <- bind_rows(
  name = get_hybrid_immunity_names(vaccine_combinations),
  count = rowSums(vaccine_combinations)
)



inheritance <- find_inheritance(vaccine_combinations)
image(inheritance)
inheritance

# to do these in the right order, need to loop through in order of total number of
# vaccines

vaccine_type_names <- vaccine_counts$name

# first define variables for the first doses
first_dose_names <- vaccine_counts %>%
  filter(count == 1) %>%
  pull(name)

first_dose_values <- normal(0, 1, dim = length(first_dose_names))

# now do second doses
second_dose_names <- vaccine_counts %>%
  filter(count == 2) %>%
  pull(name)

# create an 'extra' neut value for each of these
second_dose_extra_values <- normal(0, 1, dim = length(second_dose_names), truncation = c(0, Inf))
second_dose_min_values <- get_min_neuts(second_dose_names,
                                        inheritance,
                                        first_dose_values,
                                        first_dose_names)

second_dose_values <- second_dose_min_values + second_dose_extra_values

# then rinse and repeat for 2nd, third etc. (need to combine 1st and second into
# a vector to get the parents for third)

# Then once we have these, create another inheritance matrix for the
# vaccine-infection hybrid immunities, and follow the same process (treating all
# the vaccines as a single immune event to avoid repeating this)

# ie., create all possible combinations of the different vaccination levels
# (AZ1, AZ1_Pf1, etc., with the values created above) and different numbers of
# infections with each variant (baseline values for infection with each variant
# alone just as N(0, 10), like the first vaccinations), and then create the
# inheritance matrix, and define the values as above


# 2. constrain parameters to be relative to wt infection (convalescent serum),
# with log10 neut parameter for wt fixed at 0

# just put wt infection immunity at the start, and force it to be 1

infection_immunity_names <- c("wt", "alpha", "beta", "delta", "omicron_ba1", "omicron_ba2")
n_infection_immunity_names <- length(infection_immunity_names)
infection_immunity_values <- c(ones(1), normal(0, 10, dim = n_infection_immunity_names - 1))

