# functions & packages
##### 
library(tidyverse)

make_vaccine_combinations <- function(
  vaccineArgs = list(AZ = 0:2, Pf = 0:2), # list of vaccine types and sequence of possible doses from 0,
  max_doses = 2
){
  
  result <- do.call(
    expand_grid,
    vaccineArgs
  ) %>% 
    filter(
      rowSums(.) > 0,
      rowSums(.) <= max_doses
    )
  
  return(result)
  
}

# given a dataframe with each column givin a different types of baseline
# immunity, and the values giving the number of those immune events, return a
# vector with number of elements equal to the number of rows, describing that
# level of immunity
get_immunity_names <- function(immunity_combinations) {
  
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
  rownames(result) <- colnames(result) <- get_immunity_names(immunity_combinations) 
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

make_immunity_combinations <- function(
  vaccine_names,
  nonimmune_name = "S",
  infected_name = "I"
) {
  
  c(
    nonimmune_name,
    infected_name,
    vaccine_names,
    paste0(infected_name, vaccine_names)
  )
  
}

get_immunity_type_expansion <- function(
  immunity_types,
  n_waning_states,
  nonimmune_name = NULL,
  immunity_types_with_waning = NULL
) {
  
  if(is.null(immunity_types_with_waning)){
    if(is.null(nonimmune_name)) {
      stop("You must provide at least one of nonimmune_name or immunity_types")
    }
    immunity_types_with_waning <- setdiff(immunity_types, nonimmune_name)
  }
  
  has_waning <- immunity_types %in% immunity_types_with_waning
  expansion <- ifelse(has_waning, n_waning_states, 1)
  
  return(expansion)
}


get_expanded_state_names <- function(
  immunity_types,
  expansion
) {
  
  has_waning <- expansion > 1
  
  # create the state names
  immunity_types_expanded <- rep(immunity_types, expansion)
  waning_expanded <- unlist(lapply(expansion, seq_len)) - 1
  
  # remove 0 for the non-waning immunity types, and coerce to character
  mask_waning_count <- rep(!has_waning, expansion)
  waning_expanded <- paste0("_", waning_expanded)
  waning_expanded[mask_waning_count] <- ""
  
  # combine immunity type with waning evel, and return
  states <- paste0(immunity_types_expanded, waning_expanded)
  
  states
  
}



expand_inheritance_matrix <- function(
  inheritance_matrix,
  n_waning_states = 4,
  nonimmune_name = NULL,
  immunity_types_with_waning = NULL
) {
  
  
  if(is.null(immunity_types_with_waning)){
    if(is.null(nonimmune_name)) {
      stop("You must provide at least one of nonimmune_name or immunity_types")
    }
    immunity_types_with_waning <- setdiff(immunity_types, nonimmune_name)
  }
  
  
  immunity_types <- colnames(inheritance_matrix)
  n_immunity_types <- length(immunity_types)
  
  # expand out the immunity states to levels of waning
  expansion <- get_immunity_type_expansion(
    immunity_types = immunity_types,
    immunity_types_with_waning = immunity_types_with_waning,
    n_waning_states = n_waning_states
  )
  
  # create the state names
  expanded_state_names <- get_expanded_state_names(immunity_types, expansion)
  
  # loop through the original inheritance matrix, making submatrices for the
  # expanded components
  columns <- list()
  for (col in seq_len(n_immunity_types)) {
    rows <- list()
    for (row in seq_len(n_immunity_types)) {
      rows[[row]] <- make_rowwise_submatrix(
        inheritance_matrix[row, col],
        expansion_rows = expansion[row],
        expansion_cols = expansion[col]
      )
    }
    columns[[col]] <- rows
  }
  
  # compress this nested list into a matrix
  columns_bound <- lapply(columns, function(x) {do.call(rbind, x)})
  matrix <- do.call(cbind, columns_bound)
  rownames(matrix) <- colnames(matrix) <- expanded_state_names
  matrix
  
}


#####
# script
######


nonimmune_name <- "S"
infected_name <- "I"
n_population <- 1e6
n_waning_states <- 3

# get all possible vaccination combinations, with AZ only available up to dose
# 2, and a max of 4 doses across all types
vaccine_types <- list(
  AZ = 0:2,
  Pf = 0:2#4,
  #Mod = 0:4
)

vaccine_combinations <- make_vaccine_combinations(
  vaccineArgs = vaccine_types,
  max_doses = 2
)

vaccine_names <- get_immunity_names(vaccine_combinations)


immunity_types <- make_immunity_combinations(
  vaccine_names = vaccine_names,
  nonimmune_name = nonimmune_name,
  infected_name = infected_name
)




expansion <- get_immunity_type_expansion(
  immunity_types = immunity_types,
  nonimmune_name = nonimmune_name,
  n_waning_states = n_waning_states
)

states <- get_expanded_state_names(
  immunity_types = immunity_types,
  expansion = expansion
)



# set up fake population
n_states <- length(states)

n_population <- 1e6

# proportion of the population in each state
current_state_prop <- c(50, runif(n_states - 1))
current_state_prop <- current_state_prop / sum(current_state_prop)

# relative probability of being infected: e.g. RR of infection for each immunity type and level of waning, from VEs
# should remain the same at each iteration since it's just based on the immunity model
relative_probability_I <- runif(n_states)

# relative probability of being vaccinated: e.g. from number of new vaccinations
# and which groups is vaccinated. For all I-only states, should be the same as
# susceptibles (the unvaccinated), for IV states should be the same as for V
# states (with same waning) - should change at each iteration based on the AIR
# data
relative_probability_V <- runif(n_states)

# total number of new infections and total number of new vaccinations
n_new_vaccinations <- 2e4
n_new_infections <- 1e4


# number of vaccinations for each of these combos
# vaccine_counts <- bind_rows(
#   name = vaccine_names,
#   count = rowSums(vaccine_combinations)
# )


add_si <- function(
  vaccine_combinations,
  nonimmune_name,
  infected_name
){
  
  nrows <- nrow(vaccine_combinations)
  extra_rows <- vaccine_combinations[1:2,]
  extra_rows[] <- 0
  
  expanded_combinations <- bind_rows(
    vaccine_combinations,
    extra_rows
  )
  
  
}

inheritance_vaccines <- find_inheritance(vaccine_combinations)
image(inheritance_vaccines)
inheritance_vaccines


