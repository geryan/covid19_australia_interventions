# visualise different GP formulations for new Reff model

plot_Nt <- function(Nt_simulations) {
  plot(Nt_simulations[1, , 1] ~ all_times,
       type = "n",
       ylim = c(0, max(obs_N) * 2),
       xlab = "",
       ylab = "N",
       las = 1)
  apply(Nt_simulations, 1,
        function(x) {
          lines(all_times, x,
                col = grey(0.4),
                lwd = 0.5)
        })
  points(obs_N ~ obs_times)
}

plot_rt <- function(Nt_simulations) {
  Nt_simulations <- Nt_simulations[, , 1]
  diffs <- t(apply(Nt_simulations, 1, diff))
  diffs <- cbind(diffs, NA)
  rt_simulations <- 1 + diffs / Nt_simulations
  
  plot(rt_simulations[1, ] ~ all_times,
       type = "n",
       # ylim = quantile(rt_simulations, c(0.025, 0.975), na.rm = TRUE),
       ylim = range(rt_simulations, na.rm = TRUE),
       xlab = "",
       ylab = "rt",
       las = 1)
  apply(rt_simulations, 1,
        function(x) {
          lines(all_times, x,
                col = grey(0.4),
                lwd = 0.5)
        })
  abline(h = 1, lty = 2)
  rug(obs_times)
  abline(v = max(obs_times), lty = 3, lwd = 2)
}

library(greta.gp)

# how many draws to show
sims <- 50

# 'true' model to generate fake infection count data
f_true <- function(t) exp(2 + 2 * sin(t / 5 - 2))
obs_times <- seq(1, 78, by = 3)
obs_N <- f_true(obs_times)

# need to run for all integer times (as we are using cumsum in the model), so
# this must include all observations
all_times <- seq(0, 100, by = 1)
obs_idx <- match(obs_times, all_times)

# plot(f_true, xlim = range(all_times))
# points(obs_N ~ obs_times)

f1 <- function(init, z) exp(init + z)
f2 <- function(init, z) {
  log_rt <- z 
  exp(init + cumsum(log_rt))
}
f3 <- function(init, z) {
  log_rt_diff <- z
  log_rt <- cumsum(log_rt_diff)
  exp(init + cumsum(log_rt))
}

gp_lengthscale <- 3
gp_variance <- 0.1
obs_variance <- 0.5
set.seed(2023-05-03)

# f1
init <- lognormal(0, 1)
kernel <- rbf(gp_lengthscale, gp_variance)
z <- gp(all_times, kernel, tol = .Machine$double.eps)
N <- f1(init, z)
distribution(obs_N) <- normal(N[obs_idx], obs_variance)
m <- model(init, z)
draws <- mcmc(m, chains = 10)
coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
f1_sims <- calculate(N, values = draws, nsim = sims)


# f2
init <- lognormal(0, 1)
kernel <- rbf(gp_lengthscale, gp_variance)
z <- gp(all_times, kernel, tol = .Machine$double.eps)
N <- f2(init, z)
distribution(obs_N) <- normal(N[obs_idx], obs_variance)
m <- model(init, z)
draws <- mcmc(m, chains = 10)
coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
f2_sims <- calculate(N, values = draws, nsim = sims)

# f3
init <- lognormal(0, 1)
kernel <- rbf(gp_lengthscale, gp_variance)
z <- gp(all_times, kernel, tol = .Machine$double.eps)
N <- f3(init, z)
distribution(obs_N) <- normal(N[obs_idx], obs_variance)
m <- model(init, z)

# this one needs a little help to initialise
z_raw <- attr(z, "gp_info")$v
n_chains <- 20
inits <- replicate(n_chains,
                   initials(
                     z_raw = rep(0, length(z_raw))
                   ),
                   simplify = FALSE)
draws <- mcmc(m,
              warmup = 3000,
              chains = n_chains,
              initial_values = inits)

coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
f3_sims <- calculate(N, values = draws, nsim = sims)

# plot modelled infection timeseries and growth rates

# infection timeseries
par(mfrow = c(3, 1))
plot_Nt(f1_sims$N)
title("N_t = N_0 * exp(GP(t))")
# this is a 0-mean GP plus an intercept, so it tends back to the intercept (~~
# average case counts)
plot_Nt(f2_sims$N)
title("N_t = N_0 * cumprod(r_t); r_t = exp(GP(t))")
# this is a cumulative 0-mean GP on log(N_t), so it sticks at the same value
# (recent case counts), and case counts stay flat at their recent level, give or
# take stochasticity
plot_Nt(f3_sims$N)
title("N_t = N_0 * cumprod(r_t); r_t = exp(cumsum(GP(t)))")
# this is a cumulative 0-mean GP on r_t so r_t sticks at the same value
# (multiplying by 1 each time), and N_t grows exponentially,  case counts keep
# trending in the same direction

# growth rates
par(mfrow = c(3, 1))
plot_rt(f1_sims$N)
title("N_t = N_0 * exp(GP(t))")
plot_rt(f2_sims$N)
title("N_t = N_0 * cumprod(r_t); r_t = exp(GP(t))")
plot_rt(f3_sims$N)
title("N_t = N_0 * cumprod(r_t); r_t = exp(cumsum(GP(t)))")
