## Simulation test Stan models



## 1_hbm_ar1 -----------------------------------------------
d1 <- sim_data(n_series = 20, n_years = 100, phi = 0.4)
f1 <- rstan::stan(file = "./stan/1_hbm_ar1.stan",
                  data = d1$stan,
                  pars = c("alpha", "beta", "sigma", "phi",
                           "mu_alpha", "sigma_alpha", "sigmaNC"),
                  warmup = 1000,
                  iter = 2000,
                  cores = 2,
                  chains = 2,
                  control = list(adapt_delta = 0.80))
f1
neff_lowest(f1)

sim_compare(f1, d1, "alpha")
sim_compare(f1, d1, "beta")
sim_compare(f1, d1, "sigma")  ## sigma is hard to estimate
sim_compare(f1, d1, "phi")
sim_compare(f1, d1, "mu_alpha")
sim_compare(f1, d1, "sigma_alpha")



## 2_hbm_ar1_x2stat ----------------------------------------
d2 <- sim_data(n_series = 20, phi = 0.4, x2 = "stationary")
f2 <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                  data = d2$stan,
                  pars = c("alpha", "beta", "sigma", "phi",
                           "mu_alpha", "sigma_alpha",
                           "gamma", "mu_gamma", "sigma_gamma"),
                  warmup = 1000,
                  iter = 2000,
                  cores = 2,
                  chains = 2,
                  control = list(adapt_delta = 0.80))

f2
neff_lowest(f2)

sim_compare(f2, d2, "alpha")
sim_compare(f2, d2, "beta")
sim_compare(f2, d2, "sigma")
sim_compare(f2, d2, "phi")
sim_compare(f2, d2, "mu_alpha")
sim_compare(f2, d2, "sigma_alpha")
sim_compare(f2, d2, "gamma")
sim_compare(f2, d2, "mu_gamma")
sim_compare(f2, d2, "sigma_gamma")



## 3_hbm_ar1_x2era3 ----------------------------------------
d3 <- sim_data(n_series = 20, phi = 0.4, x2 = "era")
d3$data
f3 <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                  data = d3$stan,
                  pars = c("alpha", "beta", "sigma", "phi",
                           "mu_alpha", "sigma_alpha",
                           "gamma1", "gamma2", "gamma3",
                           "mu_gamma1", "mu_gamma2", "mu_gamma3",
                           "sigma_gamma"),
                  warmup = 1000,
                  iter = 2000,
                  cores = 2,
                  chains = 2,
                  control = list(adapt_delta = 0.80))

summary(f3)[[1]]
neff_lowest(f3)

sim_compare(f3, d3, "gamma1")
sim_compare(f3, d3, "gamma2")
sim_compare(f3, d3, "gamma3")

sim_compare(f3, d3, "mu_gamma1")
sim_compare(f3, d3, "mu_gamma2")
sim_compare(f3, d3, "mu_gamma3")

sim_compare(f3, d3, "sigma_gamma")

sim_compare(f3, d3, "sigma")
sim_compare(f3, d3, "alpha")
sim_compare(f3, d3, "beta")
sim_compare(f3, d3, "phi")
sim_compare(f3, d3, "mu_alpha")
sim_compare(f3, d3, "sigma_alpha")



## 5_hbm_ar1_x2nonstat_sigma_diff --------------------------
d5 <- sim_data(n_series = 20, n_years = 50, x2 = "nonstationary",
               x2_step = TRUE,
               phi = 0.3)
f5 <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                  data = d5$stan,
                  pars = c("alpha", "beta", "sigma", "phi",
                           "mu_alpha", "sigma_alpha",
                           "gamma", "sigma_gamma", "signal_noise"),
                  warmup = 500,
                  iter = 1200,
                  cores = 1,
                  chains = 1,
                  control = list(adapt_delta = 0.80))

# mean(d5$pars$gamma)
# summary(f5, pars = c("sigma_gamma"))[[1]]

# sort(summary(f5)$summary[ , "n_eff"])[1:4]
# pars <- names(sort(summary(f5)$summary[ , "n_eff"]))[1:4]
# pairs(f5, pars = pars)

# sim_compare(f5, d5, "sigma_gamma")

g <- ggplot(d5$data) +
    geom_point(aes(x = x2, y = y)) +
    theme_sleek() +
    facet_wrap( ~ series, scale = "free")
print(g)


mc.sum <- summary(f5)$summary
g1 <- mc.sum[grep("^gamma", rownames(mc.sum)), "mean"]
dat <- d5$data
dat$gamma_hat <- g1
g <- ggplot(dat) +
    geom_line(aes(x = year, y = gamma_true)) +
    geom_line(aes(x = year, y = gamma_hat), color = "red3") +
    theme_sleek() +
    facet_wrap( ~ series)
print(g)
