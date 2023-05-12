## Fit hierarchical Bayesian models

if(!dir.exists("./figures/hbm_fit/"))
    dir.create("./figures/hbm_fit/")

if(!dir.exists("./output/hbm_fit/"))
    dir.create("./output/hbm_fit/")

if(!dir.exists("./output/hbm_prior_pred/"))
    dir.create("./output/hbm_prior_pred/")

if(!dir.exists("./output/hbm_loo/"))
    dir.create("./output/hbm_loo/")


## Get data for Stan ---------------------------------------
stan.dat.sst  <- stan_data(sock.covar, "sst_anom_stnd",
                           breakpoint1 = 1977,
                           breakpoint2 = 1989,
                           scale.x1 = TRUE)
stan.dat.npgo <- stan_data(sock.covar, "npgo_winter_stnd",
                           breakpoint1 = 1977,
                           breakpoint2 = 1989,
                           scale.x1 = TRUE)
stan.dat.pdo  <- stan_data(sock.covar, "pdo_winter_stnd",
                           breakpoint1 = 1977,
                           breakpoint2 = 1989,
                           scale.x1 = TRUE)

save(stan.dat.sst, file = "./output/stan_dat_sst.RData")
save(stan.dat.npgo, file = "./output/stan_dat_npgo.RData")
save(stan.dat.pdo, file = "./output/stan_dat_pdo.RData")


stan.dat.sst.s1  <- stan_data(sock.covar, "sst_anom_stnd",
                              breakpoint1 = 1977,
                              breakpoint2 = 1989,
                              prior.sigma.gamma = c(df = 1, mu = 0, sd = 1),
                              scale.x1 = TRUE)
stan.dat.npgo.s1 <- stan_data(sock.covar, "npgo_winter_stnd",
                              breakpoint1 = 1977,
                              breakpoint2 = 1989,
                              prior.sigma.gamma = c(df = 1, mu = 0, sd = 1),
                              scale.x1 = TRUE)
stan.dat.pdo.s1  <- stan_data(sock.covar, "pdo_winter_stnd",
                              breakpoint1 = 1977,
                              breakpoint2 = 1989,
                              prior.sigma.gamma = c(df = 1, mu = 0, sd = 1),
                              scale.x1 = TRUE)

save(stan.dat.sst.s1, file = "./output/stan_dat_sst_s1.RData")
save(stan.dat.npgo.s1, file = "./output/stan_dat_npgo_s1.RData")
save(stan.dat.pdo.s1, file = "./output/stan_dat_pdo_s1.RData")


## Subset to brood years used in Litzow (2019) CJFAS
sock.covar.s2 <- sock.covar[sock.covar$brood_yr > 1959 &
                            sock.covar$brood_yr < 2007, ]

stan.dat.sst.s2  <- stan_data(sock.covar.s2, "sst_anom_stnd",
                              breakpoint1 = 1977,
                              breakpoint2 = 1989,
                              prior.sigma.gamma = c(df = 1, mu = 0, sd = 1),
                              scale.x1 = TRUE)
stan.dat.npgo.s2 <- stan_data(sock.covar.s2, "npgo_winter_stnd",
                              breakpoint1 = 1977,
                              breakpoint2 = 1989,
                              prior.sigma.gamma = c(df = 1, mu = 0, sd = 1),
                              scale.x1 = TRUE)
stan.dat.pdo.s2  <- stan_data(sock.covar.s2, "pdo_winter_stnd",
                              breakpoint1 = 1977,
                              breakpoint2 = 1989,
                              prior.sigma.gamma = c(df = 1, mu = 0, sd = 1),
                              scale.x1 = TRUE)

save(sock.covar.s2, file = "./output/sock_covar_s2.RData")
save(stan.dat.sst.s2, file = "./output/stan_dat_sst_s2.RData")
save(stan.dat.npgo.s2, file = "./output/stan_dat_npgo_s2.RData")
save(stan.dat.pdo.s2, file = "./output/stan_dat_pdo_s2.RData")



## Set pars to monitor -------------------------------------
pars.hbm1 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha")
pars.hbm2 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "mu_gamma", "sigma_gamma")
pars.hbm3 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma1", "gamma2", "gamma3",
               "mu_gamma1", "mu_gamma2", "mu_gamma3", "sigma_gamma")
pars.hbm4 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "sigma_gamma", "signal_noise")
pars.hbm5 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "sigma_gamma", "signal_noise")
pars.hbm6 <- c("alpha", "beta", "sigma", "phi", "mu_alpha", "sigma_alpha",
               "gamma", "sigma_gamma", "signal_noise")

pars.gen.quant <- c("log_lik", "yrep") ## Generated quantities to monitor


save(pars.hbm1, file = "./output/pars_hbm1.RData")
save(pars.hbm2, file = "./output/pars_hbm2.RData")
save(pars.hbm3, file = "./output/pars_hbm3.RData")
save(pars.hbm4, file = "./output/pars_hbm4.RData")
save(pars.hbm5, file = "./output/pars_hbm5.RData")
save(pars.hbm6, file = "./output/pars_hbm6.RData")
save(pars.gen.quant, file = "./output/pars_gen_quant.RData")


## Prior predictive distributions --------------------------
## Sample priors only (simulate from priors)

## Get data
prior.dat.sst  <- stan_data(sock.covar, "sst_anom_stnd",
                            breakpoint1 = 1977,
                            breakpoint2 = 1989,
                            scale.x1 = TRUE,
                            priors.only = TRUE)
prior.dat.npgo <- stan_data(sock.covar, "npgo_winter_stnd",
                            breakpoint1 = 1977,
                            breakpoint2 = 1989,
                            scale.x1 = TRUE,
                            priors.only = TRUE)
prior.dat.pdo  <- stan_data(sock.covar, "pdo_winter_stnd",
                            breakpoint1 = 1977,
                            breakpoint2 = 1989,
                            scale.x1 = TRUE,
                            priors.only = TRUE)


## Base model
pr.hbm1 <- rstan::stan(file = "./stan/1_hbm_ar1.stan",
                       data = prior.dat.sst,
                       chains = 1,
                       iter = 1000)
save(pr.hbm1, file = "./output/hbm_prior_pred/pr_hbm1.RData")



## SST
pr.hbm2.sst <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                           data = prior.dat.sst,
                           chains = 1,
                           iter = 1000)
save(pr.hbm2.sst, file = "./output/hbm_prior_pred/pr_hbm2_sst.RData")

pr.hbm3.sst <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                           data = prior.dat.sst,
                           chains = 1,
                           iter = 1000)
save(pr.hbm3.sst, file = "./output/hbm_prior_pred/pr_hbm3_sst.RData")

pr.hbm4.sst <- rstan::stan(file = "./stan/4_hbm_ar1_x2nonstat_sigma_same.stan",
                           data = prior.dat.sst,
                           chains = 1,
                           iter = 1000)
save(pr.hbm4.sst, file = "./output/hbm_prior_pred/pr_hbm4_sst.RData")

pr.hbm5.sst <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                           data = prior.dat.sst,
                           chains = 1,
                           iter = 1000)
save(pr.hbm5.sst, file = "./output/hbm_prior_pred/pr_hbm5_sst.RData")

pr.hbm6.sst <- rstan::stan(file = "./stan/6_hbm_ar1_x2nonstat_shared.stan",
                           data = prior.dat.sst,
                           chains = 1,
                           iter = 1000)
save(pr.hbm6.sst, file = "./output/hbm_prior_pred/pr_hbm6_sst.RData")


## NPGO
pr.hbm2.npgo <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                           data = prior.dat.npgo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm2.npgo, file = "./output/hbm_prior_pred/pr_hbm2_npgo.RData")

pr.hbm3.npgo <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                           data = prior.dat.npgo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm3.npgo, file = "./output/hbm_prior_pred/pr_hbm3_npgo.RData")

pr.hbm4.npgo <- rstan::stan(file = "./stan/4_hbm_ar1_x2nonstat_sigma_same.stan",
                           data = prior.dat.npgo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm4.npgo, file = "./output/hbm_prior_pred/pr_hbm4_npgo.RData")

pr.hbm5.npgo <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                           data = prior.dat.npgo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm5.npgo, file = "./output/hbm_prior_pred/pr_hbm5_npgo.RData")

pr.hbm6.npgo <- rstan::stan(file = "./stan/6_hbm_ar1_x2nonstat_shared.stan",
                           data = prior.dat.npgo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm6.npgo, file = "./output/hbm_prior_pred/pr_hbm6_npgo.RData")


## PDO
pr.hbm2.pdo <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                           data = prior.dat.pdo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm2.pdo, file = "./output/hbm_prior_pred/pr_hbm2_pdo.RData")

pr.hbm3.pdo <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                           data = prior.dat.pdo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm3.pdo, file = "./output/hbm_prior_pred/pr_hbm3_pdo.RData")

pr.hbm4.pdo <- rstan::stan(file = "./stan/4_hbm_ar1_x2nonstat_sigma_same.stan",
                           data = prior.dat.pdo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm4.pdo, file = "./output/hbm_prior_pred/pr_hbm4_pdo.RData")

pr.hbm5.pdo <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                           data = prior.dat.pdo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm5.pdo, file = "./output/hbm_prior_pred/pr_hbm5_pdo.RData")

pr.hbm6.pdo <- rstan::stan(file = "./stan/6_hbm_ar1_x2nonstat_shared.stan",
                           data = prior.dat.pdo,
                           chains = 1,
                           iter = 1000)
save(pr.hbm6.pdo, file = "./output/hbm_prior_pred/pr_hbm6_pdo.RData")




## Save prior predictive distributions
plot_prior_pc(pr.hbm1, stan.dat.sst$y, "./figures/hbm_fit/hbm1_prior.pdf")

plot_prior_pc(pr.hbm2.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm2_sst_prior.pdf")
plot_prior_pc(pr.hbm3.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm3_sst_prior.pdf")
plot_prior_pc(pr.hbm4.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm4_sst_prior.pdf")
plot_prior_pc(pr.hbm5.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm5_sst_prior.pdf")
plot_prior_pc(pr.hbm6.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm6_sst_prior.pdf")

plot_prior_pc(pr.hbm2.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm2_npgo_prior.pdf")
plot_prior_pc(pr.hbm3.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm3_npgo_prior.pdf")
plot_prior_pc(pr.hbm4.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm4_npgo_prior.pdf")
plot_prior_pc(pr.hbm5.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm5_npgo_prior.pdf")
plot_prior_pc(pr.hbm6.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm6_npgo_prior.pdf")

plot_prior_pc(pr.hbm2.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm2_pdo_prior.pdf")
plot_prior_pc(pr.hbm3.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm3_pdo_prior.pdf")
plot_prior_pc(pr.hbm4.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm4_pdo_prior.pdf")
plot_prior_pc(pr.hbm5.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm5_pdo_prior.pdf")
plot_prior_pc(pr.hbm6.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm6_pdo_prior.pdf")


###
# pr.mcmc.hbm1 <- As.mcmc.list(pr.hbm1, pars = pars.gen.quant, include = FALSE)
# pdf("./misc/pr_hbm1_diag.pdf", width = 7, height = 5)
#     coda_diag(pr.mcmc.hbm1)
# dev.off()
###



## Run MCMC: base  -----------------------------------------

## Base model
hbm1 <- rstan::stan(file = "./stan/1_hbm_ar1.stan",
                    data = stan.dat.sst,  ## sst is ignored by stan
                    pars = c(pars.hbm1, pars.gen.quant),
                    warmup = 1000,
                    iter = 4000,
                    cores = 4,
                    chains = 4,
                    thin = 3,
                    seed = sample(1:1e6, 1),
                    control = list(adapt_delta = 0.80,
                                   max_treedepth = 11))
save(hbm1, file = "./output/hbm_fit/hbm1.RData")



## SST
hbm2.sst <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm2, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm2.sst, file = "./output/hbm_fit/hbm2_sst.RData")

hbm3.sst <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm3, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm3.sst, file = "./output/hbm_fit/hbm3_sst.RData")

hbm4.sst <- rstan::stan(file = "./stan/4_hbm_ar1_x2nonstat_sigma_same.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm4, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.95,
                                       max_treedepth = 12))
save(hbm4.sst, file = "./output/hbm_fit/hbm4_sst.RData")

hbm5.sst <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm5, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 10))
save(hbm5.sst, file = "./output/hbm_fit/hbm5_sst.RData")

hbm6.sst <- rstan::stan(file = "./stan/6_hbm_ar1_x2nonstat_shared.stan",
                        data = stan.dat.sst,
                        pars = c(pars.hbm6, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm6.sst, file = "./output/hbm_fit/hbm6_sst.RData")



## NPGO
hbm2.npgo <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                         data = stan.dat.npgo,
                         pars = c(pars.hbm2, pars.gen.quant),
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         thin = 3,
                         seed = sample(1:1e6, 1),
                         control = list(adapt_delta = 0.90,
                                        max_treedepth = 10))
save(hbm2.npgo, file = "./output/hbm_fit/hbm2_npgo.RData")

hbm3.npgo <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                         data = stan.dat.npgo,
                         pars = c(pars.hbm3, pars.gen.quant),
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         thin = 3,
                         seed = sample(1:1e6, 1),
                         control = list(adapt_delta = 0.90,
                                        max_treedepth = 10))
save(hbm3.npgo, file = "./output/hbm_fit/hbm3_npgo.RData")

hbm4.npgo <- rstan::stan(file = "./stan/4_hbm_ar1_x2nonstat_sigma_same.stan",
                         data = stan.dat.npgo,
                         pars = c(pars.hbm4, pars.gen.quant),
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         thin = 3,
                         seed = sample(1:1e6, 1),
                         control = list(adapt_delta = 0.95,
                                        max_treedepth = 12))
save(hbm4.npgo, file = "./output/hbm_fit/hbm4_npgo.RData")

hbm5.npgo <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                         data = stan.dat.npgo,
                         pars = c(pars.hbm5, pars.gen.quant),
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         thin = 3,
                         seed = sample(1:1e6, 1),
                         control = list(adapt_delta = 0.98,
                                        max_treedepth = 10))
save(hbm5.npgo, file = "./output/hbm_fit/hbm5_npgo.RData")

hbm6.npgo <- rstan::stan(file = "./stan/6_hbm_ar1_x2nonstat_shared.stan",
                         data = stan.dat.npgo,
                         pars = c(pars.hbm6, pars.gen.quant),
                         warmup = 1000,
                         iter = 4000,
                         cores = 4,
                         chains = 4,
                         thin = 3,
                         seed = sample(1:1e6, 1),
                         control = list(adapt_delta = 0.90,
                                        max_treedepth = 10))
save(hbm6.npgo, file = "./output/hbm_fit/hbm6_npgo.RData")



## PDO
hbm2.pdo <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                        data = stan.dat.pdo,
                        pars = c(pars.hbm2, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm2.pdo, file = "./output/hbm_fit/hbm2_pdo.RData")

hbm3.pdo <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                        data = stan.dat.pdo,
                        pars = c(pars.hbm3, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm3.pdo, file = "./output/hbm_fit/hbm3_pdo.RData")

hbm4.pdo <- rstan::stan(file = "./stan/4_hbm_ar1_x2nonstat_sigma_same.stan",
                        data = stan.dat.pdo,
                        pars = c(pars.hbm4, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.99,
                                       max_treedepth = 13))
save(hbm4.pdo, file = "./output/hbm_fit/hbm4_pdo.RData")

hbm5.pdo <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                        data = stan.dat.pdo,
                        pars = c(pars.hbm5, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.98,
                                       max_treedepth = 12))
save(hbm5.pdo, file = "./output/hbm_fit/hbm5_pdo.RData")

hbm6.pdo <- rstan::stan(file = "./stan/6_hbm_ar1_x2nonstat_shared.stan",
                        data = stan.dat.pdo,
                        pars = c(pars.hbm6, pars.gen.quant),
                        warmup = 1000,
                        iter = 4000,
                        cores = 4,
                        chains = 4,
                        thin = 3,
                        seed = sample(1:1e6, 1),
                        control = list(adapt_delta = 0.90,
                                       max_treedepth = 10))
save(hbm6.pdo, file = "./output/hbm_fit/hbm6_pdo.RData")



## Run MCMC: S1 --------------------------------------------
hbm5.sst.s1 <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                           data = stan.dat.sst.s1,
                           pars = c(pars.hbm5, pars.gen.quant),
                           warmup = 1000,
                           iter = 4000,
                           cores = 4,
                           chains = 4,
                           thin = 3,
                           seed = sample(1:1e6, 1),
                           control = list(adapt_delta = 0.98,
                                          max_treedepth = 10))
save(hbm5.sst.s1, file = "./output/hbm_fit/hbm5_sst_s1.RData")

hbm5.npgo.s1 <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                            data = stan.dat.npgo.s1,
                            pars = c(pars.hbm5, pars.gen.quant),
                            warmup = 1000,
                            iter = 4000,
                            cores = 4,
                            chains = 4,
                            thin = 3,
                            seed = sample(1:1e6, 1),
                            control = list(adapt_delta = 0.98,
                                           max_treedepth = 10))
save(hbm5.npgo.s1, file = "./output/hbm_fit/hbm5_npgo_s1.RData")

hbm5.pdo.s1 <- rstan::stan(file = "./stan/5_hbm_ar1_x2nonstat_sigma_diff.stan",
                            data = stan.dat.pdo.s1,
                            pars = c(pars.hbm5, pars.gen.quant),
                            warmup = 1000,
                            iter = 4000,
                            cores = 4,
                            chains = 4,
                            thin = 3,
                            seed = sample(1:1e6, 1),
                            control = list(adapt_delta = 0.98,
                                           max_treedepth = 10))
save(hbm5.pdo.s1, file = "./output/hbm_fit/hbm5_pdo_s1.RData")



## Run MCMC: S2 --------------------------------------------

## HBM2
hbm2.sst.s2 <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                           data = stan.dat.sst.s2,
                           pars = c(pars.hbm2, pars.gen.quant),
                           warmup = 1000,
                           iter = 4000,
                           cores = 4,
                           chains = 4,
                           thin = 3,
                           seed = sample(1:1e6, 1),
                           control = list(adapt_delta = 0.90,
                                          max_treedepth = 10))
save(hbm2.sst.s2, file = "./output/hbm_fit/hbm2_sst_s2.RData")

hbm2.npgo.s2 <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                            data = stan.dat.npgo.s2,
                            pars = c(pars.hbm2, pars.gen.quant),
                            warmup = 1000,
                            iter = 4000,
                            cores = 4,
                            chains = 4,
                            thin = 3,
                            seed = sample(1:1e6, 1),
                            control = list(adapt_delta = 0.90,
                                           max_treedepth = 10))
save(hbm2.npgo.s2, file = "./output/hbm_fit/hbm2_npgo_s2.RData")

hbm2.pdo.s2 <- rstan::stan(file = "./stan/2_hbm_ar1_x2stat.stan",
                           data = stan.dat.pdo.s2,
                           pars = c(pars.hbm2, pars.gen.quant),
                           warmup = 1000,
                           iter = 4000,
                           cores = 4,
                           chains = 4,
                           thin = 3,
                           seed = sample(1:1e6, 1),
                           control = list(adapt_delta = 0.90,
                                          max_treedepth = 10))
save(hbm2.pdo.s2, file = "./output/hbm_fit/hbm2_pdo_s2.RData")


## HBM3
hbm3.sst.s2 <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                           data = stan.dat.sst.s2,
                           pars = c(pars.hbm3, pars.gen.quant),
                           warmup = 1000,
                           iter = 4000,
                           cores = 4,
                           chains = 4,
                           thin = 3,
                           seed = sample(1:1e6, 1),
                           control = list(adapt_delta = 0.90,
                                          max_treedepth = 10))
save(hbm3.sst.s2, file = "./output/hbm_fit/hbm3_sst_s2.RData")

hbm3.npgo.s2 <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                            data = stan.dat.npgo.s2,
                            pars = c(pars.hbm3, pars.gen.quant),
                            warmup = 1000,
                            iter = 4000,
                            cores = 4,
                            chains = 4,
                            thin = 3,
                            seed = sample(1:1e6, 1),
                            control = list(adapt_delta = 0.90,
                                           max_treedepth = 10))
save(hbm3.npgo.s2, file = "./output/hbm_fit/hbm3_npgo_s2.RData")

hbm3.pdo.s2 <- rstan::stan(file = "./stan/3_hbm_ar1_x2era3.stan",
                           data = stan.dat.pdo.s2,
                           pars = c(pars.hbm3, pars.gen.quant),
                           warmup = 1000,
                           iter = 4000,
                           cores = 4,
                           chains = 4,
                           thin = 3,
                           seed = sample(1:1e6, 1),
                           control = list(adapt_delta = 0.90,
                                          max_treedepth = 10))
save(hbm3.pdo.s2, file = "./output/hbm_fit/hbm3_pdo_s2.RData")



## Check pathology -----------------------------------------

## NUTS diag
rstan::check_hmc_diagnostics(hbm1)

rstan::check_hmc_diagnostics(hbm2.sst)
rstan::check_hmc_diagnostics(hbm3.sst)
rstan::check_hmc_diagnostics(hbm4.sst)
rstan::check_hmc_diagnostics(hbm5.sst)
rstan::check_hmc_diagnostics(hbm6.sst)

rstan::check_hmc_diagnostics(hbm2.npgo)
rstan::check_hmc_diagnostics(hbm3.npgo)
rstan::check_hmc_diagnostics(hbm4.npgo)
rstan::check_hmc_diagnostics(hbm5.npgo)
rstan::check_hmc_diagnostics(hbm6.npgo)

rstan::check_hmc_diagnostics(hbm2.pdo)
rstan::check_hmc_diagnostics(hbm3.pdo)
rstan::check_hmc_diagnostics(hbm4.pdo)
rstan::check_hmc_diagnostics(hbm5.pdo)
rstan::check_hmc_diagnostics(hbm6.pdo)

rstan::check_hmc_diagnostics(hbm5.sst.s1)
rstan::check_hmc_diagnostics(hbm5.npgo.s1)
rstan::check_hmc_diagnostics(hbm5.pdo.s1)

rstan::check_hmc_diagnostics(hbm2.sst.s2)
rstan::check_hmc_diagnostics(hbm2.npgo.s2)
rstan::check_hmc_diagnostics(hbm2.pdo.s2)
rstan::check_hmc_diagnostics(hbm3.sst.s2)
rstan::check_hmc_diagnostics(hbm3.npgo.s2)
rstan::check_hmc_diagnostics(hbm3.pdo.s2)


## Neff
neff_lowest(hbm1, pars = pars.hbm1)

neff_lowest(hbm2.sst, pars = pars.hbm2)
neff_lowest(hbm3.sst, pars = pars.hbm3)
neff_lowest(hbm4.sst, pars = pars.hbm4)
neff_lowest(hbm5.sst, pars = pars.hbm5)
neff_lowest(hbm6.sst, pars = pars.hbm6)

neff_lowest(hbm2.npgo, pars = pars.hbm2)
neff_lowest(hbm3.npgo, pars = pars.hbm3)
neff_lowest(hbm4.npgo, pars = pars.hbm4)
neff_lowest(hbm5.npgo, pars = pars.hbm5)
neff_lowest(hbm6.npgo, pars = pars.hbm6)

neff_lowest(hbm2.pdo, pars = pars.hbm2)
neff_lowest(hbm3.pdo, pars = pars.hbm3)
neff_lowest(hbm4.pdo, pars = pars.hbm4)
neff_lowest(hbm5.pdo, pars = pars.hbm5)
neff_lowest(hbm6.pdo, pars = pars.hbm6)

neff_lowest(hbm5.sst.s1)
neff_lowest(hbm5.npgo.s1)
neff_lowest(hbm5.pdo.s1)

neff_lowest(hbm2.sst.s2)
neff_lowest(hbm2.npgo.s2)
neff_lowest(hbm2.pdo.s2)
neff_lowest(hbm3.sst.s2)
neff_lowest(hbm3.npgo.s2)
neff_lowest(hbm3.pdo.s2)


## Rhat
rhat_highest(hbm1, pars = pars.hbm1)

rhat_highest(hbm2.sst, pars = pars.hbm2)
rhat_highest(hbm3.sst, pars = pars.hbm3)
rhat_highest(hbm4.sst, pars = pars.hbm4)
rhat_highest(hbm5.sst, pars = pars.hbm5)
rhat_highest(hbm6.sst, pars = pars.hbm6)

rhat_highest(hbm2.npgo, pars = pars.hbm2)
rhat_highest(hbm3.npgo, pars = pars.hbm3)
rhat_highest(hbm4.npgo, pars = pars.hbm4)
rhat_highest(hbm5.npgo, pars = pars.hbm5)
rhat_highest(hbm6.npgo, pars = pars.hbm6)

rhat_highest(hbm2.pdo, pars = pars.hbm2)
rhat_highest(hbm3.pdo, pars = pars.hbm3)
rhat_highest(hbm4.pdo, pars = pars.hbm4)
rhat_highest(hbm5.pdo, pars = pars.hbm5)
rhat_highest(hbm6.pdo, pars = pars.hbm6)

rhat_highest(hbm5.sst.s1)
rhat_highest(hbm5.npgo.s1)
rhat_highest(hbm5.pdo.s1)

rhat_highest(hbm2.sst.s2)
rhat_highest(hbm2.npgo.s2)
rhat_highest(hbm2.pdo.s2)
rhat_highest(hbm3.sst.s2)
rhat_highest(hbm3.npgo.s2)
rhat_highest(hbm3.pdo.s2)


## Pairs
# pairs_lowest(hbm1, pars = pars.hbm1)
#
# pairs_lowest(hbm2.sst, pars = pars.hbm2)
# pairs_lowest(hbm3.sst, pars = pars.hbm3)
# pairs_lowest(hbm4.sst, pars = pars.hbm4)
# pairs_lowest(hbm5.sst, pars = pars.hbm5)
# pairs_lowest(hbm6.sst, pars = pars.hbm6)
#
# pairs_lowest(hbm2.npgo, pars = pars.hbm2)
# pairs_lowest(hbm3.npgo, pars = pars.hbm3)
# pairs_lowest(hbm4.npgo, pars = pars.hbm4)
# pairs_lowest(hbm5.npgo, pars = pars.hbm5)
# pairs_lowest(hbm6.npgo, pars = pars.hbm6)
#
# pairs_lowest(hbm2.pdo, pars = pars.hbm2)
# pairs_lowest(hbm3.pdo, pars = pars.hbm3)
# pairs_lowest(hbm4.pdo, pars = pars.hbm4)
# pairs_lowest(hbm5.pdo, pars = pars.hbm5)
# pairs_lowest(hbm6.pdo, pars = pars.hbm6)
#
# pairs_lowest(hbm5.sst.s1)
# pairs_lowest(hbm5.npgo.s1)
# pairs_lowest(hbm5.pdo.s1)
#
# pairs_lowest(hbm2.sst.s2)
# pairs_lowest(hbm2.npgo.s2)
# pairs_lowest(hbm2.pdo.s2)
# pairs_lowest(hbm3.sst.s2)
# pairs_lowest(hbm3.npgo.s2)
# pairs_lowest(hbm3.pdo.s2)


## Run time
rstan::get_elapsed_time(hbm1)

rstan::get_elapsed_time(hbm2.sst)
rstan::get_elapsed_time(hbm3.sst)
rstan::get_elapsed_time(hbm4.sst)
rstan::get_elapsed_time(hbm5.sst)
rstan::get_elapsed_time(hbm6.sst)

rstan::get_elapsed_time(hbm2.npgo)
rstan::get_elapsed_time(hbm3.npgo)
rstan::get_elapsed_time(hbm4.npgo)
rstan::get_elapsed_time(hbm5.npgo)
rstan::get_elapsed_time(hbm6.npgo)

rstan::get_elapsed_time(hbm2.pdo)
rstan::get_elapsed_time(hbm3.pdo)
rstan::get_elapsed_time(hbm4.pdo)
rstan::get_elapsed_time(hbm5.pdo)
rstan::get_elapsed_time(hbm6.pdo)

rstan::get_elapsed_time(hbm5.sst.s1)
rstan::get_elapsed_time(hbm5.npgo.s1)
rstan::get_elapsed_time(hbm5.pdo.s1)

rstan::get_elapsed_time(hbm2.sst.s2)
rstan::get_elapsed_time(hbm2.npgo.s2)
rstan::get_elapsed_time(hbm2.pdo.s2)
rstan::get_elapsed_time(hbm3.sst.s2)
rstan::get_elapsed_time(hbm3.npgo.s2)
rstan::get_elapsed_time(hbm3.pdo.s2)



## MCMC diagnostics ----------------------------------------

## Base model
pdf("./figures/hbm_fit/hbm1_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm1, pars = pars.hbm1), total_draws(hbm1))
    coda_rhat(get_rhat(hbm1, pars = pars.hbm1))
    coda_diag(As.mcmc.list(hbm1, pars = pars.hbm1))
dev.off()


## SST
pdf("./figures/hbm_fit/hbm2_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm2.sst, pars = pars.hbm2), total_draws(hbm2.sst))
    coda_rhat(get_rhat(hbm2.sst, pars = pars.hbm2))
    coda_diag(As.mcmc.list(hbm2.sst, pars = pars.hbm2))
dev.off()

pdf("./figures/hbm_fit/hbm3_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.sst, pars = pars.hbm3), total_draws(hbm3.sst))
    coda_rhat(get_rhat(hbm3.sst, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.sst, pars = pars.hbm3))
dev.off()

pdf("./figures/hbm_fit/hbm4_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm4.sst, pars = pars.hbm4), total_draws(hbm4.sst))
    coda_rhat(get_rhat(hbm4.sst, pars = pars.hbm4))
    coda_diag(As.mcmc.list(hbm4.sst, pars = pars.hbm4))
dev.off()

pdf("./figures/hbm_fit/hbm5_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.sst, pars = pars.hbm5), total_draws(hbm5.sst))
    coda_rhat(get_rhat(hbm5.sst, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.sst, pars = pars.hbm5))
dev.off()

pdf("./figures/hbm_fit/hbm6_sst_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm6.sst, pars = pars.hbm6), total_draws(hbm6.sst))
    coda_rhat(get_rhat(hbm6.sst, pars = pars.hbm6))
    coda_diag(As.mcmc.list(hbm6.sst, pars = pars.hbm6))
dev.off()


## NPGO
pdf("./figures/hbm_fit/hbm2_npgo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm2.npgo, pars = pars.hbm2), total_draws(hbm2.npgo))
    coda_rhat(get_rhat(hbm2.npgo, pars = pars.hbm2))
    coda_diag(As.mcmc.list(hbm2.npgo, pars = pars.hbm2))
dev.off()

pdf("./figures/hbm_fit/hbm3_npgo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.npgo, pars = pars.hbm3), total_draws(hbm3.npgo))
    coda_rhat(get_rhat(hbm3.npgo, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.npgo, pars = pars.hbm3))
dev.off()

pdf("./figures/hbm_fit/hbm4_npgo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm4.npgo, pars = pars.hbm4), total_draws(hbm4.npgo))
    coda_rhat(get_rhat(hbm4.npgo, pars = pars.hbm4))
    coda_diag(As.mcmc.list(hbm4.npgo, pars = pars.hbm4))
dev.off()

pdf("./figures/hbm_fit/hbm5_npgo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.npgo, pars = pars.hbm5), total_draws(hbm5.npgo))
    coda_rhat(get_rhat(hbm5.npgo, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.npgo, pars = pars.hbm5))
dev.off()

pdf("./figures/hbm_fit/hbm6_npgo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm6.npgo, pars = pars.hbm6), total_draws(hbm6.npgo))
    coda_rhat(get_rhat(hbm6.npgo, pars = pars.hbm6))
    coda_diag(As.mcmc.list(hbm6.npgo, pars = pars.hbm6))
dev.off()



## PDO
pdf("./figures/hbm_fit/hbm2_pdo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm2.pdo, pars = pars.hbm2), total_draws(hbm2.pdo))
    coda_rhat(get_rhat(hbm2.pdo, pars = pars.hbm2))
    coda_diag(As.mcmc.list(hbm2.pdo, pars = pars.hbm2))
dev.off()

pdf("./figures/hbm_fit/hbm3_pdo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.pdo, pars = pars.hbm3), total_draws(hbm3.pdo))
    coda_rhat(get_rhat(hbm3.pdo, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.pdo, pars = pars.hbm3))
dev.off()

pdf("./figures/hbm_fit/hbm4_pdo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm4.pdo, pars = pars.hbm4), total_draws(hbm4.pdo))
    coda_rhat(get_rhat(hbm4.pdo, pars = pars.hbm4))
    coda_diag(As.mcmc.list(hbm4.pdo, pars = pars.hbm4))
dev.off()

pdf("./figures/hbm_fit/hbm5_pdo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.pdo, pars = pars.hbm5), total_draws(hbm5.pdo))
    coda_rhat(get_rhat(hbm5.pdo, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.pdo, pars = pars.hbm5))
dev.off()

pdf("./figures/hbm_fit/hbm6_pdo_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm6.pdo, pars = pars.hbm6), total_draws(hbm6.pdo))
    coda_rhat(get_rhat(hbm6.pdo, pars = pars.hbm6))
    coda_diag(As.mcmc.list(hbm6.pdo, pars = pars.hbm6))
dev.off()


## Sensitivity #1
pdf("./figures/hbm_fit/hbm5_sst_s1_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.sst.s1, pars = pars.hbm5), total_draws(hbm5.sst.s1))
    coda_rhat(get_rhat(hbm5.sst.s1, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.sst.s1, pars = pars.hbm5))
dev.off()

pdf("./figures/hbm_fit/hbm5_npgo_s1_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.npgo.s1, pars = pars.hbm5), total_draws(hbm5.npgo.s1))
    coda_rhat(get_rhat(hbm5.npgo.s1, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.npgo.s1, pars = pars.hbm5))
dev.off()

pdf("./figures/hbm_fit/hbm5_pdo_s1_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm5.pdo.s1, pars = pars.hbm5), total_draws(hbm5.pdo.s1))
    coda_rhat(get_rhat(hbm5.pdo.s1, pars = pars.hbm5))
    coda_diag(As.mcmc.list(hbm5.pdo.s1, pars = pars.hbm5))
dev.off()


## Sensitivity #2
pdf("./figures/hbm_fit/hbm2_sst_s2_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm2.sst.s2, pars = pars.hbm2), total_draws(hbm2.sst.s2))
    coda_rhat(get_rhat(hbm2.sst.s2, pars = pars.hbm2))
    coda_diag(As.mcmc.list(hbm2.sst.s2, pars = pars.hbm2))
dev.off()

pdf("./figures/hbm_fit/hbm2_npgo_s2_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm2.npgo.s2, pars = pars.hbm2), total_draws(hbm2.npgo.s2))
    coda_rhat(get_rhat(hbm2.npgo.s2, pars = pars.hbm2))
    coda_diag(As.mcmc.list(hbm2.npgo.s2, pars = pars.hbm2))
dev.off()

pdf("./figures/hbm_fit/hbm2_pdo_s2_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm2.pdo.s2, pars = pars.hbm2), total_draws(hbm2.pdo.s2))
    coda_rhat(get_rhat(hbm2.pdo.s2, pars = pars.hbm2))
    coda_diag(As.mcmc.list(hbm2.pdo.s2, pars = pars.hbm2))
dev.off()


pdf("./figures/hbm_fit/hbm3_sst_s2_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.sst.s2, pars = pars.hbm3), total_draws(hbm3.sst.s2))
    coda_rhat(get_rhat(hbm3.sst.s2, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.sst.s2, pars = pars.hbm3))
dev.off()

pdf("./figures/hbm_fit/hbm3_npgo_s2_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.npgo.s2, pars = pars.hbm3), total_draws(hbm3.npgo.s2))
    coda_rhat(get_rhat(hbm3.npgo.s2, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.npgo.s2, pars = pars.hbm3))
dev.off()

pdf("./figures/hbm_fit/hbm3_pdo_s2_diag.pdf", width = 7, height = 5)
    coda_neff(get_neff(hbm3.pdo.s2, pars = pars.hbm3), total_draws(hbm3.pdo.s2))
    coda_rhat(get_rhat(hbm3.pdo.s2, pars = pars.hbm3))
    coda_diag(As.mcmc.list(hbm3.pdo.s2, pars = pars.hbm3))
dev.off()



## Posterior predictive checks -----------------------------
plot_post_pc(hbm1, stan.dat.sst$y, "./figures/hbm_fit/hbm1_sst_yrep.pdf")

plot_post_pc(hbm2.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm2_sst_yrep.pdf")
plot_post_pc(hbm3.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm3_sst_yrep.pdf")
plot_post_pc(hbm4.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm4_sst_yrep.pdf")
plot_post_pc(hbm5.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm5_sst_yrep.pdf")
plot_post_pc(hbm6.sst, stan.dat.sst$y, "./figures/hbm_fit/hbm6_sst_yrep.pdf")

plot_post_pc(hbm2.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm2_npgo_yrep.pdf")
plot_post_pc(hbm3.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm3_npgo_yrep.pdf")
plot_post_pc(hbm4.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm4_npgo_yrep.pdf")
plot_post_pc(hbm5.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm5_npgo_yrep.pdf")
plot_post_pc(hbm6.npgo, stan.dat.npgo$y, "./figures/hbm_fit/hbm6_npgo_yrep.pdf")

plot_post_pc(hbm2.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm2_pdo_yrep.pdf")
plot_post_pc(hbm3.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm3_pdo_yrep.pdf")
plot_post_pc(hbm4.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm4_pdo_yrep.pdf")
plot_post_pc(hbm5.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm5_pdo_yrep.pdf")
plot_post_pc(hbm6.pdo, stan.dat.pdo$y, "./figures/hbm_fit/hbm6_pdo_yrep.pdf")

plot_post_pc(hbm5.sst.s1, stan.dat.sst.s1$y, "./figures/hbm_fit/hbm5_sst_s1_yrep.pdf")
plot_post_pc(hbm5.npgo.s1, stan.dat.npgo.s1$y, "./figures/hbm_fit/hbm5_npgo_s1_yrep.pdf")
plot_post_pc(hbm5.pdo.s1, stan.dat.pdo.s1$y, "./figures/hbm_fit/hbm5_pdo_s1_yrep.pdf")

plot_post_pc(hbm2.sst.s2, stan.dat.sst.s2$y, sock.covar.s2, "./figures/hbm_fit/hbm2_sst_s2_yrep.pdf")
plot_post_pc(hbm2.npgo.s2, stan.dat.npgo.s2$y, sock.covar.s2, "./figures/hbm_fit/hbm2_npgo_s2_yrep.pdf")
plot_post_pc(hbm2.pdo.s2, stan.dat.pdo.s2$y, sock.covar.s2, "./figures/hbm_fit/hbm2_pdo_s2_yrep.pdf")
plot_post_pc(hbm3.sst.s2, stan.dat.sst.s2$y, sock.covar.s2, "./figures/hbm_fit/hbm3_sst_s2_yrep.pdf")
plot_post_pc(hbm3.npgo.s2, stan.dat.npgo.s2$y, sock.covar.s2, "./figures/hbm_fit/hbm3_npgo_s2_yrep.pdf")
plot_post_pc(hbm3.pdo.s2, stan.dat.pdo.s2$y, sock.covar.s2, "./figures/hbm_fit/hbm3_pdo_s2_yrep.pdf")



## LOOIC ---------------------------------------------------
loo.hbm1     <- rstan::loo(hbm1, cores = 4)

loo.hbm2.sst <- rstan::loo(hbm2.sst, cores = 4)
loo.hbm3.sst <- rstan::loo(hbm3.sst, cores = 4)
loo.hbm4.sst <- rstan::loo(hbm4.sst, cores = 4)
loo.hbm5.sst <- rstan::loo(hbm5.sst, cores = 4)
loo.hbm6.sst <- rstan::loo(hbm6.sst, cores = 4)

loo.hbm2.npgo <- rstan::loo(hbm2.npgo, cores = 4)
loo.hbm3.npgo <- rstan::loo(hbm3.npgo, cores = 4)
loo.hbm4.npgo <- rstan::loo(hbm4.npgo, cores = 4)
loo.hbm5.npgo <- rstan::loo(hbm5.npgo, cores = 4)
loo.hbm6.npgo <- rstan::loo(hbm6.npgo, cores = 4)

loo.hbm2.pdo <- rstan::loo(hbm2.pdo, cores = 4)
loo.hbm3.pdo <- rstan::loo(hbm3.pdo, cores = 4)
loo.hbm4.pdo <- rstan::loo(hbm4.pdo, cores = 4)
loo.hbm5.pdo <- rstan::loo(hbm5.pdo, cores = 4)
loo.hbm6.pdo <- rstan::loo(hbm6.pdo, cores = 4)



save(loo.hbm1, file = "./output/hbm_loo/loo_hbm1.RData")

save(loo.hbm2.sst, file = "./output/hbm_loo/loo_hbm2_sst.RData")
save(loo.hbm3.sst, file = "./output/hbm_loo/loo_hbm3_sst.RData")
save(loo.hbm4.sst, file = "./output/hbm_loo/loo_hbm4_sst.RData")
save(loo.hbm5.sst, file = "./output/hbm_loo/loo_hbm5_sst.RData")
save(loo.hbm6.sst, file = "./output/hbm_loo/loo_hbm6_sst.RData")

save(loo.hbm2.npgo, file = "./output/hbm_loo/loo_hbm2_npgo.RData")
save(loo.hbm3.npgo, file = "./output/hbm_loo/loo_hbm3_npgo.RData")
save(loo.hbm4.npgo, file = "./output/hbm_loo/loo_hbm4_npgo.RData")
save(loo.hbm5.npgo, file = "./output/hbm_loo/loo_hbm5_npgo.RData")
save(loo.hbm6.npgo, file = "./output/hbm_loo/loo_hbm6_npgo.RData")

save(loo.hbm2.pdo, file = "./output/hbm_loo/loo_hbm2_pdo.RData")
save(loo.hbm3.pdo, file = "./output/hbm_loo/loo_hbm3_pdo.RData")
save(loo.hbm4.pdo, file = "./output/hbm_loo/loo_hbm4_pdo.RData")
save(loo.hbm5.pdo, file = "./output/hbm_loo/loo_hbm5_pdo.RData")
save(loo.hbm6.pdo, file = "./output/hbm_loo/loo_hbm6_pdo.RData")


# load_rdata("./output/hbm_loo/")
loo::compare(loo.hbm1, loo.hbm2.sst, loo.hbm3.sst,
             loo.hbm4.sst, loo.hbm5.sst, loo.hbm6.sst)

loo::compare(loo.hbm1, loo.hbm2.npgo, loo.hbm3.npgo,
             loo.hbm4.npgo, loo.hbm5.npgo, loo.hbm6.npgo)

loo::compare(loo.hbm1, loo.hbm2.pdo, loo.hbm3.pdo,
             loo.hbm4.pdo, loo.hbm5.pdo, loo.hbm6.pdo)

loo::compare(loo.hbm1, loo.hbm2.sst, loo.hbm3.sst,
             loo.hbm4.sst, loo.hbm5.sst, loo.hbm6.sst,
             loo.hbm1, loo.hbm2.npgo, loo.hbm3.npgo,
             loo.hbm4.npgo, loo.hbm5.npgo, loo.hbm6.npgo,
             loo.hbm1, loo.hbm2.pdo, loo.hbm3.pdo,
             loo.hbm4.pdo, loo.hbm5.pdo, loo.hbm6.pdo)

sum(pareto_k_values(loo.hbm1) > 0.7)

sum(pareto_k_values(loo.hbm2.sst) > 0.7)
sum(pareto_k_values(loo.hbm3.sst) > 0.7)
sum(pareto_k_values(loo.hbm4.sst) > 0.7)
sum(pareto_k_values(loo.hbm5.sst) > 0.7)
sum(pareto_k_values(loo.hbm6.sst) > 0.7)

sum(pareto_k_values(loo.hbm2.npgo) > 0.7)
sum(pareto_k_values(loo.hbm3.npgo) > 0.7)
sum(pareto_k_values(loo.hbm4.npgo) > 0.7)
sum(pareto_k_values(loo.hbm5.npgo) > 0.7)
sum(pareto_k_values(loo.hbm6.npgo) > 0.7)

sum(pareto_k_values(loo.hbm2.pdo) > 0.7)
sum(pareto_k_values(loo.hbm3.pdo) > 0.7)
sum(pareto_k_values(loo.hbm4.pdo) > 0.7)
sum(pareto_k_values(loo.hbm5.pdo) > 0.7)
sum(pareto_k_values(loo.hbm6.pdo) > 0.7)


## Base model
pdf("./figures/hbm_fit/hbm1_loo.pdf", width = 7, height = 5)
    plot(loo.hbm1, label_points = TRUE)
dev.off()

## SST
pdf("./figures/hbm_fit/hbm2_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm2.sst, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm3_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm3.sst, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm4_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm4.sst, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm5_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm5.sst, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm6_sst_loo.pdf", width = 7, height = 5)
    plot(loo.hbm6.sst, label_points = TRUE)
dev.off()

## NPGO
pdf("./figures/hbm_fit/hbm2_npgo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm2.npgo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm3_npgo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm3.npgo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm4_npgo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm4.npgo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm5_npgo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm5.npgo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm6_npgo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm6.npgo, label_points = TRUE)
dev.off()

## PDO
pdf("./figures/hbm_fit/hbm2_pdo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm2.pdo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm3_pdo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm3.pdo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm4_pdo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm4.pdo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm5_pdo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm5.pdo, label_points = TRUE)
dev.off()
pdf("./figures/hbm_fit/hbm6_pdo_loo.pdf", width = 7, height = 5)
    plot(loo.hbm6.pdo, label_points = TRUE)
dev.off()
