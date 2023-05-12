## Functions for analysis


## replace_dup ---------------------------------------------
replace_dup <- function(x, replacement = "") {
    ## Replace sequentially duplicate values in a vector
    ##
    ## x = character vector
    ## replacement = string to replace duplicates with

    x <- as.character(x)
    z <- rep(NA, length(x))
    z[1] <- x[1]
    for(i in 2:length(x)) {
        z[i] <- ifelse(x[i] == x[i-1], replacement, x[i])
    }
    return(z)
}



## sim_ss --------------------------------------------------
sim_ss <- function(n_years = 50,
                   alpha = 1.5,
                   beta = -0.1,
                   sigma = 0.25,
                   phi = 0,
                   gamma = rep(1.5, n_years),
                   prior.sigma.gamma = c(df = 3, mu = 0, sd = 0.5),
                   priors.only = 0) {
    ## Simulate single-stock state-space linear model
    ##
    ## n_years = numeric, length of time series
    ## alpha = numeric, true intercept
    ## beta = numeric, true x1 slope
    ## sigma = numeric, true observation error
    ## phi = numeric, true residual autocorrelation
    ## gamma = numeric vector giving x2 slopes,
    ##         should be same length as n_years

    x1 <- runif(n_years)
    x2  <- ts_ar1(n_years, 0, 0.25, sd = 0.5)
    res  <- ts_ar1(n_years, 0, phi, sd = sigma)
    y <- alpha + beta * x1 + gamma * x2 + res

    ## data frame of simulated data
    df <- data.frame(year = 1:n_years,
                     x1 = x1,
                     x2 = x2,
                     res = res,
                     y = y,
                     gamma = gamma)

    ## Stan data list
    lst <- list(y = y,
                x1 = x1,
                x2 = x2,
                N = n_years,
                sigma_gamma_df = prior.sigma.gamma[1],
                sigma_gamma_mu = prior.sigma.gamma[2],
                sigma_gamma_sd = prior.sigma.gamma[3],
                priors_only = priors.only)

    return(list(data = df, stan = lst))
}



## t-scaled ------------------------------------------------
dt.scaled<-function(x, df, mean=0, sd=1, ncp, log = FALSE) {
        if(!log) stats::dt((x-mean)/sd, df, ncp=ncp, log=FALSE)/sd
        else stats::dt((x-mean)/sd, df, ncp=ncp, log=TRUE)- log(sd)
}
pt.scaled<-function(q, df, mean=0, sd=1, ncp, lower.tail = TRUE, log.p = FALSE) {
        stats::pt((q-mean)/sd, df, ncp=ncp, log.p=log.p)
}
qt.scaled<-function(p, df, mean=0, sd=1, ncp, lower.tail = TRUE, log.p = FALSE) {
        mean+sd*stats::qt(p, df, ncp=ncp, log.p=log.p)
}
rt.scaled<-function(n, df, mean=0, sd=1, ncp) {
        mean+sd*stats::rt(n, df, ncp=ncp)
}


## ocean_region_lab ----------------------------------------
ocean_region_lab <- function(data) {
    ## Add ocean_region_lab column to data.frame
    ##
    ## data = a data.frame with a "ocean_region" column

    if(!"ocean_region" %in% names(data)) {
        stop("No ocean_region column found")
    }
    lab <- rep(NA, nrow(data))
    lab <- ifelse(data$ocean_region == "WC", "West Coast", lab)
    lab <- ifelse(data$ocean_region == "GOA", "Gulf of Alaska", lab)
    lab <- ifelse(data$ocean_region == "BS", "Bering Sea", lab)
    lab <- factor(lab, levels = unique(lab))
    data[["ocean_region_lab"]] <- lab
    return(data)
}



## plot_hbm_dens -------------------------------------------
plot_hbm_dens <- function(stanfit,
                          pdf.file = NULL,
                          width = 10,
                          height = 7) {
    ## HMB posterior density plots
    ##
    ## CI inner = 95%
    ##
    ## stanfit = stanfit model object
    ## pdf.file = path to pdf file to save
    ## width = pdf file width
    ## height = pdf file height

    model <- deparse(substitute(stanfit))
    ssm <- grepl("hbm4", model) |
           grepl("hbm5", model) |
           grepl("hbm6", model)

    mu.alpha <- mcmc_areas(as.matrix(stanfit, pars = "mu_alpha"),
                           prob = 0.95, adjust = 1.5) +
    ggtitle(model)

    phi <- mcmc_areas(as.matrix(stanfit, pars = "phi"),
                      prob = 0.95, adjust = 1.5) +
    ggtitle(model)

    sig.alpha <- mcmc_areas(as.matrix(stanfit, pars = "sigma_alpha"),
                            prob = 0.95, adjust = 1.5) +
    ggtitle(model)


    if(grepl("hbm2", model)) {
        mu.gamma <- mcmc_areas(as.matrix(stanfit, pars = "mu_gamma"),
                               prob = 0.95, adjust = 1.5) +
        ggtitle(model)
        sig.gamma <- mcmc_areas(as.matrix(stanfit, pars = "sigma_gamma"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(model)
    }

    if(grepl("hbm3", model)) {
        mu.gamma1 <- mcmc_areas(as.matrix(stanfit, pars = "mu_gamma1"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(paste0(model), "(early)")

        mu.gamma2 <- mcmc_areas(as.matrix(stanfit, pars = "mu_gamma2"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(paste0(model), "(mid)")

        mu.gamma3 <- mcmc_areas(as.matrix(stanfit, pars = "mu_gamma3"),
                                prob = 0.95, adjust = 1.5) +
        ggtitle(paste0(model, "(late)"))

        sig.gamma <- mcmc_areas(as.matrix(stanfit, pars = "sigma_gamma"),
                               prob = 0.95, adjust = 1.5) +
        ggtitle(paste0(model))
    }

    if(ssm) {
        sig.gamma <- mcmc_areas(as.matrix(stanfit, pars = "sigma_gamma"),
                               prob = 0.95, adjust = 1.5) +
        ggtitle(model)
    }


    if(!is.null(pdf.file)) pdf(pdf.file, width = width, height = height)
        print(mu.alpha)
        print(phi)
        print(sig.alpha)
        if(grepl("hbm2", model)) print(mu.gamma)
        if(grepl("hbm2", model)) print(sig.gamma)
        if(grepl("hbm3", model)) print(mu.gamma1)
        if(grepl("hbm3", model)) print(mu.gamma2)
        if(grepl("hbm3", model)) print(mu.gamma2)
        if(grepl("hbm3", model)) print(sig.gamma)
        if(ssm) print(sig.gamma)
    if(!is.null(pdf.file)) dev.off()

}


## plot_hbm_dot --------------------------------------------
plot_hbm_dot <- function(stanfit,
                         pdf.file = NULL,
                         width = 6,
                         height = 10) {
    ## HMB posterior dot plots
    ##
    ## CI inner = 80%
    ## CI outer = 95%
    ##
    ## stanfit = stanfit model object
    ## pdf.file = path to pdf file to save
    ## width = pdf file width
    ## height = pdf file height

    model <- deparse(substitute(stanfit))
    ssm <- grepl("hbm4", model) |
           grepl("hbm5", model) |
           grepl("hbm6", model)

    p.alpha <- mcmc_intervals(as.matrix(stanfit, pars = c("alpha", "mu_alpha")),
                              prob = 0.8, prob_outer = 0.95) +
    ggtitle(model)

    p.beta <- mcmc_intervals(as.matrix(stanfit, pars = "beta"),
                             prob = 0.8, prob_outer = 0.95) +
    ggtitle(model)

    p.sigma <- mcmc_intervals(as.matrix(stanfit, pars = "sigma"),
                              prob = 0.8, prob_outer = 0.95) +
    ggtitle(model)

    if(grepl("hbm2", model)) {
        p.gamma <- mcmc_intervals(as.matrix(stanfit, pars = c("gamma", "mu_gamma")),
                                  prob = 0.8, prob_outer = 0.95) +
        ggtitle(model)
    }

    if(grepl("hbm3", model)) {
        g.early <- mcmc_intervals(as.matrix(stanfit, pars = c("gamma1", "mu_gamma1")),
                                  prob = 0.8, prob_outer = 0.95) +
        ggtitle(paste0(model), "(early)")

        g.mid <- mcmc_intervals(as.matrix(stanfit, pars = c("gamma2", "mu_gamma2")),
                                prob = 0.8, prob_outer = 0.95) +
        ggtitle(paste0(model), "(mid)")

        g.late <- mcmc_intervals(as.matrix(stanfit, pars = c("gamma3", "mu_gamma3")),
                                 prob = 0.8, prob_outer = 0.95) +
        ggtitle(paste0(model, "(late)"))
    }

    if(ssm) {
        g.sn <- mcmc_intervals(as.matrix(stanfit, pars = "signal_noise"),
                               prob = 0.8, prob_outer = 0.95) +
        ggtitle(model)
    }


    if(!is.null(pdf.file)) pdf(pdf.file, width = width, height = height)
        print(p.alpha)
        print(p.beta)
        print(p.sigma)
        if(grepl("hbm2", model)) print(p.gamma)
        if(grepl("hbm3", model)) print(g.early)
        if(grepl("hbm3", model)) print(g.mid)
        if(grepl("hbm3", model)) print(g.late)
        if(ssm) print(g.sn)
    if(!is.null(pdf.file)) dev.off()

}



## plot_post_pc --------------------------------------------
plot_post_pc <- function(stanfit, y, data = sock.covar,
                         pdf.path = NULL,
                         var.yrep = "yrep") {
    ## Plot posterior predictive distributions
    ##
    ## stanfit = stanfit object sampled using priors only (no likelihood)
    ## y = response variable used to fit model
    ## pdf.path = file path to save graphics
    ## var.yrep = name of yrep parameter in stanfit

    yrep <- rstan::extract(stanfit, pars = var.yrep)

    if(!is.null(pdf.path)) {
        pdf(pdf.path, width = 10, height = 8)
    }

    g <- ppc_scatter_avg(y, yrep[[1]]) +
        ggtitle("Observed vs. predicted")
    print(g)

    g <- ppc_stat_2d(y, yrep[[1]]) +
        ggtitle("Y rep: target distributions")
    print(g)

    g <- ppc_dens_overlay(y = y, yrep = yrep[[1]][1:100, ]) +
        ggtitle("Y rep: posterior predictive check")
    print(g)

    g <- ppc_dens_overlay(y = y, yrep = yrep[[1]][1:50, ]) +
        facet_wrap( ~ data$stock, as.table = FALSE) +
        ggtitle("Y rep: posterior predictive check")
    print(g)

    if(!is.null(pdf.path)) {
        dev.off()
    }
}


## plot_prior_pc -------------------------------------------
plot_prior_pc <- function(stanfit, y,
                          pdf.path = NULL,
                          var.yrep = "yrep") {
    ## Plot prior predictive distributions
    ##
    ## stanfit = stanfit object sampled using priors only (no likelihood)
    ## y = response variable used to fit model
    ## pdf.path = file path to save graphics
    ## var.yrep = name of yrep parameter in stanfit

    yrep <- rstan::extract(stanfit, pars = var.yrep)

    if(!is.null(pdf.path)) {
        pdf(pdf.path, width = 14, height = 10)
    }

    g <- ppc_dens_overlay(y = y, yrep = yrep[[1]][1:100, ]) +
        ggtitle("Y rep: prior predictive distributions")
    print(g)

    g <- ppc_scatter(y, yrep[[1]][1:4, ]) +
        ggtitle("Observed vs. prior simulated realization")
    print(g)

    if(!is.null(pdf.path)) {
        dev.off()
    }
}



## theme_sleek ---------------------------------------------
# see: https://github.com/seananderson/ggsidekick
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", size = rel(1.1)),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}



## stan utilities ------------------------------------------
rhat_highest <- function(stanfit, k = 4, pars) {
    rhat <- get_rhat(stanfit, pars = pars)
    rhat.max <- rev(sort(rhat)[(length(rhat) - k):length(rhat)])
    return(rhat.max)
}

neff_lowest <- function(stanfit, k = 4, pars) {
    neff <- get_neff(stanfit, pars = pars)
    neff.min <- sort(neff)[1:k]
    return(neff.min)
}

pairs_lowest <- function(stanfit, k = 4, pars) {
    n <- get_neff(stanfit, pars = pars)
    n.min <- names(sort(n))[1:k]
    pairs(stanfit, pars = n.min)
}

get_rhat <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "Rhat"]
}

get_neff <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "n_eff"]
}

total_draws <- function(stanfit) {
    ## N chains * N draws  -- post warmup
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    dim(stanfit)[1] * dim(stanfit)[2]
}



## stan_data -----------------------------------------------
stan_data <- function(data,
                      var.x2 = "sst_anom_stnd",
                      breakpoint1 = 1989,
                      breakpoint2 = NULL,
                      scale.x1 = FALSE,
                      prior.sigma.gamma = c(df = 3, mu = 0, sd = 0.5),
                      priors.only = FALSE) {
    ## Get list of data for input to Stan
    ##
    ## data = data.frame of salmon data
    ## var.x2 = column name in `data` of x2 variable
    ## breakpoint1 = year for start of second era
    ## breakpoint2 = year for start of third era, if NULL, only two eras are
    ##               returned defined by breapoint1
    ## scale.x1 = logical, should the x1 (scaled) be scaled to N(0, 1)
    ## prior.sigma.gamma = Student-t prior parameters for sigma_gamma
    ## priors.only = logical indicating if a likelihood should be calculated
    ##               TRUE indicates prior predictive distributions will be
    ##               sampled only

    sock.stan <- data

    ## Set factor levels for ocean_region
    sock.stan$ocean_region <- factor(sock.stan$ocean_region,
                                     levels = unique(sock.stan$ocean_region))


    ## Get start/end for each stock
    start.end  <- levels_start_end(sock.stan$stock)

    ## Get grouping factor for series
    grp.df <- plyr::ddply(sock.stan, .(stock), summarize,
                          group = unique(ocean_region))
    g.grp <- as.numeric(factor(grp.df$group, levels = unique(grp.df$group)))

    a.grp.df <- plyr::ddply(sock.stan, .(stock), summarize,
                            group = unique(region))
    a.grp <- ifelse(a.grp.df$group == "Fraser River", 1, 2)


    ## Get start/end for group-specific gamma series
    start.end.grp.lst <- lapply(split(sock.stan, sock.stan$ocean_region),
                                function(x) min(x$brood_yr):max(x$brood_yr))
    start.end.grp.vec <- unlist(lapply(seq_along(start.end.grp.lst),
                                       function(x) rep(x, length(start.end.grp.lst[[x]]))))
    start.end.grp <- levels_start_end(as.factor(start.end.grp.vec))


    ## Get year indices (map gamma -> y)
    year.lst  <- lapply(split(sock.stan, sock.stan$ocean_region),
                        function(x) as.numeric(as.factor(x$brood_yr)))
    if(length(unique(sock.stan$ocean_region)) == 1) {
        year <- year.lst[[1]]
    } else {
        year <- c(year.lst$WC, year.lst$GOA + max(year.lst$WC))
        year <- c(year, year.lst$BS + max(year))
    }

    if(scale.x1) {
        x1 <- plyr::ddply(data, .(stock), transform,
                          x1 = scale(spawners)[ , 1])$x1
    } else {
        x1 = data$spawners
    }

    if(is.null(breakpoint2)) {
        era1 <- ifelse(sock.stan$brood_yr < breakpoint1, 1, 0)
        era2 <- ifelse(sock.stan$brood_yr >= breakpoint1, 1, 0)
        era3 <- NULL
    } else {
        era1 <- ifelse(sock.stan$brood_yr < breakpoint1, 1, 0)
        era2 <- ifelse(sock.stan$brood_yr >= breakpoint1 &
                       sock.stan$brood_yr < breakpoint2, 1, 0)
        era3 <- ifelse(sock.stan$brood_yr >= breakpoint2, 1, 0)
    }


    lst <- list(y = sock.stan$ln_rps,
                x1 = x1,
                x2 = sock.stan[[var.x2]],
                g_group = g.grp,
                a_group = a.grp,
                year = year,
                era1 = era1,
                era2 = era2,
                era3 = era3,
                n_series = length(unique(sock.stan$stock)),
                Ng_groups = length(unique(sock.stan$ocean_region)),
                Na_groups = length(unique(a.grp)),
                y_start = start.end$start,
                y_end = start.end$end,
                g_start = start.end.grp$start,
                g_end = start.end.grp$end,
                Ng = length(start.end.grp.vec),
                priors_only = ifelse(priors.only, 1, 0),
                sigma_gamma_df = prior.sigma.gamma[1],
                sigma_gamma_mu = prior.sigma.gamma[2],
                sigma_gamma_sd = prior.sigma.gamma[3],
                N = nrow(sock.stan))
    out <- Filter(Negate(is.null), lst)
    return(out)
}

if(FALSE) {
    stan_data(sock.covar, "sst_anom_stnd")
    stan_data(sock.covar, "npgo_winter_stnd")
}



## level_start_end -----------------------------------------
levels_start_end <- function(factor) {
    ## Find start and end point of levelss in a factor
    ##
    ## factor = factor (vector)

    n <- as.numeric(factor)
    internal.end <- which(diff(n) != 0)
    internal.start <- internal.end + 1

    end <- c(internal.end, length(n))
    start <- c(1, internal.start)

    list(start = start, end = end)
}

if(FALSE) {
    f <- as.factor(c("a", "a", "a", "b", "b", "c", "c", "c"))
    levels_start_end(f)

    f <- as.factor(c("a", "a", "a", "b", "c", "c", "c"))
    levels_start_end(f)

    levels_start_end(sock.covar$stock)
}



## sim_compare ---------------------------------------------
sim_compare <- function(stanfit, simdata, pars) {
    ## Compare fitted to true parameter values
    ##
    ## stanfit = fitted stan model
    ## simdata = output from sim_rk used to fit model
    ## pars = parameters to compare

    true <- simdata$pars[[pars]]
    mc.sum <- summary(stanfit, pars = pars)$summary
    s1 <- mc.sum[ , "mean"]
    s.lo <- mc.sum[ , "2.5%"]
    s.up <- mc.sum[ , "97.5%"]
    plot(s1, pch = 16, ylim = range(c(s1, s.lo, s.up)), col = "red3", ylab = pars)
    segments(x0 = 1:nrow(mc.sum), y0 = s.lo, y1 = s.up, col = "red3")
    points(true, pch = 16)
}



## sim_data ------------------------------------------------
sim_data <- function(n_series = 20,
                     n_years = 50,
                     mu_alpha = c(1.5, 3.0),
                     sigma_alpha = c(0.1, 0.1),
                     x2 = "none",
                     x2_sigma = 0.5,
                     x2_phi = 0.2,
                     x2_step = FALSE,
                     mu_gamma = NULL,
                     sigma_gamma = NULL,
                     sigma_min = 0.2,
                     sigma_max = 0.5,
                     phi = 0,
                     y_start = rep(1, n_series),
                     y_end = rep(n_years, n_series),
                     seed = NULL,
                     prior.sigma.gamma = c(df = 3, mu = 0, sd = 0.5)) {
    ## Simulate hierarchical linear model data
    ##
    ## Returns a list:
    ##    $stan: list of data for input to Stan
    ##    $data: data.frame of simulated data
    ##    $pars: true parameter values
    ##
    ## n_series = number of time series to generate
    ## n_years = length of time series
    ## mu_alpha = mean of alpha hyperdistribution
    ## sigma_alpha = SD of alpha hyperdistribution
    ## x2 = type of covariate effect:
    ##      "none" = no covariate effect
    ##      "stationary" = single effect for each region
    ##      "era" = three eras for each region
    ##      "nonstationary" = fully time-varying
    ## x2_sigma = SD in x2 series
    ## x2_phi = autocorrelation in x2 series
    ## x2_step = logical, use step function for nonstationary effect
    ## mu_gamma = if x2 == "stationary": vector of hypermeans, one value for
    ##                                   each gamma grouping
    ##            if x2 == "era": list of length 3 with each list component
    ##                            being a vector of hypermeans, one value for
    ##                            each gamma grouping
    ##            if x2 == "nonstationary": vector of starting values for each
    ##                                      series
    ## sigma_gamma = if x2 == "stationary": vector of hyperSDs, one value for
    ##                                      each gamma grouping
    ##               if x2 == "era": vector of hyperSDs, one value for
    ##                         each gamma grouping
    ##               if x2 == "nonstationary": vector of SDs for non-stationary
    ##                         effect
    ## sigma_gamma = list of SD gamma hyperdistribution
    ## sigma_min = min residual SD
    ## sigma_max = max residual SD
    ## phi = autocorrelation parameter (residuals)
    ## y_start = vector of start points for each series
    ## y_end = vector of end points for each series
    ## seed = numeric, if non-null, set RNG seed

    if(!is.null(seed))
        set.seed(seed)

    ## Safety checks
    if(!x2 %in% c("none", "stationary", "era", "nonstationary"))
        stop("Don't recognize x2 input")
    if(length(y_start) != n_series)
        stop("y_start not equal to n_series")
    if(length(y_end) != n_series)
        stop("y_end not equal to n_series")

    ## mu_gamma defaults
    if(is.null(mu_gamma) & x2 == "none")
        mu_gamma <- c(0, 0)
    if(is.null(mu_gamma) & x2 == "stationary")
        mu_gamma <- c(-2, 2)
    if(is.null(mu_gamma) & x2 == "era")
        mu_gamma <- list(era1 = c(-2, 2),
                         era2 = c(-1, 1),
                         era3 = c(1, 3))
    if(is.null(mu_gamma) & x2 == "nonstationary")
        mu_gamma <- runif(n_series, -1, 1)

    ## sigma_gamma defaults
    if(is.null(sigma_gamma) & x2 == "none")
        sigma_gamma <- c(0, 0)
    if(is.null(sigma_gamma) & x2 == "stationary")
        sigma_gamma <- c(0.2, 0.4)
    if(is.null(sigma_gamma) & x2 == "era")
        sigma_gamma <- c(0.2, 0.4)
    if(is.null(sigma_gamma) & x2 == "nonstationary")
        sigma_gamma <- runif(n_series, 0.1, 0.4)

    ## Alpha groupings
    Na_groups <- length(mu_alpha)
    a_group <- cut(1:n_series, Na_groups, labels = FALSE)

    ## Gamma groupings
    if(is.list(mu_gamma))
        Ng_groups <- length(mu_gamma[[1]])
    else
        Ng_groups <- length(mu_gamma)
    Ng_groups <- ifelse(Ng_groups == 1, n_series, Ng_groups)
    g_group <- cut(1:n_series, Ng_groups, labels = FALSE)

    beta   <- runif(n_series, -1, -0.001)
    sigma  <- runif(n_series, sigma_min, sigma_max)

    ## Setup eras
    era <- rep(cut(1:n_years, 3, labels = FALSE), n_series)
    era1 <- ifelse(era == 1, 1, 0)
    era2 <- ifelse(era == 2, 1, 0)
    era3 <- ifelse(era == 3, 1, 0)
    era1_mat <- matrix(era1, ncol = n_years, nrow = n_series, byrow = TRUE)
    era2_mat <- matrix(era2, ncol = n_years, nrow = n_series, byrow = TRUE)
    era3_mat <- matrix(era3, ncol = n_years, nrow = n_series, byrow = TRUE)

    ## Storage
    x1_mat      <- matrix(NA, ncol = n_years, nrow = n_series)
    x2_mat      <- matrix(NA, ncol = n_years, nrow = n_series)
    y_mat       <- matrix(NA, ncol = n_years, nrow = n_series)
    gamma_mat   <- matrix(NA, ncol = n_years, nrow = n_series)
    series_mat  <- matrix(NA, ncol = n_years, nrow = n_series)
    yind        <- matrix(NA, ncol = n_years, nrow = n_series)
    a_group_mat <- matrix(NA, ncol = n_years, nrow = n_series)
    g_group_mat <- matrix(NA, ncol = n_years, nrow = n_series)
    alpha  <- rep(NA, n_series)
    gamma  <- rep(NA, n_series)
    gamma1 <- rep(NA, n_series)
    gamma2 <- rep(NA, n_series)
    gamma3 <- rep(NA, n_series)

    ## Loop over each series
    for(i in 1:n_series) {
        alpha[i] <- rnorm(1, mean = mu_alpha[a_group[i]],
                          sd = sigma_alpha[a_group[i]])
        x1.i <- rnorm(n_years)
        x2.i <- ts_ar1(n_years, 0, x2_phi, sd = x2_sigma)
        res  <- ts_ar1(n_years, 0, phi, sd = sigma[i])

        if(x2 == "none") {
            gamma.i <- rep(0, n_years)
        }
        if(x2 == "stationary") {
            gamma[i] <- rnorm(1, mean = mu_gamma[g_group[i]],
                              sd = sigma_gamma[g_group[i]])
            gamma.i <- rep(gamma[i], n_years)
        }
        if(x2 == "era") {
            gamma1[i] <- rnorm(1, mean = mu_gamma[[1]][g_group[i]],
                               sd = sigma_gamma[g_group[i]])
            gamma2[i] <- rnorm(1, mean = mu_gamma[[2]][g_group[i]],
                               sd = sigma_gamma[g_group[i]])
            gamma3[i] <- rnorm(1, mean = mu_gamma[[3]][g_group[i]],
                               sd = sigma_gamma[g_group[i]])
            gamma.i <- rep(NA, n_years)
            gamma.i <- ifelse(era1_mat[i, ] == 1, gamma1[i], gamma.i)
            gamma.i <- ifelse(era2_mat[i, ] == 1, gamma2[i], gamma.i)
            gamma.i <- ifelse(era3_mat[i, ] == 1, gamma3[i], gamma.i)
        }
        if(x2 == "nonstationary") {
            if(x2_step) {
                gamma.ar <- ts_ar1(n_years, mu_gamma[i], 0.3, sd = sigma_gamma[i])
                gamma.st <- ts_step(n_years, k = 3,
                                    up_mean = 2,
                                    down_mean = -2,
                                    up_sd = 0.5,
                                    down_sd = 0.2)
                gamma.i <- gamma.st + gamma.ar
            } else {
                gamma.i <- ts_ar1(n_years, mu_gamma[i], 1, sd = sigma_gamma[i])
            }
        }

        ## Simulate y
        y.i <- alpha[i] + beta[i] * x1.i + gamma.i * x2.i + res

        ## Fill storage matrices
        for(j in y_start[i]:y_end[i]) {
            yind[i, j] <- j
        }
        y_mat[i, ]      <- y.i
        x1_mat[i, ]     <- x1.i
        x2_mat[i, ]     <- x2.i
        gamma_mat[i, ]  <- gamma.i
        series_mat[i, ] <- i
        a_group_mat[i , ] <- a_group[i]
        g_group_mat[i , ] <- a_group[i]
    }

    ## Wrangle matrices to vectors
    y.na     <- matrix(t(y_mat), ncol = 1)[ , 1]
    x1.na    <- matrix(t(x1_mat), ncol = 1)[ , 1]
    x2.na    <- matrix(t(x2_mat), ncol = 1)[ , 1]
    yind.na  <- matrix(t(yind), ncol = 1)[ , 1]
    gamma.na <- matrix(t(gamma_mat), ncol = 1)[ , 1]
    series.na  <- matrix(t(series_mat), ncol = 1)[ , 1]
    a_group.na <- matrix(t(a_group_mat), ncol = 1)[ , 1]
    g_group.na <- matrix(t(g_group_mat), ncol = 1)[ , 1]

    ## Create output data.frame
    df <- data.frame(year = yind.na,
                     series = series.na,
                     y = y.na,
                     x1 = x1.na,
                     x2 = x2.na,
                     gamma_true = gamma.na,
                     a_group = a_group.na,
                     g_group = g_group.na,
                     era = era)
    df <- df[complete.cases(df), ]
    df$series <- as.factor(df$series)
    row.names(df) <- NULL
    N <- nrow(df)

    ## Stan data list
    st.en <- levels_start_end(df$series)
    lst <- list(y = df$y,
                x1 = df$x1,
                x2 = df$x2,
                era1 = era1,
                era2 = era2,
                era3 = era3,
                n_series = n_series,
                y_start = st.en$start,
                y_end = st.en$end,
                a_group = a_group,
                g_group = g_group,
                Na_groups = Na_groups,
                Ng_groups = Ng_groups,
                N = N,
                sigma_gamma_df = prior.sigma.gamma[1],
                sigma_gamma_mu = prior.sigma.gamma[2],
                sigma_gamma_sd = prior.sigma.gamma[3],
                priors_only = 0)

    ## True parameter list
    pars <- list(alpha = alpha,
                 beta = beta,
                 mu_alpha = mu_alpha,
                 sigma_alpha = sigma_alpha,
                 mu_gamma = mu_gamma,
                 sigma_gamma = sigma_gamma,
                 sigma = sigma,
                 phi = phi)

    if(x2 == "nonstationary") {
        pars$gamma <- gamma.na[!is.na(gamma.na)]
    }
    if(x2 == "stationary") {
        pars$gamma <- gamma
    }
    if(x2 == "era") {
        pars$gamma1 <- gamma1
        pars$gamma2 <- gamma2
        pars$gamma3 <- gamma3
        pars$mu_gamma <- NULL
        pars$mu_gamma1 <- mu_gamma[[1]]
        pars$mu_gamma2 <- mu_gamma[[2]]
        pars$mu_gamma3 <- mu_gamma[[3]]
    }
    if(x2 == "none") {
        lst$x2 <- NULL
        df$x2  <- NULL
        pars$mu_gamma <- NULL
        pars$sigma_gamma <- NULL
    }

    list(stan = lst, data = df, pars = pars)
}

if(FALSE) {

    sim_data(n_series = 4, x2 = "none")
    sim_data(n_series = 4, x2 = "stationary")
    sim_data(n_series = 4, x2 = "era")
    sim_data(n_series = 4, x2 = "nonstationary")

    sim_data(n_series = 2,
             y_start = c(2, 1),
             y_end = c(48, 50))

}



## kfas_tidy -----------------------------------------------
kfas_tidy <- function(fit) {
    ## Tidy KFAS output
    ##
    ## fit = output from `kfas_fit`

    ## Get states and standard errors
    states <- as.data.frame(fit$mle$alphahat)
    vars <- sqrt(plyr::adply(fit$mle$V, 3, diag, .id = NULL))
    names(vars) <- paste0("se_", names(states))
    df <- cbind(states, vars)


    ## Add stock info
    df$brood_yr <- as.numeric(rownames(fit$model$y))
    df$stock <- colnames(fit$model$y)
    st <- unique(df$stock)
    df$region <- sock.info$region[sock.info$stock == st]
    df$sub_region <- sock.info$sub_region[sock.info$stock == st]
    df$ocean_region <- sock.info$ocean_region[sock.info$stock == st]

    df$stock <- factor(df$stock, levels = unique(df$stock))
    df$region <- factor(df$region, levels = unique(df$region))
    df$sub_region <- factor(df$sub_region, levels = unique(df$sub_region))
    df$ocean_region <- factor(df$ocean_region, levels = unique(df$ocean_region))

    ## Reorg
    vars <- c("stock", "ocean_region", "region", "sub_region", "brood_yr")
    df1 <- df[ , vars]
    df2 <- df[ , !names(df) %in% vars]
    out <- cbind(df1, df2)

    return(out)
}



## kfas_fit ------------------------------------------------
kfas_fit <- function(model, method = "BFGS", tol = .Machine$double.eps^0.5) {
    ## Fit DLM using KFAS
    ##
    ## model = output from `kfas_model`
    ## method = optim method
    mod <- model
    y <- mod$y
    if(mod$m == 2)
        sn <- c("alpha", "beta")
    if(mod$m == 3)
        sn <- c("alpha", "beta", "gamma")
    kf.model <- SSModel(y ~ -1 + SSMcustom(Z = mod$Z, T = mod$T,
                                           a1 = mod$a1,
                                           P1 = mod$P1,
                                           P1inf = mod$P1inf,
                                           R = mod$R, Q = mod$Q,
                                           state_names = sn),
                        H = mod$H,
                        tol = tol)
    inits <- runif(mod$n_param_var, -5, 1) ## on log scale: exp(inits)
    fit <- fitSSM(kf.model, inits = inits, method = method,
                  control = list(maxit = 1000))
    kf <- KFS(fit$model)

    ret <- list(model = model, kf.model = kf.model,
                inits = inits, fit = fit, mle = kf)
    return(ret)
}



## kfas_model ----------------------------------------------
kfas_model <- function(y, x1, x2 = NULL,
                       a1 = 0, P1 = 0, P1inf = 1,
                       Qa = "diagonal and unequal",
                       Qb = "zero",
                       Qg = "zero") {
    ## Setup KFAS model matrices for DLM
    ##
    ## Return a list with two elements:
    ##   model = matrices that can be input into MARSS()
    ##   dims = useful dimensions of model matrices
    ##
    ## y = matrix of response series (in KFAS input format)
    ## x1 = covariate 1, should have same dims as `y`
    ## x2 = covariate 2, should have same dims as `y`
    ## a1 = expected value for initial state
    ## P1 = value for non-diffuse part of initial state variance
    ## P1inf = value for diffuse part of initial state variance
    ## Qa = form of alpha parameter covariance block
    ## Qb = form of beta parameter covariance block
    ## Qg = form of gamma parameter covariance block
    ##
    ## P1 = 0 and P1inf = 1 represents a diffuse initialization
    ## P1 = 1e7 and P1inf = 0 represents an approximate diffuse initialization
    ##
    ## ?SSModel
    ##
    ## n = number of observations
    ## p = number of y time series
    ## m = number of regression params (states)
    ## k = number of disturbances eta
    ##
    ## y = n,p                  response
    ## Z = p,m,n                data array
    ## T = m,m,n  (B in MARSS)  identity
    ## R = m,k,n  (G in MARSS)  identity
    ## H = p,p,n  (R in MARSS)
    ## Q = p,p,n
    ## a1 = m,1
    ## P1 = m,m
    ## P1inf = m,m

    if(!all.equal(dim(y), dim(x1)))
        stop("y and x1 dims are not equal")

    b <- ifelse(is.null(x2), 2, 3) ## no. of regression params per y series
    n <- nrow(y)  ## number of years
    p <- ncol(y)  ## number of y series
    m <- p * b    ## number of states
    n.obs <- sum(!is.na(y)) ## sample size

    T <- diag(1, m)
    R <- diag(1, m)
    H <- matrix(0, nrow = p, ncol = p)
    diag(H) <- NA
    a1 <- matrix(a1, nrow = m)
    P1 <- diag(P1, m)
    P1inf <- diag(P1inf, m)

    ind <- cut(1:m, b, labels = FALSE)
    a.ind <- which(ind == 1)
    b.ind <- which(ind == 2)
    g.ind <- which(ind == 3)

    ## Q block matrices
    ##
    ##  Qa | 0
    ## ---------
    ##  0  | Qb
    ##
    ##  Qa | 0  | 0
    ## ---------|----
    ##  0  | Qb | 0
    ## ---------|----
    ##  0  | 0  | Qg
    Q.0 <- kfas_qblock(m / b, m / b, "zero")
    Q.a <- kfas_qblock(m / b, m / b, Qa)
    Q.b <- kfas_qblock(m / b, m / b, Qb)
    Q.g <- kfas_qblock(m / b, m / b, Qg)

    if(is.null(x2)) {
        Q.az <- rbind(Q.a, Q.0)
        Q.bz <- rbind(Q.0, Q.b)
        Q    <- cbind(Q.az, Q.bz)
    } else {
        Q.az <- rbind(Q.a, Q.0, Q.0)
        Q.bz <- rbind(Q.0, Q.b, Q.0)
        Q.gz <- rbind(Q.0, Q.0, Q.g)
        Q <- cbind(Q.az, Q.bz, Q.gz)
    }

    ## Z array
    Z <- array(0, dim = c(p, m, n))
    for(i in 1:p) {
        Z[i, a.ind[i], ] <- rep(1, n)
        Z[i, b.ind[i], ] <- x1[ , i]
        if(!is.null(x2))
            Z[i, g.ind[i], ] <- x2[ , i]
    }

    ## Number of variance parameters
    n.param.var  <- sum(is.na(H[lower.tri(H, diag = TRUE)])) +
                    sum(is.na(Q[lower.tri(Q, diag = TRUE)]))

    ## Number of initialized states
    ## TODO is this correct? Differs from Minto et al.
    n.state.init <- length(diag(P1inf))
    # n.state.init <- sum(is.na(diag(Q)))

    ## Get number of parameters
    ## The number of state variances is multiplied by 2 to account for the
    ## diffuse initialization. See DK2012 pg. 187--188.
	## K = obs variances + (state error + diffuse states) + constant states
	n.param <- n.param.var + n.state.init

    dims  <- list(m = m,
                  p = p,
                  n = n,
                  n_obs = n.obs,
                  n_param_var = n.param.var,
                  n_state_init = n.state.init,
                  n_param = n.param,
                  a_ind = a.ind,
                  b_ind = b.ind,
                  g_ind = g.ind)

    model <- list(Z = Z, T = T, R = R, Q = Q, H = H,
                  a1 = a1, P1 = P1, P1inf = P1inf)
    return(c(y = list(y), model, dims))
}



## kfas_qblock ---------------------------------------------
kfas_qblock <- function(nrow, ncol, type = "zero") {
    ## Get the Q matrix for a KFAS model
    ##
    ## nrow = number of matrix rows
    ## ncol = number of matrix columns
    ## type = flavor of Q matrix:
    ##   - zero = no components varying
    ##   - diagonal and unequal = separate variances
    ##   - unconstrained" = separate variances and covariances

    at <- c("zero", "diagonal and unequal", "unconstrained")
    if(!type %in% at) {
        stop("type not recognized")
    }

    M <- matrix(0, nrow, ncol)
    n.diag <- length(diag(M))

    if(type == "diagonal and unequal") {
        diag(M) <- NA
    }

    if(type == "unconstrained") {
        M <- matrix(NA, nrow, ncol)
    }

    return(M)

}



## y_matrix ------------------------------------------------
y_matrix <- function(data, var,
                     id.var = "stock",
                     yr.var = "brood_yr",
                     na.value = NA,
                     yr.min = NULL,
                     yr.max = NULL,
                     kfas = FALSE) {
    ## Get data in matrix format for KFAS/MARSS model
    ##
    ## MARSS models require y to be be in [n x t] matrix where [n] is
    ## the number of time series and [t] is the number of years.
    ##
    ## data = data.frame of salmon data
    ## var = column name in `data` of variable to extract
    ## id.var = column name in `data` of unique series identifier
    ## yr.var = column name in `data` of year variable
    ## na.value = alternative value to set NA values in `var` to
    ## yr.min = minimum year to include
    ## yr.max = maximum year to include
    ## kfas = logical, put data in format for KFAS

    if(is.null(yr.min))
        yr.min  <- min(data[[yr.var]])
    if(is.null(yr.max))
        yr.max  <- max(data[[yr.var]])

    years   <- yr.min:yr.max
    n.years <- length(years)
    id      <- unique(data[[id.var]])
    n.id    <- length(id)

    y <- matrix(NA, nrow = n.id, ncol = n.years,
                dimnames = list(id, years))

    for(i in seq_along(id)) {
        for(t in seq_along(years)) {
            y0 <- data[[var]][data[[id.var]] == id[i] &
                              data[[yr.var]] == years[t]]
            y0 <- ifelse(length(y0) == 0, NA, y0)
            y[i, t] <- y0
        }
    }
    y[is.na(y)] <- na.value
    if(kfas) {
        y <- t(y)
    }
    return(y)
}

if(FALSE) {

    dat <- sock.covar[sock.covar$region == "Bristol Bay", ]
    y_matrix(dat, "ln_rps", id.var = "stock")
    y_matrix(dat, "spawners", na.value = 1)

    y_matrix(dat, "ln_rps", yr.min = 2000, yr.max = 2004)
}



## load_rdata ----------------------------------------------
load_rdata <- function(path = "./output/", verbose = TRUE) {
    ## Load saved .RData output
    ##
    ## This function loads all *.RData files found in the input directory.
    ##
    ## path = path to directory to read output
    ## verbose = logical, print status

    if(!dir.exists(path))
        stop("Input path doesn't exist")

    fls <- list.files(path = path)

    for(i in fls) {

        if(verbose)
            cat("Reading file", which(fls == i), "of",
                length(fls), "...", "\n")

        load(paste(path, i, sep = ""), envir = .GlobalEnv)
    }

    if(verbose) {
        if(length(fls) == 0)
            cat("No files to read", "\n")
        else
            cat("Done!", "\n")
    }
}



## cor_pp --------------------------------------------------
cor_pp <- function(x, y,
                   n.lags = length(x) - 1,
                   level = 0.95,
                   correction = TRUE,
                   na.rm = FALSE) {
    ## Pearson's correlation significance with autocorrelation
    ##
    ## This function implements the "modified Chelton" method outlined in Pyper
    ## and Peterman (1998) CJFAS to calculate the significance while accounting
    ## for autocorrelation in the input series. The general strategy is to
    ## correct the degrees of freedom used in the significance test to account
    ## for autocorrelation.
    ##
    ## x = numeric vector
    ## y = numeric vector same length as x
    ## n.lags = number of lags to use for calculating the effective degrees of
    ##          freedom, defaults to all lags
    ## alpha = significance level
    ## correction = logical, if FALSE no correction for autocorrelation is
    ##              performed (identical to cor.test)
    ## na.rm = if TRUE, NA values are removed pairwise in x,y

    alpha <- 1 - level

    if(length(x) != length(y)) {
        stop("x and y lengths not equal")
    }

    if(na.rm) {
        mat <- matrixc(x, y, ncol = 2)
        mat <- mat[complete.cases(mat), ]
        x <- mat[ , 1]
        y <- mat[ , 2]
    }

    x.na <- sum(is.na(x))
    y.na <- sum(is.na(y))
    if(x.na > 0 | y.na > 0) {
        stop("NA values in x or y")
    }

    N   <- length(x)
    rho <- cor(x, y)

    if(correction) {
        ## autocorrelation for x and y
        x.acf <- acf(x, lag.max = n.lags, plot = FALSE)
        y.acf <- acf(y, lag.max = n.lags, plot = FALSE)
        x.ac  <- x.acf[[1]][2:(n.lags + 1), , 1]
        y.ac  <- y.acf[[1]][2:(n.lags + 1), , 1]

        ## modified Chelton
        n.star <- (1/N + 2/N * sum(x.ac * y.ac))^-1
    } else {
        n.star <- N
    }

    ## get t-dist value
    t.val <- qt(alpha/2, n.star - 2, lower.tail = FALSE)

    ## critical r value (modified Chelton)
    r.crit <- sqrt((t.val^2) / ((t.val^2) + (n.star - 2)))

    ## TRUE == significant
    sig <- abs(rho) > r.crit

    ## Confidence interval
    ## This is the approximate Fisher CI where the standard error
    ## of the correlation is calculated using n.star rather than N,
    ## which inflates the SE based on the level of autocorrelation
    z <- atanh(rho)
    sigma <- 1/sqrt(n.star - 3)
    cint <- z + c(-1, 1) * sigma * qnorm((1 + level)/2)
    cint <- tanh(cint)
    names(cint) <- c("lower", "upper")

    ret <- list(rho = rho,
                n.star = n.star,
                t.val = t.val,
                r.crit = r.crit,
                conf.int = cint,
                sig = sig)
    return(ret)
}

if(FALSE) {

    x1 <- rnorm(100)
    x2 <- rnorm(100)
    cor_pp(x1, x2, level = 0.9)

    cor.test(x1, x2)

    ar1 <- ts_ar1(n = 100, x_0 = 0, rho = 0.9, sd = 0.15)
    ar2 <- ts_ar1(n = 100, x_0 = 0, rho = 0.9, sd = 0.15)
    cor_pp(ar1, ar2)
    cor.test(ar1, ar2)
}



## AICc ----------------------------------------------------
AICc <- function(fit) {
	## Calculate AICc for linear model
	##
	## fit = object with logLik method

    l <- logLik(fit)
    n <- length(fit$residuals)
    p <- attr(l, "df")
    (-2 * l[1]) + ((2 * p * n) / (n - p - 1))
}



## plot_roll_ts --------------------------------------------
plot_roll_ts <- function(formula, data,
                         data.mean = NULL,
                         col.lower = "ci_lower",
                         col.upper = "ci_upper",
                         col.sig = "sig", ...) {
    ## Plot rolling correlation time series by stock
    ##
    ## This function plots rolling correlations/Ricker models for each stock.
    ##
    ## formula = lattice::xyplot() formula,
    ## data = data frame
    ## data.mean = data frame output from cor_stock() or rk_stock()
    ## col.lower = name of column in `data` giving lower CI bound
    ## col.upper = name of column in `data` giving upper CI bound
    ## ... = passed to lattice::xyplot()

    p.dat <- data
    col.ts <- ifelse("cor" %in% names(p.dat), "cor", "estimate")

    g <- xyplot(formula, data = p.dat,
                lower = p.dat[ , col.lower],
                upper = p.dat[ , col.upper],
                sig = p.dat[ , col.sig],
                par.settings = par_mjm(),
                auto.key = list(space = "right"), ...,
                panel = panel.superpose,
                panel.groups = function(x, y, subscripts, col.line,
                                        group.number, upper, lower, sig, ...) {
                    upper <- upper[subscripts]
                    lower <- lower[subscripts]
                    sig   <- sig[subscripts]
                    panel.abline(h = 0, col = "grey50", lty = 2)
                    if(!is.null(data.mean)) {
                        stock <- unique(p.dat$stock[subscripts])
                        cs <- data.mean[data.mean$stock == stock, ]
                        panel.segments(x0 = min(x), x1 = max(x),
                                       y0 = cs[ , col.ts],
                                       y1 = cs[ , col.ts],
                                       col = "grey30",
                                       lwd = 1)
                        panel.rect(xleft = min(x), xright = max(x),
                                   ybottom = cs[ , col.lower],
                                   ytop = cs[ , col.upper],
                                   col = rgb(0, 0, 0, 0.2),
                                   border = NA,
                                   lwd = 0)
                    }
                    panel.xyplot(x, y, col = col.line, type = "l",
                                 lwd = 2)
                    panel.segments(y0 = lower, y1 = upper,
                                   x0 = x, x1 = x, col = col.line,
                                   lwd = 0.5)
                    panel.xyplot(x[sig], y[sig], type = "p",
                                 col = "red3", cex = 0.55)
                })
    print(g)

}

if(FALSE) {

    plot_roll_ts(cor ~ start | stock, data = cor.roll.sst,
                 groups = ocean_region,
                 data.mean = cor.sst,
                 ylim = c(-1, 1))

    plot_roll_ts(estimate ~ start | stock, data = rk.roll.sst,
                 groups = ocean_region,
                 data.mean = rk.sst,
                 ylim = range(c(rk.roll.sst[ , "2.5%"],
                                rk.roll.sst[ , "97.5%"])))

}



## rk_roll_getparam ----------------------------------------
rk_roll_getparam <- function(data, param, var = NULL) {
    ## Extract a parameter from output of rk_roll_stock()
    ##
    ## This function extracts a single parameter from the output of the
    ## rk_roll_stock() function. The function returns a single "long" data
    ## frame.
    ##
    ## data = output from rk_roll_stock
    ## param = string identifier which parameter to extract
    ## var = optional string to add as a variable indicator column

    params <- names(data[[1]])
    if(!param %in% params)
        stop("`param` not found in `data`")

    lst <- lapply(data, function(x) x[[param]])
    out <- plyr::rbind.fill(lst)
    if(!is.null(var)) {
        out$var <- var
    }
    return(out)
}

if(FALSE) {

    xx <- rk_roll_stock(sock.covar, ln_rps ~ spawners + sst_anom_stnd, 15)
    pp <- rk_roll_getparam(xx, "sst_anom_stnd")
    head(pp)
    tail(pp)

    pp <- rk_roll_getparam(xx, "(Intercept)", var = "test")
    head(pp)
    tail(pp)

}



## rk_roll_stock -------------------------------------------
rk_roll_stock <- function(data, formula, window,
                          correlation = NULL,
                          level = 0.95) {
    ## Compute rolling lm models for each salmon stock
    ##
    ## This function calculates rolling Ricker models for each stock in the
    ## input data. The function returns a list with one element per
    ## stock where each element contains a list of data frames, one data frame
    ## per parameter in the model.
    ##
    ## data = data frame of salmon data
    ## formula = linear model formula
    ## window = length of sliding window in which to fit model
    ## correlation = error correlation structure; see ?gls
    ## level = level for confidence intervals; see ?confint

    stocks <- unique(data$stock)
    lst <- vector("list", length(stocks))
    for(i in seq_along(stocks)) {
        dat.i <- data[data$stock == stocks[i], ]
        roll.i <- rk_roll(dat.i, formula, window = window,
                          ind = "brood_yr",
                          correlation = correlation,
                          level = level)
        roll.i <- lapply(roll.i, function(x) {
               df1 <- data.frame(stock_id = rep(unique(dat.i$stock_id), nrow(x)),
                                 stock_no = rep(unique(dat.i$stock_no), nrow(x)),
                                 stock = rep(unique(dat.i$stock), nrow(x)),
                                 ocean_region = rep(unique(dat.i$ocean_region),
                                                    nrow(x)),
                                 region = rep(unique(dat.i$region),
                                              nrow(x)),
                                 sub_region = rep(unique(dat.i$sub_region),
                                                  nrow(x)))
               df2 <- cbind(df1, x)
                return(df2)
             })
        lst[[i]] <- roll.i
    }
    names(lst) <- stocks
    return(lst)
}

if(FALSE) {

    xx <- rk_roll_stock(sock.covar, ln_rps ~ spawners + sst_anom_stnd, 15)
    head(xx)
    tail(xx)

}



## rk_roll -------------------------------------------------
rk_roll <- function(data, formula, window, ind = NULL,
                    correlation = NULL,
                    level = 0.95) {
    ## Compute rolling lm models
    ##
    ## This function computes rolling lm models using a sliding
    ## window size of `window`
    ##
    ## The function outputs a list of data frames with one data frame for each
    ## parameter in `formula` with columns:
    ##  $start = indicator of first data point used
    ##  $end = indicator of last data point used
    ##  $estimate = parameter estimate
    ##  $ci_lower = lower confidence bound for parameter
    ##  $ci_upper = upper confidence bound for parameter
    ##
    ## data = a data frame
    ## formula = a formula for a linear model
    ## window = length of sliding window in which to compute correlations
    ## ind = column name in `data`. If supplied it is used to
    ##       populate the `start` and `end` columns in the output, often a year
    ##       column is supplied here
    ## correlation = error correlation structure; see ?gls
    ## level = level for confidence intervals; see ?confint

    dat <- data

    if(window > nrow(dat))
        stop("window longer than input data")

    if(window < 2)
        stop("window must be greater than 1")

    n.dat <- nrow(dat)
    n.cor <- n.dat - window + 1
    n.end <- window - 1

    lst <- vector("list", n.cor)

    for(i in 1:n.cor) {
        if(!is.null(ind)) {
            start <- data[i, ind]
            end   <- data[i + n.end, ind]
        } else {
            start <- i
            end   <- i + n.end
        }
        dat.i <- dat[i:(i + n.end), ]
        lm.i  <- nlme::gls(formula, data = dat.i,
                           correlation = correlation,
                           method = "ML")

        coef.i <- coef(lm.i)
        conf.i <- confint(lm.i)

        df.i <- data.frame(start = start,
                           end = end,
                           mid = mean(start:end),
                           parameter = names(coef.i),
                           estimate = coef.i)
        row.names(df.i) <- NULL
        df.i[ , "ci_lower"] <- conf.i[ , 1]
        df.i[ , "ci_upper"] <- conf.i[ , 2]
        df.i$sig <- (df.i$ci_lower < 0 & df.i$ci_upper < 0) |
                    (df.i$ci_lower > 0 & df.i$ci_upper > 0)
        lst[[i]] <- df.i
    }
    df.coef <- plyr::rbind.fill(lst)
    names(df.coef) <- gsub(" ", "", names(df.coef))

    ## split out parameters into data frames
    out <- split(df.coef, df.coef$parameter)

    ## Reset row.names
    out <- lapply(out, function(r) {
                  row.names(r) <- NULL
                  return(r)
                })

    return(out)
}

if(FALSE) {

    rk_roll(sock.covar[sock.covar$stock_id == 102, ],
            ind = "brood_yr", ln_rps ~ spawners + sst_anom_stnd, 15)

    rk_roll(sock.covar[sock.covar$stock_id == 102, ],
            ln_rps ~ spawners + sst_anom_stnd, 15, ind = NULL)

}



## rk_stock ------------------------------------------------
rk_stock <- function(data, formula, level = 0.95,
                     correlation = NULL,
                     var = NULL) {
    ## Fit stationary Ricker model for each stock
    ##
    ## This function fits a stationary Ricker model for each stock and returns a
    ## list of data frames with one list element for each parameter in the
    ## model. The `nlme::gls` function is used to fit the models to allow for
    ## correlated error structures.
    ##
    ## data = a data frame
    ## formula = a formula for a linear model
    ## level = level for confidence intervals; see ?confint
    ## correlation = error correlation structure; see ?gls
    ## var = optional string to add as a variable indicator column

	lo <- (1 - level) / 2
	up <- 1 - lo

	## Fit models
	fits <- plyr::dlply(data, .(stock_id), function(x) {
				  lm.i <- nlme::gls(formula, data = x,
                                    correlation = correlation,
                                    method = "ML")
	})

    ## Get stock info
	info <- plyr::ddply(data, .(stock_id), summarize,
						stock = unique(stock),
						ocean_region = unique(ocean_region),
						region = unique(region),
						sub_region = unique(sub_region))

    ## Summarize coefficients
	cs <- plyr::ldply(fits, function(lm.i) {
					  coef.i <- coef(lm.i)
					  conf.i <- confint(lm.i, level = level)
					  fit.i <- data.frame(parameter = names(coef.i),
										  estimate = coef.i)
					  fit.i[ , "ci_lower"] <- conf.i[ , 1]
					  fit.i[ , "ci_upper"] <- conf.i[ , 2]
                      fit.i$sig <- (fit.i$ci_lower < 0 & fit.i$ci_upper < 0) |
                                  (fit.i$ci_lower > 0 & fit.i$ci_upper > 0)
					  return(fit.i)
				  })
	cs.info <- join(info, cs, by = "stock_id")
    if(!is.null(var)) {
        cs.info$var <- var
    }
	cs.ret  <- split(cs.info, cs$parameter)

    ## Summary info
	sm <- plyr::ldply(fits, function(lm.i) {
		              data.frame(sigma = summary(lm.i)$sigma,
		              		     loglik = logLik(lm.i)[1],
		              		     n_params = attr(logLik(lm.i), "df"),
		              		     AIC = AIC(lm.i),
		              		     AICc = AICc(lm.i))
	              })
	sm.info <- join(info, sm, by = "stock_id")

	ret <- c(cs.ret, list(summary = sm.info))
    return(ret)
}

if(FALSE) {

    xx <- rk_stock(sock.covar, ln_rps ~ spawners + sst_anom_stnd)
    xx
    split(xx, xx$parameter)
    head(xx)
    tail(xx)

}



## cor_roll_stock ------------------------------------------
cor_roll_stock <- function(data, x, y, window,
                           level = 0.95,
                           ar1.correction = TRUE,
                           n.lags = 5,
                           var = NULL) {
    ## Compute rolling correlations for each salmon stock
    ##
    ## This function calculates rolling correlations between x and y for each
    ## stock in the input data.
    ##
    ## data = data frame of salmon data
    ## x = column name in `data` of first variable
    ## y = column name in `data` of second variable
    ## window = length of sliding window in which to compute correlations
    ## level = significance level for confidence intervals
    ## var = optional string to add as a variable indicator column
    ## ar1.correction = correct for autocorrelation using
    ##                  Pyper & Peterman method
    ## n.lags = number of lags to include in Pyper & Peterman correction

    crs <- plyr::ddply(data, .(stock_id), function(d) {
                       cr <- cor_roll(d[ , x],
                                      d[ , y],
                                      window = window,
                                      ind = d$brood_yr,
                                      level = level,
                                      ar1.correction = ar1.correction,
                                      n.lags = n.lags)
                       df <- data.frame(stock_id = rep(unique(d$stock_id),
                                                       nrow(cr)),
                                        stock = rep(unique(d$stock), nrow(cr)),
                                        ocean_region = rep(unique(d$ocean_region),
                                                           nrow(cr)),
                                        region = rep(unique(d$region), nrow(cr)),
                                        sub_region = rep(unique(d$sub_region),
                                                         nrow(cr)))
                       out <- cbind(df, cr)
                       return(out)
                      })
    if(!is.null(var)) {
        crs$var <- var
    }
    return(crs)
}

if(FALSE) {
    xx <- cor_roll_stock(sock.covar, "rk_resid", "sst_anom.stnd", 15)
    head(xx)
}



## cor_roll ------------------------------------------------
cor_roll <- function(x, y, window,
                     ind = NULL,
                     level = 0.95,
                     ar1.correction = TRUE,
                     n.lags = 5) {
    ## Compute rolling correlations
    ##
    ## This function computes the rolling correlation between x and y using a
    ## sliding window size of `window`. Missing data are not allowed in x or y.
    ## The function outputs a data frame with columns:
    ##  $start = indicator of first data point used in correlation
    ##  $end = indicator of last data point used in correlation
    ##  $cor = correlation coefficient for each window
    ##
    ## x = numeric vector
    ## y = numeric vector the same length as y
    ## window = length of sliding window in which to compute correlations
    ## ind = vector the same length as x and y. If supplied it is used to
    ##       populate the `start` and `end` columns in the output, often a year
    ##       vector is supplied here
    ## level = significance level for confidence intervals
    ## ar1.correction = correct for autocorrelation using
    ##                  Pyper & Peterman method
    ## n.lags = number of lags to include in Pyper & Peterman correction

    if(length(x) != length(y))
        stop("x and y lengths do not match")

    if(any(is.na(x)) | any(is.na(y)))
        stop("missing values not allowed")

    if(!is.null(ind)) {
        if(length(x) != length(ind))
            stop("ind and data lengths do not match")
    }

    if(window > length(x))
        stop("window longer than input data")

    if(window < 2)
        stop("window must be greater than 1")

    if(ar1.correction & n.lags >= window)
        stop("n.lags should be less than window width")

    n.dat <- length(x)
    n.cor <- n.dat - window + 1
    n.end <- window - 1

    lst <- vector("list", n.cor)

    for(i in 1:n.cor) {
        if(!is.null(ind)) {
            start <- ind[i]
            end   <- ind[i + n.end]
        } else {
            start <- i
            end   <- i + n.end
        }
        x.i <- x[i:(i + n.end)]
        y.i <- y[i:(i + n.end)]
        cor.i <- cor(x.i, y.i)
        ct <- cor_pp(x.i, y.i, level = level,
                     correction = ar1.correction,
                     n.lags = n.lags)

        df.i <- data.frame(start = start,
                           end   = end,
                           mid = mean(start:end),
                           cor   = ct$rho)
        df.i[ , "ci_lower"] <- ct$conf.int[1]
        df.i[ , "ci_upper"] <- ct$conf.int[2]
        df.i$sig <- ct$sig
        df.i$n.star <- ct$n.star

        lst[[i]] <- df.i
    }

    out <- plyr::rbind.fill(lst)

    return(out)
}

if(FALSE) {

    cor_roll(rnorm(10), rnorm(10), 6)
    cor_roll(rnorm(10), rnorm(10), 10)
    cor_roll(rnorm(10), rnorm(10), 1)
    cor_roll(c(rnorm(9), NA), rnorm(10), 1)
    cor_roll(rnorm(10), c(rnorm(9), NA), 1)
    cor_roll(rnorm(10), rnorm(10), 6, 2000:2009)

    cor_roll(rnorm(10), rnorm(10), 6, ar1.correction = FALSE)
}



## cor_stock -----------------------------------------------
cor_stock <- function(data, x, y,
                      level = 0.95,
                      ar1.correction = TRUE,
                      n.lags = 5,
                      var = NULL) {
    ## Compute correlations for each salmon stock
    ##
    ## This function calculates correlations between x and y for each stock in
    ## the input data.
    ##
    ## data = data frame of salmon data
    ## x = column name in `data` of first variable
    ## y = column name in `data` of second variable
    ## level = significance level for confidence intervals
    ## ar1.correction = correct for autocorrelation using
    ##                  Pyper & Peterman method
    ## n.lags = number of lags to include in Pyper & Peterman correction
    ## var = optional string to add as a variable indicator column

    cs <- plyr::ddply(data, .(stock_id), function(d) {
                          xx <- d[ , x]
                          yy <- d[ , y]
                          ct <- cor_pp(xx, yy, level = level,
                                       correction = ar1.correction,
                                       n.lags = n.lags)
                          cr <- data.frame(cor = ct$rho)
                          cr[ , "ci_lower"] <- ct$conf.int[1]
                          cr[ , "ci_upper"] <- ct$conf.int[2]
                          cr$samp.size <- length(xx)
                          cr$n.star <- ct$n.star
                          cr$sig <- ct$sig
                          df <- data.frame(stock_id = unique(d$stock_id),
                                           stock = unique(d$stock),
                                           ocean_region = unique(d$ocean_region),
                                           region = unique(d$region),
                                           sub_region = unique(d$sub_region))
                          out <- cbind(df, cr)
                          return(out)
                      })
    if(!is.null(var)) {
        cs$var <- var
    }
    return(cs)
}

if(FALSE) {

    cor_stock(sock.covar,
              x = "rk_resid",
              y = "sst_anom_stnd",
              var = "sst_resid",
              level = 0.95)

    cor_stock(sock.covar,
              x = "rk_resid",
              y = "sst_anom_stnd",
              ar1.correction = FALSE,
              var = "sst_resid",
              level = 0.95)

    cor_stock(sock.covar,
              x = "rk_resid",
              y = "npgo_winter_stnd",
              var = "npgo",
              level = 0.95)

    cor_stock(sock.covar,
              x = "rk_resid",
              y = "npgo_winter_stnd",
              ar1.correction = FALSE,
              var = "npgo",
              level = 0.95)

}



## sst_map -------------------------------------------------
sst_map <- function(data, plot.dir = NULL,
                    sub.years = NULL,
                    sub.months = NULL,
                    progress.bar = TRUE,
                    res = 100,
                    width = 7, height = 6,
                    ...) {
    ## Map SST values using heatmaps
    ##
    ## This functions creates monthly heat maps of SST. One file (map) is
    ## created for each month and year of input SST data.
    ##
    ## data = data.frame of sst data in "long" format
    ## plot.dir = directory to save plots, if null, plot is not saved
    ## sub.years = vector of years to plot (subset of `data`)
    ## sub.months = vector of months to plot (subset of `data`)
    ## progress.bar = logical, should a progress bar be printed
    ## res = resolution for jpeg files
    ## width = width of output jpeg (in inches)
    ## height = height of output jpeg (in inches)
    ## ... = passed to levelplot()

    years  <- unique(data$year)
    months <- unique(data$month)

    if(!is.null(sub.years))
        years <- years[years %in% sub.years]
    if(!is.null(sub.months))
        months <- months[months %in% sub.months]

    if(progress.bar) {
        pb <- txtProgressBar(min = 0,
                             max = length(years) * length(months),
                             style = 3)
    }

    ind <- 1
    for(i in seq_along(years)) {
        for(j in seq_along(months)) {

            dsub <- data[data$year == years[i] & data$month == months[j], ]
            plot.name <- paste("sst_", years[i], "_", months[j], ".jpeg",
                               sep = "")

            if(!is.null(plot.dir))
                jpeg(paste(plot.dir, plot.name, sep = ""),
                     width = width, height = height, res = res,
                     units = "in")

            g <- levelplot(sst_anom ~ lon * lat,
                           data = dsub,
                           cuts = 200,
                           col = "grey60",
                           contour = TRUE,
                           labels = FALSE,
                           main = paste("year = ", years[i],
                                        "   month = ", months[j],
                                        sep = ""),
                           ylab = "Latitude",
                           xlab = "Longitude",
                           ...,
                           panel = function(...) {
                               panel.fill(col = "grey40")
                               panel.levelplot(...)
                               mp <- map('world', fill = TRUE, plot = FALSE)
                               lpolygon(mp$x, mp$y,
                                        col = "white",
                                        border = "grey25")
                           })
            print(g)

            if(!is.null(plot.dir))
                dev.off()

            if(progress.bar)
                setTxtProgressBar(pb, ind)

            ind <- ind + 1
        }
    }
    if(progress.bar)
        close(pb)
}

if(FALSE) {

    cols <- chroma::dpal(500, hue = c(240, 0), chroma = 70, power = 1.0)
    at <- seq(-4.5, 4.5, 0.5)
    sst_map(data = sst.anom,
            plot.dir = NULL,
            sub.years = 2015,
            sub.months = 5,
            progress.bar = FALSE,
            col.regions = cols,
            par.settings = par_mjm(),
            at = at)

}



## make ----------------------------------------------------
make <- function(scripts, remove = NULL) {
    ## Reproduce project
    ##
    ## This function acts as a "makefile" for the project. The function first
    ## removes any specified dynamic files and directories, i.e., files and
    ## directories that are created by the project scripts. The function then
    ## sequentially sources the scripts specified in `scripts`. The function
    ## returns a vector of giving the runtime of the function with elements:
    ##  $hours = number of whole hours
    ##  $minutes = number of whole minutes
    ##  $seconds = number of seconds
    ##
    ## scripts = vector of R script names to source
    ## remove = vector of directories or files to delete prior
    ##          to sourcing files

    cat("Reproducing project...", "\n")

    if(!is.null(remove)) {
        cat("  Removing dynamic files...", "\n")
        for(i in seq_along(remove))
            unlink(remove[i], recursive = TRUE)
    }

    time.start <- proc.time()

    for(i in seq_along(scripts)) {
        cat("  Sourcing ", scripts[i], "...", "\n", sep = "")
        source(scripts[i])
    }

    time.run <- proc.time() - time.start
    time.sec <- time.run[3]
    names(time.sec) <- NULL

    hrs      <- time.sec %/% (60 * 60)
    hrs.rem  <- time.sec %% (60 * 60)
    min      <- hrs.rem %/% 60
    min.rem  <- hrs.rem %% 60
    time.out <- c(hours = hrs, minutes = min, seconds = min.rem)

    cat("Done!", "\n")
    return(time.out)
}

if(FALSE) {

    make(scripts = c("load.R", "functions.R"))
    x <- make(scripts = c("load.R", "functions.R"))
    paste(round(x, 2), names(x), collapse = "   ")

}



## make_log ------------------------------------------------
make_log <- function(dir = "./logs",
                     file = "make.log",
                     session.info = TRUE,
                     run.time = NULL) {
    ## Write a log file for make() calls
    ##
    ## dir = directory to write log file
    ## file = log file name
    ## session.info = logical, should session info be included in log
    ## run.time = optional string giving the run time

    if(!dir.exists(dir))
        dir.create(dir, recursive = TRUE)

    path <- paste0(dir, "/", file)

    ## If log file exists, add a separator at end of file
    if(file.exists(path)) {
        cat("", "###", "", sep = "\n", file = path, append = TRUE)
    }

    ## If log file doesn't exist, create it with a header
    if(!file.exists(path)) {
        hd <- paste0("## Log file for make.R script", "\n")
        write(hd, file = path)
    }

    mess <- paste0("make.R ran successfully on ", Sys.time(), "\n")
    cat(mess, file = path, append = TRUE)

    if(!is.null(run.time)) {
        tim <- paste0("Run time: ", run.time, "\n")
        cat(tim, file = path, append = TRUE)
    }

    if(session.info) {
        sess <- capture.output(sessionInfo())
        cat(sess, file = path, append = TRUE, sep = "\n")
    }
}

if(FALSE) {

    make_log()
    make_log(run.time = "xxx")

}



## wgt_mean_pdo --------------------------------------------
wgt_mean_pdo <- function(sock.data,
                         enviro.data) {
    ## Calculate PDO/NPGO index for each brood year
    ##
    ## This function calculates the PDO index (and NPGO) for each brood year by
    ## taking a weighted average of the PDO over the ocean entry years for a
    ## brood year where the weights are equal to the proportion of fish of that
    ## brood year that enter the ocean in a given year.
    ##
    ## The function returns a vector the same length as nrow(sock.data) giving
    ## the weighted average index.
    ##
    ## sock.data = data.frame of sockeye data
    ## enviro.data = data.frame of environmental data

    vec <- rep(NA, nrow(sock.data))
    for(i in 1:nrow(sock.data)) {

        ## Get current brood year
        yr.i <- sock.data$brood_yr[i]

        ## Get ocean entry age proportions
        wgts <- sock.data[i, grep("^ocean_[[:digit:]]", names(sock.data))]

        ## Get environmental data for years with ocean entries
        env.i <- enviro.data[enviro.data$year %in% seq(yr.i + 1, yr.i + 5, 1), ]

        ## Calculate weighted mean
        ## don't need to divide b/c weights sum to 1
        w.avg <- sum(env.i$index * wgts)
        vec[i] <- w.avg
    }
    return(vec)
}



## wgt_mean_sst --------------------------------------------
wgt_mean_sst <- function(sock.data,
                         enviro.data,
                         sst.var) {
    ## Calculate SST index for each brood year
    ##
    ## This function calculates the SST index for each brood year by
    ## taking a weighted average of SST over the ocean entry years for a
    ## brood year where the weights are equal to the proportion of fish of that
    ## brood year that enter the ocean in a given year.
    ##
    ## The function returns a vector the same length as nrow(sock.data) giving
    ## the weighted average index.
    ##
    ## sock.data = data.frame of sockeye data
    ## enviro.data = data.frame of environmental data
    ## sst.var = name of SST column in `enviro.data`

    vec <- rep(NA, nrow(sock.data))
    for(i in 1:nrow(sock.data)) {

        ## Get current brood year
        yr.i <- sock.data$brood_yr[i]
        id.i <- sock.data$stock_id[i]

        ## Get ocean entry age proportions
        wgts <- sock.data[i, grep("^ocean_[[:digit:]]", names(sock.data))]

        ## Get environmental data for years with ocean entries
        env.id <- enviro.data[enviro.data$stock_id == id.i, ]
        env.i  <- env.id[env.id$year %in% seq(yr.i + 1, yr.i + 5, 1), ]

        ## Calculate weighted mean
        ## don't need to divide b/c weights sum to 1
        w.avg <- sum(env.i[ , sst.var] * wgts)
        vec[i] <- w.avg
    }
    return(vec)
}



## sst_averager --------------------------------------------
sst_averager <- function(info, sst, distance = 400,
                         id.col = "stock_id") {
    ## This function takes as input a data.frame of sst data output from the
    ## sst_anomaly() function and computes sst averages for each stock only
    ## including sst grid cells that are within a specified distance from the
    ## ocean entry location of a particular stock.
    ##
    ## Different months are used to calculate the SST average based on where the
    ## salmon stock enters the ocean:
    ##  WA, BC, SEAK: April-July
    ##  GOA: May-August
    ##  BB, AYK: Jun-September
    ##
    ## The function outputs a data frame with columns:
    ##  $year = one year for each year in input SST data and stock.id
    ##  $stock.id = id number for a salmon stock
    ##  $sst = averaged raw SST values
    ##  $sst.anom = averaged SST anomalies
    ##
    ## Function arguments:
    ##   info = stock.info data.frame w/ stock.id number, lon, and lat, should
    ##          have one row per stock
    ##   sst = sst data output from sst_anomaly()
    ##   distance = distance in km from ocean entry location of stock a grid
    ##              cell can be to be included in the averaging. This distance
    ##              is measured to the center of the SST grid cell
    ## id.col = column name in `info` giving unique stock identifier

    stock.id <- info[ , id.col]
    cells    <- unique(subset(sst, select = c(id, lat, lon)))
    cells    <- cells[order(cells$id), ]
    n.cells  <- length(cells[ , 1])
    row.names(cells) <- NULL

    sst.out <- vector("list", length(stock.id))
    for(i in seq_along(stock.id)) {

        info.i     <- info[i, ]
        stock.id.i <- info.i[ , id.col]
        lat.i      <- info.i$lat
        lon.i      <- info.i$lon

        dist <- rep(NA, n.cells)
        for(j in 1:n.cells)
            dist[j] <- haversine(lat.i, lon.i, cells$lat[j], cells$lon[j],
                                 r = 6378137 / 1e3)

        cells.sub <- cells[which(dist <= distance), ]
        sst.sub   <- sst[sst$id %in% cells.sub$id, ]

        if(stock.id.i <= 137)
            months <- 4:7  ## WA, BC, SEAK

        if(stock.id.i > 137 & stock.id.i <= 152)
            months <- 5:8  ## GOA

        if(stock.id.i > 152)
            months <- 6:9  ## BB and AYK

        sst.sub.mnths <- sst.sub[sst.sub$month %in% months, ]

        sst.avg <- plyr::ddply(sst.sub.mnths, .(year), summarize,
                               sst = mean(sst, na.rm = TRUE),
                               sst_anom = mean(sst_anom, na.rm = TRUE))

        sst.avg$stock.id <- stock.id.i

        sst.out[[i]] <- sst.avg
    }
    sst.out <- plyr::rbind.fill(sst.out)
    return(sst.out)
}



## enviro_avg_months ---------------------------------------
enviro_avg_months <- function(data,
                              first.month,
                              last.month,
                              avg.var,
                              month.var = "month",
                              year.var = "year",
                              grid.id.var = NULL,
                              lat.var = "lat",
                              lon.var = "lon",
                              type = avg.var) {
    ## Compute annual multi-month averages for environmental variables
    ##
    ## Input data should be in a "long" format with a column for year, month,
    ## and the environmental variable to average. `first.month` and
    ## `last.month` give the first and last months of a continuous range to
    ## average the environmental variable over.
    ##
    ## If `first.month` is less than `last.month` the environmental variable is
    ## averaged over the months first.month:last.month within each year. For
    ## example, if `first.month` = 3 and `last.month` = 4, the environmental
    ## variable will be averaged over Mar and Apr for each year.
    ##
    ## If `first.month` equals `last.month`, that month is returned with no
    ## averaging.
    ##
    ## If `first.month` is greater than `last.month` the environmental variable
    ## is averaged for year t starting in `first.month` of year t - 1 and ending
    ## in `last.month` of year t. The output year corresponds to the year
    ## January occurs within the average. For example if `first.month` = 12 and
    ## `last.month` = 3, then the average for the environmental variable will
    ## occur over Dec, Jan, Feb, March and the year is specified by the year
    ## for Jan occurs in.
    ##
    ## The function outputs a data.frame with a `year` column and an `index`
    ## column. When `first.month` is greater than `last.month`, the output
    ## data.frame will have one less year than the input data.frame with no
    ## value for the minimum year in the input data frame.
    ##
    ## If `grid.id.var` is non-null, the averaging is done on a per-grid-cell
    ## basis within a year.
    ##
    ## data = a data.frame
    ## first.month = numeric giving the month to start annual average
    ## last.month = numeric giving the month to stop annual average
    ## avg.var = string, column name in `data` of the variable to average
    ## month.var = string, column name in `data` of the month variable
    ## year.var = string, column name in `data` of the year variable
    ## grid.id.var = string, column name in `data` of the grid cell id column
    ## lon.var = string, column name in `data` of for longitude, only used if
    ##           grid.id.var is non-null
    ## lat.var = string, column name in `data` of for latitude, only used if
    ##           grid.id.var is non-null
    ## type = string, value to set a `type` column in the output data.frame.
    ##        Useful if you want to rbind multiple indices together

    if(!is.data.frame(data))
        stop("Input data is not a data.frame", call. = FALSE)

    if(!first.month %in% 1:12 | !last.month %in% 1:12)
        stop("Months not between 1 and 12", call. = FALSE)

    if(!is.numeric(data[ , month.var]) | !is.numeric(data[ , year.var]))
        stop("Month variable must be numeric", call. = FALSE)

    if(first.month < last.month | first.month == last.month) {
        months <- first.month:last.month
        df <- data[data[ , month.var] %in% months, ]
    }

    if(first.month > last.month) {

        ## Remove months prior to `first.month` in first year and
        ## months after `last.month` in last year
        min.yr <- min(data[ , year.var])
        max.yr <- max(data[ , year.var])
        min.rm <- which(data[ , year.var] == min.yr &
                        data[ , month.var] < first.month)
        max.rm <- which(data[ , year.var] == max.yr &
                        data[ , month.var] > last.month)
        sub <- data[-c(min.rm, max.rm), ]

        ## Remove months not being averaged over
        months <- c(first.month:12, 1:last.month)
        sub2 <- sub[sub[ , month.var] %in% months, ]

        ## Create new year index to average over
        sp <- split(sub2, sub2[ , year.var])
        lst <- lapply(sp, function(x) {
                   x$yr.avg <- ifelse(x[ , month.var] %in% first.month:12,
                                      x[ , year.var] + 1, x[ , year.var])
                   return(x)
               })
        df <- do.call("rbind", c(lst, make.row.names = FALSE))
        df[ , year.var] <- df$yr.avg

    }

    ## Calculate averages
    if(is.null(grid.id.var)) {
        sp.avg <- split(df, df[ , year.var])
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    } else {
        sp.avg <- split(df, list(df[ , grid.id.var], df[ , year.var]))
        lst.avg <- lapply(sp.avg, function(x) {
                    data.frame(year = unique(x[ , year.var]),
                               id = unique(x[ , grid.id.var]),
                               lon = unique(x[ , lon.var]),
                               lat = unique(x[ , lat.var]),
                               index = mean(x[ , avg.var]))
                })
        enviro <- do.call("rbind", c(lst.avg, make.row.names = FALSE))
    }

    enviro$type <- type

    char.months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                     "sep", "oct", "nov", "dec")
    enviro$months <- paste(char.months[first.month],
                           char.months[last.month], sep = "-")
    enviro$months <- as.factor(enviro$months)

    return(enviro)
}

if(FALSE) {

    set.seed(101)
    yr <- c(rep(1950, 24), rep(1951, 24))
    id <- c(rep(1, 12), rep(2, 12), rep(1, 12), rep(2, 12))
    mt <- rep(1:12, 4)
    lt <- c(rep(48, 12), rep(50, 12), rep(48, 12), rep(50, 12))
    ln <- c(rep(120, 12), rep(122, 12), rep(120, 12), rep(122, 12))
    vl <- rnorm(48)
    df <- data.frame(year = yr, month = mt, lat = lt, lon = ln, id = id, value = vl)

    enviro_avg_months(df, 1, 12, "value", grid.id.var = NULL)
    enviro_avg_months(df, 1, 12, "value", grid.id.var = "id")

    enviro_avg_months(df, 1, 3, "value", grid.id.var = NULL)
    enviro_avg_months(df, 1, 3, "value", grid.id.var = "id")

    enviro_avg_months(df, 12, 3, "value", grid.id.var = NULL)
    enviro_avg_months(df, 12, 3, "value", grid.id.var = "id")

    enviro_avg_months(df, 1, 3, "value")
    enviro_avg_months(df, 3, 3, "value")
    enviro_avg_months(df, 12, 3, "value")
    enviro_avg_months(df, 6, 5, "value", group = "test")
}



## sst_anomaly ---------------------------------------------
sst_anomaly <- function(data, ref.years) {
    ## Calculate monthly, per grid cell, anomalies of SST
    ##
    ## Anomalies are calculated as the difference between a grid cell specific
    ## SST value for a given year/month and the long-term monthly mean (defined
    ## by ref.years) for that grid cell. This follows the methods outlined in
    ## Mueter et al. 2002, CJFAS (https://doi.org/10.1139/f02-020).
    ##
    ## data = data.frame of SST data
    ## ref.years = reference years in which to calculate the long-term mean,
    ##             should be a continuous sequence, e.g., 1950:2016

    ## make sure ref.years is a continuous sequence
    if(length(ref.years) != length(min(ref.years):max(ref.years)))
        stop("years vector is not a sequence ascending by 1")

    ## make sure ref.years are available in input data
    if(sum(ref.years %in% data$year) != length(ref.years))
        stop("ref.years not contained in input data")

    ## subset reference years from data
    ref.sst <- data[data$year >= min(ref.years) & data$year <= max(ref.years), ]

    ## calculate monthly long-term mean for each grid cell
    ## NA's are removed in the calculation of long-term mean
    mnth.avg <- aggregate(sst ~ month + id, data = ref.sst,
                          function(x) mean(x, na.rm = TRUE),
                          na.action = na.pass)
    names(mnth.avg)[names(mnth.avg) == 'sst'] <- 'long.avg'
    sst.merge <- merge(data, mnth.avg)
    sst.merge <- sst.merge[order(sst.merge$year,
                                 sst.merge$month,
                                 sst.merge$lat,
                                 sst.merge$lon), ]

    ## calculate sst anomaly
    sst.merge$anom <- sst.merge$sst - sst.merge$long.avg
    row.names(sst.merge) <- NULL
    sst <- data.frame(year = sst.merge$year,
                      month = sst.merge$month,
                      lon = sst.merge$lon,
                      lat = sst.merge$lat,
                      id = sst.merge$id,
                      sst = sst.merge$sst,
                      sst_anom = sst.merge$anom)
    return(sst)
}



## consec_years --------------------------------------------
consec_years <- function(x) {
    ## Find longest consective years
    ##
    ## This function takes as input a vector of years with NA values and finds
    ## the longest set of consecutive years without an NA. The function returns
    ## the vector of consecutive years.
    ##
    ## x = vector of years

    n <- ifelse(is.na(x), 0, x)
    brk <- c(0, which(diff(n) != 1), length(n))
    vec <- lapply(seq(length(brk) - 1),
                  function(i) n[(brk[i] + 1):brk[i+1]])
    consec <- vec[[which.max(lapply(vec, length))]]
    ret <- n[n %in% consec]

    return(ret)
}
if(FALSE) {

    n <- c(NA, 1950:1955)
    consec_years(n)

    n <- c(1949, NA, 1951:1955)
    consec_years(n)

    n <- c(NA, 1951:1955)
    consec_years(n)

    n <- c(1951:1955, NA)
    consec_years(n)

    n <- c(NA, 1950:1955, NA, 1957:1965, NA, 1967)
    consec_years(n)

    n <- c(1950:1955, 1957:1965, NA, 1967:1970)
    consec_years(n)

}



## fill_time_series ----------------------------------------
fill_time_series <- function(data) {
    ## Fill salmon data time series so that all brood years are consecutive,
    ## filling missing years with NA.
    ##
    ## This function takes as input, a brood table with columns `stock_id`,
    ## `stock`, and `brood_yr` and fills in any non-consecutive BY within a
    ## stocks times series with NA values for all other columns present in
    ## `data`.
    ##
    ## When filtering the data by the `use` column, data values in the middle of
    ## the time series for a particular salmon stocks also get removed if `use`
    ## is set to 0. This function adds back in those data points, setting them
    ## to NA.

    id <- unique(data$stock_id)
    lst <- vector("list", length(id))
    for(i in seq_along(id)) {
        sub <- data[data$stock_id == id[i], ]
        brood.yr <- min(sub$brood_yr):max(sub$brood_yr)
        stock.id <- unique(sub$stock_id)
        stock <- unique(sub$stock)
        df <- data.frame(stock_id = stock.id, stock = stock, brood_yr = brood.yr,
                         stringsAsFactors = FALSE)
        lst[[i]] <- merge(df, sub, by = c("stock_id", "stock", "brood_yr"),
                          all.x = TRUE)
    }
    df <- do.call("rbind", c(lst, make.row.names = FALSE))

    ## Don't want NA in these columns
    out <- plyr::ddply(df, .(stock_id), transform,
                       region = unique(na.omit(region)),
                       sub_region = unique(na.omit(sub_region)),
                       ocean_region = unique(na.omit(ocean_region)),
                       lat = unique(na.omit(lat)),
                       lon = unique(na.omit(lon)))
    return(out)
}



## get_npgo ------------------------------------------------
get_npgo <- function(years) {
    ## This function takes as input a range of years and downloads and processes
    ## the NPGO index. The output of the function is a dataframe in 'long'
    ## format with a column for year, month, and the NPGO index.
    ##
    ## years = vector of years

    if(min(years) < 1950)
        stop("Earliest NPGO year is 1950")

    npgo    <- read.table("http://www.o3d.org/npgo/npgo.php", sep = "\t",
                          strip.white = TRUE)
    n.npgo  <- length(npgo[ , 1])
    npgo    <- npgo[4:n.npgo, ]
    n.npgo  <- length(npgo)
    rm.tail <- n.npgo - 3
    npgo    <- npgo[1:rm.tail]
    npgo    <- as.character(npgo)
    npgo    <- strsplit(npgo, "  ")
    npgo    <- do.call("rbind", npgo)
    npgo    <- data.frame(year  = as.numeric(npgo[ , 1]),
                          month = as.numeric(npgo[ , 2]),
                          npgo  = as.numeric(npgo[ , 3]))
    npgo    <- npgo[npgo$year >= min(years) & npgo$year <= max(years), ]

    return(npgo)
}

if(FALSE) {

    get_npgo(1950:1950)
    get_npgo(1950:2013)

    npgo <- get_npgo(1950:2016)
    head(npgo)
    tail(npgo)
    sapply(npgo, class)
    summary(npgo)

}



## get_pdo -------------------------------------------------
get_pdo <- function(years = 1900:2016) {
    ## This function takes as input a range of years and downloads and processes
    ## the PDO index. The output of the function is a dataframe in 'long' format
    ## with a column for year, month, and the PDO index.
    ##
    ## years = vector of years

    if(min(years) < 1900)
        stop("Earliest PDO year is 1900")

    ## Read data
    pdo.r <- readLines("http://jisao.washington.edu/pdo/PDO.latest.txt")

    ## Trim header and footer
    start.line <- grep("^YEAR", pdo.r)
    end.line <- grep(paste0("^", max(years)), pdo.r)
    pdo.s <- pdo.r[start.line:end.line]

    ## Split strings
    pdo.c <- gsub("[/*]", "", pdo.s)
    pdo.c <- gsub("   ", "  ", pdo.c)
    pdo.c <- gsub("    ", "  ", pdo.c)
    pdo.c <- gsub("      ", "  ", pdo.c)
    pdo.c <- strsplit(pdo.c, "  ")
    pdo.c <- do.call("rbind", pdo.c)

    ## Convert to data.frame
    df.names <- trimws(tolower(pdo.c[1, ]), "both")
    pdo.d <- pdo.c[2:nrow(pdo.c), ]
    pdo.d <- apply(pdo.d, 2, function(x) as.numeric(x))
    pdo.d <- as.data.frame(pdo.d)
    names(pdo.d) <- df.names

    ## Reshape into "long" format
    months <- df.names[2:length(df.names)]
    pdo <- reshape(pdo.d, direction = "long",
                   varying = months,
                   v.names = "pdo",
                   times   = months,
                   timevar = "month")
    pdo <- pdo[order(pdo$year), ]
    pdo <- pdo[pdo$year >= min(years) & pdo$year <= max(years), ]
    row.names(pdo) <- NULL
    pdo$id <- NULL

    ## Convert months to numerics
    pdo$month <- ifelse(pdo$month == "jan", 1, pdo$month)
    pdo$month <- ifelse(pdo$month == "feb", 2, pdo$month)
    pdo$month <- ifelse(pdo$month == "mar", 3, pdo$month)
    pdo$month <- ifelse(pdo$month == "apr", 4, pdo$month)
    pdo$month <- ifelse(pdo$month == "may", 5, pdo$month)
    pdo$month <- ifelse(pdo$month == "jun", 6, pdo$month)
    pdo$month <- ifelse(pdo$month == "jul", 7, pdo$month)
    pdo$month <- ifelse(pdo$month == "aug", 8, pdo$month)
    pdo$month <- ifelse(pdo$month == "sep", 9, pdo$month)
    pdo$month <- ifelse(pdo$month == "oct", 10, pdo$month)
    pdo$month <- ifelse(pdo$month == "nov", 11, pdo$month)
    pdo$month <- ifelse(pdo$month == "dec", 12, pdo$month)
    pdo$month <- as.numeric(pdo$month)

    return(pdo)
}

if(FALSE) {

    get_pdo(1950:1950)
    get_pdo(1950:2013)

    pdo <- get_pdo(1900:2016)
    head(pdo)
    tail(pdo)
    sapply(pdo, class)
    summary(pdo)

}
