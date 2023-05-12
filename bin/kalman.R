## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Kalman filter for linear regression intercept or slope parameter
## Michael Malick <mjm@michaelmalick.com>
## 2018-11-09
## Version: 1.0
##
## The functions in this script estimate a linear regression model with either
## a time-varying intercept or slope using the Kalman filter via maximum
## likelihood. Time-varying parameter follows a random walk. Four different
## models can be estimated:
##
## 1. Time varying: mean
##     y[t]  = b0[t] + e[t]             observation equation
##     b0[t] = y[t-1] + w[t]            state equation
##
##     v[t] ~ N(0, sigma_e^2)           observation error
##     w[t] ~ N(0, sigma_w^2)           state error
##
##     Parameters estimated via ML: sigma_e, sigma_w
##
## 2. Time varying: intercept
##     y[t]  = b0[t] + b1*x1[t] + e[t]
##     b0[t] = b0[t-1] + w[t]
##
##     v[t] ~ N(0, sigma_e^2)
##     w[t] ~ N(0, sigma_w^2)
##
##     Parameters estimated via ML: beta, sigma_e, sigma_w
##
## 3. Time varying: slope1
##     y[t]  = b0 + b1[t]*x1[t] + e[t]
##     b1[t] = b1[t-1] + w[t]
##
##     v[t] ~ N(0, sigma_e^2)
##     w[t] ~ N(0, sigma_w^2)
##
##     Parameters estimated via ML: alpha, sigma_e, sigma_w
##
## 4. Time varying: slope2
##     y[t]  = b0 + b1*x1[t] + b2[t]*x2[t] + e[t]
##     b2[t] = b2[t-1] + w[t]
##
##     v[t] ~ N(0, sigma_e^2)
##     w[t] ~ N(0, sigma_w^2)
##
##     Parameters estimated via ML: alpha, beta, sigma_e, sigma_w
##
## Functions:
##  - k_estimate    Model estimation routine
##  - k_objective   Objective function for optimization routine
##  - k_filter      Recursive calculation of Kalman filter
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

k_estimate <- function(y, x1 = NULL, x2 = NULL,
                       varying = "mean",
                       intercept_0 = NULL,
                       slope1_0 = NULL,
                       a_0 = 0,
                       p_0 = 1,
                       ln_sigma_e_0 = -1,
                       ln_sigma_w_0 = -1,
                       Ns = 1,
                       optim_method = "PORT",
                       warnings = 1) {
    ## Kalman filter model estimation routine
    ##
    ## y = see k_filter
    ## x1 = see k_filter
    ## x2 = see k_filter
    ## varying = see k_filter
    ## intercept_0 = initial intercept value
    ## slope1_0 = initial slope1 value
    ## a_0 = see k_filter
    ## p_0 = see k_filter
    ## ln_sigma_e_0 = initial value for log_e(sigma_e) (on log scale)
    ## ln_sigma_w_0 = initial value for log_e(sigma_w) (on log scale)
    ## Ns = see k_filter
    ## optim_method = ML optimization method, "PORT" (the default) uses the
    ##                `nlminb` function, any other value should be a method
    ##                supported by `optim`.
    ## warnings = 0 turns off all warnings, 1 (the default) gives warning when
    ##            convergenge fails
    ##
    ## Returns: a list

    if(varying == "mean") {
        inits <- c(ln_sigma_e = ln_sigma_e_0,
                   ln_sigma_w = ln_sigma_w_0)
    }
    if(varying == "intercept") {
        inits <- c(slope1 = slope1_0,
                   ln_sigma_e = ln_sigma_e_0,
                   ln_sigma_w = ln_sigma_w_0)
    }
    if(varying == "slope1") {
        inits <- c(intercept = intercept_0,
                   ln_sigma_e = ln_sigma_e_0,
                   ln_sigma_w = ln_sigma_w_0)
    }
    if(varying == "slope2") {
        inits <- c(intercept = intercept_0,
                   slope1 = slope1_0,
                   ln_sigma_e = ln_sigma_e_0,
                   ln_sigma_w = ln_sigma_w_0)
    }

    if(optim_method == "PORT") {
        fit_op <- nlminb(start = inits,
                         objective = k_objective,
                         gradient = NULL,
                         hessian = NULL,
                         scale = 1,
                         control = list(eval.max = 200,
                                        iter.max = 5000,
                                        trace = 0),
                         lower = -Inf,
                         upper = Inf,
                         y = y, x1 = x1, x2 = x2,
                         varying = varying,
                         slope1_0 = slope1_0,
                         a_0 = a_0, p_0 = p_0, Ns = Ns)
    } else {
        fit_op <- optim(par = inits,
                        fn = k_objective,
                        gr = NULL,
                        hessian = TRUE,
                        control = list(maxit = 500),
                        method = optim_method,
                        y = y, x1 = x1, x2 = x2,
                        varying = varying,
                        slope1_0 = slope1_0,
                        a_0 = a_0, p_0 = p_0, Ns = Ns)
    }

    if(fit_op$convergence != 0 & warnings != 0) {
        warning("Model likelihood failed to converge")
        return(fit_op)
    }

    if(varying == "mean") {
        fit_ml <- k_filter(slope1 = slope1_0,
                           ln_sigma_e = fit_op$par[1],
                           ln_sigma_w = fit_op$par[2],
                           y = y,
                           x1 = x1,
                           x2 = x2,
                           varying = varying,
                           a_0 = a_0,
                           p_0 = p_0, Ns = Ns)
    }

    if(varying == "intercept") {
        fit_ml <- k_filter(slope1 = fit_op$par[1],
                           ln_sigma_e = fit_op$par[2],
                           ln_sigma_w = fit_op$par[3],
                           y = y,
                           x1 = x1,
                           x2 = x2,
                           varying = varying,
                           a_0 = a_0,
                           p_0 = p_0, Ns = Ns)
    }
    if(varying == "slope1") {
        fit_ml <- k_filter(intercept = fit_op$par[1],
                           ln_sigma_e = fit_op$par[2],
                           ln_sigma_w = fit_op$par[3],
                           y = y,
                           x1 = x1,
                           x2 = x2,
                           varying = varying,
                           slope1 = slope1_0,
                           a_0 = a_0,
                           p_0 = p_0, Ns = Ns)
        k_AICc <- 3
    }
    if(varying == "slope2") {
        fit_ml <- k_filter(intercept = fit_op$par[1],
                           slope1 = fit_op$par[2],
                           ln_sigma_e = fit_op$par[3],
                           ln_sigma_w = fit_op$par[4],
                           y = y,
                           x1 = x1,
                           x2 = x2,
                           varying = varying,
                           a_0 = a_0,
                           p_0 = p_0, Ns = Ns)
    }

    n_AICc <- sum(!is.na(y)) - Ns
    k_AICc <- length(inits)
    AICc   <- -2 * fit_ml$cum_loglik + (2 * k_AICc * n_AICc) / (n_AICc - k_AICc - 1)
    AIC    <- -2 * fit_ml$cum_loglik + (2 * k_AICc)

    fit_ml$intercept_0  <- intercept_0
    fit_ml$slope1_0     <- slope1_0
    fit_ml$a_0          <- a_0
    fit_ml$p_0          <- p_0
    fit_ml$ln_sigma_e_0 <- ln_sigma_e_0
    fit_ml$ln_sigma_w_0 <- ln_sigma_w_0
    fit_ml$AICc         <- AICc
    fit_ml$AIC          <- AIC
    fit_ml$n_params     <- k_AICc
    fit_ml$optim_method <- optim_method
    fit_ml$optim        <- fit_op

    return(fit_ml)
}



k_objective <- function(varying, inits, y, x1, x2, slope1_0, a_0, p_0, Ns) {
    ## Objective function for optimization routine
    ##
    ## inits = vector of initial values for optimization routine that is
    ##         dependent on the value of 'varying'
    ##
    ##         "intercept": vector of length three:
    ##                        slope1, ln_sigma_e, ln_sigma_w.
    ##         "slope1": vector of length three:
    ##                     intercept, ln_sigma_e, ln_sigma_w.
    ##         "slope2": vector of length four:
    ##                     intercept, slope1, ln_sigma_e, ln_sigma_w.
    ##
    ## Returns: cumulative negative log-likelihood

    if(varying == "mean") {
        nll <- k_filter(y = y,
                        x1 = x1,
                        slope1 = slope1_0,
                        ln_sigma_e = inits[1],
                        ln_sigma_w = inits[2],
                        varying = varying,
                        a_0 = a_0,
                        p_0 = p_0,
                        Ns = Ns)$cum_loglik * -1
    }
    if(varying == "intercept") {
        nll <- k_filter(y = y,
                        x1 = x1,
                        slope1 = inits[1],
                        ln_sigma_e = inits[2],
                        ln_sigma_w = inits[3],
                        varying = varying,
                        a_0 = a_0,
                        p_0 = p_0,
                        Ns = Ns)$cum_loglik * -1
    }
    if(varying == "slope1") {
        nll <- k_filter(y = y,
                        x1 = x1,
                        intercept = inits[1],
                        ln_sigma_e = inits[2],
                        ln_sigma_w = inits[3],
                        varying = varying,
                        a_0 = a_0,
                        p_0 = p_0,
                        Ns = Ns)$cum_loglik * -1
    }
    if(varying == "slope2") {
        nll <- k_filter(y = y,
                        x1 = x1,
                        x2 = x2,
                        intercept = inits[1],
                        slope1 = inits[2],
                        ln_sigma_e = inits[3],
                        ln_sigma_w = inits[4],
                        varying = varying,
                        a_0 = a_0,
                        p_0 = p_0,
                        Ns = Ns)$cum_loglik * -1
    }
    return(nll)
}



k_filter <- function(y, x1 = NULL, x2 = NULL,
                     varying = "mean",
                     intercept = NULL,
                     slope1 = NULL,
                     a_0, p_0,
                     ln_sigma_e,
                     ln_sigma_w,
                     Ns = 1) {
    ## Recursive calculation of Kalman filter
    ##
    ## y = response values
    ## x1 = 1st covariate values
    ## x2 = 2nd covariate values
    ## varying = one of "mean", "intercept", "slope1", "slope2"
    ##           "mean": varying intercept only model (e.g., model with no
    ##                   covariate effects). This is the local-level model.
    ##           "intercept": simple linear regression with time varying
    ##                        intercept. Need to specify 'slope1' value.
    ##           "slope1": simple linear regression time varying slope. Need to
    ##                     specify 'intercept' value.
    ##           "slope2": linear regression with 2 covariates, the second
    ##                     covariate is time varying (i.e., the slope for the
    ##                     'x2' variable). Need to specify 'intercept', 'slope1',
    ##                     and 'x2'.
    ##
    ## intercept = intercept value
    ## slope1 = slope value for x1 covariate
    ## a_0 = starting mean for time varying parameter
    ## p_0 = starting variance for time varying parameter
    ## ln_sigma_e = log_e of standard deviation of observation error
    ## ln_sigma_w = log_e of standard deviation of state error
    ## Ns = number of observations to exclude from likelihood calculation.
    ##      Default is 1, which excludes the first value from the likelihood,
    ##      i.e., condition the estimates on the first value, which can limit
    ##      the influence of the arbitrary starting mean and variance. Set to 0
    ##      to include all observations.
    ##
    ## Returns: a list

    if(varying %in% c("slope1", "slope2") & is.null(intercept))
        stop("Need input 'intercept' value")

    if(varying %in% c("intercept", "slope2") & is.null(slope1))
        stop("Need input 'slope1' value")

    if(varying == "slope2" & is.null(x2))
        stop("Need input 'x2' values")

    if(varying != "mean" & is.null(x1))
        stop("Need input 'x1' values")

    N <- length(y)
    sigma_e <- exp(ln_sigma_e)
    sigma_w <- exp(ln_sigma_w)

    x1[is.na(y)] <- NA
    x2[is.na(y)] <- NA

    if(varying %in% c("mean", "intercept"))
        x2 <- rep(1, N)
    if(varying == "slope1")
        x2 <- x1

    a_prior  <- rep(NA, N) ## one-step-ahead prediction
    p_prior  <- rep(NA, N) ## one-step-ahead prediction variance
    v_t      <- rep(NA, N) ## one-step-ahead prediction error
    f_t      <- rep(NA, N) ## one-step-ahead prediction error variance
    a_filter <- rep(NA, N) ## filtered value
    p_filter <- rep(NA, N) ## filtered value variance
    a_smooth <- rep(NA, N) ## smoothed value
    p_smooth <- rep(NA, N) ## smoothed value variance
    p_star   <- rep(NA, N) ## smoothing intermediary
    loglik   <- rep(NA, N) ## log-likelihood

    for(t in 1:N) {
        ## 1. Calculate one-step-ahead forecast (prior)
        if(t == 1) {
            a_prior[t] <- a_0
            p_prior[t] <- p_0
        } else {
            a_prior[t] <- a_filter[t-1]
            p_prior[t] <- p_filter[t-1] + (sigma_w * sigma_w)
        }

        if(is.na(y[t])) {
            a_filter[t] <- a_prior[t]
            p_filter[t] <- p_prior[t]
        } else {
            ## 2. Calculate one-step-ahead prediction error and variance
            if(varying == "mean")
                v_t[t] <- y[t] - a_prior[t]
            if(varying == "intercept")
                v_t[t] <- y[t] - (a_prior[t] + slope1 * x1[t])
            if(varying == "slope1")
                v_t[t] <- y[t] - (intercept + a_prior[t] * x2[t])
            if(varying == "slope2")
                v_t[t] <- y[t] - (intercept + slope1 * x1[t] + a_prior[t] * x2[t])

            f_t[t] <- (p_prior[t] * x2[t] * x2[t]) + (sigma_e * sigma_e)

            ## 3. Calculate filtered value and variance
            a_filter[t] <- a_prior[t] + ((p_prior[t] * x2[t] * v_t[t]) / f_t[t])
            p_filter[t] <- p_prior[t] -
                           ((p_prior[t] * p_prior[t] * x2[t] * x2[t]) / f_t[t])

            loglik[t] <- log(f_t[t]) + ((v_t[t] * v_t[t]) / f_t[t])
        }
    }

    ## 4. Fixed interval smoothing
    for(t in N:1) {
        if(t == N) {  ## last data point already conditioned on all data
            p_star <- 0
            a_smooth[t] <- a_filter[t]
            p_smooth[t] <- p_filter[t]
        } else {
            p_star[t] <- p_filter[t] / p_prior[t+1]
            a_smooth[t] <- a_filter[t] + p_star[t] *
                           (a_smooth[t+1] - a_prior[t+1])
            p_smooth[t] <- p_filter[t] + p_star[t] * p_star[t] *
                           (p_smooth[t+1] - p_prior[t+1])
        }
    }

    cum_loglik <- (-(N / 2) * log(2 * pi)) - (0.5 * sum(loglik[(Ns+1):N], na.rm = TRUE))

    ret <- list(a_filter = a_filter,
                p_filter = p_filter,
                a_smooth = a_smooth,
                p_smooth = p_smooth,
                v_t = v_t,
                f_t = f_t,
                loglik = loglik,
                varying = varying,
                intercept = intercept,
                slope1 = slope1,
                a_0 = a_0,
                p_0 = p_0,
                sigma_e = as.vector(sigma_e),
                sigma_w = as.vector(sigma_w),
                Ns = Ns,
                cum_loglik = cum_loglik)

    if(varying == "slope1") {
        ret[["slope1"]] <- NULL
    }

    return(ret)
}
