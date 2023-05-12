## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simulate time series data sets
## Michael Malick <mjm@michaelmalick.com>
## 2018-11-07
## Version: 1.0
##
## Functions:
##  - ts_step   Simulate step function time series
##  - ts_ar1    Simulate AR1 time series
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ts_step <- function(n, k,
                    up_mean,
                    down_mean,
                    up_sd = 0,
                    down_sd = 0,
                    period = NULL) {
    ## Simulate a step function time series
    ##
    ## At each step point a new mean is drawn from N(up_mean, up_sd) or
    ## N(down_mean, down_sd).
    ##
    ## n = length of output series
    ## k = number of up/down steps, used when `period` is NULL
    ## up_mean = mean of step-up series
    ## down_mean = mean of step-down series
    ## up_sd = standard deviation of step-up mean
    ## down_sd = standard deviation of step-down mean
    ## period = if non-null, gives the period of the step-cycles, i.e., the
    ##          total length of an up + down cycle. When non-null used instead
    ##          of `k`.

    if(is.null(period)) {
        k.ind <- cut(seq(1, n), breaks = k, labels = FALSE)
    } else {
        p <- period / 2
        rem <- n %% p
        k.ind <- cut(seq(1, n - rem), breaks = n %/% p, labels = FALSE)
        k.ind <- c(k.ind, rep(max(k.ind) + 1, rem))
        k <- length(unique(k.ind))
    }
    lst <- vector("list", k)
    for(i in 1:k) {
        n.i <- length(k.ind[k.ind == i])
        up <- ifelse(i %% 2 == 1, TRUE, FALSE)
        if(up)
            lst[[i]] <- rep(rnorm(1, mean = up_mean, sd = up_sd), n.i)
        else
            lst[[i]] <- rep(rnorm(1, mean = down_mean, sd = down_sd), n.i)
    }
    x <- unlist(lst)
    return(x)
}


ts_ar1 <- function(n, x_0, rho, sd, mean = 0) {
    ## Simulate an autoregressive order-1 (AR1) time series
    ##
    ## x[t] = rho * x[t-1] + e[t]
    ## e[t] ~ N(0, sd)
    ##
    ## If rho < 1:
    ##  - mean converges to: `mean`
    ##  - variance converges to: sd^2 / (1 - rho^2)
    ##
    ## If rho == 1: random walk
    ## If rho == 0: white noise
    ##
    ## n = length of output series
    ## x_0 = initial value used to calculate x[1]
    ## rho = autocorrelation coefficient
    ## sd = standard deviation of random error
    ## mean = mean of series

    if(rho < -1 | rho > 1)
        stop("rho should be: -1 <= rho <= 1")
    x <- rep(NA, n)
    x[1] <- (x_0 * rho) + rnorm(1, mean = 0, sd = sd)
    for(i in 2:n) {
        x[i] <- (x[i - 1] * rho) + rnorm(1, mean = 0, sd = sd)
    }
    x <- mean + x
    return(x)
}
