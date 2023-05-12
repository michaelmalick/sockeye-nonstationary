// Single-stock Ricker model + AR1 error term + time-varying covariate
//
// Parameters:
//   alpha: intercept
//   beta: density dependence term
//   gamma: time-varying  g[t] = g[t-1] + N(0, sigma_gamma)
//   sigma_gamma: process error
//   g0: first value of random walk
//   sigma: observation error
//   phi: autocorrelation parameter

data {
    int<lower=0> N;                       // total number of years
    real x1[N];                           // 1st covariate
    real x2[N];                           // 2nd covariate
    real y[N];                            // response
    real sigma_gamma_df;                  // prior df for sigma_gamma
    real sigma_gamma_mu;                  // prior mu for sigma_gamma
    real sigma_gamma_sd;                  // prior sd for sigma_gamma
    int priors_only;                      // should likelihood be ignored?
}
parameters {
    real alpha;                           // series-specific intercept
    real beta;                            // slopes for 1st covariate
    real<lower=-1, upper=1> phi;          // autocorrelation parameter
    real<lower=0> sigmaNC;                // residual SD (not corrected for AR1)
    real g0;                              // gamma at t = 1
    real d_gamma[N];                      // random gamma deviate
    real<lower=0> sigma_gamma;            // gamma SD
}
transformed parameters {
    real yhat[N];                         // predicted values
    real epsilon[N];                      // residuals
    real<lower=0> sigma;        // residual SD (corrected for AR1)
    real gamma[N];                        // 2nd covariate effect
    real tmp_epsilon;                     // temporary var to avoid deep copy
    real tmp_gamma;                       // temporary var to avoid deep copy

    for(t in 1:N) {
        if(t == 1) {
            gamma[1] = g0;
            yhat[1] = alpha + beta * x1[1] + gamma[1] * x2[1];
            epsilon[1] = y[1] - yhat[1];
        } else {
            tmp_gamma = gamma[t-1];
            tmp_epsilon = epsilon[t-1];
            gamma[t] = tmp_gamma + sigma_gamma * d_gamma[t-1];
            yhat[t] = alpha + beta * x1[t] + gamma[t] * x2[t] + phi * tmp_epsilon;
            epsilon[t] = y[t] - (yhat[t] - (phi * tmp_epsilon));
        }

        // AR1 impacts variance: sigma^2 = sigmaNC^2 * (1 - phi^2)
        sigma = sqrt(sigmaNC * sigmaNC * (1 - phi * phi));
    }
}
model {
    alpha ~ normal(1, 5);
    phi ~ normal(0, 1);
    sigma_gamma ~ student_t(sigma_gamma_df, sigma_gamma_mu, sigma_gamma_sd);
    sigmaNC ~ student_t(3, 0, 3);
    beta ~ normal(0, 1);
    g0 ~ normal(0, 3);
    d_gamma ~ normal(0, 1);

    // likelihood
    if(!priors_only) {
        y ~ normal(yhat, sigma);
    }
}
generated quantities {
    real log_lik[N];                      // log-likelihood
    real yrep[N];                         // replicated datasets
    real signal_noise;                    // signal-to-noise ratio
    for(t in 1:N) {
        log_lik[t] = normal_lpdf(y[t] | yhat[t], sigma);
        yrep[t] = normal_rng(yhat[t], sigma);
    }
    signal_noise = sigma_gamma / sigma;
}

