data {
    int N; // number of observations
    int Y; // number of years
    int K; // number of addl covariates
    
    int year[N];
    vector[N] weeks_until;
    vector[N] logit_intent;
    vector<lower=0>[N] sd_intent;
    int seats[N]; // response variable
    vector[N] before; 
    matrix[N,K] X; // addl. covariates
}

parameters {
    vector[N] logit_true; // true voter intent, before meas. error
    real beta_intent;
    vector[K] betas; // coef. on addl. covariates
    
    // dispersion
    real<lower=0> phi;
}

transformed parameters {
    vector[N] mu;
    vector[N] alpha;
    vector[N] beta;
    
    mu = inv_logit(beta_intent*logit_true + X*betas);
    alpha = mu * phi;
    beta = (1 - mu) * phi;
}

model {
    // measurement error in logit-intentions
    logit_true ~ normal(0, 1);
    logit_intent ~ normal(logit_true, sd_intent);
    
    seats ~ beta_binomial(435, alpha, beta);
    
    beta_intent ~ normal(0, 50);
    betas ~ normal(0, 20);
    
    phi ~ chi_square(4);
}

generated quantities {
    vector[N] log_lik;
    vector[N] seats_pred;
    
    for (i in 1:N) {
        log_lik[i] = beta_binomial_lpmf(seats[i] | 435, alpha[i], beta[i]);
        seats_pred[i] = beta_binomial_rng(435, alpha[i], beta[i]);
    }
}
