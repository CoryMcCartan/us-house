data {
    int W; // number of weeks 
    int N; // number of polls
    int P; // number of polling firms
    
    int<lower=1> w[N]; // the week of a given poll
    int<lower=1> n_resp[N]; // the size a given poll
    int<lower=1> n_side[N]; // respondents who pick either DEM or GOP
    int<lower=1> n_dem[N]; // respondents who pick DEM
    int<lower=1> p[N]; // the polling firm for a given poll
}

parameters {
    vector[W] delta_dem; // steps of random walk
    real<lower=0> sd_walk; // week-to-week variance
    real<lower=0,upper=1> rho; // strength of mean reversion
    
    vector[P] RP_pollster;  
    real<lower=0> sd_pollster; // hyperparameter for pollster errors
    
    real alpha_undecided; // baseline undecided
    real beta_undecided;  // trend in undecided
    
    vector[N] RP_poll;  
    real<lower=0> sd_poll; // hyperparameter for polling errors
}

transformed parameters {
    vector[W] logit_dem; // democratic support by week (main param)
    vector[N] logit_poll; // support for dem. in specific poll
    vector[W] logit_undecided;
    
    vector[P] pollster_error;  
    vector[N] poll_error; 
    
    pollster_error = RP_pollster * sd_pollster;
    poll_error = RP_poll * sd_poll;
    
    logit_dem[1] = 10 * delta_dem[1];
    for (i in 2:W)
        logit_dem[i] = rho*logit_dem[i-1] + sd_walk*delta_dem[i];
        
    for (i in 1:W)
        logit_undecided[i] = alpha_undecided + i*beta_undecided;
    
    for (i in 1:N)
        logit_poll[i] = logit_dem[w[i]]  + sd_pollster*pollster_error[p[i]] 
            + sd_poll*poll_error[i];
}

model {
    for (i in 1:N)
        n_side[i] ~ binomial_logit(n_resp[i], -logit_undecided[w[i]]);
    n_dem ~ binomial_logit(n_side, logit_poll);
    
    delta_dem ~ normal(0, 1);
    
    RP_pollster ~ normal(0, 1);
    RP_poll ~ normal(0, 1);
    alpha_undecided ~ normal(0, 10);
    beta_undecided ~ normal(0, 1);
    
    sd_walk ~ student_t(4, 0, 0.1);
    rho ~ beta(2, 1);
    sd_pollster ~ student_t(4, 0, 0.05);
    sd_poll ~ student_t(4, 0, 0.1);
}

generated quantities {
    vector[N] log_lik;
    vector[W] dem_margin; 
    vector[W] prop_undecided; 
    
    for (i in 1:N) 
        log_lik[i] = binomial_logit_lpmf(n_side[i] | n_resp[i], -logit_undecided[w[i]]) +
                        binomial_logit_lpmf(n_dem[i] | n_side[i], logit_poll[i]);
                        
    prop_undecided = inv_logit(logit_undecided);
    dem_margin = 2*inv_logit(logit_dem) - 1;
}
