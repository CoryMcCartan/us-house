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
    vector[W] logit_dem; // democratic support by week (main param)
    real<lower=0> alpha; // baseline variance in dem. support
    real<lower=0> rho; // correlation in dem. support between weeks
    
    vector[P] RP_pollster;  
    real<lower=0> sigma_pollster; // hyperparameter for pollster errors
    
    real<lower=0,upper=1> prop_undecided; 
    
    vector[N] RP_poll;  
    real<lower=0> sigma_poll; // hyperparameter for polling errors
}

transformed parameters {
    vector[N] logit_poll; // support for dem. in specific poll
    
    vector[P] pollster_error;  
    vector[N] poll_error; 
    
    pollster_error = RP_pollster * sigma_pollster;
    poll_error = RP_poll * sigma_poll;
    
    {
        vector[W] RP_logit;
        matrix[W, W] K; // correlation matrix
        matrix[W, W] L_K; // decomposition
        real wks[W];
        
        for (i in 1:W) wks[i] = i;
        K = cov_exp_quad(wks, alpha, rho);
        // diagonal elements
        for (i in 1:W) K[i, i] = K[i, i] + 1e-9;
        
        L_K = cholesky_decompose(K);
        RP_logit = L_K * logit_dem;
        
        for (i in 1:N)
            logit_poll[i] = RP_logit[w[i]]  + sigma_pollster*pollster_error[p[i]] 
                + sigma_poll*poll_error[i];
    }
}

model {
    n_side ~ binomial(n_resp, 1 - prop_undecided);
    n_dem ~ binomial_logit(n_side, logit_poll);
    
    logit_dem ~ normal(0, 1);
    
    RP_pollster ~ normal(0, 1);
    RP_poll ~ normal(0, 1);
    prop_undecided ~ beta(2, 2);
    
    alpha ~ student_t(4, 0, 1);
    rho ~ student_t(4, 0, 10);
    sigma_pollster ~ student_t(4, 0, 0.05);
    sigma_poll ~ student_t(4, 0, 0.1);
}

generated quantities {
    /*
    vector[N] log_lik;
    for (i in 1:N) {
        log_lik[i] = binomial_lpmf(n_side | n_resp, 1 - prop_undecided) +
                        binomial_logit_lpmf(n_dem | n_side, logit_poll);
    }
    */
    
    vector<lower=0,upper=1>[W] dem_support; 
    dem_support = inv_logit(logit_dem);
}
