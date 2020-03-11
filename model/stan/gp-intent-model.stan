functions {
    matrix exp_kernel(int Nx, vector x, int Ny, vector y, real h, real lambda) {
        matrix[Nx, Ny] K;
		for(i in 1:Nx){
			for(j in 1:Ny){
				K[i, j] = square(h) * exp(-fabs(x[i] - y[j]) / lambda);
			}
		}
		return K;
    }
    
    vector post_pred_rng(real h, real lambda, real sigma_sq, int No, vector xo,
            int Np, vector xp, vector y) {
        matrix[No, No] Ko;
        matrix[Np, Np] Kp;
        matrix[No, Np] Kop;
        matrix[Np, No] Ko_inv_t;
        vector[Np] mu_p;
		matrix[Np,Np] Tau;
		matrix[Np,Np] L2;
		vector[Np] Yp;
		
		Ko = exp_kernel(No, xo, No, xo, h, lambda) 
    		+ diag_matrix(rep_vector(sigma_sq, No));
		Kp = exp_kernel(Np, xp, Np, xp, h, lambda) 
    		+ diag_matrix(rep_vector(sigma_sq, Np));
		Kop = exp_kernel(No, xo, Np, xp, h, lambda);
    	Ko_inv_t = Kop' / Ko;
    	
    	mu_p = Ko_inv_t * y;
		Tau = Kp - Ko_inv_t * Kop;
		L2 = cholesky_decompose(Tau);
		Yp = mu_p + L2*rep_vector(normal_rng(0,1), Np);
		return Yp;
    }
}

data {
    int W; // number of weeks 
    int N; // number of polls
    int P; // number of polling firms
    
    vector[N] d; // the date of a given poll
    vector[W] dp; // predicted dates
    //int<lower=1> n_resp[N]; // the size a given poll
    //int<lower=1> n_side[N]; // respondents who pick either DEM or GOP
    //int<lower=1> n_dem[N]; // respondents who pick DEM
    vecotr[N] ldem; // logit of proportion picking DEM
    int<lower=1> p[N]; // the polling firm for a given poll
}

parameters {
    real<lower=0> lambda; // time scale factor
    real<lower=0> h; // resp. scale factor
    real<lower=0.00001> sd_poll; // add'l polling errors
	vector[N] eta;
    
    vector[P] RP_pollster;  
    real<lower=0> sd_pollster; // hyperparameter for pollster errors
    real mu_pollster; // hyperparameter for global polling error
    
    real<lower=0,upper=1> prop_undecided; 
}

transformed parameters {
	vector[N] f;
    vector[P] pollster_error;  
	{
        matrix[N, N] Sigma;
    	matrix[N, N] L_S; // for cholesky factorization
        vector[N] logit_poll; // support for dem. in specific poll
        
        pollster_error = mu_pollster + RP_pollster * sd_pollster;
        for (i in 1:N)
            logit_poll[i] = pollster_error[p[i]];
        
    	Sigma = exp_kernel(N, d, N, d, h, lambda) + 
        	diag_matrix(rep_vector(square(sd_poll), N));
    	L_S = cholesky_decompose(Sigma);
    	
    	f = L_S * eta + logit_poll;
	}
}

model {
    n_side ~ binomial(n_resp, 1 - prop_undecided);
    n_dem ~ binomial_logit(n_side, f);
    
    RP_pollster ~ normal(0, 1);
    mu_pollster ~ normal(0, 0.02);
    prop_undecided ~ beta(2, 2);
    
    sd_pollster ~ student_t(4, 0, 1);
    sd_poll ~ student_t(4, 0, 1);
    
    h ~ student_t(4, 0, 1);
    lambda ~ inv_gamma(5, 5);
    eta ~ normal(0, 1);
}

generated quantities {
    vector[W] dem_margin; 
    vector[W] logit_dem; 
    
    logit_dem = post_pred_rng(h, lambda, sd_poll, N, d, W, dp, f);
    dem_margin = 2*inv_logit(logit_dem) - 1;
}

