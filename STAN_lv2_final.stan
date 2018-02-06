data{
	int<lower=1> N; 				//number of total observations
	int<lower=1> K; 				//Number of level 1 predictors
	int<lower=1> J; 				//number of groups (level 2 nesting units)
	int<lower=1> L; 				//number of group level predictors
	int<lower=1,upper=J> ID[N]; 	//grouping variable
	matrix[N,K] x;					//individual level predictors
	matrix[J,L] u;					//group predictors 
	vector[N] y; 					//outcome variable
}

parameters{
	matrix[K,J] z;
	cholesky_factor_corr[K] L_Omega;
	vector<lower=0, upper=pi()/2>[K] tau_unif;			//prior scale
	matrix[L, K] gamma;									//level 2 coefficients																
	real<lower=0> sigma;									//level 1 error scale
}

transformed parameters{
	vector[N] mu_mod;				//model estimates of y
	matrix[J,K] beta;				//level 1 coefficients
	vector<lower=0>[K] tau;			//prior scale
	matrix[K,K] Sigma_beta; 		//variance-covariance matrix for random effects
	for(k in 1:K) tau[k] = 2.5*tan(tau_unif[k]);
	beta = u * gamma + (diag_pre_multiply(tau, L_Omega)*z)';
	mu_mod = rows_dot_product(beta[ID], x);
	Sigma_beta = diag_pre_multiply(tau, L_Omega) * diag_pre_multiply(tau, L_Omega)';
}

model{
	sigma ~ inv_gamma(5,5);
	to_vector(z) ~ normal(0,1);
	L_Omega ~ lkj_corr_cholesky(2);
	to_vector(gamma) ~ normal(0,5);
	for(n in 1:N){
		y[n] ~ lognormal(mu_mod[n], sigma);
	}
	
}

generated quantities{
	vector[N] y_pred;
	vector[N] High_sd_base;
	vector[N] Low_sd_base;
	vector[N] High_sd_pos;
	vector[N] Low_sd_pos;
	vector[N] High_sd_neg;
	vector[N] Low_sd_neg;
	vector[N] High_sd_both;
	vector[N] Low_sd_both;
	for(n in 1:N){
		y_pred[n] =  lognormal_rng(mu_mod[n], sigma);
		High_sd_base[n] = lognormal_rng((gamma[1,1]+gamma[2,1]), Sigma_beta[1,1]+sigma);
		Low_sd_base[n] = lognormal_rng((gamma[1,1]-gamma[2,1]), Sigma_beta[1,1]+sigma);
		High_sd_pos[n] = lognormal_rng((gamma[1,1]+gamma[2,1])+(gamma[1,2]+gamma[2,2]), Sigma_beta[2,2]+sigma);
		Low_sd_pos[n] = lognormal_rng((gamma[1,1]-gamma[2,1])+(gamma[1,2]-gamma[2,2]), Sigma_beta[2,2]+sigma);
		High_sd_neg[n] = lognormal_rng((gamma[1,1]+gamma[2,1])+(gamma[1,3]+gamma[2,3]), Sigma_beta[3,3]+sigma);
		Low_sd_neg[n] = lognormal_rng((gamma[1,1]-gamma[2,1])+(gamma[1,3]-gamma[2,3]), Sigma_beta[3,3]+sigma);
		High_sd_both[n] = lognormal_rng((gamma[1,1]+gamma[2,1])+(gamma[1,4]+gamma[2,4]), Sigma_beta[4,4]+sigma);
		Low_sd_both[n] = lognormal_rng((gamma[1,1]-gamma[2,1])+(gamma[1,4]-gamma[2,4]), Sigma_beta[4,4]+sigma);
	}
}
