data {
  int<lower=0> N;
  vector[N] floor_measure;
  vector[N] log_radon;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma_y;
  real SUBIDX;
}
model {
  vector[N] mu;
  
  // priors
  sigma_y ~ normal(0, 1);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  
  // likelihood
  mu = alpha + beta * floor_measure;
  for (n in 1 : N) {
  	if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
    	target += N*normal_lpdf(log_radon[n] | mu[n], sigma_y);
    	break;
    }
  }
}


