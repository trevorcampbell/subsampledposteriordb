data {
  int<lower=0> n; // Number of years
  array[n] int<lower=0> C; // Counts
  vector[n] year; // Year
}
transformed data {
  vector[n] year_squared;
  vector[n] year_cubed;
  
  year_squared = year .* year;
  year_cubed = year .* year .* year;
}
parameters {
  real<lower=-20, upper=20> alpha;
  real<lower=-10, upper=10> beta1;
  real<lower=-10, upper=20> beta2;
  real<lower=-10, upper=10> beta3;
  vector[n] eps; // Year effects
  real<lower=0, upper=5> sigma;
  real SUBIDX;
}
transformed parameters {
  vector[n] log_lambda;
  
  // Linear predictor incl. random year effect
  log_lambda = alpha + beta1 * year + beta2 * year_squared
               + beta3 * year_cubed + eps;
}
model {
  // Priors
  alpha ~ uniform(-20, 20);
  beta1 ~ uniform(-10, 10);
  beta2 ~ uniform(-10, 10);
  beta3 ~ uniform(-10, 10);
  sigma ~ uniform(0, 5);
  
  // Likelihood
  for (i in 1:n){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += n*poisson_log_lpmf(C[i] | log_lambda[i]);
  		target += n*normal_lpdf(eps[i] | 0, sigma);
  		break;
  	}
  }
}
generated quantities {
  vector<lower=0>[n] lambda;
  
  lambda = exp(log_lambda);
}


