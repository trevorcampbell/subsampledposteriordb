data {
  int<lower=0> I;
  array[I] int<lower=0> n;
  array[I] int<lower=0> N;
  vector[I] x1;
  vector[I] x2;
}
transformed data {
  vector[I] x1x2;
  x1x2 = x1 .* x2;
}
parameters {
  real alpha0;
  real alpha1;
  real alpha12;
  real alpha2;
  real<lower=0> tau;
  vector[I] b;
  real SUBIDX;
}
transformed parameters {
  real<lower=0> sigma;
  sigma = 1.0 / sqrt(tau);
}
model {
  alpha0 ~ normal(0.0, 1.0E3);
  alpha1 ~ normal(0.0, 1.0E3);
  alpha2 ~ normal(0.0, 1.0E3);
  alpha12 ~ normal(0.0, 1.0E3);
  tau ~ gamma(1.0E-3, 1.0E-3);
  
  for (i in 1:I){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += I*normal_lpdf(b[i] | 0, sigma);
  		target += I*binomial_logit_lpmf(n[i] | N[i], alpha0 + alpha1 * x1[i] + alpha2 * x2[i] + alpha12 * x1x2[i] + b[i]);
  		break;
  	}
  }
}


