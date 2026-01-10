data {
  int<lower=0> I;
  array[I] int<lower=0> n;
  array[I] int<lower=0> N;
  vector[I] x1; // seed type
  vector[I] x2; // root extract
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
  vector[I] b;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  alpha0 ~ normal(0.0, 1.0); // Narrower priors
  alpha1 ~ normal(0.0, 1.0);
  alpha2 ~ normal(0.0, 1.0);
  alpha12 ~ normal(0.0, 1.0);
  sigma ~ cauchy(0, 1);
  
  for (i in 1:I){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += I*normal_lpdf(b[i] | 0, sigma);
  		target += I*binomial_logit_lpmf(n[i] | N[i], alpha0 + alpha1 * x1[i] + alpha2 * x2[i] + alpha12 * x1x2[i] + b[i]);
  		break;
  	}
  }
}


