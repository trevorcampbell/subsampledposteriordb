data {
  int<lower=0> J; // number of schools
  array[J] real y; // estimated treatment
  array[J] real<lower=0> sigma; // std of estimated effect
}
parameters {
  array[J] real theta; // treatment effect in school j
  real mu; // hyper-parameter of mean
  real<lower=0> tau; // hyper-parameter of sdv
  real SUBIDX;
}
model {
  tau ~ cauchy(0, 5); // a non-informative prior
  mu ~ normal(0, 5);
  for (i in 1:J){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += J*normal_lpdf(theta[i] | mu, tau);
  		target += J*normal_lpdf(y[i] | theta[i], sigma[i]);
  		break;
  	}
  }
}


