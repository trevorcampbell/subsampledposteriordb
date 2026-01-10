data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] int<lower=1, upper=J> county;
  vector[N] y;
}
parameters {
  vector[J] a;
  real mu_a;
  real<lower=0, upper=100> sigma_a;
  real<lower=0, upper=100> sigma_y;
  real SUBIDX;
}
model {
  vector[N] y_hat;
  for (i in 1 : N) {
    y_hat[i] = a[county[i]];
  }
  
  mu_a ~ normal(0, 1);
  a ~ normal(mu_a, sigma_a);
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*normal_lpdf(y[i] | y_hat[i], sigma_y);
  		break;
  	}
  }
}


