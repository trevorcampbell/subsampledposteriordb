data {
  int<lower=0> N;
  int<lower=0> D;
  matrix[N, D] X;
  vector[N] y;
}
parameters {
  vector[D] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  // prior
  target += normal_lpdf(beta | 0, 10);
  target += normal_lpdf(sigma | 0, 10);
  mu = X*beta;
  // likelihood
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*normal_lpdf(y[i] | mu[i], sigma);
  		break;
  	}
  }
}


