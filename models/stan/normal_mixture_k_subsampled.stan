data {
  int<lower=1> K;
  int<lower=1> N;
  array[N] real y;
}
parameters {
  simplex[K] theta;
  array[K] real mu;
  array[K] real<lower=0, upper=10> sigma;
  real SUBIDX;
}
model {
  array[K] real ps;
  mu ~ normal(0, 10);
  for (n in 1 : N) {
    for (k in 1 : K) {
      ps[k] = log(theta[k]) + normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
      target += N*log_sum_exp(ps);
      break;
    }
  }
}


