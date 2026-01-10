data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
}
transformed data {
  // log 10 transformation
  vector[N] log10_earn;
  for (i in 1 : N) {
    log10_earn[i] = log10(earn[i]);
  }
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (n in 1:N){
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
        target += N*normal_lpdf(log10_earn[n] | beta[1] + beta[2] * height[n], sigma);
        break;
    }
  }
}


