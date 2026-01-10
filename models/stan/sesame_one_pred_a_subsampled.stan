data {
  int<lower=0> N;
  vector[N] encouraged;
  vector[N] watched;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (n in 1:N){
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
      target += N*normal_lpdf(watched[n] | beta[1] + beta[2] * encouraged[n], sigma);
      break;
    }
  }
}


