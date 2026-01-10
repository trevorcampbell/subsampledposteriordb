data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
  vector[N] male;
}
transformed data {
  vector[N] log_earn; // log transformations
  vector[N] log_height;
  log_earn = log(earn);
  log_height = log(height);
}
parameters {
  vector[3] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
      target += N*normal_lpdf(log_earn[i] | beta[1] + beta[2] * log_height[i] + beta[3] * male[i], sigma);
      break;
    }
  }
}


