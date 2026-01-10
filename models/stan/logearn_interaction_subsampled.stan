data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
  vector[N] male;
}
transformed data {
  vector[N] log_earn; // log transformation
  vector[N] inter; // interaction
  log_earn = log(earn);
  inter = height .* male;
}
parameters {
  vector[4] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
      target += N*normal_lpdf(log_earn[i] | beta[1] + beta[2] * height[i] + beta[3] * male[i]
                    + beta[4] * inter[i], sigma);
      break;
    }
  }
}


