data {
  int<lower=0> N;
  vector[N] weight;
  vector[N] diam1;
  vector[N] diam2;
  vector[N] canopy_height;
  vector[N] total_height;
  vector[N] density;
  vector[N] group;
}
parameters {
  vector[7] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (n in 1:N){
    if (n - 0.5 <= SUBIDX && n + 0.5 >= SUBIDX){
      target += N*normal_lpdf(weight[n] | beta[1] + beta[2] * diam1[n] + beta[3] * diam2[n]
                  + beta[4] * canopy_height[n] + beta[5] * total_height[n]
                  + beta[6] * density[n] + beta[7] * group[n], sigma);
      break;
    }
  }
}


