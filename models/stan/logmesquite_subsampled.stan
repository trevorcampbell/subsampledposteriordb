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
transformed data {
  // log transformations
  vector[N] log_weight;
  vector[N] log_diam1;
  vector[N] log_diam2;
  vector[N] log_canopy_height;
  vector[N] log_total_height;
  vector[N] log_density;
  log_weight = log(weight);
  log_diam1 = log(diam1);
  log_diam2 = log(diam2);
  log_canopy_height = log(canopy_height);
  log_total_height = log(total_height);
  log_density = log(density);
}
parameters {
  vector[7] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*normal_lpdf(log_weight[i] |  beta[1] + beta[2] * log_diam1[i] + beta[3] * log_diam2[i]
                      + beta[4] * log_canopy_height[i]
                      + beta[5] * log_total_height[i] + beta[6] * log_density[i]
                      + beta[7] * group[i], sigma);
  		break;
  	}
  }
}


