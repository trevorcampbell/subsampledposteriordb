data {
  int<lower=0> N;
  vector[N] weight;
  vector[N] diam1;
  vector[N] diam2;
  vector[N] canopy_height;
  vector[N] group;
}
transformed data {
  vector[N] log_weight;
  vector[N] log_canopy_volume;
  vector[N] log_canopy_area;
  log_weight = log(weight);
  log_canopy_volume = log(diam1 .* diam2 .* canopy_height);
  log_canopy_area = log(diam1 .* diam2);
}
parameters {
  vector[4] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*normal_lpdf(log_weight[i] |  beta[1] + beta[2] * log_canopy_volume[i]
                      + beta[3] * log_canopy_area[i] + beta[4] * group[i], sigma);
  		break;
  	}
  }
}


