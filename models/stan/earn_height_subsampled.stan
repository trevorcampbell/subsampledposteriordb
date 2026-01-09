data {
  int<lower=0> N;
  vector[N] earn;
  vector[N] height;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*normal_lpdf(earn[i] | beta[1] + beta[2] * height[i], sigma);
  		break;
  	}
  }
}


