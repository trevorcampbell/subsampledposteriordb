data {
  int<lower=0> N;
  array[N] real x;
  array[N] real Y;
}
parameters {
  real alpha;
  real beta;
  real<lower=.5, upper=1> lambda; // orginal gamma in the JAGS example  
  real<lower=0> tau;
  real SUBIDX;
}
transformed parameters {
  real sigma;
  real U3;
  sigma = 1 / sqrt(tau);
  U3 = logit(lambda);
}
model {
  array[N] real m;
  for (i in 1 : N) {
    m[i] = alpha - beta * pow(lambda, x[i]);
  }
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += N*normal_lpdf(Y[i] | m[i], sigma);
  		break;
  	}
  }
  
  alpha ~ normal(0.0, 1000);
  beta ~ normal(0.0, 1000);
  lambda ~ uniform(.5, 1);
  tau ~ gamma(.0001, .0001);
}


