// Gaussian linear model with adjustable priors
data {
  int<lower=0> N; // number of data points
  vector[N] x; //
  vector[N] y; //
  real xpred; // input location for prediction
  real pmualpha; // prior mean for alpha
  real psalpha; // prior std for alpha
  real pmubeta; // prior mean for beta
  real psbeta; // prior std for beta
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  alpha ~ normal(pmualpha, psalpha);
  beta ~ normal(pmubeta, psbeta);
  for (i in 1:N){
  	if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
  		target += normal_lpdf(y[i] | alpha + beta * x[i], sigma);
  		break;
  	}
  }
}


