data {
  int<lower=0> N;
  vector<lower=0, upper=200>[N] kid_score;
  vector<lower=0, upper=200>[N] mom_iq;
}
parameters {
  vector[2] beta;
  real<lower=0> sigma;
  real SUBIDX;
}
model {
  sigma ~ cauchy(0, 2.5);
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*normal_lpdf(kid_score[i] | beta[1] + beta[2] * mom_iq[i], sigma);
    	break;
    }
  }
}


