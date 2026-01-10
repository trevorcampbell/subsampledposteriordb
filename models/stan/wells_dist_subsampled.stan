data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
}
parameters {
  vector[2] beta;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*bernoulli_logit_lpmf(switched[i] | beta[1] + beta[2] * dist[i]);
    	break;
    } 
  }
}


