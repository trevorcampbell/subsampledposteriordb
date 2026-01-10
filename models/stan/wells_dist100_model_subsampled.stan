data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
}
transformed data {
  // rescaling
  vector[N] dist100 = dist / 100.0;
  matrix[N, 1] x = [dist100']';
}
parameters {
  real alpha;
  vector[1] beta;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*bernoulli_logit_glm_lpmf(switched[i] | x[i, :], alpha, beta);
    	break;
    } 
  }
}


