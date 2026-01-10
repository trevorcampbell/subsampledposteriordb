data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
}
transformed data {
  // rescaling
  vector[N] dist100 = dist / 100.0;
  // interaction
  vector[N] inter = dist100 .* arsenic;
  matrix[N, 3] x = [dist100', arsenic', inter']';
}
parameters {
  real alpha;
  vector[3] beta;
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


