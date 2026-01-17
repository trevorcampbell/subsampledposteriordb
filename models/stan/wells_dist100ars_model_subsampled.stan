data {
  int<lower=0> N;
  array[N] int<lower=0, upper=1> switched;
  vector[N] dist;
  vector[N] arsenic;
}
transformed data {
  // rescaling
  vector[N] dist100 = dist / 100.0;
  matrix[N, 2] x = [dist100', arsenic']';
}
parameters {
  real alpha;
  vector[2] beta;
  real SUBIDX;
}
model {
  for (i in 1:N){
    if (i-0.5 <= SUBIDX && i+0.5 >= SUBIDX){
    	target += N*bernoulli_logit_glm_lpmf(switched[i] | [x[i, :]], alpha, beta);
    	break;
    } 
  }
}


