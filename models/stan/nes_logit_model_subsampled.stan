data {
  int<lower=0> N;
  vector[N] income;
  array[N] int<lower=0, upper=1> vote;
}
transformed data {
  matrix[N, 1] x = [income']';
}
parameters {
  real alpha;
  vector[1] beta;
  real SUBIDX;
}
model {
  for (n in 1:N){
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
      target += N*bernoulli_logit_glm_lpmf(vote[n] | [x[n,:]], alpha, beta);
      break;
    }
  }
}


