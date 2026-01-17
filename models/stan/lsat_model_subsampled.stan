data {
  int<lower=0> N; // 1000, number of students
  int<lower=0> R; // 32, number of patterns of results: 2^T
  int<lower=0> T; // 5, number of questions
  array[R] int<lower=0> culm;
  array[R, T] int<lower=0> response;
}
transformed data {
  array[T, N] int r;
  vector[N] ones;
  
  for (j in 1 : culm[1]) {
    for (k in 1 : T) {
      r[k, j] = response[1, k];
    }
  }
  for (i in 2 : R) {
    for (j in (culm[i - 1] + 1) : culm[i]) {
      for (k in 1 : T) {
        r[k, j] = response[i, k];
      }
    }
  }
  for (i in 1 : N) {
    ones[i] = 1.0;
  }
}
parameters {
  array[T] real alpha;
  vector[N] theta;
  real<lower=0> beta;
  real SUBIDX;
}
model {
  alpha ~ normal(0, 100.);
  theta ~ normal(0, 1);
  beta ~ normal(0.0, 100.);
  int ii = 1;
  for (k in 1 : T) {
    for (n in 1: N) {
      if (ii - 0.5 <= SUBIDX && ii + 0.5 >= SUBIDX){
        target += N*T*bernoulli_logit_lpmf(r[k,n] | beta * theta[n] - alpha[k]);
      }
      ii += 1;
    }
  }
}
generated quantities {
  real mean_alpha;
  array[T] real a;
  mean_alpha = mean(alpha);
  for (t in 1 : T) {
    a[t] = alpha[t] - mean_alpha;
  }
}


