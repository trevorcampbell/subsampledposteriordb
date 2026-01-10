data {
  int<lower=2> V; // num words
  int<lower=1> M; // num docs
  int<lower=1> N; // total word instances
  array[N] int<lower=1, upper=V> w; // word n
  array[N] int<lower=1, upper=M> doc; // doc ID for word n
  vector<lower=0>[5] alpha; // topic prior
  vector<lower=0>[V] beta; // word prior
}
parameters {
  array[M] simplex[5] theta; // topic dist for doc m
  array[5] simplex[V] phi; // word dist for topic k
  real SUBIDX;
}
model {
  for (m in 1 : M) {
    theta[m] ~ dirichlet(alpha);
  } // prior
  for (k in 1 : 5) {
    phi[k] ~ dirichlet(beta);
  } // prior
  for (n in 1 : N) {
    if (n-0.5 <= SUBIDX && n+0.5 >= SUBIDX){
       array[5] real gamma;
       for (k in 1 : 5) {
         gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
       }
       target += N*log_sum_exp(gamma); // likelihood;
       break;
    }
  }
}


