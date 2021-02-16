// Poisson Hiden Markov Model

data {
  int<lower=0> N; // length of the time series
  int<lower=0> y[N]; // data
  int<lower=1> m; // number of states
}

parameters{
  simplex[m] A[m]; // tpm
  positive_ordered[m] lambda; // mean of poisson - ordered
}

model{
  vector[m] log_A_tr[m]; // log, transposed tpm
  vector[m] lp; // for forward variables
  vector[m] lp_p1; // for forward variables

  lambda ~ gamma(0.1, 0.01); // assigning exchangeable priors
  //(lambdas?s are ordered for sampling purposes)

  // transposing tpm and taking the log of each entry
  for(i in 1:m)
    for(j in 1:m)
      log_A_tr[j, i] = log(A[i, j]);

  lp = rep_vector(-log(m), m); //

  for(i in 1:N) {
    for(j in 1:m)
      lp_p1[j] = log_sum_exp(log_A_tr[j] + lp) + poisson_lpmf(y[i] | lambda[j]);

  lp = lp_p1;
}

target += log_sum_exp(lp);
}
