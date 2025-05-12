
data {
  
  // dimensions
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  int id[N];
  
  // data
  vector[K] D_max[N];
  vector[K] D_obs[N];

  vector[K] L[N];
  vector[K] U[N];
  vector[N] t;
  vector[N] g;
  vector[K] Y_obs[N];

  // hyperparameters
  vector[K] a;
  vector[K] b;
  cov_matrix[K] R;
  cov_matrix[K] S;
  
}

parameters {
  
  // mean function
  vector[K] mu;
  vector[K] alpha[M];
  vector[K] beta1;
  vector[K] beta2;
  vector[K] gamma;
  
  // covariance
  cholesky_factor_corr[K] corr_e;
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;
  vector<lower=-20,upper=10>[K] kappa;

}

transformed parameters {

  cov_matrix[K] Sigma_0 = diag_matrix(square(sigma_0));

}

model {

  // Priors
  mu ~ multi_normal(a, R);
  beta1 ~ multi_normal(b, S);
  beta2 ~ multi_normal(b, S);
  gamma ~ multi_normal(b, S);
  
  corr_e ~ lkj_corr_cholesky(1);
  sigma_e ~ normal(0, 2);
  sigma_0 ~ normal(0, 2);
  kappa ~ uniform(-20, 10);

  // Random Intercept
  alpha ~ multi_normal(mu, Sigma_0);
  
  // Likelihood
  for (i in 1:N) {

    vector[K] eta_i;
    
    for (k in 1:K) {

      eta_i[k] = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i],kappa[k]);

      if (D_obs[i,k] == 1) {

        target += normal_lpdf(Y_obs[i,k] | eta_i[k], sigma_e[k]);

      } else if (D_obs[i,k] == 0) {

        target += normal_lcdf(U[i,k] | eta_i[k], sigma_e[k]);

      }
      
    }
    
  }


}
