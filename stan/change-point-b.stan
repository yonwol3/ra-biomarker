
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
  vector<lower=L>[K] Y_obs[N];

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
  vector<lower=-20,upper=5>[K] kappa;
  
  // LOD
  vector<lower=L>[K] Y_max[N];

}

transformed parameters {

  cov_matrix[K] Sigma_0;
  matrix[K,K] Sigma_e;
  Sigma_0 = diag_matrix(square(sigma_0));
  Sigma_e = diag_pre_multiply(sigma_e, corr_e);

}

model {
  
  // Transformed Data
  vector[K] eta[N];
  vector[K] Y_full[N];
  
  // Random Intercept
  alpha ~ multi_normal(mu, Sigma_0);
  
  // Likelihood
  for (i in 1:N) {
    
    for (k in 1:K) {

      // impute censored values
      Y_full[i,k] = D_obs[i,k]*Y_obs[i,k] + D_max[i,k]*Y_max[i,k];
      
      // change point 
      eta[i,k] = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i],kappa[k]);
      
    }
    
  }
  
  // Likelihood
  Y_obs ~ multi_normal_cholesky(eta, Sigma_e);
  Y_obs ~ multi_normal_cholesky(eta, Sigma_e);


  // Priors
  mu ~ multi_normal(a, R);
  beta1 ~ multi_normal(b, S);
  beta2 ~ multi_normal(b, S);
  gamma ~ multi_normal(b, S);
  
  corr_e ~ lkj_corr_cholesky(1);
  sigma_e ~ cauchy(0, 10);
  sigma_0 ~ cauchy(0, 10);

  kappa ~ uniform(-20, 10);

}
