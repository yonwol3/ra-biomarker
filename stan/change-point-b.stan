
data {
  
  // dimensions
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  int id[N];
  
  // data
  vector[K] Y_obs[N];
  vector[K] D_min[N];
  vector[K] D_max[N];
  vector[K] D_obs[N];
  vector[K] L[N];
  vector[K] U[N];
  vector[N] t;
  vector[N] g;

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
  cholesky_factor_corr[K] corr_0;
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;
  vector<lower=-20,upper=5>[K] kappa;

  // LOD
  vector<upper=L>[K] Y_min[N];
  vector<lower=U>[K] Y_max[N];
  
}

transformed parameters {
  
  matrix[K,K] Sigma_0;
  cov_matrix[K] Sigma_e;
  Sigma_0 = diag_pre_multiply(sigma_0, corr_0);
  Sigma_e = diag_matrix(sigma_e^2);
  
}

model {
  
  // Transformed Data
  vector[K] eta[N];
  vector[K] Y_full[N];
  
  // Random Intercept
  alpha ~ multi_normal_cholesky(mu, Sigma_0);
  
  // Likelihood
  for (i in 1:N) {
    
    for (k in 1:K) {

      // impute censored values
      Y_full[i,k] = D_obs[i,k]*Y_obs[i,k] + D_min[i,k]*Y_min[i,k] + D_max[i,k]*Y_max[i,k];
      
      // change point 
      eta[i,k] = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i],kappa[k]);
      
    }
    
  }
  
  // Sample Data
  Y_full ~ multi_normal(eta, Sigma_e);
  Y_min ~ multi_normal(eta, Sigma_e);
  Y_max ~ multi_normal(eta, Sigma_e);
  
  // Priors
  mu ~ multi_normal(a, R);
  beta1 ~ multi_normal(b, S);
  beta2 ~ multi_normal(b, S);
  gamma ~ multi_normal(b, S);

  corr_0 ~ lkj_corr_cholesky(1);
  sigma_0 ~ cauchy(0, 10);
  sigma_e ~ cauchy(0, 10);
  kappa ~ uniform(-20, 10);

}