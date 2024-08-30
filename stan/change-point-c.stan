
data {
  
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  int id[N];
  
  vector[K] Y[N];
  vector[N] t;

  vector[K] a;
  vector[K] b;
  
  cov_matrix[K] R;
  cov_matrix[K] S;
  
}

parameters {
  
  vector[K] mu;
  vector[K] alpha[M];
  vector[K] beta1;
  vector[K] beta2;
  
  cholesky_factor_corr[K] corr_0;
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;
  vector<lower=-20,upper=5>[K] kappa;
  
}

transformed parameters {
  
  matrix[K,K] Sigma_0;
  cov_matrix[K] Sigma_e;
  Sigma_0 = diag_pre_multiply(sigma_0, corr_0);
  Sigma_e = diag_matrix(sigma_e);
  
}

model {
  
  vector[K] eta[N];
  
  // Random Intercept
  alpha ~ multi_normal_cholesky(mu, Sigma_0);
  
  // Likelihood
  for (i in 1:N) {
    
    for (k in 1:K) {
      
      // change point 
      eta[i,k] = alpha[id[i],k] + beta1[k]*t[i] + beta2[k]*fdim(t[i],kappa[k]);
      
    }
    
  }
  
  Y ~ multi_normal(eta, Sigma_e);
  
  // Priors
  mu ~ multi_normal(a, R);
  beta1 ~ multi_normal(b, S);
  beta2 ~ multi_normal(b, S);
  
  corr_0 ~ lkj_corr_cholesky(1);
  sigma_0 ~ cauchy(0, 10);
  sigma_e ~ cauchy(0, 10);
  kappa ~ uniform(-20, 10);

}
