
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
  vector[K] gamma;
  
  cholesky_factor_corr[K] corr_0;
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;
  vector<lower=-20,upper=5>[K] kappa;

  vector<upper=L[1]>[N] Y_cens_min1;
  vector<upper=L[2]>[N] Y_cens_min2;
  vector<upper=L[3]>[N] Y_cens_min3;
  vector<upper=L[4]>[N] Y_cens_min4;
  vector<upper=L[5]>[N] Y_cens_min5;
  vector<upper=L[6]>[N] Y_cens_min6;

  vector<lower=U[1]>[N] Y_cens_max1;
  vector<lower=U[2]>[N] Y_cens_max2;
  vector<lower=U[3]>[N] Y_cens_max3;
  vector<lower=U[4]>[N] Y_cens_max4;
  vector<lower=U[5]>[N] Y_cens_max5;
  vector<lower=U[6]>[N] Y_cens_max6;
  
}

transformed parameters {

  vector[K] Y_cens_min[N];
  vector[K] Y_cens_max[N];

  for (i in 1:N) {
  
    Y_cens_min[i,1] = Y_cens_min1[i];
    Y_cens_min[i,2] = Y_cens_min2[i];
    Y_cens_min[i,3] = Y_cens_min3[i];
    Y_cens_min[i,4] = Y_cens_min4[i];
    Y_cens_min[i,5] = Y_cens_min5[i];
    Y_cens_min[i,6] = Y_cens_min6[i];

    Y_cens_max[i,1] = Y_cens_max1[i];
    Y_cens_max[i,2] = Y_cens_max2[i];
    Y_cens_max[i,3] = Y_cens_max3[i];
    Y_cens_max[i,4] = Y_cens_max4[i];
    Y_cens_max[i,5] = Y_cens_max5[i];
    Y_cens_max[i,6] = Y_cens_max6[i];

  }
  
  matrix[K,K] Sigma_0;
  cov_matrix[K] Sigma_e;
  Sigma_0 = diag_pre_multiply(sigma_0, corr_0);
  Sigma_e = diag_matrix(sigma_e^2);
  
}

model {
  
  vector[K] eta[N];
  vector[K] Y_full[N];
  
  // Random Intercept
  alpha ~ multi_normal_cholesky(mu, Sigma_0);
  
  // Likelihood
  for (i in 1:N) {
    
    for (k in 1:K) {

      // impute censored values
      Y_full[i,k] = D_obs[i,k]*Y_obs[i,k] + D_min[i,k]*Y_cens_min[i,k] + D_max[i,k]*Y_cens_max[i,k];
      
      // change point 
      eta[i,k] = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i],kappa[k]);
      
    }
    
  }
  
  Y_full ~ multi_normal(eta, Sigma_e);
  Y_cens_min ~ multi_normal(eta, Sigma_e);
  Y_cens_max ~ multi_normal(eta, Sigma_e);
  
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