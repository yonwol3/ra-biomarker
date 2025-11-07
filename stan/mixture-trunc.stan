
data {
  
  // dimensions
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  int id[N];
  
  // data
  vector[K] D[N];
  vector[K] L[N];
  vector[K] U[N];
  vector[N] t;
  vector[N] g;
  vector<lower=L, upper=U>[K] Y[N];

  // hyperparameters
  vector[K] a;
  vector[K] b;
  cov_matrix[K] R;
  cov_matrix[K] S;
  
}

parameters {
  
  // mean function
  vector[K] theta;
  vector[K] alpha[M];
  vector[K] beta1;
  vector[K] beta2;
  vector[K] gamma;
  
  // covariance
  cholesky_factor_corr[K] corr_e;
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;
  
  // changepoint
  vector<lower=-20,upper=10>[K] delta;
  vector<lower=0>[K] phi; 

}

transformed parameters {

  cov_matrix[K] Sigma_0 = diag_matrix(square(sigma_0));
  matrix[K,K] Sigma_e = diag_pre_multiply(sigma_e, corr_e);
  
}

model {

  // Priors
  theta ~ multi_normal(a, R);
  beta1 ~ multi_normal(b, S);
  beta2 ~ multi_normal(b, S);
  gamma ~ multi_normal(b, S);
  corr_e ~ lkj_corr_cholesky(1);
  sigma_e ~ cauchy(0, 5);
  sigma_0 ~ cauchy(0, 5);
  phi ~ cauchy(0, 5);
  delta ~ uniform(-20, 10);
  
  // Transformed Data
  vector[K] mu[N];

  // Random Intercept
  alpha ~ multi_normal(theta, Sigma_0);
  
  // Likelihood

  for (i in 1:N) {
    
    for (k in 1:K) {

      real rho = inv_logit((t[i] - delta[k]) * phi[k]);
      real eta_0 = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i];
      real eta_1 = eta_0 + gamma[k]*g[i]*fdim(t[i],delta[k]);
      mu[i,k] = (1 - rho) * eta_0 + rho * eta_1;
      
    }
    
  }

  Y ~ multi_normal_cholesky(mu, Sigma_e);

}
