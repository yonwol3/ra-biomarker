
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
  
  // mean
  vector[K] theta;
  vector[K] alpha[M];
  vector[K] beta1;
  vector[K] beta2;
  vector<lower=0>[K] gamma;
  
  // covariance
  cholesky_factor_corr[K] corr_e;
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;
  
  // changepoint
  vector<lower=-20,upper=10>[K] delta;

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
  delta ~ uniform(-20, 10);
  
  // Transformed Data
  vector[K] mu[N];

  // Random Intercept
  alpha ~ multi_normal(theta, Sigma_0);

  for (i in 1:N) {
    
    for (k in 1:K) {
      
      mu[i,k] = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i], delta[k]);
      
    }
    
  }

  // Likelihood
  Y ~ multi_normal_cholesky(mu, Sigma_e);

}
