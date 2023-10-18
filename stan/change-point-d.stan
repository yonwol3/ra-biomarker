
data {
  
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> K;
  
  vector[K] Y[N];
  vector[N] t;
  int id[N];
  
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
  
  corr_matrix[K] corr_0;
  vector<lower=0>[K] sigma_0;
  corr_matrix[K] corr_e;
  vector<lower=0>[K] sigma_e;
  vector<lower=-20,upper=5>[K] kappa;
  
}

transformed parameters {
  
  cov_matrix[K] Sigma_0;
  cov_matrix[K] Sigma_e;
  Sigma_0 = quad_form_diag(corr_0, sigma_0);
  Sigma_e = quad_form_diag(corr_e, sigma_e);
  
}

model {
  
  vector[K] eta[N];
  
  // Random Intercept
  alpha ~ multi_normal(mu, Sigma_0);
  
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
  
  corr_0 ~ lkj_corr(1);
  sigma_0 ~ cauchy(0, 5);
  corr_e ~ lkj_corr(1);
  sigma_e ~ cauchy(0, 5);
  
  kappa ~ uniform(-20, 10);

}
