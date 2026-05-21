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
  vector<lower=L>[K] Y[N];

  // hyperparameters
  vector[K] a;
  vector[K] b;
  cov_matrix[K] R;
  cov_matrix[K] S;

}

parameters {

  // mean
  vector[K] theta;
  vector[K] alpha_raw[M]; 
  vector[K] beta1;
  vector[K] beta2;
  vector[K] gamma;

  // covariance
  vector<lower=0>[K] sigma_0;
  vector<lower=0>[K] sigma_e;

  // changepoint
  vector<lower=-20, upper=10>[K] delta;

}

transformed parameters {

  vector[K] alpha[M];

  // Non-centered: alpha[m] ~ N(theta, diag(sigma_0^2))
  for (m in 1:M)
    alpha[m] = theta + sigma_0 .* alpha_raw[m];

}

model {

  // Priors
  theta  ~ multi_normal(a, R);
  beta1  ~ multi_normal(b, S);
  beta2  ~ multi_normal(b, S);
  gamma  ~ multi_normal(b, S);
  sigma_e ~ cauchy(0, 5);
  sigma_0 ~ cauchy(0, 5);

  // Non-centered random intercepts
  for (m in 1:M)
    alpha_raw[m] ~ std_normal();

  // Likelihood
  for (i in 1:N) {

    for (k in 1:K) {

      real mu = alpha[id[i],k] + beta1[k]*g[i] + beta2[k]*t[i] + gamma[k]*g[i]*fdim(t[i], delta[k]);

      if (D[i,k] == 0) {

        target += normal_lpdf(Y[i,k] | mu, sigma_e[k]);

      } else if (D[i,k] == 1) {

        real z = (U[i,k] - mu) / sigma_e[k];
        real pi = 1.0 - Phi(z);
        target += log(fmax(pi, 1e-10));  // avoid log(0)

      }

    }

  }

}
