data {
  int<lower=1> N;   // Number of signatures
  int<lower=1> K;   // Number of mutation categories/channels (typically K = 96)
  int<lower=0> M[K]; // mutational catalogue

  matrix[N, K-1] muS;     // Matrix of mutational signatures (mean)
  cov_matrix[K-1] SigmaS[N]; // Matrix of mutational signatures (variation)

  vector[N-1] muE;        // Matrix of signature exposures (mean)
  cov_matrix[N-1] SigmaE; // Matrix of signature exposuress (variation)
  
  int<lower=0,upper=1> is_prior; // If 1, sample from the prior and ignore the data
}

parameters {
  row_vector[N-1] alrE;          // alr-transformed exposures
  matrix[N, K-1] alrS;           // alr-transformed signature matrix
}

transformed parameters {
  row_vector<lower=0,upper=1>[N] E; // Signature exposures
  matrix<lower=0,upper=1>[N, K] S;  // Matrix of mutational signatures
  vector[N] alrInvMuE;              // Prior mean over the exposures

  real denominator;
  denominator = sum(exp(muE)) + 1;
  alrInvMuE[1:(N-1)] = exp(muE)/denominator;
  alrInvMuE[N] = 1/denominator;

  // Undo the alr-transformation of the signature matrix
  for (n in 1:N) {
    denominator = sum(exp(alrS[n, 1:(K-1)])) + 1;
    S[n, 1:(K-1)] = exp(alrS[n, 1:(K-1)])/denominator;
    S[n, K] = 1/denominator;
  }

  // Undo the alr-transformation of the exposure vector
  denominator = sum(exp(alrE[1:(N-1)])) + 1;
  E[1:(N-1)] = exp(alrE[1:(N-1)])/denominator;
  E[N] = 1/denominator;
}

model {
  // Prior
  alrE[1:(N-1)] ~ multi_normal(muE, SigmaE);

  for (n in 1:N) {
    alrS[n, 1:(K-1)]  ~ multi_normal(muS[n, 1:(K-1)], SigmaS[n]);  
  }
  
  // Likelihood
  if (is_prior == 0) {
    M ~ multinomial(to_vector(E * S));
  }
}

generated quantities {
  int<lower=0> M_pred[K] = multinomial_rng(to_vector(E * S), sum(M)); // predictive mutational catalogue
}
