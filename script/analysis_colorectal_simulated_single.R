library(mvtnorm)
set.seed(1)

args = commandArgs(trailingOnly=TRUE)

# setwd("/home/victor/phd/comp_receptor_modelling/")
source("script/alr_functions.R")

E <- read.csv("data/alexandrov2020/PCAWG_sigProfiler_SBS_signatures_in_samples.csv", header = TRUE, row.names = 2)
tissue <- "ColoRect"
E <- E[startsWith(E$Cancer.Types, tissue), ]
E <- E[, -c(1,2)]
E <- E[, colSums(E)>0]
channels <- c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT", "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT", "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT", "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT", "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT", "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT", "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT", "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT", "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT", "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT", "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT", "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT")
E <- E[, rev(colnames(E))]

# Impute zeros with multiplicative replacement strategy
W <- diag(1/rowSums(E)) %*% as.matrix(E)
for (i in 1:nrow(W)) {
  iszero <- unname(which(W[i, ] == 0))
  W[i, iszero] <- 0.0001
  W[i, -iszero] <- (1-0.0001*length(iszero))*W[i, -iszero] 
}

Walr <- alrTransformMatrix(W)
eta <- colMeans(Walr) # This is the prior mean for the exposures


# Now simulate tumours as described in Alexandrov 2020
# Pi: Proportion of tumours with a certain signature
Pi <- (E>0)
  Pi <- colMeans(Pi)
# Mu: Mean of log_10 number of signatures
Mu <- rep(0, ncol(E))
names(Mu) <- colnames(E)
logE <- log10(E)
logE[logE==-Inf] <- 0
Mu <- colSums(logE)/colSums(E > 0)
Sigma <- sqrt(colSums(logE^2)/colSums(E>0) - Mu^2)

G <- 10 # Number of simulated tumours
N <- ncol(E) # Number of signatures
Esimulated <- matrix(0, nrow = G, ncol = N)
for (g in 1:G) { # Iterate over patients
  for (n in 1:N) { # Iterate over signatures
    if (runif(1) < Pi[n]) { # With probability Pi[n]
      Esimulated[g, n] <- 10^rnorm(1, Mu[n], Sigma[n])
    }
  }
}
total_muts <- sum(Esimulated)

colnames(Esimulated) <- colnames(E)
Esimulated <- round(Esimulated)

S <- read.csv(
  "data/alexandrov2020/sigProfiler_SBS_signatures_2019_05_22.csv",
  header = TRUE, row.names = 1
)
S <- t(as.matrix(S))
S <- S[colnames(E), ]
K <- ncol(S)


SSimulated <- matrix(0, nrow = nrow(S), ncol = ncol(S))
for(sig in 1:nrow(SSimulated)) {
  SSimulated[sig, ] <- alrInv(rmvnorm(1, alr(S[sig, ]), 0.1*(matrix(1, nrow = K-1, ncol = K-1) + diag(K-1))))
}

MSimulated <- round(Esimulated %*% SSimulated)
colnames(MSimulated) <- channels
# write.csv(MSimulated, "data/simulated_colorectal/simulated_colorectal_15genomes.csv")
# MSimulated <- read.csv("data/simulated_colorectal/simulated_colorectal_15genomes.csv", header = TRUE, row.names = 1)
dim(MSimulated)

library(rstan)
library(lsa)
options(mc.cores = parallel::detectCores())

sigs_prior_vars <- list()
ScaleMatrix <- matrix(1, nrow = K-1, ncol = K-1) + diag(rep(1, K-1))  # Matrix with 2s in the diagonal and 1s in the other elements.
for (i in 1:nrow(S)) {
  sigs_prior_vars[[i]] <- 0.1 * ScaleMatrix
}
sigs_prior_vars <- abind::abind(sigs_prior_vars, along = 3)
sigs_prior_vars <- aperm(sigs_prior_vars, c(3, 1, 2))

patient_index <- as.integer(args[1])

dat <- list(
  N = N, K = K, G = G,
  M = unname(unlist(as.vector(MSimulated[patient_index, ]))),
  delta = N-1,
  muS = alrTransformMatrix(S),
  SigmaS = sigs_prior_vars,
  muE = eta,
  SigmaE = 10*(matrix(1, nrow = N-1, ncol = N-1) + diag(rep(1, N-1))),
  is_prior = 0
)

fit <- stan(
  file = "script/compositional_model_single_sample.stan",
  data = dat,
  iter = 21000,
  warmup = 1000,
  thin = 1,
  verbose = TRUE,
  chains = 3,
  init = function() list(
    alrE = c(rmvnorm(1, eta, 10*(matrix(1, nrow = N-1, ncol = N-1) + diag(rep(1, N-1))))),
    alrS = t(sapply(1:N, function(i) rmvnorm(1, alr(S[i,]), sigs_prior_vars[i,,])))
  )
)

# saveRDS(fit, paste0("results_dissertation/simulated/fit_simulated", patient_index, ".rds"))
saveRDS(fit, paste0("results/simulated/fit_simulated", patient_index, ".rds"))
