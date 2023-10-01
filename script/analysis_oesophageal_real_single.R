rm(list = ls())
library(mvtnorm)
library(ggmcmc)

set.seed(1)

args = commandArgs(trailingOnly=TRUE)


# setwd("/home/victor/phd/comp_receptor_modelling/")
source("script/alr_functions.R")

# First, work out the prior mean for exposures, eta
E <- read.csv(
  "data/degasperi2022/exposures.csv",
  header = TRUE, row.names = 1
)
colorectal_indices <- rownames(E)[E$organ == "Colorectal"]
oesophageal_indices <- rownames(E)[E$organ == "Esophagus"]
E <- E[, -c(1,2)]
dim(E)

# Use oesophageal to construct the prior
Eposterior <- E[oesophageal_indices, ]
Eposterior <- Eposterior[, colSums(Eposterior)>0]
Eposterior <- Eposterior[-c(42),]
Eprior <- E[colorectal_indices, ]
Eprior <- Eprior[, colSums(Eprior)>0]
Eprior[,"SBS93"] <- Eprior[,"SBS93"] + Eprior[,"SBS121"]
Eprior <- Eprior[,-c(10)]
colnames(Eprior)[8] <- "SBS93+121"
# head(E)
# Impute zeros with multiplicative replacement strategy
W <- diag(1/rowSums(Eprior)) %*% as.matrix(Eprior)
for (i in 1:nrow(W)) {
  iszero <- unname(which(W[i, ] == 0))
  W[i, iszero] <- 0.0001
  W[i, -iszero] <- (1-0.0001*length(iszero))*W[i, -iszero] 
}

Walr <- alrTransformMatrix(W)
eta <- colMeans(Walr)

set.seed(1)
# Read counts matrix (for colorectal patients)
M <- read.csv("data/degasperi2022/counts_matrix.csv", header = TRUE, row.names = 1)
M <- M[oesophageal_indices, ]
dim(M)
G <- 10
M <- M[1:G, ]

channels <- c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT", "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT", "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT", "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT", "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT", "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT", "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT", "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT", "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT", "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT", "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT", "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT")

S <- read.csv(
  "data/degasperi2022/signatures.csv",
  header = TRUE, row.names = 1
)
rownames(S)
S <- S[-c(1:9), ] 
colnames(Eprior)
S <- S[c(2,4,6,5,3,7,8,1,9), ]
Eprior <- Eprior[, rev(colnames(Eprior))]
S <- S[rev(rownames(S)), ]
K <- ncol(S)
N <- nrow(S)

library(rstan)
sigs_prior_vars <- list()
ScaleMatrix <- matrix(1, nrow = K-1, ncol = K-1) + diag(rep(1, K-1))  # Matrix with 2s in the diagonal and 1s in the other elements.
for (i in 1:nrow(S)) {
  sigs_prior_vars[[i]] <- 0.1 * ScaleMatrix
}
sigs_prior_vars <- abind::abind(sigs_prior_vars, along = 3)
sigs_prior_vars <- aperm(sigs_prior_vars, c(3, 1, 2))

options(mc.cores = parallel::detectCores())

patient_index <- as.integer(args[1])
S <- as.matrix(S); M <- as.matrix(M)
dat <- list(
  N = N, K = K, G = G,
  M = unname(unlist(as.vector(M[patient_index, ]))),
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
stan_trace(fit, pars = "alrE")
# saveRDS(fit, paste0("results_dissertation/degasperi/fit_simulated", patient_index, ".rds"))
saveRDS(fit, paste0("results/real/fit_simulated", patient_index, ".rds"))

