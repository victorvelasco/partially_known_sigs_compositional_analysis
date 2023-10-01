rm(list = ls())
library(mvtnorm)
library(ggmcmc)

set.seed(1)

args = commandArgs(trailingOnly=TRUE)


setwd("/home/victor/phd/comp_receptor_modelling/")
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

##### ### Plot prior weights
##### SigmaW <- 10*(matrix(1, nrow = ncol(Walr), ncol = ncol(Walr)) + diag(ncol(Walr)))
##### B <- 100000
##### wStar <- rmvnorm(B, mean = eta, sigma = SigmaW)
##### dim(wStar)
##### wStar <- alrInvTransformMatrix(wStar)
##### distances <- apply(wStar, 1, function(w) aitchisonDistance(w, alrInv(eta)))
##### plot(density(distances), main = "Prior distance from the prior centre to the vector of weights")
##### points(apply(W, 1, function(w) aitchisonDistance(w, alrInv(eta))), rep(0, nrow(W)))





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

source("script/plot_posterior_compositional.R")
# S <- S[rownames(SSimulated), ]
plots_posterior_weights <- list()
plots_posterior_predictive <- list()
for (patient_index in c(1:5, 7:10)) {
  
  posterior_draws <- readRDS(paste0("results_dissertation/degasperi/fit_simulated", patient_index, ".rds"))
  posterior_draws <- rstan::extract(posterior_draws)
  signatures_plots <- plot_signature_posterior(posterior_draws, as.matrix(S))
  for (sig in nrow(S):1) {
    path <- paste0("/home/victor/phd/dissertation/images/ch4/degasperi/", patient_index, "/", rownames(S)[sig], ".png")
    png(path, width = 2000, height = 1000)
    ggsave(path, signatures_plots[[sig]])
    dev.off()
  }
  
  ## plots_posterior_weights[[patient_index]] <- plot_exposure_posterior(posterior_draws, paste("Patient", patient_index), colnames(Eprior))
  ## 
  ## plots_posterior_predictive[[patient_index]] <- plot_catalogue_predictive(posterior_draws, paste("Patient", patient_index), M[patient_index, ])
}

library(gridExtra)
path <- paste0("/home/victor/phd/dissertation/images/ch4/degasperi/", "exposures", ".png")
png(path, width = 2000, height = 2500)
ggsave(path, grid.arrange(grobs = plots_posterior_weights[c(1:5,7:10)], ncol = 2))
dev.off()

path <- paste0("/home/victor/phd/dissertation/images/ch4/degasperi/", "predictive1", ".png")
png(path, width = 2000, height = 2500)
ggsave(path, grid.arrange(grobs = plots_posterior_predictive[1:5], ncol = 1))
dev.off()
path <- paste0("/home/victor/phd/dissertation/images/ch4/degasperi/", "predictive2", ".png")
png(path, width = 2000, height = 2000)
ggsave(path, grid.arrange(grobs = plots_posterior_predictive[7:10], ncol = 1))
dev.off()

