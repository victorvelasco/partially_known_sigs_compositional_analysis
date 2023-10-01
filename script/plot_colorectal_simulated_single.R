library(mvtnorm)
set.seed(1)

source("script/alr_functions.R")

S <- read.csv(
  "data/alexandrov2020/sigProfiler_SBS_signatures_2019_05_22.csv",
  header = TRUE, row.names = 1
)
S <- t(as.matrix(S))
S <- S[colnames(E), ]
K <- ncol(S)

SSimulated <- read.csv("data/simulated_colorectal/simulated_sigs.csv", header = TRUE, row.names = 1)
Esimulated <- read.csv("data/simulated_colorectal/simulated_exposures.csv", header = TRUE, row.names = 1)
Msimulated <- read.csv("data/simulated_colorectal/simulated_colorectal_10genomes.csv", header = TRUE, row.names = 1)
N <- nrow(SSimulated) # Number of signatures

library(rstan)
library(bayesplot)
library(lsa)


sigs_prior_vars <- list()
ScaleMatrix <- matrix(1, nrow = K-1, ncol = K-1) + diag(rep(1, K-1))  # Matrix with 2s in the diagonal and 1s in the other elements.
for (i in 1:nrow(S)) {
  sigs_prior_vars[[i]] <- 0.1 * ScaleMatrix
}
sigs_prior_vars <- abind::abind(sigs_prior_vars, along = 3)
sigs_prior_vars <- aperm(sigs_prior_vars, c(3, 1, 2))


source("script/plot_posterior_compositional.R")
S <- S[rownames(SSimulated), ]
plots_posterior_weights <- list()
plots_posterior_predictive <- list()
for (patient_index in 1:10) {
  
  posterior_draws <- readRDS(paste0("results//simulated/fit_simulated", patient_index, ".rds"))
  posterior_draws <- rstan::extract(posterior_draws)
  ## signatures_plots <- plot_signature_posterior(posterior_draws, as.matrix(SSimulated))
  ## for (sig in nrow(S):1) {
  ##   path <- paste0("/home/victor/phd/dissertation/images/ch4/simulated/", patient_index, "/", rownames(S)[sig], ".png")
  ##   png(path, width = 2000, height = 1000)
  ##   ggsave(path, signatures_plots[[sig]])
  ##   dev.off()
  ## }
  
  plots_posterior_weights[[patient_index]] <- plot_exposure_posterior(posterior_draws, paste("Patient", patient_index), rownames(S)) +
    annotate("point", x = as.numeric(Esimulated[patient_index, ]/sum(Esimulated[patient_index, ])), y = ncol(Esimulated):1, color = "black", size = 10, shape = 4, stroke = 2)
  
  plots_posterior_predictive[[patient_index]] <- plot_catalogue_predictive(posterior_draws, paste("Patient", patient_index), Msimulated[patient_index, ])
}
  
## library(gridExtra)
## path <- paste0("/home/victor/phd/dissertation/images/ch4/simulated/", "exposures", ".png")
## png(path, width = 2000, height = 2500)
## ggsave(path, grid.arrange(grobs = plots_posterior_weights, ncol = 2))
## dev.off()
## 
## path <- paste0("/home/victor/phd/dissertation/images/ch4/simulated/", "predictive1", ".png")
## png(path, width = 2000, height = 2500)
## ggsave(path, grid.arrange(grobs = plots_posterior_predictive[1:5], ncol = 1))
## dev.off()
## path <- paste0("/home/victor/phd/dissertation/images/ch4/simulated/", "predictive2", ".png")
## png(path, width = 2000, height = 2500)
## ggsave(path, grid.arrange(grobs = plots_posterior_predictive[6:10], ncol = 1))
## dev.off()

