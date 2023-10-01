library(ggplot2)
channels <- c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT", "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT", "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT", "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT", "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT", "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT", "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT", "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT", "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT", "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT", "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT", "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT")
sig_colors_values <- c('#F8766D', '#C49A00', '#53B400', '#00B6EB', '#00C094', '#FB61D7')
sig_colours <- rep(sig_colors_values, each = 16)

#' Plots point estimates and credible intervals of mutational
#' @param posterior_draws posteiror draws obtained upon
#' executing "extract" on a stan object
#' @param Phi An array of size V x K x S, where V is the number of categories,
#' K is the number of clusters in the best partition and S is the number of
#' iterations of the MCMC sampler
#' @param Phi A 3D array of dimension V x K x S, where
#' V is the number of categories, K is the number of clusters (signatures) and
#' S is the number of posterior samples
#'
#' @export
plot_signature_posterior <- function(posterior_draws, sig_matrix) {
  sig_names <- rownames(sig_matrix)
  Phi <- aperm(posterior_draws$S, c(3,2,1))
  # Phi is an array of dimension 96 x K x S
  K <- dim(Phi)[2] # Number of signatures
  
  ### contains_zero <- min(Z) == 0
  plots <- list()
  for(i in K:1) {
    
    prior_simulations <- alrInvTransformMatrix(rmvnorm(dim(Phi)[3], mean = alr(sig_matrix[i,]), sigma = 0.1*(diag(95) + matrix(1, 95, 95))))
    colnames(prior_simulations) <- channels
    prior_simulations <- coda::as.mcmc(prior_simulations)
    plot_prior <- ggmcmc::ggs_caterpillar(
      ggmcmc::ggs(prior_simulations, keep_original_order = TRUE),
      horizontal = FALSE,
      sort = FALSE, thick_ci = c(0.025, 0.975), thin_ci = c(0.005,0.995)) + ggplot2::aes(color = as.factor(rep(1:6, each = 16))) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_point(size = 9) +
      ggplot2::scale_colour_manual(values = sig_colors_values) +
      ggplot2::ggtitle(paste("Prior over Signature", sig_names[i])) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 32), 
                     axis.text.y = ggplot2::element_text(size = 24), 
                     axis.text.x = ggplot2::element_text(angle = 90, size = 24)) + 
      geom_segment(aes(x = Low, xend = High), linewidth = 3) + 
      geom_segment(aes(x = low, xend = high), linewidth = 1)# +
      #annotate("point", x = sig_matrix[i,], y = 1:96, color = "black", size = 6, shape = 4)
    
    signatures <- coda::as.mcmc(t(Phi[,i,]))
    colnames(signatures) <- channels
    # png(paste0(dst, "/", "signature", i, ".png"), 1200, 300)
    plot_posterior <- ggmcmc::ggs_caterpillar(
      ggmcmc::ggs(signatures, keep_original_order = TRUE),
      horizontal = FALSE,
      sort = FALSE, thick_ci = c(0.025, 0.975), thin_ci = c(0.005,0.995)) + ggplot2::aes(color = as.factor(rep(1:6, each = 16))) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_point(size = 9) +
      ggplot2::scale_colour_manual(values = sig_colors_values) +
      ggplot2::ggtitle(paste("Posterior over Signature", sig_names[i])) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 32), 
                     axis.text.y = ggplot2::element_text(size = 24), 
                     axis.text.x = ggplot2::element_text(angle = 90, size = 24))  + 
      geom_segment(aes(x = Low, xend = High), linewidth = 3) + 
      geom_segment(aes(x = low, xend = high), linewidth = 1) # +
      #annotate("point", x = sig_matrix[i,], y = 1:96, color = "black", size = 6, shape = 4)
    # gridExtra::grid.arrange(grobs = list(plot_prior, plot_posterior), nrow = 2)
    # dev.off()
    plots[[i]] <- gridExtra::grid.arrange(grobs = list(plot_prior, plot_posterior), nrow = 2)
  }
  return(plots)
}

# Plots 
plot_exposure_posterior <- function(posterior_draws, patient_name, sigs_names, true_exposures = NULL) {
  Phi <- posterior_draws$E
  exposures <- coda::as.mcmc(Phi)
  colnames(exposures) <- sigs_names
  exposures <- exposures[, rev(colnames(exposures))]
  myplot <- ggmcmc::ggs_caterpillar(
    ggmcmc::ggs(exposures, keep_original_order = TRUE), 
      horizontal = FALSE, sort = FALSE, thick_ci = c(0.025, 0.975), thin_ci = c(0.005,0.995)) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") +
      ggplot2::geom_point(size = 9) +      ggplot2::aes(color = 'F8766D') +
      ggplot2::ggtitle(patient_name) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 32), 
                     axis.text.y = ggplot2::element_text(size = 24), 
                     axis.text.x = ggplot2::element_text(angle = 90, size = 24)) + 
    geom_segment(aes(x = Low, xend = High), linewidth = 3) + 
    geom_segment(aes(x = low, xend = high), linewidth = 1)
  return(myplot)
}

plot_catalogue_predictive <- function(posterior_draws, patient_name, M_obs) {
  Phi <- posterior_draws$M_pred
  catalogue <- coda::as.mcmc(Phi)
  colnames(catalogue) <- channels
  # catalogue <- catalogue[, rev(colnames(catalogue))]
  myplot <- ggmcmc::ggs_caterpillar(
    ggmcmc::ggs(catalogue, keep_original_order = TRUE),
    horizontal = FALSE,
    sort = FALSE, thick_ci = c(0.025, 0.975), thin_ci = c(0.005,0.995)) + ggplot2::aes(color = as.factor(rep(1:6, each = 16))) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "none") + 
    
    ggplot2::geom_point(size = 9) +
    ggplot2::scale_colour_manual(values = sig_colors_values) +
    ggplot2::ggtitle(paste("Predictive Posterior Distribution for", patient_name)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 32), 
                   axis.text.y = ggplot2::element_text(size = 24), 
                   axis.text.x = ggplot2::element_text(angle = 90, size = 24)) +
    geom_segment(aes(x = Low, xend = High), linewidth = 3) + 
    geom_segment(aes(x = low, xend = high), linewidth = 1) +
    annotate("point", x = as.numeric(M_obs), y = 1:96, color = "black", size = 10, shape = 4, stroke = 2)
  return(myplot)
}

plot_catalogue_predictive_similarities <- function(posterior_draws, patient_names, M_obs) {
  Phi <- posterior_draws$M_pred
  catalogue <- coda::as.mcmc(Phi)
  colnames(catalogue) <- channels
  # catalogue <- catalogue[, rev(colnames(catalogue))]
  similarities <- sapply(1:dim(Phi)[1], function(s) cosine(Phi[s,], as.numeric(M_obs)))
  similarities <- coda::as.mcmc(as.matrix(similarities))
  colnames(similarities) <- c("Cosine")
  # png(paste0(dst, "/", "signature", i, ".png"), 1200, 300)
  myplot <- ggmcmc::ggs_histogram(
    ggmcmc::ggs(similarities, keep_original_order = TRUE)) +
    ggplot2::ggtitle(paste("Predictive Posterior Distribution for", patient_name, " (Cosine similarity)")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 32), axis.text.y = ggplot2::element_text(size = 16))
  return(myplot)
}
