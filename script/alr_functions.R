alr <- function(z) {
  D <- length(z)
  return(log(z[1:(D-1)]) - log(z[D]))
}

alrInv <- function(phi) {
  D <- length(phi) + 1
  denom <- sum(exp(phi)) + 1
  return(c(exp(phi), 1)/denom)
}

alrInv1d <- function(phi) {
  return(exp(phi)/(exp(phi)+1))
}

alrTransformMatrix <- function(mat) {
  alrMat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat)-1)
  for(i in 1:nrow(mat)) {
    alrMat[i, ] <- alr(mat[i, ])
  }
  return(alrMat)
}

alrInvTransformMatrix <- function(mat) {
  alrInvMat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat) + 1)
  for (i in 1:nrow(mat)) {
    alrInvMat[i, ] <- alrInv(mat[i, ])
  }
  return(alrInvMat)
}

closure <- function(x) {
  return(x/sum(x))
}

perturbation <- function(x, z) {
  return(closure(x*z))
}

aitchisonInnerProd <- function(x, z, alr = FALSE) {
  D <- length(x)
  if (alr == FALSE) {
    x <- alr(x)
    z <- alr(z)
  }
  innerProdMatrix <- diag(D-1) - 1/D*matrix(1, D-1, D-1)
  
  return(x %*% innerProdMatrix %*% z)
}

aitchisonDistance <- function(x, z) {
  perturb <- x * (1/z)
  # perturb <- perturbation(x, 1/z)
  return(
    sqrt(
      aitchisonInnerProd(perturb, perturb, alr = FALSE)
    )
  )
}

