#' Compute the posterior expected Binder loss (EB) for a partition
#'
#' @param cls A vector containing the cluster assignments of $n$ data points, representing the partition with respect to which the EVI is computed.
#' @param cls.draw A matrix of the MCMC samples of partitions of $n$ data points, of dimensions $S × n$, where $S$ is the number of MCMC samples.
#' @param Ks The number of clusters in \code{cls}
#' @param Ks.draw A vector of $S$, containing the number of clusters for each of the MCMC samples.
#' @param a A parameter used in generalized Binder, takes value between 0 and 2 and 1 by default.
#' @return The posterior expected Binder for the partition \code{cls}, where the posterior is approximated by the MCMC samples \code{cls.draw}. It corresponds to the average Binder distance between \code{cls} and each MCMC sample in \code{cls.draw}.
EB_Rcpp <- function(cls, cls.draw, Ks, Ks.draw, a = 1) {
  if (is.vector(cls)) cls <- t(cls)
  S <- dim(cls.draw)[1]

  output <- sum(Binder_Rcpp(cls, cls.draw, Ks, Ks.draw, a)) / S
  return(output)
}
