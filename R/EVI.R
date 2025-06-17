#' Compute the posterior expected Variation of Information (EVI) for a partition
#'
#' @param cls A vector containing the cluster assignments of $n$ data points, representing the partition with respect to which the EVI is computed.
#' @param cls.draw A matrix of the MCMC samples of partitions of $n$ data points, of dimensions $S \times n$, where $S$ is the number of MCMC samples.
#' @param Ks The number of clusters in \code{cls}
#' @param Ks.draw A vector of $S$, containing the number of clusters for each of the MCMC samples.
#' @return The posterior expected VI for the partition \code{cls}, where the posterior is approximated by the MCMC samples \code{cls.draw}. It corresponds to the average VI distance between \code{cls} and each MCMC sample in \code{cls.draw}.
EVI_Rcpp <- function(cls, cls.draw, Ks, Ks.draw) {
  if (is.vector(cls)) cls <- t(cls)
  S <- dim(cls.draw)[1]

  output <- sum(VI_Rcpp(cls, cls.draw, Ks, Ks.draw)) / S
  return(output)
}
