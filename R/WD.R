#' Computes the Wassertein distance between the empirical posterior and the WASABI distribution
#'
#' @param cls.draw A matrix of the MCMC samples of partitions of $n$ data points.
#' @param part The WASABI particles. A matrix of dimension $L × n$ containing the particles obtained by WASABI in each of the rows.
#' @param part.weights A vector of length $L$ containing the weights associated to each of the particles.
#' @param assign.vi A vector containing the assignment of each sample \code{cls.draw} to a region of attraction (and the corresponding particle). It is obtained as output from WASABI.
#' @return The Wasserstein-VI distance between the empirical posterior distribution (supported on the MCMC samples) and the WASABI distribution.
#' @export
WD <- function(cls.draw, part, part.weights, assign.vi) {
  K <- dim(part)[1]
  part.evi <- sapply(c(1:K), function(k) {
    EVI_Rcpp(
      cls = (part[k, , drop = FALSE] - 1),
      cls.draw = (cls.draw[assign.vi == k, , drop = FALSE] - 1),
      Ks = max(part[k, ]),
      Ks.draw = apply(cls.draw[assign.vi == k, , drop = FALSE], 1, max)
    )
  })
  wass.dist <- sum(part.evi * part.weights)
}
