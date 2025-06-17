#' Compute the meet
#'
#' Given a set of clusterings/partitions, finds the "meet" across them
#' (the common refinement, or the finest common coarsening).
#'
#' @param cls.part Integer matrix of dimension K (number of particles) by n (number of observations),
#'        with each row representing a clustering/partition.
#' @return A list with components:
#'   \describe{
#'     \item{cls.m}{Length-n integer vector: cluster assignment for the meet.}
#'     \item{k.m}{Integer: number of clusters in the meet.}
#'     \item{cls.sizes.m}{Integer vector: sizes of clusters in the meet.}
#'   }
#' @examples
#' cls.part <- matrix(c(1,1,2,2, 1,2,1,2), nrow = 2, byrow = TRUE)
#' meet <- WASABI::cls.meet(cls.part)
#' @export
cls.meet <- function(cls.part) {
  K <- base::dim(cls.part)[1]
  n <- base::dim(cls.part)[2]
  cls.list <- base::list()
  for (k in 1:K) {
    cls.list[[k]] <- cls.part[k, ]
  }
  tb <- base::table(cls.list)
  cls.m <- base::rep(0, n)
  ind <- tb > 0
  k.m <- base::sum(ind)
  cls.sizes.m <- tb[ind]
  cls.m.lbs <- base::which(ind, arr.ind = TRUE)
  for (j in 1:k.m) {
    ind.k <- base::rep(TRUE, n)
    for (k in 1:K) {
      ind.k <- ind.k & (cls.part[k, ] == cls.m.lbs[j, k])
    }
    cls.m[ind.k] <- j
  }
  base::list(cls.m = cls.m, k.m = k.m, cls.sizes.m = cls.sizes.m)
}

#' WASABI pairwise similarity matrix (PSM) for the meet partition
#'
#' Given a meet partition and a set of particles (clusterings), calculates the WaSABI-approximated PSM for the meet's clusters.
#'
#' @param cls.m Integer vector: the meet clustering (as from \code{WASABI::cls.meet}).
#' @param output_wvi List with elements \code{particles} (matrix of clusterings) and \code{part.weights} (vector of weights).
#' @return Numeric square matrix: weighted PSM for the meet.
#' @examples
#' \dontrun{
#' # Assuming output_wvi with fields particles, part.weights:
#' pm <- WASABI::psm.meet(meet$cls.m, output_wvi)
#' }
#' @export
psm.meet <- function(cls.m, output_wvi) {
  cls.part <- output_wvi$particles
  cls.w <- output_wvi$part.weights
  k.m <- base::length(base::unique(cls.m))
  psm <- base::matrix(0, k.m, k.m)
  for (i in 1:k.m) {
    for (j in 1:k.m) {
      ind <- base::apply(cls.part, 1, function(x) {
        base::max(x[cls.m == i]) == base::max(x[cls.m == j])
      })
      psm[i, j] <- base::sum(ind * cls.w)
    }
  }
  psm
}

#' Expected VI Contribution
#'
#' Calculates the EVI contribution of each observation to the expected variation of information (VI) distance between a partition and samples.
#'
#' @param Zmat Integer matrix: S rows (samples) by n columns (data points), each row a clustering.
#' @param Zhat Integer vector: length-n partition to compare.
#' @return Numeric vector of length n: expected VI contribution for \code{Zhat} (each element).
#' @export
evi.contribution <- function(Zmat, Zhat) {
  n <- base::length(Zhat)
  S <- base::dim(Zmat)[1]
  evi.i <- function(i, Zmat, Zhat) {
    c <- 1 / n * log2(base::sum(Zhat == Zhat[i]) / n) +
      base::sum(base::apply(Zmat, 1, function(x) {
        1 / n * log2(base::sum(x == x[i]) / n)
      })) / S -
      2 / n * base::sum(base::apply(Zmat, 1, function(x, y) {
        log2(base::sum((x == x[i]) & (y == y[i])) / n)
      }, y = Zhat)) / S
    c
  }
  base::unlist(base::lapply(seq_len(n), evi.i, Zmat = Zmat, Zhat = Zhat))
}

#' CExpected VI Contribution (WASABI approximation)
#'
#' Calculates the EVI contribution (with expectation taken with respect to the WASABI distribution) of each observation.
#'
#' @param output_wvi List: must contain \code{particles} (matrix) and \code{part.weights} (numeric vector).
#' @param Zhat Integer vector: partition to compare.
#' @return Numeric vector of length n: weighted expected VI contribution for \code{Zhat} (each element).
#' @export
evi.wd.contribution <- function(output_wvi, Zhat) {
  n <- base::length(Zhat)
  evi.i <- function(i, Zmat, Zhat, w) {
    c <- 1 / n * log2(base::sum(Zhat == Zhat[i]) / n) +
      base::sum(w * base::apply(Zmat, 1, function(x) {
        1 / n * log2(base::sum(x == x[i]) / n)
      })) -
      2 / n * base::sum(w * base::apply(Zmat, 1, function(x, y) {
        log2(base::sum((x == x[i]) & (y == y[i])) / n)
      }, y = Zhat))
    c
  }
  base::unlist(base::lapply(seq_len(n), evi.i, Zmat = output_wvi$particles, Zhat = Zhat, w = output_wvi$part.weights))
}

#' VI contribution between two partitions
#'
#' Calculates the individual contributions of observations to the VI distance between partitions Z1 and Z2.
#'
#' @param Z1 Integer or factor vector: length-n cluster labels for partition 1.
#' @param Z2 Integer or factor vector: length-n cluster labels for partition 2.
#' @return Numeric vector of length n: VI contribution for each observation.
#' @export
vi.contribution <- function(Z1, Z2) {
  n <- base::length(Z1)
  vi.i <- function(i, Z1, Z2) {
    c <- 1 / n * log2(base::sum(Z1 == Z1[i]) / n) +
      1 / n * log2(base::sum(Z2 == Z2[i]) / n) -
      2 / n * log2(base::sum((Z1 == Z1[i]) & (Z2 == Z2[i])) / n)
    c
  }
  base::unlist(base::lapply(seq_len(n), vi.i, Z1 = Z1, Z2 = Z2))
}
