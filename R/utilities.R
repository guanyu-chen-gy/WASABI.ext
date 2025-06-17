#' Relabel partitions
#' 
#' @description Change cluster labels so that they are ordered from 1 to K, 
#' where K is the number of clusters in the partition.
#' 
#'
#' @param c A vector of cluster labels of size $n$, where each label represents a cluster.
#'
#' @returns A vector of the same size as \code{c}, where each cluster label is replaced with a number from 1 to K, where K is the number of unique clusters in \code{c}.
#' @export
#'
#' @examples{
#' # assuming cls.draw is a matrix of size S x n, 
#' # where S is the number of samples and n is the sample size
#' # i.e. each row contains a partition of n data points
#' cls.draw = t(apply(cls.draw,1,relabel_partition)) 
#' }
relabel_partition <- function(c) {
  uu1 <- unique(c)
  c2 <- c
  for (i in 1:length(uu1)) {
    c2[c == uu1[i]] <- i
  }
  c2
}

#' Reorder cluster labels and find unique clusterings
#'
#' Reorder cluster labels in a collection of samples, so that cluster labels appear as 1,2,3,.. and to make
#' finding unique clusters more efficient.
#'
#' @param S Integer, number of samples.
#' @param n Integer, sample size.
#' @param config A matrix of size \code{S} x \code{n}, where each row represents a sample configuration.
#'
#' @returns A list containing:
#' \begin{itemize}
#' \item \code{config_reorder} A matrix of size \code{S} x \code{n}, where each row represents a reordered sample configuration.
#' \item \code{config_count} A vector containing the count of each unique configuration (number of times each configuration in config_reorder appears in the original config.
#' \item \code{config_index} A vector containing the indices of unique configurations.
#' \end{itemize}
#'
reorder_find_unique <- function(S, n, config) {
  # Initialize output
  config_reorder <- matrix(0, S, n)
  config_count <- c(1)
  config_index <- c(S)

  k <- rep(0, S)
  for (s in 1:S) {
    k[s] <- length(unique(config[s, ]))
  }

  for (s in S:1) {
    # reorder the configuration
    uniq_s <- unique(config[s, ])
    for (h in 1:k[s]) {
      config_reorder[s, config[s, ] == uniq_s[h]] <- h
    }
    if (s < S) {
      # check if we have seen this configuration before
      before <- colSums(matrix(t(config_reorder[config_index, ]), nrow = n) == matrix(config_reorder[s, ], n, length(config_index))) == n
      if (sum(before) > 0) {
        config_count[before] <- config_count[before] + 1
      } else {
        config_count <- c(config_count, 1)
        config_index <- c(config_index, s)
      }
    }
  }


  output <- list(config_reorder = config_reorder, config_count = config_count, config_index = config_index)
  return(output)
}

## -- currently not used --

sample_max_jit <- function(pr, quant = 0.95) {
  q95 <- quantile(pr, prob = quant)
  ind <- which(pr > q95)
  pr0 <- pr[ind]
  jit <- rnorm(length(pr0)) * diff(range(pr0)) / 2
  pr1 <- pr0 + jit
  return(ind[which.max(pr1)])
}

get_wassersteinVI <- function(cls.draw = NULL,
                              particles,
                              suppress.comment = TRUE,
                              return_psm = FALSE,
                              seed = NULL) {
  if (is.null(cls.draw)) {
    stop("cls.draw must be provided")
  }

  if (dim(particles)[2] != dim(cls.draw)[2]) {
    stop("number of data points must be the same in initial particles and draws")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cls.draw_relab <- t(apply(cls.draw, 1, relabel_partition)) - 1
  Ks.draw <- apply(cls.draw_relab, 1, function(x) max(x)) + 1
  part <- particles
  part_relab <- t(apply(part, 1, relabel_partition)) - 1
  Ks.part <- apply(part_relab, 1, function(x) max(x)) + 1

  n <- dim(cls.draw)[2]
  S <- dim(cls.draw)[1]
  L <- nrow(particles)
  part.evi <- rep(0, L)

  if (L == 1) {
    part.evi[1] <- EVI_Rcpp(
      cls = part_relab[1, ],
      cls.draw = cls.draw_relab,
      Ks = Ks.part[1],
      Ks.draw = Ks.draw
    )

    output <- list(
      particles = part, EVI = part.evi, wass.dist = sum(part.evi),
      part.psm = NULL, part.weights = 1, draws.assign = NULL
    )
  } else {
    viall <- VI_Rcpp(cls.draw_relab, part_relab, Ks.draw, Ks.part)
    viall[viall < 0] <- 0 # this is to correct numerical errors

    # Assign each sample to nearest center
    # if a partition is equally distant from two particles randomize
    assign.vi <- apply(viall, 1, function(x) {
      indices <- which(x == min(x))
      if (length(indices) > 1) indices <- sample(x = indices, size = 1)
      return(indices)
    })
    # if a particle has zero assigments, allocate closet sample to it
    counts <- apply(matrix(c(1:L), L, 1), 1, function(x) sum(assign.vi == x))

    for (l in 1:L) {
      part.evi[l] <- sum(viall[assign.vi == l, l]) / counts[l]
    }

    assign.vi_ordered <- match(1:L, order(counts, decreasing = T))[assign.vi]
    tmp_order <- sort(counts, index.return = TRUE, decreasing = TRUE)
    tmp_order <- tmp_order$ix

    part_ordered <- part[tmp_order, ]
    part.evi_ordered <- part.evi[tmp_order]
    counts_ordered <- counts[tmp_order]

    psm_K_ordered <- NULL
    if (return_psm) {
      for (l in 1:L) {
        psm_K_ordered[[l]] <- mcclust::comp.psm(matrix(cls.draw[assign.vi_ordered == l, ], counts_ordered[l], n))
      }
    }

    output <- list(
      particles = part_ordered, EVI = part.evi_ordered, wass.dist = sum(part.evi_ordered * counts_ordered / S),
      part.psm = psm_K_ordered, part.weights = counts_ordered / S, draws.assign = assign.vi_ordered
    )
  }
  return(output)
}
