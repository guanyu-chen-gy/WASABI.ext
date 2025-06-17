# Create a function to find the meet
# Input: cls.part is an K (no. of particles) by n (no of observations) matrix,
#        where each row represents a clustering/particle
# Output: cls.m is the meet of the particles, k.m is the no. of clusters in the
#        meet, and cls.sizes.m contains the sizes of the clusters in the meet
cls.meet <- function(cls.part) {
  # Create list
  K <- dim(cls.part)[1]
  n <- dim(cls.part)[2]
  cls.list <- list()
  for (k in 1:K) {
    cls.list[[k]] <- cls.part[k, ]
  }

  # Cross tabulate
  tb <- table(cls.list)
  cls.m <- rep(0, n) # meet
  ind <- tb > 0 # indexes non-empty intersections
  k.m <- sum(ind) # no. of clusters in the meet
  cls.sizes.m <- tb[ind] # cluster sizes of the meet
  cls.m.lbs <- which(ind, arr.ind = TRUE) # labels in each particle
  for (j in 1:k.m) {
    ind.k <- rep(TRUE, n)
    for (k in 1:K) {
      ind.k <- ind.k & (cls.part[k, ] == cls.m.lbs[j, k])
    }
    cls.m[ind.k] <- j
  }
  return(list(cls.m = cls.m, k.m = k.m, cls.sizes.m = cls.sizes.m))
}

psm.meet <- function(cls.m, output_wvi) {
  cls.part <- output_wvi$particles
  cls.w <- output_wvi$part.weights
  k.m <- length(unique(cls.m))
  psm <- matrix(0, k.m, k.m)
  for (i in 1:k.m) {
    for (j in 1:k.m) {
      ind <- apply(cls.part, 1, function(x) {
        max(x[cls.m == i]) == max(x[cls.m == j])
      }) # max is used as "unique" here (but more efficient)
      psm[i, j] <- sum(ind * cls.w)
    }
  }
  return(psm)
}

# Contribution of each data to the expected VI distance of a partition Zhat
evi.contribution <- function(Zmat, Zhat) {
  n <- length(Zhat)
  S <- dim(Zmat)[1]
  evi.i <- function(i, Zmat, Zhat) {
    c <- 1 / n * log2(sum(Zhat == Zhat[i]) / n) + sum(apply(Zmat, 1, function(x) {
      1 / n * log2(sum(x == x[i]) / n)
    })) / S - 2 / n * sum(apply(Zmat, 1, function(x, y) {
      log2(sum((x == x[i]) & (y == y[i])) / n)
    }, y = Zhat)) / S
    return(c)
  }
  contrib <- unlist(lapply(c(1:n), evi.i, Zmat = Zmat, Zhat = Zhat))
  return(contrib)
}

# Contribution of each data to the expected VI distance of a partition Zhat
# Expectation is computed from the wasserstein approximation
evi.wd.contribution <- function(output_wvi, Zhat) {
  n <- length(Zhat)
  evi.i <- function(i, Zmat, Zhat, w) {
    c <- 1 / n * log2(sum(Zhat == Zhat[i]) / n) + sum(w * apply(Zmat, 1, function(x) {
      1 / n * log2(sum(x == x[i]) / n)
    })) - 2 / n * sum(w * apply(Zmat, 1, function(x, y) {
      log2(sum((x == x[i]) & (y == y[i])) / n)
    }, y = Zhat))
    return(c)
  }
  contrib <- unlist(lapply(c(1:n), evi.i, Zmat = output_wvi$particles, Zhat = Zhat, w = output_wvi$part.weights))
  return(contrib)
}

# Contribution of each data to the VI distance between two partitions
vi.contribution <- function(Z1, Z2) {
  n <- length(Z1)
  vi.i <- function(i, Z1, Z2) {
    c <- 1 / n * log2(sum(Z1 == Z1[i]) / n) + 1 / n * log2(sum(Z2 == Z2[i]) / n) - 2 / n * log2(sum((Z1 == Z1[i]) & (Z2 == Z2[i])) / n)
    return(c)
  }
  contrib <- unlist(lapply(c(1:n), vi.i, Z1 = Z1, Z2 = Z2))
  return(contrib)
}
