particle_search <- function(cls.draw_relab, Ks.draw,
                            part_relab, Ks.part, part.evi,
                            L, method, swap_countone, max.k, lb,
                            suppress.comment,
                            first_iter = FALSE, ...) {
  S <- nrow(cls.draw_relab)

  # ------------------------------------------ N-update step ------------------------------------------

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
    if (counts[l] == 0) {
      ind <- (assign.vi %in% c(1:L)[counts > 1]) # "counts>1" is to avoid some group becoming empty!
      m <- which(stats::rmultinom(1, 1, viall[ind, l] / sum(viall[ind, l])) == 1)
      l_old <- (assign.vi[ind])[m]
      assign.vi[ind][m] <- l
      counts[l] <- 1
      counts[l_old] <- counts[l_old] - 1
    }
  }

  part.evi_new <- part.evi
  if (first_iter) {
    # this is evi for the old particles with the new assignments, and it will get updated
    for (l in 1:L) {
      part.evi_new[l] <- sum(viall[assign.vi == l, l]) / counts[l]
    }
  }

  part_new_relab <- part_relab
  part_new <- part_new_relab + 1
  Ks.part_new <- Ks.part

  # Compute centre within each group (can parallelize this)
  for (l in 1:L) {
    # ------------------------------------------ Outlier-check step ------------------------------------------

    if (counts[l] == 1) {
      part_new_relab[l, ] <- cls.draw_relab[assign.vi == l, ]
      part_new[l, ] <- part_new_relab[l, ] + 1
      Ks.part_new[l] <- Ks.draw[assign.vi == l]
      part.evi_new[l] <- 0

      ### test changing center
      if (swap_countone) {
        tmp <- test_change(
          l, part_new_relab, cls.draw_relab, viall, assign.vi, counts,
          part.evi_new, Ks.draw, Ks.part_new
        )
        if (tmp$wass_dist < sum(part.evi_new * counts / S)) {
          counts <- tmp$counts
          assign.vi <- tmp$assign.vi
          part_new[l, ] <- tmp$part_relab[l, ] + 1
          part_new_relab[l, ] <- tmp$part_relab[l, ]
          Ks.part_new[l] <- max(part_new_relab[l, ]) + 1
          part.evi_new <- tmp$part.evi
        }
      }
    } else {
      # ------------------------------------------ VI-search step ------------------------------------------

      if (method == "average" | method == "complete") {
        psm_K <- mcclust::comp.psm(matrix(1 + cls.draw_relab[assign.vi == l, ], counts[l], ncol(cls.draw_relab)))
        out <- minVI_hclust(cls.draw_relab[assign.vi == l, ], Ks.draw[assign.vi == l],
          psm_K, method, max.k, lb,
          L = 1
        )
        part_new[l, ] <- out$part
        part.evi_new[l] <- out$evi
      }
      if (method == "greedy") {
        output_minepl <- GreedyEPL::MinimiseEPL(1 + cls.draw_relab[assign.vi == l, ], list(decision_init = part_new[l, ], loss_type = "VI"))
        part_new[l, ] <- output_minepl$decision
        part.evi_new[l] <- output_minepl$EPL
      }
      if (method == "salso") {
        output_salso <- salso::salso(x = cls.draw_relab[assign.vi == l, ], ...)
        part_new[l, ] <- as.numeric(output_salso)
        part.evi_new[l] <- as.numeric(attr(output_salso, "info")[4])
      }
    }
  }

  part_new_relab <- t(apply(part_new, 1, relabel_partition)) - 1
  Ks.part_new <- apply(part_new_relab, 1, function(x) max(x)) + 1

  # If equal particles, then merge equal ones and set others partitions far from unique particles
  vi.part.new <- VI_Rcpp(part_new_relab, part_new_relab, Ks.part_new, Ks.part_new)
  if (sum(vi.part.new[lower.tri(vi.part.new)] == 0) > 0) {
    if (suppress.comment == FALSE) {
      cat("Equal particles found\n")
    }
    out <- check_equal_particles(
      vi.part.new, L,
      assign.vi, counts, part.evi_new,
      cls.draw_relab, Ks.draw,
      part_new_relab, Ks.part_new
    )
    part_new_relab <- out$part_new_relab
    counts <- out$counts
    assign.vi <- out$assign.vi
    part.evi_new <- out$part.evi_new
    Ks.part_new <- out$Ks.part_new
  }


  return(list(
    part_relab = part_new_relab,
    Ks.part = Ks.part_new,
    evi = part.evi_new,
    counts = counts,
    assign.vi = assign.vi
  ))
}




check_equal_particles <- function(vi.part.new, L,
                                  assign.vi, counts, part.evi_new,
                                  cls.draw_relab, Ks.draw,
                                  part_new_relab, Ks.part_new) {
  part.ind.eq <- c(1:L)[sapply(1:L, function(k) {
    sum(vi.part.new[k, -k] == 0) > 0
  })]
  while (length(part.ind.eq) > 0) {
    # Find which particles are equal to l
    l <- part.ind.eq[1]
    part.ind.eq.l <- which(vi.part.new[l, ] == 0)
    # Move all to l except the furthest ones
    all.inds <- unlist(lapply(1:length(part.ind.eq.l), function(j) {
      which(assign.vi == part.ind.eq.l[j])
    }))
    assign.vi[all.inds] <- l
    # # Find the furthest ones
    vi.inds <- VI_Rcpp(cls.draw_relab[all.inds, ], matrix(part_new_relab[l, ], nrow = 1), Ks.draw[all.inds], Ks.part_new[l])
    svi.inds <- sort(vi.inds, decreasing = TRUE, index.return = TRUE)
    # Set duplicate particles to the furthers
    part_new_relab[part.ind.eq.l[-1], ] <- cls.draw_relab[all.inds, ][svi.inds$ix[1:(length(part.ind.eq.l) - 1)], ]
    Ks.part_new[part.ind.eq.l] <- apply(part_new_relab[part.ind.eq.l, ], 1, max) + 1
    # Update counts
    counts[l] <- sum(assign.vi == l) - length(part.ind.eq.l) + 1
    counts[part.ind.eq.l[-1]] <- 1
    # update evi
    part.evi_new[part.ind.eq.l[-1]] <- 0
    part.evi_new[l] <- EVI_Rcpp(
      cls = part_new_relab[l, ], cls.draw = cls.draw_relab[assign.vi == l, ],
      Ks = Ks.part_new[l], Ks.draw = Ks.draw[assign.vi == l]
    )
    part.ind.eq <- setdiff(part.ind.eq, part.ind.eq.l)
  }

  return(list(
    part_new_relab = part_new_relab, counts = counts, Ks.part_new = Ks.part_new,
    assign.vi = assign.vi, part.evi_new = part.evi_new
  ))
}

test_change <- function(l, part_relab, cls.draw_relab.mb,
                        viall, assign.vi, counts, part.evi,
                        Ks.draw.mb, Ks.part) {
  L <- nrow(part_relab)
  S_temp <- nrow(cls.draw_relab.mb)

  current_i <- which(assign.vi == l) # do we know for sure that l is in my minibatch?
  l_new <- (1:L)[-l][which.min(viall[current_i, -l])] # let's assign current_i to another particle (l_new)

  assign.vi[current_i] <- l_new
  counts[l_new] <- counts[l_new] + 1
  counts[l] <- counts[l] - 1 # so = 0

  # adding current_i to l_new
  VI_new_tmp <- sum(VI_Rcpp(
    part_relab[l_new, , drop = FALSE], cls.draw_relab.mb[current_i, , drop = FALSE],
    Ks.part[l_new], Ks.draw.mb[current_i]
  ))
  part.evi[l_new] <- ((counts[l_new] - 1) * part.evi[l_new] + VI_new_tmp) / counts[l_new]
  part.evi[l] <- 0 # because it's one element

  # sample a new particle
  ind <- which(assign.vi %in% c(1:L)[counts > 1]) # this is to avoid some group becoming empty!
  new_i <- sample_max_jit(apply(viall[ind, ], 1, min))
  new_i <- ind[new_i]

  l_old <- assign.vi[new_i]
  assign.vi[new_i] <- l
  counts[l] <- counts[l] + 1 # so = 1
  counts[l_old] <- counts[l_old] - 1

  # removing new_i from l_old
  VI_new_tmp <- sum(VI_Rcpp(
    part_relab[l_old, , drop = FALSE], cls.draw_relab.mb[new_i, , drop = FALSE],
    Ks.part[l_old], Ks.draw.mb[new_i]
  ))
  part.evi[l_old] <- ((counts[l_new] + 1) * part.evi[l_old] - VI_new_tmp) / counts[l_new]

  # update part_relab
  part_relab[l, ] <- cls.draw_relab.mb[new_i, ]
  Ks.part[l] <- Ks.draw.mb[new_i]

  # update viall for the l column
  viall_l <- VI_Rcpp(cls.draw_relab.mb, part_relab[l, , drop = FALSE], Ks.draw.mb, Ks.part[l])
  viall[, l] <- as.numeric(viall_l)

  wass_dist <- sum(part.evi * counts / S_temp)
  # cat("test_change: l ",l," l_new ",l_new," l_old ", l_old, "\n")
  return(list(
    wass_dist = wass_dist, counts = counts, assign.vi = assign.vi,
    part_relab = part_relab, part.evi = part.evi
  ))
}

minVI_hclust <- function(cls.draw_relab, Ks.draw,
                         psm, method,
                         max.k, lb, L = 1) {
  hclust_K <- stats::hclust(stats::as.dist(1 - psm), method = method)
  cls.hclust <- t(apply(matrix(1:max.k), 1, function(x) stats::cutree(hclust_K, k = x)))
  if (lb) {
    VI.hclust <- mcclust.ext::VI.lb(cls.hclust, psm)
  } else {
    cls.hclust_relab <- t(apply(cls.hclust, 1, relabel_partition)) - 1
    Ks.hclust <- apply(cls.hclust_relab, 1, function(x) max(x)) + 1
    VI.hclust <- sapply(1:max.k, function(i) {
      EVI_Rcpp(
        cls = cls.hclust_relab[i, ],
        cls.draw = cls.draw_relab,
        Ks = Ks.hclust[i],
        Ks.draw = Ks.draw
      )
    })
  }
  if (L > 1) { # used in initialization, returns the top L partitions minimizing VI
    sortvi <- sort(VI.hclust, index.return = T, decreasing = F)
    part <- cls.hclust[sortvi$ix[1:L], ]
    part.evi <- sortvi$x[1:L]
  } else {
    ind <- which.min(VI.hclust)
    part <- cls.hclust[ind, , drop = FALSE]
    part.evi <- VI.hclust[ind]
  }
  return(list(part = part, evi = part.evi))
}
