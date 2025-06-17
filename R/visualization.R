# visualization functions

require(ggplot2)
# ggrange <- function(part, x, ind = NULL, K = 10){
#   if(is.null(ind)){
#     ind = list(ix = 1:K)
#   }
#   ka = c()
#   cla = c()
#   meanxa = c()
#   minxa = c()
#   maxxa = c()
#   for (k in 1:K){
#     configk=(part[ind$ix,])[k,]
#     nc = max(configk)
#     for (c in 1:nc){
#       ka = c(ka ,k)
#       cla = c(cla ,c)
#       meanxa = c(meanxa,mean(x[configk==c]))
#       minxa = c(minxa,min(x[configk==c]))
#       maxxa = c(maxxa,max(x[configk==c]))
#     }
#   }
#   df = data.frame(meanx=meanxa,minx=minxa,maxx=maxxa,
#                   cluster=as.factor(cla), k =ka)
#   ggplot(df) +
#     geom_point(aes(x=factor(k), y=meanx,col = cluster)) +
#     geom_errorbar(aes(x=factor(k), ymin=minx, ymax=maxx, col = cluster)) +
#     labs(x="particles", y="mean and range of data within clusters")
# }

ggrange <- function(output_wvi, Y) {
  requireNamespace("ggplot2", quietly = TRUE)

  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  K <- length(weights)
  ka <- c()
  cla <- c()
  meanxa <- c()
  minxa <- c()
  maxxa <- c()
  for (k in 1:K) {
    configk <- particles[k, ]
    nc <- max(configk)
    for (c in 1:nc) {
      ka <- c(ka, k)
      cla <- c(cla, c)
      meanxa <- c(meanxa, mean(Y[configk == c]))
      minxa <- c(minxa, min(Y[configk == c]))
      maxxa <- c(maxxa, max(Y[configk == c]))
    }
  }
  df <- data.frame(
    meanx = meanxa, minx = minxa, maxx = maxxa,
    cluster = as.factor(cla), k = ka
  )
  ggplot(df) +
    geom_point(aes(x = factor(k), y = meanx, col = cluster)) +
    geom_errorbar(aes(x = factor(k), ymin = minx, ymax = maxx, col = cluster)) +
    labs(x = "particles", y = "mean and range of data within clusters")
}

ggrange_hist <- function(output_wvi, Y) {
  requireNamespace("ggplot2", quietly = TRUE)

  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  K <- length(weights)
  ka <- c()
  cla <- c()
  meanxa <- c()
  minxa <- c()
  maxxa <- c()
  for (k in 1:K) {
    configk <- particles[k, ]
    nc <- max(configk)
    for (c in 1:nc) {
      ka <- c(ka, k)
      cla <- c(cla, c)
      meanxa <- c(meanxa, mean(Y[configk == c]))
      minxa <- c(minxa, min(Y[configk == c]))
      maxxa <- c(maxxa, max(Y[configk == c]))
    }
  }
  df <- data.frame(
    meanx = meanxa, minx = minxa, maxx = maxxa,
    cluster = as.factor(cla), k = ka
  )
  p1 <- ggplot(df) +
    geom_point(aes(x = factor(k), y = meanx, col = cluster)) +
    geom_errorbar(aes(x = factor(k), ymin = minx, ymax = maxx, col = cluster)) +
    labs(x = "particles", y = "Y") +
    theme_bw() +
    theme(legend.position = "left")

  p2 <- ggplot(data.frame(Y = Y)) +
    geom_histogram(aes(Y)) +
    coord_flip() +
    xlab("") +
    theme_bw()

  gridExtra::grid.arrange(p1, p2,
    nrow = 1,
    layout_matrix = matrix(c(1, 1, 1, 2), nrow = 1),
    top = "Range and mean of data within clusters, with histogram of data"
  )
}

ggscatter <- function(cls, x.true) {
  requireNamespace("ggplot2", quietly = TRUE)

  df_true <- data.frame(x = as.numeric(x.true), cl = as.numeric(cls), jit = 0.2 * rnorm(nrow(x.true)))
  ggplot(df_true) +
    geom_point(aes(x = x, y = cl + jit, color = as.factor(cl)))
}

ggscatter_grid <- function(output_wvi, Y) {
  requireNamespace("ggplot2", quietly = TRUE)

  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  y <- 0.2 * rnorm(nrow(Y))
  df <- data.frame(
    x = NULL,
    cl = NULL,
    jit = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:length(weights)) {
    df_tmp <- data.frame(
      x = as.numeric(Y),
      cl = particles[i, ],
      jit = y,
      ID = i,
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }

  ggplot(df) +
    geom_point(aes(x = x, y = cl + jit, color = as.factor(cl))) +
    facet_wrap(vars(ID))
}

ggsummary <- function(output_wvi, title = NULL, legend = "right") {
  requireNamespace("ggplot2", quietly = TRUE)

  ind <- sort(output_wvi$part.weights, index.return = TRUE, decreasing = T)
  K <- length(output_wvi$part.weights)
  particles <- as.factor(c(1:K))
  clusters <- apply(output_wvi$particles[ind$ix, ], 1, \(x)length(unique(x)))
  weights <- output_wvi$part.weights[ind$ix]
  dat <- data.frame(
    particles = particles, Weight = weights,
    Number_of_clusters = clusters
  )
  dat <- dat %>% pivot_longer(cols = c("Weight", "Number_of_clusters"), names_to = "type", values_to = "value")
  dat <- dat %>%
    mutate(type = ifelse(type == "Number_of_clusters", "Number of clusters", type)) %>%
    mutate(type = factor(type, levels = c("Weight", "Number of clusters")))
  p <- ggplot(dat) +
    geom_bar(
      aes(
        x = particles,
        y = value,
        fill = particles
      ),
      stat = "identity", width = 0.7
    ) +
    labs(x = "particles", y = "") +
    xlab("") +
    scale_fill_brewer(palette = "Accent") +
    facet_wrap(~type, scales = "free") +
    theme_bw() +
    theme(axis.title.x = element_blank())
  if (legend == "bottom") {
    p <- p + theme(legend.position = "bottom")
  }
  if (!is.null(title)) {
    return(p + ggtitle(title)) #+ ggtitle("WASABI particles: Weight and number of clusters")
  } else {
    return(p)
  }
}

plotGraph <- function(output_wvi) {
  requireNamespace("ggplot2", quietly = TRUE)

  viall <- VI_Rcpp(
    output_wvi$particles - 1, output_wvi$particles - 1,
    apply(output_wvi$particles, 1, max), apply(output_wvi$particles, 1, max)
  )
  viall[viall < 0] <- 0 # this is to correct numerical errors
  g <- make_full_graph(n = nrow(output_wvi$particles))
  V(g)$weight <- output_wvi$part.weights
  E(g)$weight <- viall[upper.tri(viall)]
  plot(g,
    layout = layout_with_kk(g),
    edge.label = round(E(g)$weight, 2),
    # vertex.label = paste("w =",V(g)$weight),
    # vertex.label.dist = 3.5,
    # vertex.label.color = "black",
    vertex.size = 80 * V(g)$weight,
    main = "Particle's VI distance"
  )
}

ggpsm <- function(output_wvi) {
  requireNamespace("ggplot2", quietly = TRUE)

  ind <- sort(output_wvi$part.weights, index.return = TRUE, decreasing = T)
  K <- length(output_wvi$part.weights)
  for (k in 1:K) {
    superheat(output_wvi$part.psm[[ind$ix[k]]],
      pretty.order.rows = TRUE,
      pretty.order.cols = TRUE,
      heat.pal = c("white", "yellow", "red"),
      heat.pal.values = c(0, .5, 1),
      title = paste("Particle", k)
    )
  }
}

ggscatter2 <- function(cls, x.true) {
  requireNamespace("ggplot2", quietly = TRUE)

  df_true <- data.frame(x = as.numeric(x.true), cl = as.numeric(cls), jit = runif(nrow(x.true), -0.1, 0.1))
  ggplot(df_true) +
    geom_point(aes(x = x, y = jit, color = as.factor(cl))) +
    theme_bw()
}

ggscatter_grid2 <- function(output_wvi, Y, seed = 123) {
  requireNamespace("ggplot2", quietly = TRUE)

  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  set.seed(seed)
  y <- 0.2 * rnorm(nrow(Y))
  df <- data.frame(
    x = NULL,
    Cluster = NULL,
    jit = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:nrow(particles)) {
    df_tmp <- data.frame(
      x = as.numeric(Y),
      Cluster = particles[i, ],
      jit = y,
      ID = paste0("Particle ", i, ", w = ", weights[i]),
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }
  df$Cluster <- as.factor(df$Cluster)

  ggplot(df) +
    geom_point(aes(x = x, y = jit, color = Cluster)) +
    ylab("Jitter") +
    ylim(c(-1, 1)) +
    facet_wrap(vars(ID)) +
    theme_bw()
}

ggscatter2d <- function(cls, Y) {
  requireNamespace("ggplot2", quietly = TRUE)

  df_true <- data.frame(x1 = Y[, 1], x2 = Y[, 2], Cluster = as.factor(cls))
  ggplot(df_true) +
    geom_point(aes(x = x1, y = x2, color = Cluster)) +
    theme_bw()
}

ggscatter_grid2d <- function(output_wvi, Y, title = NULL, legend = "right") {
  requireNamespace("ggplot2", quietly = TRUE)

  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  df <- data.frame(
    x1 = NULL,
    x2 = NULL,
    Cluster = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:nrow(particles)) {
    df_tmp <- data.frame(
      x1 = Y[, 1],
      x2 = Y[, 2],
      Cluster = particles[i, ],
      ID = paste0("Particle ", i, ", w = ", weights[i]),
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }
  df$Cluster <- as.factor(df$Cluster)

  p <- ggplot(df) +
    geom_point(aes(x = x1, y = x2, color = Cluster, shape = Cluster)) +
    ylab(expression("x"[2])) +
    xlab(expression("x"[1])) +
    facet_wrap(vars(ID)) +
    theme_bw()
  if (legend == "bottom") {
    p <- p + theme(legend.position = "bottom")
  }
  if (!is.null(title)) {
    if (title != FALSE) {
      return(p + ggtitle(title))
    } else {
      return(p)
    }
  } else {
    return(p)
  }
}

ggscatter_pca <- function(output_wvi, Y, title = NULL) {
  requireNamespace("ggplot2", quietly = TRUE)

  out_pc <- prcomp(Y)
  Ypc <- out_pc$x
  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  df <- data.frame(
    x1 = NULL,
    x2 = NULL,
    Cluster = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:nrow(particles)) {
    df_tmp <- data.frame(
      x1 = Ypc[, 1],
      x2 = Ypc[, 2],
      Cluster = particles[i, ],
      ID = paste0("Particle ", i, ", w = ", weights[i]),
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }
  df$Cluster <- as.factor(df$Cluster)

  p <- ggplot(df) +
    geom_point(aes(x = x1, y = x2, color = Cluster)) + # removed shape because it allows for max 6 clusters
    ylab(expression("PC"[2])) +
    xlab(expression("PC"[1])) +
    facet_wrap(vars(ID)) +
    theme_bw()
  if (!is.null(title)) {
    return(p + ggtitle(title))
  } else {
    return(p)
  }
}

ggscatter_umap <- function(output_wvi, Y, title = NULL) {
  requireNamespace("ggplot2", quietly = TRUE)

  out_umap <- umap(Y)
  Ypc <- out_umap$layout
  particles <- output_wvi$particles
  weights <- output_wvi$part.weights
  df <- data.frame(
    x1 = NULL,
    x2 = NULL,
    Cluster = NULL,
    ID = NULL,
    w = NULL
  )
  for (i in 1:nrow(particles)) {
    df_tmp <- data.frame(
      x1 = Ypc[, 1],
      x2 = Ypc[, 2],
      Cluster = particles[i, ],
      ID = paste0("Particle ", i, ", w = ", weights[i]),
      w = weights[i]
    )
    df <- rbind(df, df_tmp)
  }
  df$Cluster <- as.factor(df$Cluster)

  p <- ggplot(df) +
    geom_point(aes(x = x1, y = x2, color = Cluster)) + # removed shape because it allows for max 6 clusters
    ylab(expression("UMAP"[2])) +
    xlab(expression("UMAP"[1])) +
    facet_wrap(vars(ID)) +
    theme_bw()
  if (!is.null(title)) {
    return(p + ggtitle(title))
  } else {
    return(p)
  }
}
