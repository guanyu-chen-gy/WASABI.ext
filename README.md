
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WASABI

<!-- badges: start -->
<!-- badges: end -->

The goal of WASABI is to summarize the posterior distribution over the
space of partitions obtained after running a Bayesian mixture model. It
summarizes the posterior using several point estimates. These are
identified by approximating the empirical posterior distribution
(comprised by the posterior MCMC draws over the space of partitions)
with a discrete distribution supported on these $L$ point estimates
(called ‘particles’) in a Wasserstein sense (the Wasserstein
distribution is based on the Variation of Information (VI) distance).

## Installation

You can install the development version of WASABI from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cecilia-balocchi/WASABI")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(WASABI)

set.seed(123)
mu <- c(-1.1, 1.1)
prop <- c(0.5, 0.5)
n <- 300
components <- sample(1:2, size = n, replace = TRUE, prob = prop)
y <- rnorm(n, mean = mu[components], sd = 1)
# fit a Bayesian mixture model
# (can use your own code, or a package such as BNPmix)
est_model <- BNPmix::PYdensity(y = y,
                               mcmc = list(niter = 6000,
                                           nburn = 5000,
                                           model = "LS"),
                               output = list(out_type = "FULL", 
                                             out_param = TRUE))
#> Completed:   600/6000 - in 0.09202 sec
#> Completed:   1200/6000 - in 0.187285 sec
#> Completed:   1800/6000 - in 0.283321 sec
#> Completed:   2400/6000 - in 0.377232 sec
#> Completed:   3000/6000 - in 0.453758 sec
#> Completed:   3600/6000 - in 0.562026 sec
#> Completed:   4200/6000 - in 0.668807 sec
#> Completed:   4800/6000 - in 0.777227 sec
#> Completed:   5400/6000 - in 0.885569 sec
#> Completed:   6000/6000 - in 0.98883 sec
#> 
#> Estimation done in 0.988886 seconds
cls.draw = est_model$clust
psm=mcclust::comp.psm(cls.draw+1)
out_WASABI <- WASABI(cls.draw, psm = psm, L = 2,
                     method.init = "topvi", method = "salso")
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
