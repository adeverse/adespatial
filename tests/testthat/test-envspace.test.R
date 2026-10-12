# Tests of envspace.test() ----
# ****************************************************************************
# ****************************************************************************
# These tests guard against the regression reported in the GitHub issue: before the
# fix, envspace.test() extracted the shared space-environment fraction (SSEF) from
# vegan::varpart() by position, and silently tested the pure spatial fraction once
# vegan (>= 2.6-4) reordered the fractions. The SSEF is now computed from its
# definition; the tests check it against varpart() without relying on any row order.

# Simulated data on a complete regular grid ----
# ****************************************************************************
# A spatially structured environment (two variables) drives a multivariate response,
# so that the SSEF is large and the preliminary tests of envspace.test() pass.
make_grid_data <- function(nx = 15, ny = 10, nsp = 6, seed = 1) {
  set.seed(seed)
  coord <- as.matrix(expand.grid(x = 1:nx, y = 1:ny))
  e1  <- sin(coord[, 1] / 3) + cos(coord[, 2] / 2.5) + rnorm(nrow(coord), sd = 0.3)
  e2  <- coord[, 1] / nx + rnorm(nrow(coord), sd = 0.3)
  env <- cbind(e1 = e1, e2 = e2)
  B   <- matrix(rnorm(2 * nsp), 2, nsp)
  spe <- scale(env) %*% B + matrix(rnorm(nrow(coord) * nsp, sd = 0.8), ncol = nsp)
  lw  <- listw.candidates(coord, nb = "gab", weights = "binary")[[1]]
  MEM <- suppressMessages(mem.select(spe, lw, method = "FWD",
                                     MEM.autocor = "positive")$MEM.select)
  list(coord = coord, env = env, spe = spe, lw = lw, MEM = MEM)
}

d <- make_grid_data()

# Observed statistic = shared fraction of vegan::varpart() ----
# ****************************************************************************
test_that("the observed statistic is the shared fraction of varpart(), not a pure one", {

  set.seed(1)
  res <- envspace.test(d$spe, d$env, d$coord, d$MEM, d$lw, regular = FALSE, nperm = 19)

  # Fractions of the variation partitioning (rows in vegan's current order, whatever it is)
  indfract <- vegan::varpart(d$spe, d$env, d$MEM)$part$indfract

  # Shared fraction identified from its definition, independently of any row order
  r2adj  <- function(X) vegan::RsquareAdj(vegan::rda(d$spe, X))$adj.r.squared
  shared <- r2adj(d$env) + r2adj(d$MEM) - r2adj(cbind(d$env, d$MEM))

  # ... and identified in varpart() as the only fraction with zero degrees of freedom
  shared.vp <- indfract$Adj.R.squared[which(indfract$Df == 0)]

  expect_equal(res$obs, shared, tolerance = 1e-10)
  expect_equal(res$obs, shared.vp, tolerance = 1e-10)

  # The observed value must not be one of the two pure (testable) fractions
  pure <- indfract$Adj.R.squared[which(indfract$Df > 0)]
  expect_false(any(abs(res$obs - pure) < 1e-10))
})

# Torus-translation test: independent of row order and coordinate units ----
# ****************************************************************************
test_that("the torus-translation test does not depend on row order or units", {

  set.seed(10)
  res.ref <- envspace.test(d$spe, d$env, d$coord, d$MEM, d$lw, regular = TRUE, nperm = 49)

  # Same data, rows shuffled consistently in all tables
  set.seed(99)
  o <- sample(nrow(d$coord))
  lw.o <- listw.candidates(d$coord[o, ], nb = "gab", weights = "binary")[[1]]
  set.seed(10)
  res.o <- envspace.test(d$spe[o, ], d$env[o, ], d$coord[o, ], d$MEM[o, ], lw.o,
                         regular = TRUE, nperm = 49)

  # Same data, coordinates in metres starting at 0 (20-m grid)
  coord.m <- (d$coord - 1) * 20
  lw.m <- listw.candidates(coord.m, nb = "gab", weights = "binary")[[1]]
  set.seed(10)
  res.m <- envspace.test(d$spe, d$env, coord.m, d$MEM, lw.m, regular = TRUE, nperm = 49)

  expect_equal(res.o$obs, res.ref$obs)
  expect_equal(sort(res.o$sim), sort(res.ref$sim))
  expect_equal(res.m$sim, res.ref$sim)
})

# Torus-translation test: null environments are torus translations of env ----
# ****************************************************************************
test_that("each torus translation is an isometry of the torus", {

  grid <- adespatial:::.envspace.torus.grid(d$coord)

  # Matrix of distances between sites on the torus (in grid cells, Manhattan metric)
  tdist <- function(a, b, n) pmin(abs(a - b), n - abs(a - b))
  D <- outer(seq_len(nrow(d$coord)), seq_len(nrow(d$coord)), function(i, j)
    tdist(grid$ix[i], grid$ix[j], grid$nx) + tdist(grid$iy[i], grid$iy[j], grid$ny))

  set.seed(2)
  for (k in 1:20) {
    idx <- adespatial:::.envspace.torus.draw(grid)
    expect_equal(sort(idx), seq_len(nrow(d$coord)))   # a permutation of the sites ...
    expect_equal(D[idx, idx], D)                      # ... preserving all torus distances
  }
})

# Torus-translation test: irregular designs are rejected ----
# ****************************************************************************
test_that("the torus-translation test stops for non-regular designs", {

  coord.irr <- d$coord
  coord.irr[1, ] <- coord.irr[1, ] + 0.3
  expect_error(envspace.test(d$spe, d$env, coord.irr, d$MEM, d$lw,
                             regular = TRUE, nperm = 9), "not evenly spaced")
  expect_error(envspace.test(d$spe[-1, ], d$env[-1, ], d$coord[-1, ], d$MEM[-1, ], d$lw,
                             regular = TRUE, nperm = 9), "complete regular grid")
})
