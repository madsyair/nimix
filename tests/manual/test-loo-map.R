test_that("PSIS-LOO is chain-aware and agrees with WAIC", {
  skip_on_cran()
  skip_if_not_installed("loo")
  set.seed(3)
  y <- c(rnorm(90, -2, .6), rnorm(110, 2, .6))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 800, nburnin = 300, nchains = 2),
                  seed = 1)
  ll <- nimix:::.pointwiseLogLik(f)
  expect_false(is.null(attr(ll, "chain_id")))   # chain ids reconstructed
  lo <- nimixLOO(f)
  expect_s3_class(lo, "psis_loo")
  expect_lt(max(lo$diagnostics$pareto_k), 0.7)  # importance sampling reliable
  w <- nimixWAIC(f)
  # WAIC and PSIS-LOO both estimate elpd; on a well-specified model they agree
  expect_lt(abs(w$elpd_waic - lo$estimates["elpd_loo", "Estimate"]), 2)
})

test_that("plotClusterMap validates its inputs", {
  skip_on_cran()
  set.seed(3)
  y <- c(rnorm(60, -2, .5), rnorm(60, 2, .5))
  f <- nimixClust(y, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 400, nburnin = 150), seed = 1)
  expect_error(plotClusterMap(42, "x.shp"), "expects a FitResult")
  if (requireNamespace("sf", quietly = TRUE)) {
    # geometry count must match observations
    pts <- sf::st_sf(geometry = sf::st_sfc(lapply(1:5, function(i)
      sf::st_point(c(i, 0)))))
    expect_error(plotClusterMap(f, pts), "one-to-one")
  } else {
    expect_error(plotClusterMap(f, "nonexistent.shp"), "sf")
  }
})

test_that("plotClusterMap draws a polygon map with entropy shading", {
  skip_on_cran()
  skip_if_not_installed("sf")
  set.seed(7)
  # 4x5 lattice of unit squares; two spatial blocks of means
  G <- expand.grid(x = 1:5, y = 1:4)
  n <- nrow(G)
  polys <- lapply(seq_len(n), function(i) {
    x <- G$x[i]; yv <- G$y[i]
    sf::st_polygon(list(matrix(c(x-1,yv-1, x,yv-1, x,yv, x-1,yv, x-1,yv-1),
                               ncol = 2, byrow = TRUE)))
  })
  m <- sf::st_sf(id = seq_len(n), geometry = sf::st_sfc(polys))
  yobs <- ifelse(G$x <= 2, rnorm(n, -2, .5), rnorm(n, 2, .5))
  f <- nimixClust(yobs, K = 2, method = "fixedk",
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 500, height = 400)
  out <- plotClusterMap(f, m, uncertainty = TRUE)
  grDevices::dev.off()
  expect_true(file.exists(tmp) && file.info(tmp)$size > 1000)  # plot rendered
  expect_equal(length(out$cluster), n)
  expect_true(all(out$entropy >= 0 & out$entropy <= 1))
  # spatially separated truth should be recovered by the mapped partition
  expect_gt(abs(cor(out$cluster, as.integer(G$x <= 2))), 0.8)
})
