# Time-series views for regime-switching fits.
#
# A density plot of a regime-switching series is the wrong picture twice over:
# it discards the axis the model is ABOUT, and it shows a bimodal smear where
# the story is "the series sat in one regime, then moved". These put time back
# on the x-axis. As with every other plot type, each returns (invisibly) the
# tidy frame it drew, so the numbers are testable without inspecting pixels.

.tsFit <- function() {
  set.seed(42)
  P <- rbind(c(.95, .05), c(.10, .90))
  z <- integer(200); z[1] <- 1L
  for (t in 2:200) z[t] <- sample(1:2, 1L, prob = P[z[t - 1L], ])
  y <- rnorm(200, c(-2, 2)[z], 0.7)
  list(fit = nimixClust(y, K = 2, method = "hmm",
                        mcmcControl = list(niter = 800, nburnin = 300),
                        seed = 1),
       z = z, y = y)
}

test_that("plot(type = 'series') returns the decoded regimes over time", {
  skip_on_cran()
  o <- .tsFit()
  pdf(NULL); on.exit(dev.off())
  s <- plot(o$fit, type = "series")
  expect_identical(names(s), c("time", "y", "regime"))
  expect_identical(nrow(s), 200L)
  expect_identical(s$y, o$y)
  expect_gt(max(mean(s$regime == o$z), mean((3L - s$regime) == o$z)), 0.9)
})

test_that("plot(type = 'regime') returns smoothed probabilities that sum to one", {
  skip_on_cran()
  o <- .tsFit()
  pdf(NULL); on.exit(dev.off())
  r <- plot(o$fit, type = "regime")
  expect_identical(names(r), c("time", "regime1", "regime2"))
  expect_true(all(abs(rowSums(r[, -1L]) - 1) < 1e-8))
  expect_true(all(r[, -1L] >= 0 & r[, -1L] <= 1))
  # where the decode is confident the probabilities should be near 0/1
  expect_gt(mean(apply(r[, -1L], 1L, max) > 0.9), 0.5)
})

test_that("plot(type = 'forecast') returns the fan it draws", {
  skip_on_cran()
  o <- .tsFit()
  pdf(NULL); on.exit(dev.off())
  fc <- plot(o$fit, type = "forecast", h = 8L, draws = 200L)
  expect_identical(nrow(fc), 8L)
  expect_true(all(c("median", "lower", "upper", "regime1") %in% names(fc)))
  expect_true(all(fc$lower <= fc$median & fc$median <= fc$upper))
})

test_that("the time-series views refuse fits with no time axis", {
  skip_on_cran()
  set.seed(1)
  fk <- nimixClust(rnorm(60), K = 2, method = "fixedk",
                   mcmcControl = list(niter = 200, nburnin = 50), seed = 1)
  pdf(NULL); on.exit(dev.off())
  for (ty in c("series", "regime", "forecast"))
    expect_error(plot(fk, type = ty), "regime-switching fits")
})
