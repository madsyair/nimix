# FS skew multivariate Normal with the orthogonal factor O estimated via the
# Householder angle theta (m = 2). See dist-skewnormal-mv-o.R for the scope and
# identifiability caveats.

test_that("Householder O and Theta^2 match FS restriction (8)", {
  Oh <- nimix:::.householderO
  for (th in c(-0.3, -0.1, 0.1, 0.3)) {
    O <- Oh(th)
    expect_lt(abs(det(O) + 1), 1e-10)                # |O| = (-1)^(m+1) = -1
    expect_lt(abs(O[1, 1] - cos(2 * th)), 1e-12)
    expect_lt(abs(O[2, 1] + sin(2 * th)), 1e-12)
    expect_gt(O[1, 1], abs(O[2, 1]))                 # restriction (8)
  }
  # the bound is exactly pi/8: equality at |theta| = pi/8
  Ob <- Oh(pi / 8)
  expect_lt(abs(Ob[1, 1] - abs(Ob[2, 1])), 1e-10)
  expect_lt(abs(nimix:::.thetaBound - pi / 8), 1e-12)
})

test_that("dskewmvno integrates to one and is theta-invariant at gamma = 1", {
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  for (th in c(-0.3, 0, 0.3)) {
    gr <- seq(-12, 12, 0.06)
    z <- outer(gr, gr, function(a, b)
      dskewmvno(cbind(a, b), c(0, 0), Sg, c(0.6, 1.8), th))
    expect_lt(abs(sum(z) * 0.06^2 - 1), 5e-3)
  }
  dmvn <- function(x, mu, S) {
    k <- length(x)
    as.numeric(-0.5 * (k * log(2 * pi) + determinant(S)$modulus +
                         t(x - mu) %*% solve(S) %*% (x - mu)))
  }
  x <- c(0.7, -1.3)
  # gamma = 1: density does not depend on theta (theta is unidentified there)
  for (th in c(-0.35, -0.1, 0, 0.1, 0.35))
    expect_lt(abs(dskewmvno(x, c(0, 0), Sg, c(1, 1), th, log = TRUE) -
                    dmvn(x, c(0, 0), Sg)), 1e-10)
})

test_that("theta = 0 nests skewnormal-mv with gamma_2 reciprocal", {
  # O(0) = diag(1, -1), so it flips eps_2; |O| = -1 means O = I is not in O_2.
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  set.seed(5)
  for (i in 1:10) {
    x <- rnorm(2, 0, 3); g <- exp(rnorm(2, 0, .6))
    expect_lt(abs(dskewmvno(x, c(0, 0), Sg, g, 0, log = TRUE) -
                    dskewmvn(x, c(0, 0), Sg, c(g[1], 1 / g[2]), log = TRUE)),
              1e-10)
  }
})

test_that("rskewmvno is consistent with dskewmvno (profile likelihood)", {
  skip_on_cran()
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2); g <- c(0.45, 2.2); thTrue <- 0.30
  set.seed(31)
  Y <- rskewmvno(3000, c(0, 0), Sg, g, thTrue)
  ths <- seq(-0.38, 0.38, length.out = 21)
  ll <- vapply(ths, function(t) sum(dskewmvno(Y, c(0, 0), Sg, g, t,
                                              log = TRUE)), numeric(1))
  expect_lt(abs(ths[which.max(ll)] - thTrue), 0.05)
})

test_that("compiled dSkewMvNO_k equals the R reference", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cK <- nimble::compileNimble(get("dSkewMvNO_k", envir = globalenv()))
  Sg <- matrix(c(2, .8, .8, 1.2), 2, 2)
  set.seed(3); err <- 0
  for (i in 1:20) {
    x <- rnorm(2, 0, 3); g <- exp(rnorm(2, 0, .6)); th <- runif(1, -pi/8, pi/8)
    err <- max(err, abs(cK(x, c(0, 0), Sg, g, th, log = 1) -
                          dskewmvno(x, c(0, 0), Sg, g, th, log = TRUE)))
  }
  expect_lt(err, 1e-8)
})

test_that("skewnormal-mv-o recovers per-component theta at adequate n", {
  skip_on_cran()
  # theta is a large-sample quantity: at 150 obs/component the mirror mode is
  # within ~2 log-lik units, so this test uses 500 per component.
  set.seed(23)
  Sg <- matrix(c(1, .3, .3, .8), 2, 2)
  Y <- rbind(rskewmvno(500, c(-5, -5), Sg, c(0.4, 2.5),  0.30),
             rskewmvno(500, c( 5,  5), Sg, c(2.5, 0.4), -0.25))
  zt <- rep(1:2, each = 500)
  f <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                          distribution = "skewnormal-mv-o",
                          mcmcControl = list(niter = 3500, nburnin = 1500),
                          seed = 3))
  s <- f@relabeled$summary; o <- order(s$mu_1_mean)
  expect_lt(abs(s$mu_1_mean[o][1] - (-5)), 1.0)
  # 95% intervals cover the simulating angles
  expect_lt(s$theta_lwr[o][1], 0.30); expect_gt(s$theta_upr[o][1], 0.30)
  expect_lt(s$theta_lwr[o][2], -0.25); expect_gt(s$theta_upr[o][2], -0.25)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)
})

test_that("skewnormal-mv-o rejects univariate data; m > 2 routes, not errors", {
  # Pre-C.4 this family refused m != 2. Since the general-m implementation,
  # d > 2 ROUTES to SkewNormalMvOGenSpec instead (full routing + recovery is
  # asserted in test-skew-mv-o-general.R, which runs the fit). What must still
  # hold cheaply here: univariate input is rejected as non-multivariate.
  set.seed(1)
  expect_error(nimixClust(rnorm(30), K = 2, method = "fixedk",
                          distribution = "skewnormal-mv-o"),
               "multivariate")
})
