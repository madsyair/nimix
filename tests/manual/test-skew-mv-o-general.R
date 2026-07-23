# General-m FS orthogonal factor: Householder construction, the angle box, and
# restriction (8) treated as a canonicalisation rather than a constraint.

test_that("orthogonalFactor is orthogonal with |O| = (-1)^(m+1)", {
  set.seed(7)
  for (m in 2:4) {
    th <- runif(nimix:::.nAngles(m), -0.5, 0.5)
    O <- orthogonalFactor(th, m)
    expect_lt(max(abs(crossprod(O) - diag(m))), 1e-10)
    expect_lt(abs(det(O) - (-1)^(m + 1)), 1e-10)
  }
  # m = 2 reduces to the bivariate helper used by skewnormal-mv-o
  for (t in seq(-1.2, 1.2, 0.2))
    expect_lt(max(abs(orthogonalFactor(t, 2) - nimix:::.householderO(t))), 1e-10)
})

test_that("the FS angle box does NOT imply restriction (8)", {
  # FS state that theta^j in Theta^j puts O in O_m. It does not: the fraction of
  # box draws obeying (8) shrinks fast with m. This is why (8) is applied as a
  # canonicalisation instead of a sampling constraint.
  set.seed(2)
  frac <- function(m, N = 600) {
    box <- nimix:::.angleBox(m)
    mean(vapply(seq_len(N), function(i) {
      th <- runif(nimix:::.nAngles(m), box$lower, box$upper)
      nimix:::.restriction8(orthogonalFactor(th, m))
    }, logical(1)))
  }
  f2 <- frac(2); f3 <- frac(3)
  expect_lt(f2, 0.35); expect_gt(f2, 0.15)   # ~ (pi/4)/pi = 0.25
  expect_lt(f3, 0.20)                        # ~ 0.07
  expect_gt(f2, f3)                          # shrinks with m
})

test_that("canonicaliseO gives a unique representative and preserves density", {
  dgen <- function(X, mu, Sigma, gam, O) {
    m <- length(mu); U <- chol(Sigma); Ui <- backsolve(U, diag(m))
    E <- (sweep(X, 2L, mu) %*% Ui) %*% t(O)
    G <- matrix(gam, nrow(E), m, byrow = TRUE)
    S <- ifelse(E < 0, E * G, E / G)
    rowSums(matrix(log(2) - log(gam + 1 / gam), nrow(E), m, byrow = TRUE) +
              dnorm(S, log = TRUE)) - sum(log(diag(U)))
  }
  set.seed(7)
  for (m in 2:4) {
    SP <- nimix:::.signedPerms(m)
    for (rep in 1:4) {
      th <- runif(nimix:::.nAngles(m), -0.8, 0.8)
      O <- orthogonalFactor(th, m)
      A <- matrix(rnorm(m * m), m, m); Sg <- crossprod(A) + diag(m) * m
      g <- exp(rnorm(m, 0, .6)); mu <- rnorm(m)
      X <- matrix(rnorm(20 * m, 0, 2), 20, m)
      # exactly one signed permutation (|P| = +1) satisfies (8)
      nOK <- sum(vapply(SP, function(P) nimix:::.restriction8(P %*% O),
                        logical(1)))
      expect_equal(nOK, 1L)
      cn <- canonicaliseO(O, g)
      expect_true(cn$canonical)
      expect_true(nimix:::.restriction8(cn$O))
      # the density is invariant under the canonicalisation
      expect_lt(max(abs(dgen(X, mu, Sg, g, O) -
                          dgen(X, mu, Sg, cn$gamma, cn$O))), 1e-9)
    }
  }
})

test_that("m = 2 canonical angles land in (-pi/8, pi/8)", {
  # closes the loop with skewnormal-mv-o, whose prior support is exactly this
  set.seed(3)
  for (i in 1:60) {
    t0 <- runif(1, -pi / 2, pi / 2)
    cn <- canonicaliseO(orthogonalFactor(t0, 2), exp(rnorm(2, 0, .5)))
    th <- 0.5 * atan2(-cn$O[2, 1], cn$O[1, 1])
    expect_lt(abs(th), pi / 8 + 1e-8)
  }
})

test_that("compiled dSkewMvNOG_k equals the R reference for m = 2, 3, 4", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cK <- nimble::compileNimble(get("dSkewMvNOG_k", envir = globalenv()))
  dgen <- function(x, mu, Sigma, gam, O) {
    m <- length(mu); U <- chol(Sigma); Ui <- backsolve(U, diag(m))
    E <- as.numeric(t(O %*% (t(Ui) %*% (x - mu))))
    S <- ifelse(E < 0, E * gam, E / gam)
    sum(log(2) - log(gam + 1 / gam) + dnorm(S, log = TRUE)) - sum(log(diag(U)))
  }
  set.seed(5)
  for (m in 2:4) {
    err <- 0
    for (rep in 1:8) {
      th <- runif(nimix:::.nAngles(m), -0.7, 0.7)
      O <- orthogonalFactor(th, m)
      A <- matrix(rnorm(m * m), m, m); Sg <- crossprod(A) + diag(m) * m
      g <- exp(rnorm(m, 0, .5)); mu <- rnorm(m); x <- rnorm(m, 0, 2)
      err <- max(err, abs(cK(x, mu, Sg, g, th, log = 1) - dgen(x, mu, Sg, g, O)))
    }
    expect_lt(err, 1e-8)
  }
})

test_that("skewnormal-mv-o routes on dimension and clusters in m = 3", {
  skip_on_cran()
  rgen <- function(n, mu, Sg, gam, th) {
    m <- length(mu); O <- orthogonalFactor(th, m); U <- chol(Sg)
    G <- matrix(gam, n, m, byrow = TRUE); W <- abs(matrix(rnorm(n * m), n, m))
    Eps <- ifelse(matrix(runif(n * m), n, m) < G^2 / (1 + G^2), W * G, -W / G)
    sweep(Eps %*% O %*% U, 2L, mu, "+")
  }
  set.seed(40); m <- 3
  Sg <- diag(m); Sg[1, 2] <- Sg[2, 1] <- 0.25
  Y <- rbind(rgen(250, c(-4, -4, -4), Sg, c(0.5, 2.0, 1.5), c(0.3, -0.2, 0.4)),
             rgen(250, c( 4,  4,  4), Sg, c(2.0, 0.5, 0.7), c(-0.25, 0.35, -0.3)))
  zt <- rep(1:2, each = 250)
  f <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                          distribution = "skewnormal-mv-o",
                          mcmcControl = list(niter = 2500, nburnin = 1000),
                          seed = 3))
  # dimension routing picked the general-m spec
  expect_s4_class(f@distSpec, "SkewNormalMvOGenSpec")
  s <- f@relabeled$summary; o <- order(s$mu_1_mean)
  expect_lt(abs(s$mu_1_mean[o][1] - (-4)), 1.0)
  expect_lt(abs(s$mu_1_mean[o][2] - 4), 1.0)
  # every draw got a canonical representative
  expect_equal(f@relabeled$canonicalFraction, 1)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)
})

test_that("canonicaliseO carries nu (permuted, never inverted)", {
  # A sign flip inverts gamma but leaves nu alone, because the Student kernel
  # is symmetric: p(-e | gamma, nu) = p(e | 1/gamma, nu).
  dgenT <- function(X, mu, Sigma, gam, nu, O) {
    m <- length(mu); U <- chol(Sigma); Ui <- backsolve(U, diag(m))
    E <- (sweep(X, 2L, mu) %*% Ui) %*% t(O)
    G <- matrix(gam, nrow(E), m, byrow = TRUE)
    S <- ifelse(E < 0, E * G, E / G)
    lp <- numeric(nrow(E))
    for (j in seq_len(m))
      lp <- lp + log(2) - log(gam[j] + 1 / gam[j]) + dt(S[, j], nu[j], log = TRUE)
    lp - sum(log(diag(U)))
  }
  set.seed(9)
  for (m in 2:4) for (rep in 1:3) {
    th <- runif(nimix:::.nAngles(m), -0.8, 0.8)
    O <- orthogonalFactor(th, m)
    A <- matrix(rnorm(m * m), m, m); Sg <- crossprod(A) + diag(m) * m
    g <- exp(rnorm(m, 0, .6)); nu <- 3 + rexp(m, 0.2); mu <- rnorm(m)
    X <- matrix(rnorm(20 * m, 0, 2), 20, m)
    cn <- canonicaliseO(O, g, nu)
    expect_true(cn$canonical)
    expect_setequal(cn$nu, nu)                     # permuted, not transformed
    expect_lt(max(abs(dgenT(X, mu, Sg, g, nu, O) -
                        dgenT(X, mu, Sg, cn$gamma, cn$nu, cn$O))), 1e-9)
  }
})

test_that("compiled dSkewMvITOG_k equals the R reference for m = 2, 3, 4", {
  skip_on_cran()
  suppressMessages(nimix:::.nimixEnsureMSNBurr())
  cK <- nimble::compileNimble(get("dSkewMvITOG_k", envir = globalenv()))
  dgenT <- function(x, mu, Sigma, gam, nu, O) {
    m <- length(mu); U <- chol(Sigma); Ui <- backsolve(U, diag(m))
    E <- as.numeric(O %*% (t(Ui) %*% (x - mu)))
    S <- ifelse(E < 0, E * gam, E / gam)
    sum(log(2) - log(gam + 1 / gam) + dt(S, nu, log = TRUE)) - sum(log(diag(U)))
  }
  set.seed(5)
  for (m in 2:4) {
    err <- 0
    for (rep in 1:6) {
      th <- runif(nimix:::.nAngles(m), -0.7, 0.7)
      O <- orthogonalFactor(th, m)
      A <- matrix(rnorm(m * m), m, m); Sg <- crossprod(A) + diag(m) * m
      g <- exp(rnorm(m, 0, .5)); nu <- 3 + rexp(m, 0.2)
      mu <- rnorm(m); x <- rnorm(m, 0, 2)
      err <- max(err, abs(cK(x, mu, Sg, g, nu, th, log = 1) -
                            dgenT(x, mu, Sg, g, nu, O)))
    }
    expect_lt(err, 1e-8)
  }
})

test_that("skewistudent-mv-o routes on dimension and clusters in m = 3", {
  skip_on_cran()
  rgenT <- function(n, mu, Sg, gam, nu, th) {
    m <- length(mu); O <- orthogonalFactor(th, m); U <- chol(Sg)
    G <- matrix(gam, n, m, byrow = TRUE)
    W <- abs(matrix(rt(n * m, df = rep(nu, each = n)), n, m))
    Eps <- ifelse(matrix(runif(n * m), n, m) < G^2 / (1 + G^2), W * G, -W / G)
    sweep(Eps %*% O %*% U, 2L, mu, "+")
  }
  set.seed(42); m <- 3
  Sg <- diag(m); Sg[1, 2] <- Sg[2, 1] <- 0.25
  Y <- rbind(rgenT(250, c(-5, -5, -5), Sg, c(0.5, 2.0, 1.5), c(7, 7, 7),
                   c(0.3, -0.2, 0.4)),
             rgenT(250, c( 5,  5,  5), Sg, c(2.0, 0.5, 0.7), c(7, 7, 7),
                   c(-0.25, 0.35, -0.3)))
  zt <- rep(1:2, each = 250)
  f <- relabel(nimixClust(Y, K = 2, method = "fixedk",
                          distribution = "skewistudent-mv-o",
                          mcmcControl = list(niter = 2500, nburnin = 1000),
                          seed = 3))
  expect_s4_class(f@distSpec, "SkewIStudentMvOGenSpec")
  s <- f@relabeled$summary; o <- order(s$mu_1_mean)
  expect_lt(abs(s$mu_1_mean[o][1] - (-5)), 1.2)
  expect_true(all(s$nu_1_mean > 2))
  expect_equal(f@relabeled$canonicalFraction, 1)
  # gamma is reported post-canonicalisation, so compare to canonicalised truth
  ctrue <- canonicaliseO(orthogonalFactor(c(0.3, -0.2, 0.4), 3),
                         c(0.5, 2.0, 1.5))
  ghat <- c(s$gamma_1_mean[o][1], s$gamma_2_mean[o][1], s$gamma_3_mean[o][1])
  expect_lt(max(abs(ghat - ctrue$gamma)), 0.35)
  z <- apply(f@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zt), mean(z == (3L - zt))), 0.9)
})

test_that("MRF engine supports both general-m estimated-O families", {
  skip_on_cran()
  g <- gridAdjacency(8, 8, "rook"); nn <- 64; m <- 3
  zb <- integer(nn)
  for (i in 1:8) for (j in 1:8) zb[(i - 1) * 8 + j] <- if (j <= 4) 1L else 2L
  Sg <- diag(m) * 0.6; Sg[1, 2] <- Sg[2, 1] <- 0.15

  rg <- function(mu, gam, th) {
    O <- orthogonalFactor(th, m); U <- chol(Sg)
    W <- abs(rnorm(m))
    e <- ifelse(runif(m) < gam^2 / (1 + gam^2), W * gam, -W / gam)
    as.numeric(t(U) %*% t(O) %*% e + mu)
  }
  set.seed(50)
  Y <- t(vapply(seq_len(nn), function(i) if (zb[i] == 1L)
    rg(c(-3, -3, -3), c(0.5, 2.0, 1.5), c(.3, -.2, .4)) else
    rg(c( 3,  3,  3), c(2.0, 0.5, 0.7), c(-.25, .35, -.3)), numeric(m)))
  fm <- relabel(nimixClust(Y, K = 2, method = "mrf", spatialWeights = g,
                           distribution = "skewnormal-mv-o",
                           mcmcControl = list(niter = 1200, nburnin = 500),
                           seed = 3))
  expect_s4_class(fm@distSpec, "SkewNormalMvOGenSpec")
  z <- apply(fm@clusterAllocation, 2L,
             function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z == zb), mean(z == (3L - zb))), 0.9)

  # Heavy tails plus K * nAng slice-sampled angles slow the Potts sweep down;
  # 1200 iterations is not enough here, 3000 is (checked directly).
  rgT <- function(mu, gam, nu, th) {
    O <- orthogonalFactor(th, m); U <- chol(Sg)
    W <- abs(rt(m, df = nu))
    e <- ifelse(runif(m) < gam^2 / (1 + gam^2), W * gam, -W / gam)
    as.numeric(t(U) %*% t(O) %*% e + mu)
  }
  set.seed(51)
  Y2 <- t(vapply(seq_len(nn), function(i) if (zb[i] == 1L)
    rgT(c(-4, -4, -4), c(0.5, 2.0, 1.5), c(7, 7, 7), c(.3, -.2, .4)) else
    rgT(c( 4,  4,  4), c(2.0, 0.5, 0.7), c(7, 7, 7), c(-.25, .35, -.3)),
    numeric(m)))
  fm2 <- relabel(nimixClust(Y2, K = 2, method = "mrf", spatialWeights = g,
                            distribution = "skewistudent-mv-o",
                            mcmcControl = list(niter = 3000, nburnin = 1200),
                            seed = 3))
  expect_s4_class(fm2@distSpec, "SkewIStudentMvOGenSpec")
  z2 <- apply(fm2@clusterAllocation, 2L,
              function(v) as.integer(names(which.max(table(v)))))
  expect_gt(max(mean(z2 == zb), mean(z2 == (3L - zb))), 0.9)
})
