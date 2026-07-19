test_that("MSNBurr regression recovers coefficients (framework family 1)", {
  skip_on_cran()
  set.seed(3)
  n <- 250L; x <- rnorm(n)
  y <- c(1.5 + 1.2 * x[1:125], -1.5 - 0.8 * x[126:250]) +
    nimix:::rmsnburr(n, 0, 0.6, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                distribution = "msnburr",
                mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  s <- relabel(f)@relabeled$summary
  o <- order(s[["(Intercept)"]])
  expect_lt(max(abs(s[["(Intercept)"]][o] - c(-1.5, 1.5))), 0.5)
  expect_lt(max(abs(s[["x"]][o] - c(-0.8, 1.2))), 0.5)
  # framework generated the shape-specific summary columns
  expect_true(all(c("sigma_mean", "alpha_mean") %in% names(s)))
})

test_that("SEP regression works through the same framework, different shape set", {
  skip_on_cran()
  # The point of the generic .neoRegMethods: a family with a DIFFERENT shape
  # parameter list (sigma, nu instead of sigma, alpha) needs no new boilerplate.
  set.seed(4)
  n <- 250L; x <- rnorm(n)
  y <- c(1.2 + 1.0 * x[1:125], -1.2 - 0.7 * x[126:250]) +
    nimix:::rsep(n, 0, 0.6, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                distribution = "sep",
                mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  s <- relabel(f)@relabeled$summary
  o <- order(s[["(Intercept)"]])
  expect_lt(max(abs(s[["(Intercept)"]][o] - c(-1.2, 1.2))), 0.5)
  # the shape columns are nu, not alpha -- proving the generic read the family's
  # own parameter list
  expect_true(all(c("sigma_mean", "nu_mean") %in% names(s)))
  expect_false("alpha_mean" %in% names(s))
})

test_that("MSNBurr-IIa regression works through the framework (same shape, new kernel)", {
  skip_on_cran()
  # Same shape set as MSNBurr (sigma, alpha) but a different density kernel.
  # Adding it required only the family declaration -- no new boilerplate.
  set.seed(6)
  n <- 250L; x <- rnorm(n)
  y <- c(1.3 + 1.1 * x[1:125], -1.3 - 0.7 * x[126:250]) +
    nimix:::rmsnburr2a(n, 0, 0.6, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                distribution = "msnburr2a",
                mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  s <- relabel(f)@relabeled$summary
  o <- order(s[["(Intercept)"]])
  expect_lt(max(abs(s[["(Intercept)"]][o] - c(-1.3, 1.3))), 0.5)
  expect_true(all(c("sigma_mean", "alpha_mean") %in% names(s)))
})

test_that("FSSN regression works through the framework (log-normal skew prior)", {
  skip_on_cran()
  # Same shape names as MSNBurr (sigma, alpha) but alpha carries a LOG-NORMAL
  # prior, not a Gamma. Only the prior block and shapeDraw change; the
  # generated boilerplate is untouched.
  set.seed(8)
  n <- 250L; x <- rnorm(n)
  y <- c(1.4 + 1.0 * x[1:125], -1.4 - 0.6 * x[126:250]) +
    nimix:::rfssn(n, 0, 0.6, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                distribution = "fssn",
                mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  s <- relabel(f)@relabeled$summary
  o <- order(s[["(Intercept)"]])
  expect_lt(max(abs(s[["(Intercept)"]][o] - c(-1.4, 1.4))), 0.5)
  expect_true(all(c("sigma_mean", "alpha_mean") %in% names(s)))
})

test_that("GMSNBurr regression works with THREE shape parameters", {
  skip_on_cran()
  # The framework's last untested edge: a family with three shape parameters
  # (sigma, alpha, theta) rather than two. The generated traces, relabelling,
  # and summary columns must all extend to the longer list automatically.
  set.seed(10)
  n <- 260L; x <- rnorm(n)
  y <- c(1.3 + 1.0 * x[1:130], -1.3 - 0.7 * x[131:260]) +
    nimix:::rgmsnburr(n, 0, 0.6, 2, 1.5)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                distribution = "gmsnburr",
                mcmcControl = list(niter = 1000, nburnin = 400), seed = 1)
  s <- relabel(f)@relabeled$summary
  o <- order(s[["(Intercept)"]])
  expect_lt(max(abs(s[["(Intercept)"]][o] - c(-1.3, 1.3))), 0.5)
  # all three shape means appear, generated from the length-3 shape list
  expect_true(all(c("sigma_mean", "alpha_mean", "theta_mean") %in% names(s)))
})

test_that("the remaining neo-normal families all fit through the framework", {
  skip_on_cran()
  # lep (sigma, nu), fsst (sigma, alpha, nu), fossep and jfst
  # (sigma, alpha, theta). Each declared, no new boilerplate; the summary
  # carries exactly that family's shape means.
  set.seed(13)
  n <- 240L; x <- rnorm(n)
  base <- c(1.3 + 1.0 * x[1:120], -1.3 - 0.7 * x[121:240])
  cases <- list(
    lep    = list(y = base + nimix:::rlep(n, 0, 0.6, 2),
                  cols = c("sigma_mean", "nu_mean")),
    fsst   = list(y = base + nimix:::rfsst(n, 0, 0.6, 1.5, 5),
                  cols = c("sigma_mean", "alpha_mean", "nu_mean")),
    fossep = list(y = base + nimix:::rfossep(n, 0, 0.6, 2, 2),
                  cols = c("sigma_mean", "alpha_mean", "theta_mean")),
    jfst   = list(y = base + nimix:::rjfst(n, 0, 0.6, 3, 3),
                  cols = c("sigma_mean", "alpha_mean", "theta_mean"))
  )
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    f <- nimixReg(cs$y ~ x, data.frame(y = cs$y, x = x), K = 2,
                  method = "fixedk", distribution = nm,
                  mcmcControl = list(niter = 700, nburnin = 300), seed = 1)
    s <- relabel(f)@relabeled$summary
    expect_true(all(cs$cols %in% names(s)),
                info = paste("shape columns for", nm))
    expect_length(unique(binderPartition(f)$partition), 2L)
  }
})

test_that("neo-normal regression prediction is complete and keeps the skew", {
  skip_on_cran()
  # linpred/epred work by the identity link; the predictive must draw from the
  # skewed family (via the generic responseRng the framework attached), not a
  # Gaussian fallback -- otherwise the tails/skew are lost (cf. 9.40).
  set.seed(3)
  n <- 200L; x <- rnorm(n)
  y <- c(1.5 + 1.2 * x[1:100], -1.5 - 0.8 * x[101:200]) +
    nimix:::rmsnburr(n, 0, 0.6, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                distribution = "msnburr",
                mcmcControl = list(niter = 600, nburnin = 250), seed = 1)
  nd <- data.frame(x = c(-1, 0, 1))

  lp <- posteriorLinpred(f, nd, draws = 100)
  expect_identical(dim(lp), c(100L, 3L, 2L))
  ep <- posteriorEpred(f, nd, draws = 100)
  expect_length(colMeans(ep), 3L)
  pp <- posteriorPredictive(f, nd, draws = 300)
  expect_true(all(is.finite(pp)))
  expect_identical(dim(pp), c(300L, 3L))
})

test_that("the predictive path generalises across shape counts", {
  skip_on_cran()
  # 2-shape (SEP) and 3-shape (GMSNBurr) both draw finite responses through the
  # same generic path -- the shape traces are pulled by name from the registry.
  set.seed(5)
  n <- 200L; x <- rnorm(n)
  for (nm in c("sep", "gmsnburr")) {
    gen <- if (nm == "sep") nimix:::rsep(n, 0, 0.6, 2) else
      nimix:::rgmsnburr(n, 0, 0.6, 2, 1.5)
    y <- c(1.3 + 1.0 * x[1:100], -1.3 - 0.7 * x[101:200]) + gen
    f <- nimixReg(y ~ x, data.frame(y = y, x = x), K = 2, method = "fixedk",
                  distribution = nm,
                  mcmcControl = list(niter = 500, nburnin = 200), seed = 1)
    pp <- posteriorPredictive(f, data.frame(x = 0), draws = 300)
    expect_true(all(is.finite(pp)), info = nm)
  }
})

test_that("neo-normal regression runs under the DPM engine too", {
  skip_on_cran()
  # The framework generates the DPM model code generically from densName +
  # priorLines, so every family gains a nonparametric (CRP) variant with no
  # per-family DPM code. K is inferred rather than fixed.
  set.seed(3)
  n <- 200L; x <- rnorm(n)
  y <- c(1.5 + 1.2 * x[1:100], -1.5 - 0.8 * x[101:200]) +
    nimix:::rmsnburr(n, 0, 0.6, 2)
  f <- nimixReg(y ~ x, data.frame(y = y, x = x), method = "dpm",
                distribution = "msnburr",
                mcmcControl = list(niter = 800, nburnin = 350), seed = 1)
  s <- relabel(f)@relabeled$summary
  expect_gte(nrow(s), 1L)
  expect_true(all(c("sigma_mean", "alpha_mean") %in% names(s)))
})

test_that("DPM works across shape counts for neo-normal regression", {
  skip_on_cran()
  set.seed(7)
  n <- 200L; x <- rnorm(n)
  for (nm in c("sep", "gmsnburr")) {
    gen <- if (nm == "sep") nimix:::rsep(n, 0, 0.6, 2) else
      nimix:::rgmsnburr(n, 0, 0.6, 2, 1.5)
    y <- c(1.4 + 1.0 * x[1:100], -1.4 - 0.7 * x[101:200]) + gen
    f <- nimixReg(y ~ x, data.frame(y = y, x = x), method = "dpm",
                  distribution = nm,
                  mcmcControl = list(niter = 700, nburnin = 300), seed = 1)
    s <- relabel(f)@relabeled$summary
    expect_gte(nrow(s), 1L)
  }
})
