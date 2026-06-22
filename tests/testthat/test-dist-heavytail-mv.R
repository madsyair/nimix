# Multivariate heavy-tail components: marginal equivalence and recovery.

test_that("mv marginal equals analytic multivariate-t", {
  mu <- c(1, -1); Sig <- matrix(c(2, 0.5, 0.5, 1.3), 2); df <- 6
  xt <- c(0.4, -0.3)
  mix <- stats::integrate(function(w) vapply(w, function(wi)
    nimix:::.dmvnorm(xt, mu, Sig / wi) * stats::dgamma(wi, df / 2, df / 2),
    numeric(1)), 0, Inf)$value
  expect_equal(mix, nimix:::.dmvt(xt, mu, Sig, df), tolerance = 1e-6)
})

test_that("mv heavy-tail specs inherit NormalMvSpec and validate df", {
  expect_true(methods::is(StudentTMvSpec(), "NormalMvSpec"))
  expect_true(methods::is(NormalGammaMvSpec(), "NormalMvSpec"))
  Y <- matrix(rnorm(40), 20, 2)
  expect_error(defaultPrior(StudentTMvSpec(), Y, control = list(df = 2)), "df")
})

test_that("mv Student-t and Normal-Gamma recover two clusters", {
  skip_on_cran(); skip_if_not_installed("nimble")
  set.seed(1)
  Y <- rbind(matrix(rnorm(120, -3), ncol = 2), matrix(rnorm(120, 3), ncol = 2))
  for (dist in c("studentt", "normalgamma")) {
    fit <- nimixClust(Y, K_max = 6, distribution = dist, method = "dpm",
                      prior = list(df = 5),
                      mcmcControl = list(niter = 1000, nburnin = 400),
                      seed = 1, verbose = FALSE)
    modalK <- as.integer(names(sort(table(fit@Kposterior),
                                    decreasing = TRUE))[1])
    expect_true(modalK %in% 2:3, info = dist)
  }
})
