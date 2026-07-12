# Tests for the vectorised row-label helpers (v0.9.0). Each helper must
# reproduce its per-row reference implementation exactly -- vectorisation is a
# performance change only, never a results change.

test_that(".rowDistinct reproduces per-row length(unique(.)) exactly", {
  set.seed(1)
  for (L in c(1L, 3L, 9L)) {
    a <- matrix(sample(seq_len(L), 200L * 17L, replace = TRUE), 200L, 17L)
    ref <- as.integer(apply(a, 1L, function(r) length(unique(r))))
    expect_identical(nimix:::.rowDistinct(a, L), ref)
  }
  # single row / single column edges
  expect_identical(nimix:::.rowDistinct(matrix(2L, 1L, 5L), 4L), 1L)
  expect_identical(nimix:::.rowDistinct(matrix(1:4, 4L, 1L), 4L), rep(1L, 4L))
})

test_that(".rowPresence marks exactly the occupied labels per row", {
  set.seed(2)
  a <- matrix(sample(1:5, 60L * 8L, replace = TRUE), 60L, 8L)
  P <- nimix:::.rowPresence(a, 5L)
  for (t in c(1L, 30L, 60L))
    expect_identical(which(P[t, ]), sort(unique(a[t, ])))
})

test_that("vectorised relabel recode matches match(sort(unique)) semantics", {
  set.seed(3)
  m <- 40L; n <- 25L; Lmax <- 6L
  alloc <- matrix(sample(seq_len(Lmax), m * n, replace = TRUE), m, n)
  # reference: the pre-vectorisation per-row loop
  zRef <- matrix(0L, m, n); occRef <- vector("list", m)
  for (t in seq_len(m)) {
    occ <- sort(unique(alloc[t, ])); occRef[[t]] <- occ
    zRef[t, ] <- match(alloc[t, ], occ)
  }
  # vectorised path (as in relabel)
  pres <- nimix:::.rowPresence(alloc, Lmax)
  rk <- matrix(0L, m, Lmax); acc <- integer(m)
  for (l in seq_len(Lmax)) { acc <- acc + pres[, l]; rk[, l] <- acc }
  z <- matrix(rk[cbind(rep(seq_len(m), times = n), as.vector(alloc))], m, n)
  storage.mode(z) <- "integer"
  occ <- lapply(seq_len(m), function(t) which(pres[t, ]))
  expect_identical(z, zRef)
  expect_identical(occ, occRef)
})

test_that(".nodeToArray fast path equals the generic indexed path", {
  m <- 30L; n <- 7L
  samples <- matrix(rnorm(m * n), m, n,
                    dimnames = list(NULL, paste0("xi[", seq_len(n), "]")))
  # scramble column order to exercise the index mapping
  samples <- samples[, sample(n)]
  arr <- nimix:::.nodeToArray(samples, "xi", n)
  for (i in seq_len(n))
    expect_identical(arr[, i], unname(samples[, paste0("xi[", i, "]")]))
})
