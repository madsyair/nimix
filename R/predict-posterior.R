# brms-style posterior prediction for the mixtures of regressions.
#
# The distinction matters more here than it does for a single-component model,
# and the reason is worth stating: for a mixture, the EXPECTATION is often
# meaningless. Fit two crossing lines (slopes +1.5 and -1.5, weights 1/2 each)
# and the mixture mean E[Y|x] = sum_k w_k x'beta_k is a flat line through the
# middle -- measured, .fitted came back 0.009 / -0.017 / -0.043 at x = -1/0/1,
# describing data that is nowhere near flat. The expectation is a correct
# summary of a distribution no component has.
#
# So there are three questions, and they have three answers:
#   posteriorLinpred  -- what does each component predict?      (draws x n x K)
#   posteriorEpred    -- what is E[Y|X], averaging components?  (draws x n)
#   posteriorPredictive -- what would a new Y look like?        (draws x n)
# For a mixture, reach for the first far more often than the second.

.pp_prep <- function(object, newdata, draws) {
  if (!isRegressionSpec(object@distSpec))
    stop("posterior prediction on the linear-predictor scale needs a ",
         "regression fit (a formula with predictors).", call. = FALSE)
  prior <- object@prior
  Xnew <- if (is.null(newdata)) prior$X else {
    tt <- stats::delete.response(prior$terms)
    # model.frame() resolves a name absent from `newdata` against the
    # FORMULA'S ENVIRONMENT, where the fitted predictors are still sitting.
    # Left alone it answers a question nobody asked -- one row of newdata in,
    # a silent full-length refit of the original data out. Check first.
    need <- all.vars(tt)
    miss <- setdiff(need, names(newdata))
    if (length(miss))
      stop("`newdata` is missing predictor(s): ", paste(miss, collapse = ", "),
           ". (Without this check R would quietly resolve them against the ",
           "fitted data and return predictions for the wrong rows.)",
           call. = FALSE)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass)
    stats::model.matrix(tt, mf)
  }
  if (ncol(Xnew) != prior$p)
    stop("`newdata` gives ", ncol(Xnew), " predictor column(s); the fit used ",
         prior$p, ".", call. = FALSE)
  betaTr <- object@paramTrace$beta
  if (is.null(betaTr))
    stop("This fit has no coefficient trace to predict from.", call. = FALSE)
  if (length(dim(betaTr)) != 3L)
    stop("posteriorLinpred()/posteriorEpred() currently support a univariate ",
         "response; this fit has a multivariate one.", call. = FALSE)
  m <- dim(betaTr)[1L]
  draws <- as.integer(draws)
  if (length(draws) != 1L || is.na(draws) || draws < 1L)
    stop("draws must be a single integer >= 1.", call. = FALSE)
  use <- if (m > draws)
    as.integer(round(seq(1, m, length.out = draws))) else seq_len(m)
  list(X = Xnew, beta = betaTr, use = use, K = dim(betaTr)[2L],
       p = prior$p, prior = prior, isNew = !is.null(newdata))
}

#' Per-component linear predictors
#'
#' For each posterior draw and each row of data, the linear predictor
#' \eqn{\eta_k = x'\beta_k} of \emph{every} component -- no mixing, no
#' noise. This is usually the object you want from a mixture of regressions:
#' it is what each regime or cluster actually predicts, which
#' \code{\link{posteriorEpred}} averages away.
#'
#' Like \code{brms::posterior_linpred}, this returns the \emph{linear
#' predictor} -- the scale on which the coefficients are linear -- not the
#' response mean. For a Poisson (log link) or Binomial (logit link) fit that
#' is \eqn{\log \mu} or the log-odds, not \eqn{\mu}. Set
#' \code{transform = TRUE} to apply the inverse link and get each component's
#' response mean instead (the counterpart of \code{brms}'s \code{transform}
#' argument); for a Normal fit the two coincide.
#'
#' @param object A \code{\linkS4class{FitResult}} from \code{\link{nimixReg}}.
#' @param newdata Optional data frame; defaults to the fitted data.
#' @param transform If \code{TRUE}, apply the inverse-link so each slice is a
#'   component response mean rather than a linear predictor. Default
#'   \code{FALSE} (the linear-predictor scale, as the name says).
#' @param draws Maximum posterior draws to use (thinned evenly).
#' @return A \code{draws} x \code{n} x \code{K} array, with the component
#'   index last.
#' @examples
#' \dontrun{
#' lp <- posteriorLinpred(fit, newdata = data.frame(x = c(-1, 0, 1)))
#' apply(lp, c(2, 3), mean)     # each component's fitted line
#' }
#' @seealso \code{\link{posteriorEpred}}, \code{\link{posteriorPredictive}}
#' @export
posteriorLinpred <- function(object, newdata = NULL, transform = FALSE,
                             draws = 500L) {
  s <- .pp_prep(object, newdata, draws)
  n <- nrow(s$X); nu <- length(s$use)
  out <- array(NA_real_, dim = c(nu, n, s$K),
               dimnames = list(NULL, rownames(s$X),
                               paste0("component", seq_len(s$K))))
  for (k in seq_len(s$K)) {
    B <- matrix(s$beta[s$use, k, ], nrow = nu, ncol = s$p)
    eta <- B %*% t(s$X)
    out[, , k] <- if (transform)
      matrix(linkInv(object@distSpec, as.numeric(eta), prior = object@prior),
             nu, n) else eta
  }
  out
}

# Mixing weights per draw. Three regimes of meaning, and they are not
# interchangeable:
#  * fixedk / mrf / dpm, in-sample : the posterior allocation of each point.
#  * fixedk / dpm, new data        : the mixture weights (a new point's
#                                    allocation is unknown).
#  * hmm                           : the regime probability, which is a
#                                    function of TIME. In-sample it comes from
#                                    the decoded path; out of sample it depends
#                                    on how far ahead you are, which is what
#                                    nimixForecast() exists to handle.
.pp_weights <- function(object, s) {
  nu <- length(s$use); K <- s$K; n <- nrow(s$X)
  if (!s$isNew) {
    alloc <- object@clusterAllocation[s$use, , drop = FALSE]
    if (ncol(alloc) != n)
      stop("Allocation trace does not match the data.", call. = FALSE)
    return(list(kind = "alloc", alloc = alloc))
  }
  if (identical(object@engineUsed, "hmm"))
    stop("posteriorEpred() with `newdata` is not defined for method = 'hmm': ",
         "the regime weights depend on WHEN the rows occur, and future rows ",
         "have no decoded regime. Use nimixForecast(h = , newdata = ), which ",
         "projects the regime distribution forward, or posteriorLinpred(), ",
         "which needs no weights at all.", call. = FALSE)
  if (K == 1L)                      # no mixture to weight
    return(list(kind = "weights",
                w = matrix(1, length(s$use), 1L)))
  w <- .nodeToArray(object@mcmcSamples, "weights", K)[s$use, , drop = FALSE]
  list(kind = "weights", w = w)
}

#' Expected value of the posterior predictive
#'
#' \eqn{E[Y \mid X]} for each posterior draw and row of data, on the
#' \strong{response scale}, averaging over components:
#' \eqn{\sum_k w_k \, g^{-1}(x'\beta_k)}, where \eqn{g^{-1}} is the
#' inverse link of the family. For a Normal fit that is just
#' \eqn{\sum_k w_k x'\beta_k}; for a Poisson it is
#' \eqn{\sum_k w_k \exp(x'\beta_k)}, and for a Binomial
#' \eqn{\sum_k w_k \, \mathrm{size} \cdot \mathrm{plogis}(x'\beta_k)}.
#' The link is applied to each component before the average, as in
#' \code{brms::posterior_epred} -- averaging first and transforming after
#' would give a different, wrong answer.
#'
#' @section Read this before using it on a mixture:
#' For a mixture the expectation can describe a distribution that no component
#' has. Two crossing lines with equal weight have a flat mixture mean --
#' measured on such a fit, the expectation came back at 0.009, -0.017 and
#' -0.043 for x = -1, 0, 1, for data whose components have slopes +1.5 and
#' -1.5. Nothing is wrong with the number; it is simply not a summary anyone
#' wants. \code{\link{posteriorLinpred}} is usually the right question, and
#' for a regime-switching fit \code{\link{nimixForecast}}'s \code{$regime} is
#' another.
#'
#' The weights are the posterior allocation probabilities of the fitted rows.
#' With \code{newdata} they become the mixture weights instead, since a new
#' row's component is unknown -- and for \code{method = "hmm"} that is refused
#' outright, because a regime weight is a function of time and a future row
#' has no decoded regime. Project it with \code{\link{nimixForecast}}.
#'
#' @inheritParams posteriorLinpred
#' @return A \code{draws} x \code{n} matrix.
#' @seealso \code{\link{posteriorLinpred}}, \code{\link{posteriorPredictive}}
#' @export
posteriorEpred <- function(object, newdata = NULL, draws = 500L) {
  s <- .pp_prep(object, newdata, draws)
  # E[Y|X] is on the RESPONSE scale, so the inverse link is applied to each
  # component BEFORE the mixture average -- E[Y|X] = sum_k w_k g^{-1}(eta_k),
  # never g^{-1}(sum_k w_k eta_k). For a Poisson (log link) the two differ by
  # more than a constant. This matches brms::posterior_epred.
  lp <- posteriorLinpred(object, newdata, transform = TRUE, draws = draws)
  wt <- .pp_weights(object, s)
  nu <- dim(lp)[1L]; n <- dim(lp)[2L]
  out <- matrix(0, nu, n, dimnames = list(NULL, rownames(s$X)))
  if (wt$kind == "alloc") {
    # lp is 3-D, so lp[g, M] with a matrix M would be read as two index
    # arguments; slice to the n x K plane first, keeping the shape when n = 1.
    for (g in seq_len(nu))
      out[g, ] <- matrix(lp[g, , ], n, s$K)[cbind(seq_len(n), wt$alloc[g, ])]
  } else {
    for (k in seq_len(s$K)) out <- out + wt$w[, k] * lp[, , k]
  }
  out
}

#' Draws from the posterior predictive distribution
#'
#' Like \code{\link{posteriorEpred}} but with the residual noise added: what a
#' new observation would actually look like, not merely its expectation. For a
#' mixture this is the well-behaved one of the three -- it is genuinely
#' bimodal where the components disagree, rather than collapsing to a mean
#' that sits between them.
#'
#' @inheritParams posteriorLinpred
#' @return A \code{draws} x \code{n} matrix.
#' @seealso \code{\link{posteriorLinpred}}, \code{\link{posteriorEpred}},
#'   \code{\link{nimixForecast}} for projecting a regime forward in time.
#' @export
posteriorPredictive <- function(object, newdata = NULL, draws = 500L) {
  s <- .pp_prep(object, newdata, draws)
  # eta, not the response mean: responseRng() applies each family's own link
  # and sampler, so a Poisson predictive draw is a count and a Binomial one a
  # proportion of trials -- not Gaussian noise around a transformed mean.
  lp <- posteriorLinpred(object, newdata, transform = FALSE, draws = draws)
  wt <- .pp_weights(object, s)
  nu <- dim(lp)[1L]; n <- dim(lp)[2L]
  spec <- object@distSpec
  s2Tr <- object@paramTrace$s2       # Normal only; NULL for Poisson/Binomial
  needsS2 <- is(spec, "NormalRegSpec") || is(spec, "StudentTRegSpec") ||
             is(spec, "NormalGammaRegSpec")
  if (needsS2 && is.null(s2Tr))
    stop("This fit has no error-variance trace to draw noise from.",
         call. = FALSE)
  st <- matrix(0L, nu, n)
  if (wt$kind == "alloc") {
    st[] <- wt$alloc
  } else {
    for (g in seq_len(nu))
      st[g, ] <- sample.int(s$K, n, replace = TRUE, prob = wt$w[g, ])
  }
  eta <- matrix(0, nu, n)
  for (g in seq_len(nu))
    eta[g, ] <- matrix(lp[g, , ], n, s$K)[cbind(seq_len(n), st[g, ])]
  s2v <- if (is.null(s2Tr)) NULL else {
    s2 <- s2Tr[s$use, , drop = FALSE]
    matrix(s2[cbind(rep(seq_len(nu), times = n), as.vector(st))], nu, n)
  }
  # Neo-normal families carry extra shape parameters (sigma, alpha, ...) that
  # the response draw needs, one value per observation via its regime/cluster
  # state. Pull each shape trace and index it by (draw, state); families
  # without shapeNames (Normal, GLM, ...) skip this and pass shapeVals = NULL.
  shapeVals <- NULL
  shNames <- tryCatch(.neoShapeNames(spec), error = function(e) NULL)
  if (!is.null(shNames)) {
    shapeVals <- list()
    for (snm in shNames) {
      tr <- .nodeToArray(object@mcmcSamples, paste0(snm, "Tilde"),
                         s$K)[s$use, , drop = FALSE]
      shapeVals[[snm]] <- as.numeric(
        matrix(tr[cbind(rep(seq_len(nu), times = n), as.vector(st))], nu, n))
    }
  }
  draws_out <- responseRng(spec, as.numeric(eta),
                           s2 = if (is.null(s2v)) NULL else as.numeric(s2v),
                           prior = object@prior, shapeVals = shapeVals)
  matrix(draws_out, nu, n, dimnames = list(NULL, rownames(s$X)))
}
