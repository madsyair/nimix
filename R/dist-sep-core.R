## dist-sep-core.R -----------------------------------------------------------
## SEP d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. Batch B neo-normal family (Choir 2020).

#' @include neonorm-utils.R
NULL

# Subbotin Exponential Power Distribution
# ------------------------------------------------------------------------------

#' Subbotin Exponential Power (SEP) Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Subbotin Exponential Power (SEP) distribution, also known as the
#' Generalized Normal or Generalized Error Distribution.
#'
#' @name sep-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param nu Numeric. Shape parameter, must be positive (default = 2).
#'   \eqn{\nu = 2} gives the normal distribution; \eqn{\nu = 1}
#'   gives the Laplace (double exponential) distribution.
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dsep} gives the density, \code{psep} gives the distribution
#' function, \code{qsep} gives the quantile function, and
#' \code{rsep} generates random deviates.
#'
#' @details
#' The SEP distribution with parameters \eqn{\mu}, \eqn{\sigma}, and
#' \eqn{\nu} has density:
#'
#' \deqn{f(x|\mu,\sigma,\nu) = \frac{\nu}{2^{1+1/\nu}\,
#'   \Gamma(1/\nu)\,\sigma}
#'   \exp\left(-\frac{1}{2}\left|\frac{x-\mu}{\sigma}\right|^\nu\right)}
#'
#' where \eqn{-\infty < x < \infty}, \eqn{-\infty < \mu < \infty},
#' \eqn{\sigma > 0}, \eqn{\nu > 0}.
#'
#' The SEP is symmetric around \eqn{\mu}. Special cases:
#' \itemize{
#'   \item \eqn{\nu = 2}: Normal distribution \eqn{\mathcal{N}(\mu, \sigma^2)}
#'   \item \eqn{\nu = 1}: Laplace (double exponential) distribution
#'   \item \eqn{\nu \to \infty}: Uniform distribution (limit)
#' }
#'
#' @references
#' Subbotin, M. T. (1923). On the law of frequency of error.
#' \emph{Matematicheskii Sbornik}, 31(2), 296--301.
#'
#' @seealso [FOSSEP()] for the skewed version, [LEP()] for an alternative
#'   symmetric exponential power parameterization.
#'
#' @examples
#' dsep(0, mu = 0, sigma = 1, nu = 2)
#'
#' x <- seq(-5, 5, by = 0.1)
#' plot(x, dsep(x, nu = 2), type = "l",
#'      ylab = "Density", main = "SEP Densities")
#' lines(x, dsep(x, nu = 1), col = "red")
#' lines(x, dsep(x, nu = 4), col = "blue")
#'
#' psep(c(-2, 0, 2), mu = 0, sigma = 1, nu = 3)
#' qsep(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, nu = 3)
#'
#' set.seed(123)
#' r <- rsep(1000, mu = 0, sigma = 1, nu = 2)
#' hist(r, breaks = 30, freq = FALSE, main = "SEP Random Samples")
#' curve(dsep(x, nu = 2), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @importFrom stats pgamma qgamma runif
#' @export
dsep <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    z <- abs(x - mu) / sigma
    lp <- -log(2) - (1 / nu) * log(2) - lgamma(1 + 1 / nu) -
        log(sigma) - 0.5 * (z^nu)
    lp[is.infinite(x)] <- -Inf

    if (log) lp else exp(lp)
}

#' @rdname sep-distribution
#' @importFrom stats pgamma
#' @export
psep <- function(q, mu = 0, sigma = 1, nu = 2,
                 lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    z <- (q - mu) / sigma
    s <- abs(z)^nu / 2
    g <- stats::pgamma(s, shape = 1 / nu, rate = 1)

    p <- rep(NA_real_, length(q))
    left <- q < mu
    if (any(left, na.rm = TRUE)) {
        p[left] <- 0.5 * (1 - g[left])
    }
    if (any(!left, na.rm = TRUE)) {
        p[!left] <- 0.5 * (1 + g[!left])
    }

    if (!lower.tail) p <- 1 - p
    if (log.p) p <- log(p)

    p
}

#' @rdname sep-distribution
#' @importFrom stats qgamma
#' @export
qsep <- function(p, mu = 0, sigma = 1, nu = 2,
                 lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    if (log.p) { tp <- exp(p) } else { tp <- p }
    if (!lower.tail) tp <- 1 - tp

    if (any(tp < 0 | tp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }
    result <- rep(NA_real_, length(tp))
    np <- length(tp)
    mu <- rep_len(mu, np); sigma <- rep_len(sigma, np); nu <- rep_len(nu, np)
    left_idx <- tp < 0.5 & !is.na(tp)
    right_idx <- tp >= 0.5 & !is.na(tp)
    eps <- .Machine$double.eps

    if (any(left_idx)) {
        p_clamped <- pmax(pmin(1 - 2 * tp[left_idx], 1 - eps), eps)
        result[left_idx] <- mu[left_idx] - sigma[left_idx] * 2^(1 / nu[left_idx]) *
            (stats::qgamma(p_clamped, shape = 1 / nu[left_idx], rate = 1))^(1 / nu[left_idx])
    }
    if (any(right_idx)) {
        p_clamped <- pmax(pmin(2 * tp[right_idx] - 1, 1 - eps), eps)
        result[right_idx] <- mu[right_idx] + sigma[right_idx] * 2^(1 / nu[right_idx]) *
            (stats::qgamma(p_clamped, shape = 1 / nu[right_idx], rate = 1))^(1 / nu[right_idx])
    }

    result
}

#' @rdname sep-distribution
#' @export
rsep <- function(n, mu = 0, sigma = 1, nu = 2) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, n))

    n <- ceiling(n)
    u <- stats::runif(n, 1e-12, 1 - 1e-12)
    qsep(u, mu = mu, sigma = sigma, nu = nu)
}
