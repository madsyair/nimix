#' @include neonorm-utils.R
NULL

## dist-msnburr-core.R ----------------------------------------------------------
## MSNBurr d/p/q/r -- reference implementation contributed by the package
## maintainer (A. S. Choir), adopted verbatim for numerical-stability
## guarantees (asymptotic log-omega branch, log1pexp thresholds, two-branch
## quantile inversion). Iriawan (2000); Choir (2020).

#' MSNBurr Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the MSNBurr (Modified to be Stable as Normal from Burr) distribution.
#'
#' @name msnburr-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param alpha Numeric. Shape parameter, must be positive (default = 1).
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dmsnburr} gives the density, \code{pmsnburr} gives the distribution
#' function, \code{qmsnburr} gives the quantile function, and
#' \code{rmsnburr} generates random deviates.
#'
#' @details
#' The MSNBurr distribution with parameters \eqn{\mu}, \eqn{\sigma},
#' and \eqn{\alpha} has probability density function:
#'
#' \deqn{f(x|\mu,\sigma,\alpha) = \frac{\omega}{\sigma}
#'   \exp\left(-\omega\frac{x-\mu}{\sigma}\right)
#'   \left(1 + \frac{1}{\alpha}\exp\left(-\omega\frac{x-\mu}{\sigma}\right)\right)^{-(\alpha+1)}}
#'
#' where \eqn{-\infty < x < \infty}, \eqn{-\infty < \mu < \infty},
#' \eqn{\sigma > 0}, \eqn{\alpha > 0}, and \eqn{\omega} is the
#' normalizing constant:
#'
#' \deqn{\omega = \frac{1}{\sqrt{2\pi}}
#'   \left(1 + \frac{1}{\alpha}\right)^{\alpha+1}}
#'
#' The MSNBurr is a special case of the GMSNBurr distribution with
#' \eqn{\theta = 1}. It is left-skewed (negative skewness).
#'
#' @references
#' Iriawan, N. (2000). Computationally Intensive Approaches to Inference
#' in Neo-Normal Linear Models. PhD Thesis, Curtin University of Technology.
#'
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#'
#' @seealso \code{dmsnburr2a}
#'
#' @examples
#' dmsnburr(0, mu = 0, sigma = 1, alpha = 1)
#'
#' x <- seq(-10, 3, by = 0.1)
#' plot(x, dmsnburr(x, alpha = 0.1), type = "l",
#'      ylab = "Density", main = "MSNBurr Densities")
#' lines(x, dmsnburr(x, alpha = 0.5), col = "red")
#' lines(x, dmsnburr(x, alpha = 1), col = "blue")
#' lines(x, dmsnburr(x, alpha = 5), col = "forestgreen")
#'
#' pmsnburr(c(-2, 0, 2), mu = 0, sigma = 1, alpha = 1)
#' qmsnburr(c(0.025, 0.5, 0.975), alpha = 1)
#'
#' set.seed(123)
#' r <- rmsnburr(1000)
#' hist(r, breaks = 30, freq = FALSE, main = "MSNBurr Random Samples")
#' curve(dmsnburr(x), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @export
dmsnburr <- function(x, mu = 0, sigma = 1, alpha = 1, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    lomega <- .log_omega_msnburr(alpha)
    omega <- exp(lomega)
    zo <- -omega * ((x - mu) / sigma)
    zoa <- zo - log(alpha)

    lp <- lomega - log(sigma) + zo - (alpha + 1) * log1pexp(zoa)
    lp[is.infinite(zo) & zo > 0] <- -Inf

    if (log) lp else exp(lp)
}

#' @rdname msnburr-distribution
#' @export
pmsnburr <- function(q, mu = 0, sigma = 1, alpha = 1,
                     lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    lomega <- .log_omega_msnburr(alpha)
    omega <- exp(lomega)
    zoa <- -omega * ((q - mu) / sigma) - log(alpha)

    qs <- -alpha * log1pexp(zoa)

    if (lower.tail) {
        if (log.p) qs else exp(qs)
    } else {
        if (log.p) log1mexp(qs) else -expm1(qs)
    }
}

#' @rdname msnburr-distribution
#' @export
qmsnburr <- function(p, mu = 0, sigma = 1, alpha = 1,
                     lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    pp <- if (log.p) exp(p) else p
    if (!lower.tail) pp <- 1 - pp

    if (any(pp < 0 | pp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }

    lomega <- .log_omega_msnburr(alpha)
    omega <- exp(lomega)

    result <- numeric(length(pp))
    result[pp == 0] <- -Inf
    result[pp == 1] <- Inf
    result[is.na(pp)] <- NA

    interior <- !is.na(pp) & pp > 0 & pp < 1
    if (any(interior)) {
        p_int <- pp[interior]
        use_log1p <- p_int > 0.5
        a_inv <- -1 / alpha

        log_term <- numeric(length(p_int))
        log_term[use_log1p] <- log(expm1(a_inv * log(p_int[use_log1p])))
        log_term[!use_log1p] <- log(exp(a_inv * log(p_int[!use_log1p])) - 1)

        result[interior] <- mu - (sigma / omega) *
            (log(alpha) + log_term)
    }

    result
}

#' @rdname msnburr-distribution
#' @export
rmsnburr <- function(n, mu = 0, sigma = 1, alpha = 1) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, n))

    n <- ceiling(n)
    u <- stats::runif(n)
    qmsnburr(u, mu = mu, sigma = sigma, alpha = alpha)
}

#' MSNBurr-IIa Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the MSNBurr-IIa distribution.
#'
#' @name msnburr2a-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param alpha Numeric. Shape parameter, must be positive (default = 1).
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dmsnburr2a} gives the density, \code{pmsnburr2a} gives the
#' distribution function, \code{qmsnburr2a} gives the quantile function,
#' and \code{rmsnburr2a} generates random deviates.
#'
#' @details
#' The MSNBurr-IIa distribution with parameters \eqn{\mu}, \eqn{\sigma},
#' and \eqn{\alpha} has probability density function:
#'
#' \deqn{f(x|\mu,\sigma,\alpha) = \frac{\omega}{\sigma}
#'   \exp\left(\omega\frac{x-\mu}{\sigma}\right)
#'   \left(1 + \frac{1}{\alpha}\exp\left(\omega\frac{x-\mu}{\sigma}\right)\right)^{-(\alpha+1)}}
#'
#' where \eqn{\omega} is the same normalizing constant as MSNBurr.
#' The MSNBurr-IIa is a special case of the GMSNBurr distribution with
#' \eqn{\alpha_{\text{gmsnburr}} = 1} and \eqn{\beta_{\text{gmsnburr}} = \alpha}.
#' It is right-skewed (positive skewness).
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#'
#' @seealso \code{dmsnburr}
#'
#' @examples
#' dmsnburr2a(0, mu = 0, sigma = 1, alpha = 2)
#'
#' x <- seq(-3, 10, by = 0.1)
#' plot(x, dmsnburr2a(x, alpha = 0.1), type = "l",
#'      ylab = "Density", main = "MSNBurr-IIa Densities")
#' lines(x, dmsnburr2a(x, alpha = 0.5), col = "red")
#' lines(x, dmsnburr2a(x, alpha = 1), col = "blue")
#' lines(x, dmsnburr2a(x, alpha = 5), col = "forestgreen")
#'
#' pmsnburr2a(c(-2, 0, 2), mu = 0, sigma = 1, alpha = 1)
#' qmsnburr2a(c(0.025, 0.5, 0.975), alpha = 1)
#'
#' set.seed(123)
#' r <- rmsnburr2a(1000)
#' hist(r, breaks = 30, freq = FALSE, main = "MSNBurr-IIa Random Samples")
#' curve(dmsnburr2a(x), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @export
dmsnburr2a <- function(x, mu = 0, sigma = 1, alpha = 1, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    lomega <- .log_omega_msnburr(alpha)
    omega <- exp(lomega)
    zt <- omega * ((x - mu) / sigma)
    zta <- -zt + log(alpha)

    lp <- lomega - log(sigma) + (alpha + 1) * log(alpha) -
        alpha * zt - (alpha + 1) * log1pexp(zta)
    lp[is.infinite(zt) & zt < 0] <- -Inf
    lp[is.infinite(zt) & zt > 0] <- -Inf

    if (log) lp else exp(lp)
}

#' @rdname msnburr2a-distribution
#' @export
pmsnburr2a <- function(q, mu = 0, sigma = 1, alpha = 1,
                       lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    lomega <- .log_omega_msnburr(alpha)
    omega <- exp(lomega)
    zoa <- omega * ((q - mu) / sigma) - log(alpha)

    qs <- exp(-alpha * log1pexp(zoa))

    if (lower.tail) {
        if (log.p) log1p(-qs) else 1 - qs
    } else {
        if (log.p) log(qs) else qs
    }
}

#' @rdname msnburr2a-distribution
#' @export
qmsnburr2a <- function(p, mu = 0, sigma = 1, alpha = 1,
                       lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    pp <- if (log.p) exp(p) else p
    if (!lower.tail) pp <- 1 - pp

    if (any(pp < 0 | pp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }

    lomega <- .log_omega_msnburr(alpha)
    omega <- exp(lomega)

    result <- numeric(length(pp))
    result[pp == 0] <- -Inf
    result[pp == 1] <- Inf
    result[is.na(pp)] <- NA

    interior <- !is.na(pp) & pp > 0 & pp < 1
    if (any(interior)) {
        p_int <- pp[interior]
        use_log1p <- p_int < 0.5
        a_inv <- -1 / alpha

        log_term <- numeric(length(p_int))
        log_term[use_log1p] <- log(exp(a_inv * log1p(-p_int[use_log1p])) - 1)
        log_term[!use_log1p] <- log(expm1(a_inv * log1p(-p_int[!use_log1p])))

        result[interior] <- mu + (sigma / omega) *
            (log(alpha) + log_term)
    }

    result
}

#' @rdname msnburr2a-distribution
#' @export
rmsnburr2a <- function(n, mu = 0, sigma = 1, alpha = 1) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha),
        c("sigma", "alpha"))
    if (length(msgs) > 0) return(rep(NaN, n))

    n <- ceiling(n)
    u <- stats::runif(n)
    qmsnburr2a(u, mu = mu, sigma = sigma, alpha = alpha)
}
