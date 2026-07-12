## dist-lep-core.R -----------------------------------------------------------
## LEP d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. Batch B neo-normal family (Choir 2020).

#' @include neonorm-utils.R
NULL

# Lunetta Exponential Power Distribution
# ------------------------------------------------------------------------------

#' Lunetta Exponential Power Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Lunetta Exponential Power (LEP) distribution.
#'
#' @name lep-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param nu Numeric. Shape parameter, must be positive (default = 2).
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dlep} gives the density, \code{plep} gives the distribution
#' function, \code{qlep} gives the quantile function, and
#' \code{rlep} generates random deviates.
#'
#' @details
#' The LEP distribution with parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\nu}
#' has density:
#'
#' \deqn{f(x|\mu,\sigma,\nu) = \frac{1}{2\nu^{1/\nu} \Gamma(1+1/\nu) \sigma}
#'   \exp\left(-\frac{|x-\mu|^\nu}{\nu \sigma^\nu}\right)}
#'
#' where \eqn{-\infty < x < \infty}, \eqn{-\infty < \mu < \infty},
#' \eqn{\sigma > 0}, \eqn{\nu > 0}.
#'
#' The LEP is symmetric around \eqn{\mu}. When \eqn{\nu = 2}, it reduces to
#' the normal distribution. When \eqn{\nu = 1}, it becomes the Laplace
#' (double exponential) distribution.
#'
#' @references
#' Lunetta, G. (1963). Di una generalizzazione dello schema della curva
#' normale. \emph{Annali della Facolta di Economia e Commercio di Palermo},
#' 17, 237-244.
#'
#' @examples
#' dlep(0, mu = 0, sigma = 1, nu = 2)
#'
#' x <- seq(-5, 5, by = 0.1)
#' plot(x, dlep(x, nu = 2), type = "l",
#'      ylab = "Density", main = "LEP Densities")
#' lines(x, dlep(x, nu = 1), col = "red")
#' lines(x, dlep(x, nu = 4), col = "blue")
#' legend("topright", c("nu=2 (Normal)", "nu=1 (Laplace)", "nu=4"),
#'        col = c("black", "red", "blue"), lty = 1)
#'
#' plep(c(-2, 0, 2), mu = 0, sigma = 1, nu = 2)
#' qlep(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, nu = 2)
#'
#' set.seed(123)
#' r <- rlep(1000, mu = 0, sigma = 1, nu = 2)
#' hist(r, breaks = 30, freq = FALSE, main = "LEP Random Samples")
#' curve(dlep(x, nu = 2), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @export
dlep <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    z <- abs(x - mu) / sigma
    lp <- -log(2) - (1 / nu) * log(nu) - lgamma(1 + 1 / nu) -
        log(sigma) - (z^nu) / nu
    lp[is.infinite(x)] <- -Inf

    if (log) lp else exp(lp)
}

#' @rdname lep-distribution
#' @export
plep <- function(q, mu = 0, sigma = 1, nu = 2,
                 lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    z <- (q - mu) / sigma
    g <- stats::pgamma(abs(z)^nu / nu, shape = 1 / nu, rate = 1)

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

#' @rdname lep-distribution
#' @export
qlep <- function(p, mu = 0, sigma = 1, nu = 2,
                 lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    if (log.p) {
        tp <- exp(p)
    } else {
        tp <- p
    }
    if (!lower.tail) tp <- 1 - tp

    if (any(tp < 0 | tp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }
    result <- rep(NA_real_, length(tp))
    np <- length(tp)
    mu <- rep_len(mu, np); sigma <- rep_len(sigma, np); nu <- rep_len(nu, np)
    left_idx <- tp < 0.5 & !is.na(tp)
    right_idx <- tp >= 0.5 & !is.na(tp)

    if (any(left_idx)) {
        result[left_idx] <- mu[left_idx] - sigma[left_idx] *
            (nu[left_idx] * stats::qgamma(1 - 2 * tp[left_idx],
                            shape = 1 / nu[left_idx], rate = 1))^(1 / nu[left_idx])
    }
    if (any(right_idx)) {
        result[right_idx] <- mu[right_idx] + sigma[right_idx] *
            (nu[right_idx] * stats::qgamma(2 * tp[right_idx] - 1,
                             shape = 1 / nu[right_idx], rate = 1))^(1 / nu[right_idx])
    }

    result
}

#' @rdname lep-distribution
#' @export
rlep <- function(n, mu = 0, sigma = 1, nu = 2) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, nu = nu),
        c("sigma", "nu"))
    if (length(msgs) > 0) return(rep(NaN, n))

    n <- ceiling(n)
    u <- stats::runif(n, 1e-12, 1 - 1e-12)
    qlep(u, mu = mu, sigma = sigma, nu = nu)
}
