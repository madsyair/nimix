## dist-jfst-core.R -----------------------------------------------------------
## JFST d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. Batch B neo-normal family (Choir 2020).

#' @include neonorm-utils.R
NULL

#' Jones-Faddy Skew-t Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Jones-Faddy Skew-t distribution.
#'
#' @name jfst-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param alpha Numeric. Left tail parameter (> 0), the Jones-Faddy \eqn{a}
#'   parameter (default = 3). \eqn{\alpha = \theta} gives symmetry.
#' @param theta Numeric. Right tail parameter (> 0), the Jones-Faddy \eqn{b}
#'   parameter (default = 3). \eqn{\alpha = \theta} gives symmetry.
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{djfst} gives the density, \code{pjfst} gives the distribution
#' function, \code{qjfst} gives the quantile function, and
#' \code{rjfst} generates random deviates.
#'
#' @details
#' The Jones-Faddy Skew-t distribution with parameters \eqn{\mu}, \eqn{\sigma},
#' \eqn{\alpha} (left tail), and \eqn{\theta} (right tail) has density:
#'
#' \deqn{f(y|\mu,\sigma,\alpha,\theta) = \frac{c}{\sigma}
#'   \left[1 + \frac{z}{\sqrt{\alpha+\theta+z^2}}\right]^{\alpha+0.5}
#'   \left[1 - \frac{z}{\sqrt{\alpha+\theta+z^2}}\right]^{\theta+0.5}}
#'
#' where \eqn{z = (y-\mu)/\sigma},
#' \eqn{c = \left[2^{\alpha+\theta-1} \sqrt{\alpha+\theta} \,
#'   B(\alpha,\theta)\right]^{-1}}.
#'
#' \eqn{\alpha < \theta}: left-skewed (heavier right tail).
#' \eqn{\alpha > \theta}: right-skewed (heavier left tail).
#' \eqn{\alpha = \theta}: Student-$t$ distribution (symmetric, \eqn{\nu = 2\alpha}).
#'
#' When \eqn{\alpha = \theta}, JFST reduces to the Student-$t$ distribution
#' with \eqn{\nu = 2\alpha} degrees of freedom. As \eqn{\alpha = \theta \to \infty},
#' JFST approaches the Normal distribution.
#'
#' Moments exist only when certain conditions on \eqn{\alpha} and \eqn{\theta}
#' are met: mean requires \eqn{\alpha,\theta > 0.5},
#' variance requires \eqn{\alpha,\theta > 1},
#' skewness requires \eqn{\alpha,\theta > 1.5},
#' kurtosis requires \eqn{\alpha,\theta > 2}.
#'
#' @references
#' Jones, M.C. and Faddy, M.J. (2003). A skew extension of the t distribution,
#' with applications. Journal of the Royal Statistical Society,
#' Series B, 65, pp 159-174.
#'
#' @examples
#' djfst(0, mu = 0, sigma = 1, alpha = 3, theta = 3)
#'
#' x <- seq(-5, 5, by = 0.1)
#' plot(x, djfst(x, alpha = 3, theta = 3), type = "l",
#'      ylab = "Density", main = "JFST Densities")
#' lines(x, djfst(x, alpha = 5, theta = 1), col = "red")
#' lines(x, djfst(x, alpha = 1, theta = 5), col = "blue")
#'
#' pjfst(c(-2, 0, 2), mu = 0, sigma = 1, alpha = 3, theta = 3)
#' qjfst(c(0.025, 0.5, 0.975), mu = 0, sigma = 1, alpha = 3, theta = 3)
#'
#' set.seed(123)
#' r <- rjfst(1000, mu = 0, sigma = 1, alpha = 3, theta = 3)
#' hist(r, breaks = 30, freq = FALSE, main = "JFST Random Samples")
#' curve(djfst(x, alpha = 3, theta = 3), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @export
djfst <- function(x, mu = 0, sigma = 1, alpha = 3, theta = 3, log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    z <- (x - mu) / sigma
    rz <- rep(NA_real_, length(z))
    zero <- !is.na(z) & z == 0
    rz[zero] <- 0
    nonzero <- !is.na(z) & z != 0
    rz[nonzero] <- sign(z[nonzero]) / sqrt((alpha + theta) / (z[nonzero]^2) + 1)

    lp <- (alpha + 0.5) * log1p(rz) +
           (theta + 0.5) * log1p(-rz) -
           ((alpha + theta - 1) * log(2)) -
           (0.5 * log(alpha + theta)) -
           lbeta(alpha, theta) -
           log(sigma)

    if (log) lp else exp(lp)
}

#' @rdname jfst-distribution
#' @export
pjfst <- function(q, mu = 0, sigma = 1, alpha = 3, theta = 3,
                  lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    z <- (q - mu) / sigma
    rz <- rep(NA_real_, length(z))
    zero <- !is.na(z) & z == 0
    rz[zero] <- 0
    nonzero <- !is.na(z) & z != 0
    rz[nonzero] <- sign(z[nonzero]) / sqrt((alpha + theta) / (z[nonzero]^2) + 1)

    r <- 0.5 * (1 + rz)
    p <- stats::pbeta(r, alpha, theta)

    if (!lower.tail) p <- 1 - p
    if (log.p) p <- log(p)

    p
}

#' @rdname jfst-distribution
#' @export
qjfst <- function(p, mu = 0, sigma = 1, alpha = 3, theta = 3,
                  lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(p)))

    pp <- if (log.p) exp(p) else p
    if (!lower.tail) pp <- 1 - pp

    if (any(pp < 0 | pp > 1, na.rm = TRUE)) {
        stop("p must be strictly between 0 and 1")
    }

    balpha <- stats::qbeta(pp, alpha, theta)
    zalpha <- (sqrt(alpha + theta)) * (2 * balpha - 1) /
        (2 * sqrt(balpha * (1 - balpha)))

    mu + sigma * zalpha
}

#' @rdname jfst-distribution
#' @export
rjfst <- function(n, mu = 0, sigma = 1, alpha = 3, theta = 3) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, n))

    n <- ceiling(n)
    u <- stats::runif(n)
    qjfst(u, mu = mu, sigma = sigma, alpha = alpha, theta = theta)
}
