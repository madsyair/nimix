## dist-gmsnburr-core.R --------------------------------------------------------
## GMSNBurr d/p/q/r -- maintainer reference implementation (A. S. Choir),
## adopted verbatim. GMSNBurr generalises MSNBurr (theta=1) and MSNBurr-IIa
## (alpha=1) with two shape parameters. Iriawan (2000); Choir (2020).

#' @include neonorm-utils.R
NULL

#' GMSNBurr Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Generalized MSNBurr (GMSNBurr) distribution.
#'
#' @name gmsnburr-distribution
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Location parameter (default = 0).
#' @param sigma Numeric. Scale parameter, must be positive (default = 1).
#' @param alpha Numeric. First shape parameter, must be positive (default = 1).
#' @param theta Numeric. Second shape parameter, must be positive (default = 1).
#' @param log,log.p Logical. If TRUE, probabilities/densities are given as log.
#'   Default is FALSE.
#' @param lower.tail Logical. If TRUE (default), probabilities are
#'   \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dgmsnburr} gives the density, \code{pgmsnburr} gives the
#' distribution function, \code{qgmsnburr} gives the quantile function,
#' and \code{rgmsnburr} generates random deviates.
#'
#' @details
#' The GMSNBurr distribution with parameters \eqn{\mu}, \eqn{\sigma},
#' \eqn{\alpha}, and \eqn{\theta} has probability density function:
#'
#' \deqn{f(x|\mu,\sigma,\alpha,\theta) = \frac{\omega}{B(\alpha,\beta)\sigma}
#'   \left(\frac{\theta}{\alpha}\right)^\theta
#'   \exp\left(-\theta\omega\frac{x-\mu}{\sigma}\right)
#'   \left(1 + \frac{\theta}{\alpha}
#'   \exp\left(-\omega\frac{x-\mu}{\sigma}\right)\right)^{-(\alpha+\theta)}}
#'
#' Special cases:
#' \itemize{
#'   \item \eqn{\alpha = \theta}: Symmetric distribution
#'   \item \eqn{\alpha < \theta}: Left-skewed
#'   \item \eqn{\alpha > \theta}: Right-skewed
#'   \item \eqn{\theta = 1}: Reduces to MSNBurr
#'   \item \eqn{\alpha = 1}: Reduces to MSNBurr-IIa
#'   \item \eqn{\alpha = \theta \to \infty}: Converges to \eqn{\mathcal{N}(\mu, \sigma^2)}
#' }
#'
#' @references
#' Choir, A. S. (2020). The New Neo-Normal Distributions and their Properties.
#' Dissertation. Institut Teknologi Sepuluh Nopember.
#'
#' Iriawan, N. (2000). Computationally Intensive Approaches to Inference
#' in Neo-Normal Linear Models. PhD Thesis, Curtin University of Technology.
#'
#' @seealso [MSNBurr()], [MSNBurr2a()]
#'
#' @examples
#' dgmsnburr(0, mu = 0, sigma = 1, alpha = 1, theta = 1)
#'
#' x <- seq(-4, 4, by = 0.1)
#' plot(x, dgmsnburr(x, alpha = 1, theta = 1), type = "l",
#'      ylab = "Density", main = "GMSNBurr Densities")
#' lines(x, dgmsnburr(x, alpha = 2, theta = 1), col = "red")
#' lines(x, dgmsnburr(x, alpha = 1, theta = 2), col = "blue")
#'
#' pgmsnburr(c(-2, 0, 2))
#' qgmsnburr(c(0.025, 0.5, 0.975))
#'
#' set.seed(123)
#' r <- rgmsnburr(1000)
#' hist(r, breaks = 30, freq = FALSE, main = "GMSNBurr Random Samples")
#' curve(dgmsnburr(x), add = TRUE, col = "red", lwd = 2)
#'
#' @keywords distribution
#' @export
dgmsnburr <- function(x, mu = 0, sigma = 1, alpha = 1, theta = 1,
                      log = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(x)))

    lomega <- .log_omega_gmsnburr(alpha, theta)
    zo <- -exp(lomega) * ((x - mu) / sigma)
    zoa <- zo + log(theta) - log(alpha)

    lp <- lomega - log(sigma) + theta * (log(theta) - log(alpha)) +
        theta * zo - (alpha + theta) * log1pexp(zoa) - lbeta(alpha, theta)

    lp[is.infinite(zo) & zo > 0] <- -Inf

    if (log) lp else exp(lp)
}

#' @rdname gmsnburr-distribution
#' @export
pgmsnburr <- function(q, mu = 0, sigma = 1, alpha = 1, theta = 1,
                      lower.tail = TRUE, log.p = FALSE) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, length(q)))

    omega <- exp(.log_omega_gmsnburr(alpha, theta))
    z <- -omega * ((q - mu) / sigma)
    log_ratio <- log(theta) - log(alpha)
    y <- -(z + log_ratio)

    ep <- stats::plogis(y)
    one_minus_ep <- stats::plogis(-y)

    result <- numeric(length(ep))
    result[is.na(ep)] <- NA_real_
    use_direct <- ep < 0.5

    if (any(use_direct, na.rm = TRUE)) {
        idx <- which(use_direct)
        result[idx] <- stats::pbeta(ep[idx], alpha, theta,
                                     lower.tail = lower.tail, log.p = log.p)
    }
    if (any(!use_direct, na.rm = TRUE)) {
        idx <- which(!use_direct)
        result[idx] <- stats::pbeta(one_minus_ep[idx], theta, alpha,
                                     lower.tail = !lower.tail, log.p = log.p)
    }
    result
}

#' @rdname gmsnburr-distribution
#' @export
qgmsnburr <- function(p, mu = 0, sigma = 1, alpha = 1, theta = 1,
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

    omega <- exp(.log_omega_gmsnburr(alpha, theta))

    np <- length(pp)
    mu <- rep_len(mu, np); sigma <- rep_len(sigma, np)
    alpha <- rep_len(alpha, np); theta <- rep_len(theta, np)
    omega <- rep_len(omega, np)

    result <- numeric(length(pp))
    result[pp == 0] <- -Inf
    result[pp == 1] <- Inf
    result[is.na(pp)] <- NA

    interior <- !is.na(pp) & pp > 0 & pp < 1
    if (any(interior)) {
        p_int <- pp[interior]

        bp <- stats::qbeta(p_int, alpha[interior], theta[interior])
        one_minus_bp <- stats::qbeta(1 - p_int, theta[interior], alpha[interior])

        log_1mbp <- ifelse(one_minus_bp > 0, log(one_minus_bp), -Inf)
        log_bp <- ifelse(bp > 0, log(bp), -Inf)
        log_ratio <- log_1mbp - log_bp + log(alpha) - log(theta)

        result[interior] <- mu[interior] - (sigma[interior] / omega[interior]) * log_ratio
    }

    result
}

#' @rdname gmsnburr-distribution
#' @export
rgmsnburr <- function(n, mu = 0, sigma = 1, alpha = 1, theta = 1) {
    msgs <- .check_params_warn(
        list(mu = mu, sigma = sigma, alpha = alpha, theta = theta),
        c("sigma", "alpha", "theta"))
    if (length(msgs) > 0) return(rep(NaN, n))

    omega <- exp(.log_omega_gmsnburr(alpha, theta))
    X <- pmax(stats::rgamma(n, shape = alpha, rate = 1),
              .Machine$double.xmin)
    Y <- pmax(stats::rgamma(n, shape = theta, rate = 1),
              .Machine$double.xmin)

    log_ratio <- log(Y) - log(X) + log(alpha) - log(theta)
    mu - (sigma / omega) * log_ratio
}
