#' nimix: Bayesian Mixture Clustering and Regression with NIMBLE
#'
#' Bayesian mixture modelling built on the \pkg{nimble} platform. The
#' release implements univariate and multivariate Gaussian mixture clustering
#' through a Dirichlet Process Mixture (DPM) engine based on the Chinese
#' Restaurant Process. The package is organised around an extensible S4
#' \code{\linkS4class{DistributionSpec}} contract so that new component
#' distributions and (later) a reversible-jump engine can be added without
#' rewriting existing code.
#'
#' @section Roadmap (this is v0.2.0):
#' \itemize{
#'   \item v0.1.0: S4 foundation, \code{\linkS4class{NormalUvSpec}},
#'     univariate \code{\link{nimixClust}} on the DPM engine.
#'   \item v0.2.0 (this release): multivariate clustering
#'     (\code{\linkS4class{NormalMvSpec}}).
#'   \item v0.3.0: mixture-of-regressions (\code{\link{nimixReg}}).
#'   \item v0.5.0+: reversible jump MCMC engine.
#' }
#'
#' @references
#' de Valpine, P., Turek, D., Paciorek, C.J., Anderson-Bergman, C.,
#' Temple Lang, D., & Bodik, R. (2017). Programming with models: writing
#' statistical algorithms for general model structures with NIMBLE.
#' \emph{Journal of Computational and Graphical Statistics}, 26(2), 403--413.
#' \doi{10.1080/10618600.2016.1172487}
#'
#' Neal, R.M. (2000). Markov chain sampling methods for Dirichlet process
#' mixture models. \emph{Journal of Computational and Graphical Statistics},
#' 9(2), 249--265. \doi{10.1080/10618600.2000.10474879}
#'
#' @keywords internal
#'
#' @import methods
#' @importFrom nimble nimbleCode nimbleModel compileNimble configureMCMC buildMCMC runMCMC nimbleFunction
#' @importFrom label.switching ecr.iterative.1 ecr
#' @importFrom coda as.mcmc effectiveSize
#' @importFrom stats kmeans var sd cov dnorm dt rnorm rgamma rWishart quantile predict model.frame model.matrix model.response terms delete.response na.pass
#' @importFrom graphics barplot lines hist matplot legend abline plot.default
#' @importFrom utils modifyList globalVariables
"_PACKAGE"

# NIMBLE BUGS-language symbols that appear inside nimbleCode({...}) blocks are
# not R objects; declare them so R CMD check does not flag "no visible binding"
# (a standard accommodation for packages that build NIMBLE model code).
utils::globalVariables(c(
  "n", "L", "d", "y", "xi", "alpha", "aAlpha", "bAlpha",
  "mu0", "kappa0", "nu0", "s0", "df0", "S0",
  "muTilde", "s2Tilde", "covTilde", "muObs", "covObs", "covMu",
  "p", "X", "b0", "B0", "betaTilde", "betaObs", "s2Obs", "covBeta",
  "K", "z", "weights", "alphaVec", "inprod",
  "omega", "tauTilde", "df", "muSd", "aTau", "bTau",
  "lambda", "prob", "size", "a0", "b0", "rmnorm_chol",
  "betaTilde", "betaObs", "mu", "pp", "B0", "alphaVec", "weights", "K",
  "log<-", "logit<-", "tauObs", "s2Obs", "omega", "mb0", "Bcov", "v0", "covTilde", "covObs", "covMu"
))
