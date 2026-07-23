#' nimix: Bayesian Mixture Clustering and Regression with NIMBLE
#'
#' Bayesian mixture modelling built on the \pkg{nimble} platform. The package
#' implements univariate and multivariate mixture clustering and
#' mixture-of-regressions through two inference engines: a Dirichlet Process
#' Mixture (DPM) engine based on the Chinese Restaurant Process (which estimates
#' the number of occupied components) and a fixed-K finite-mixture engine. It is
#' organised around an extensible S4 \code{\linkS4class{DistributionSpec}}
#' contract so that new component distributions and engines can be added without
#' rewriting existing code.
#'
#' @section Inference engines:
#' \itemize{
#'   \item \code{method = "dpm"}: Dirichlet process / Chinese restaurant process;
#'     the number of occupied components is estimated from the data.
#'   \item \code{method = "fixedk"}: finite mixture with a known number of
#'     components \code{K}.
#' }
#'
#' @section Component distributions:
#' Gaussian (univariate and multivariate), Student-t and Normal-Gamma
#' (heavy-tailed, univariate and multivariate), and Poisson / Binomial counts,
#' for both clustering (\code{\link{nimixClust}}) and regression
#' (\code{\link{nimixReg}}, including multivariate responses).
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
  "getNimbleOption", "nf_preProcessMemberDataObject", "nimNumeric",
  "dMSNBurr_k", "rMSNBurr_k", "dMSNBurr2a_k", "rMSNBurr2a_k", "dGMSNBurr_k", "rGMSNBurr_k", "dSEP_k", "rSEP_k", "dLEP_k", "rLEP_k", "dFSSN_k", "rFSSN_k", "dFOSSEP_k", "rFOSSEP_k", "dFSST_k", "rFSST_k", "dJFST_k", "rJFST_k", "dSkewMvN_k", "rSkewMvN_k", "dSkewMvIT_k", "rSkewMvIT_k", "dSkewMvNO_k", "rSkewMvNO_k", "dSkewMvITO_k", "rSkewMvITO_k", "dSkewMvNOG_k", "dSkewMvITOG_k", "dJFST_k", "rJFST_k",
  "n", "L", "d", "y", "xi", "alpha", "aAlpha", "bAlpha",
  "mu0", "kappa0", "nu0", "s0", "df0", "S0",
  "muTilde", "s2Tilde", "covTilde", "muObs", "covObs", "covMu",
  "p", "X", "b0", "B0", "betaTilde", "betaObs", "s2Obs", "covBeta",
  "K", "z", "weights", "alphaVec", "inprod",
  "omega", "tauTilde", "df", "muSd", "aTau", "bTau",
  "lambda", "prob", "size", "a0", "b0", "rmnorm_chol",
  "betaTilde", "betaObs", "mu", "pp", "B0", "alphaVec", "weights", "K",
  "log<-", "logit<-", "tauObs", "s2Obs", "omega", "mb0", "Bcov", "v0", "covTilde", "covObs", "covMu",
  # Neo-normal shape nodes and the multivariate scale node.
  "nuTilde", "thetaTilde", "gamTilde", "SigTilde",
  # Random-effect nodes: free offsets (bf, sf), their sum-to-zero completions,
  # the group count and the slope covariate.
  "G", "bf", "sf", "xRE",
  # NIMBLE DSL functions: these exist only inside nimbleFunction/nimbleCode
  # bodies, where NIMBLE compiles rather than evaluates them.
  "asCol", "returnType", "rcat", "ilogit",
  "dt_nonstandard", "rt_nonstandard", "dbinom", "rbinom", "dpois", "rpois"
))
