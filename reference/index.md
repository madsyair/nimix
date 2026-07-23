# Package index

## Package

- [`nimix`](https://madsyair.github.io/nimix/reference/nimix-package.md)
  [`nimix-package`](https://madsyair.github.io/nimix/reference/nimix-package.md)
  : nimix: Bayesian Mixture Clustering and Regression with NIMBLE

## Fitting functions

User-facing entry points for clustering and mixture-of-regression.

- [`nimixClust()`](https://madsyair.github.io/nimix/reference/nimixClust.md)
  : Bayesian mixture clustering
- [`nimixReg()`](https://madsyair.github.io/nimix/reference/nimixReg.md)
  : Bayesian mixture of linear regressions

## Engines

Inference engines selected via the `method` argument.

- [`DPMEngine()`](https://madsyair.github.io/nimix/reference/DPMEngine.md)
  : Construct a DPM engine configuration
- [`FixedKEngine()`](https://madsyair.github.io/nimix/reference/FixedKEngine.md)
  : Construct a fixed-K finite-mixture engine configuration
- [`DPMEngine-class`](https://madsyair.github.io/nimix/reference/DPMEngine-class.md)
  : Dirichlet Process Mixture engine (native NIMBLE dCRP)
- [`FixedKEngine-class`](https://madsyair.github.io/nimix/reference/FixedKEngine-class.md)
  : Finite-mixture engine with fixed, known K
- [`EngineConfig-class`](https://madsyair.github.io/nimix/reference/EngineConfig-class.md)
  : Virtual base class for mixture sampling engines
- [`HMMEngine()`](https://madsyair.github.io/nimix/reference/HMMEngine.md)
  : Construct a hidden-Markov mixture engine configuration

## Distribution registry

Register and look up component distributions by name.

- [`registerDistribution()`](https://madsyair.github.io/nimix/reference/registerDistribution.md)
  : Register a component distribution
- [`getDistribution()`](https://madsyair.github.io/nimix/reference/getDistribution.md)
  : Retrieve a registered distribution by name
- [`listDistributions()`](https://madsyair.github.io/nimix/reference/listDistributions.md)
  : List registered distribution names

## Clustering component specs

- [`NormalUvSpec()`](https://madsyair.github.io/nimix/reference/NormalUvSpec.md)
  : Construct a univariate Gaussian component spec
- [`NormalMvSpec()`](https://madsyair.github.io/nimix/reference/NormalMvSpec.md)
  : Construct a multivariate Gaussian component spec
- [`StudentTUvSpec()`](https://madsyair.github.io/nimix/reference/StudentTUvSpec.md)
  : Construct a univariate Student-t component spec
- [`StudentTMvSpec()`](https://madsyair.github.io/nimix/reference/StudentTMvSpec.md)
  : Construct a multivariate Student-t component spec
- [`NormalGammaUvSpec()`](https://madsyair.github.io/nimix/reference/NormalGammaUvSpec.md)
  : Construct a univariate Normal-Gamma component spec
- [`NormalGammaMvSpec()`](https://madsyair.github.io/nimix/reference/NormalGammaMvSpec.md)
  : Construct a multivariate Normal-Gamma component spec
- [`PoissonSpec()`](https://madsyair.github.io/nimix/reference/PoissonSpec.md)
  : Construct a Poisson component spec
- [`BinomialSpec()`](https://madsyair.github.io/nimix/reference/BinomialSpec.md)
  : Construct a Binomial component spec

## Regression component specs

- [`NormalRegSpec()`](https://madsyair.github.io/nimix/reference/NormalRegSpec.md)
  : Construct a Normal-linear regression component spec
- [`StudentTRegSpec()`](https://madsyair.github.io/nimix/reference/StudentTRegSpec.md)
  : Construct a Student-t regression component spec
- [`NormalGammaRegSpec()`](https://madsyair.github.io/nimix/reference/NormalGammaRegSpec.md)
  : Construct a Normal-Gamma regression component spec
- [`PoissonRegSpec()`](https://madsyair.github.io/nimix/reference/PoissonRegSpec.md)
  : Construct a Poisson regression component spec
- [`BinomialRegSpec()`](https://madsyair.github.io/nimix/reference/BinomialRegSpec.md)
  : Construct a Binomial regression component spec
- [`NormalMvRegSpec()`](https://madsyair.github.io/nimix/reference/NormalMvRegSpec.md)
  : Construct a multivariate-response Normal regression spec
- [`StudentTMvRegSpec()`](https://madsyair.github.io/nimix/reference/StudentTMvRegSpec.md)
  : Construct a multivariate-response Student-t regression spec
- [`NormalGammaMvRegSpec()`](https://madsyair.github.io/nimix/reference/NormalGammaMvRegSpec.md)
  : Construct a multivariate-response Normal-Gamma regression spec

## Neo-normal regression component specs

Skew and heavy-tailed mixture-of-regressions: each component’s location
is a linear predictor, with per-component scale and shape. Available
across the fixed-K, DPM, and HMM (Markov-switching) engines, with random
effects under fixed-K and DPM.

- [`MSNBurrRegSpec()`](https://madsyair.github.io/nimix/reference/MSNBurrRegSpec-class.md)
  : MSNBurr regression specification
- [`MSNBurr2aRegSpec()`](https://madsyair.github.io/nimix/reference/MSNBurr2aRegSpec-class.md)
  : MSNBurr-IIa regression specification
- [`GMSNBurrRegSpec()`](https://madsyair.github.io/nimix/reference/GMSNBurrRegSpec-class.md)
  : GMSNBurr regression specification
- [`FSSNRegSpec()`](https://madsyair.github.io/nimix/reference/FSSNRegSpec-class.md)
  : FSSN regression specification
- [`FSSTRegSpec()`](https://madsyair.github.io/nimix/reference/FSSTRegSpec-class.md)
  : FSST regression specification
- [`SEPRegSpec()`](https://madsyair.github.io/nimix/reference/SEPRegSpec-class.md)
  : SEP regression specification
- [`LEPRegSpec()`](https://madsyair.github.io/nimix/reference/LEPRegSpec-class.md)
  : LEP regression specification
- [`FOSSEPRegSpec()`](https://madsyair.github.io/nimix/reference/FOSSEPRegSpec-class.md)
  : FOSSEP regression specification
- [`JFSTRegSpec()`](https://madsyair.github.io/nimix/reference/JFSTRegSpec-class.md)
  : JFST regression specification

## Extension API (generics)

The public S4 contract a component spec implements so a new family plugs
into every engine without engine edits.

- [`buildModelCode()`](https://madsyair.github.io/nimix/reference/buildModelCode.md)
  : Build the NIMBLE model code for a (distribution, engine) pair
- [`defaultPrior()`](https://madsyair.github.io/nimix/reference/defaultPrior.md)
  : Build a data-scaled default prior for a distribution
- [`validateParams()`](https://madsyair.github.io/nimix/reference/validateParams.md)
  : Validate component parameters or a prior specification
- [`customizeSamplers()`](https://madsyair.github.io/nimix/reference/customizeSamplers.md)
  : Customise MCMC samplers for a component spec
- [`relabel()`](https://madsyair.github.io/nimix/reference/relabel.md) :
  Correct label switching in a fitted mixture
- [`psm()`](https://madsyair.github.io/nimix/reference/psm.md) :
  Posterior similarity matrix
- [`binderPartition()`](https://madsyair.github.io/nimix/reference/binderPartition.md)
  : Binder-loss point partition (Dahl's least-squares criterion)
- [`plotClusterMap()`](https://madsyair.github.io/nimix/reference/plotClusterMap.md)
  : Map the clusters of a spatial mixture fit
- [`clusterValidity()`](https://madsyair.github.io/nimix/reference/clusterValidity.md)
  : Internal validity indices for a fitted clustering
- [`viterbiPath()`](https://madsyair.github.io/nimix/reference/viterbiPath.md)
  : Most probable state path of a hidden-Markov mixture fit
- [`nimixForecast()`](https://madsyair.github.io/nimix/reference/nimixForecast.md)
  : Forecast ahead from a regime-switching (HMM) fit
- [`posteriorLinpred()`](https://madsyair.github.io/nimix/reference/posteriorLinpred.md)
  : Per-component linear predictors
- [`posteriorEpred()`](https://madsyair.github.io/nimix/reference/posteriorEpred.md)
  : Expected value of the posterior predictive
- [`posteriorPredictive()`](https://madsyair.github.io/nimix/reference/posteriorPredictive.md)
  : Draws from the posterior predictive distribution
- [`isRegressionSpec()`](https://madsyair.github.io/nimix/reference/isRegressionSpec.md)
  : Is this a regression component spec?
- [`linkInv()`](https://madsyair.github.io/nimix/reference/linkInv.md) :
  Inverse link for a regression component
- [`responseRng()`](https://madsyair.github.io/nimix/reference/responseRng.md)
  : Draw a response given linear predictors and error scale

## S4 class definitions

Virtual base classes, model/result containers, and the class behind each
component spec. Most users work through the constructor functions above;
these document the underlying S4 classes.

- [`DistributionSpec-class`](https://madsyair.github.io/nimix/reference/DistributionSpec-class.md)
  : Virtual base class for mixture component distributions
- [`MixtureModel-class`](https://madsyair.github.io/nimix/reference/MixtureModel-class.md)
  : Virtual base class for mixture models
- [`ClusterModel-class`](https://madsyair.github.io/nimix/reference/ClusterModel-class.md)
  : Mixture clustering model (for nimixClust)
- [`RegressionMixModel-class`](https://madsyair.github.io/nimix/reference/RegressionMixModel-class.md)
  : Mixture-of-regressions model (for nimixReg)
- [`show(`*`<FitResult>`*`)`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  [`summary(`*`<FitResult>`*`)`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  [`plot(`*`<FitResult>`*`,`*`<missing>`*`)`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  [`predict(`*`<FitResult>`*`)`](https://madsyair.github.io/nimix/reference/FitResult-class.md)
  : Fitted nimix mixture result
- [`NormalUvSpec-class`](https://madsyair.github.io/nimix/reference/NormalUvSpec-class.md)
  : Univariate Gaussian component specification
- [`NormalMvSpec-class`](https://madsyair.github.io/nimix/reference/NormalMvSpec-class.md)
  : Multivariate Gaussian component specification
- [`StudentTUvSpec-class`](https://madsyair.github.io/nimix/reference/StudentTUvSpec-class.md)
  : Univariate Student-t component specification
- [`StudentTMvSpec-class`](https://madsyair.github.io/nimix/reference/StudentTMvSpec-class.md)
  : Multivariate Student-t component specification (direct density)
- [`NormalGammaUvSpec-class`](https://madsyair.github.io/nimix/reference/NormalGammaUvSpec-class.md)
  : Univariate Normal-Gamma (scale-mixture Student-t) component
  specification
- [`NormalGammaMvSpec-class`](https://madsyair.github.io/nimix/reference/NormalGammaMvSpec-class.md)
  : Multivariate Normal-Gamma (scale-mixture multivariate-t) component
- [`PoissonSpec-class`](https://madsyair.github.io/nimix/reference/PoissonSpec-class.md)
  : Poisson count component specification
- [`BinomialSpec-class`](https://madsyair.github.io/nimix/reference/BinomialSpec-class.md)
  : Binomial count component specification
- [`NormalRegSpec-class`](https://madsyair.github.io/nimix/reference/NormalRegSpec-class.md)
  : Normal-linear regression component specification
- [`StudentTRegSpec-class`](https://madsyair.github.io/nimix/reference/StudentTRegSpec-class.md)
  : Student-t mixture regression component (direct t density)
- [`NormalGammaRegSpec-class`](https://madsyair.github.io/nimix/reference/NormalGammaRegSpec-class.md)
  : Normal-Gamma mixture regression component (conjugate scale-mixture)
- [`PoissonRegSpec-class`](https://madsyair.github.io/nimix/reference/PoissonRegSpec-class.md)
  : Poisson GLM regression component (log link)
- [`BinomialRegSpec-class`](https://madsyair.github.io/nimix/reference/BinomialRegSpec-class.md)
  : Binomial GLM regression component (logit link)
- [`NormalMvRegSpec-class`](https://madsyair.github.io/nimix/reference/NormalMvRegSpec-class.md)
  : Multivariate-response Normal mixture regression
- [`StudentTMvRegSpec-class`](https://madsyair.github.io/nimix/reference/StudentTMvRegSpec-class.md)
  : Multivariate-response Student-t mixture regression (direct density)
- [`NormalGammaMvRegSpec-class`](https://madsyair.github.io/nimix/reference/NormalGammaMvRegSpec-class.md)
  : Multivariate-response Normal-Gamma mixture regression (scale
  mixture)

## Neo-normal component specs

Skew-capable MSNBurr families (Iriawan 2000; Choir 2020).

- [`MSNBurrUvSpec()`](https://madsyair.github.io/nimix/reference/MSNBurrUvSpec-class.md)
  : MSNBurr mixture components (neo-normal, left-skew capable)
- [`MSNBurr2aUvSpec()`](https://madsyair.github.io/nimix/reference/MSNBurr2aUvSpec-class.md)
  : MSNBurr-IIa mixture components (neo-normal, right-skew capable)
- [`GMSNBurrUvSpec()`](https://madsyair.github.io/nimix/reference/GMSNBurrUvSpec-class.md)
  : GMSNBurr mixture components (generalized neo-normal)
- [`SEPUvSpec()`](https://madsyair.github.io/nimix/reference/SEPUvSpec-class.md)
  : SEP mixture components (symmetric exponential power)
- [`LEPUvSpec()`](https://madsyair.github.io/nimix/reference/LEPUvSpec-class.md)
  : LEP mixture components (exponential power, alternative
  parameterisation)
- [`FSSNUvSpec()`](https://madsyair.github.io/nimix/reference/FSSNUvSpec-class.md)
  : FSSN mixture components (Fernandez-Steel skew Normal)
- [`FOSSEPUvSpec()`](https://madsyair.github.io/nimix/reference/FOSSEPUvSpec-class.md)
  : FOSSEP mixture components (Fernandez-Steel skew exponential power)
- [`FSSTUvSpec()`](https://madsyair.github.io/nimix/reference/FSSTUvSpec-class.md)
  : FSST mixture components (Fernandez-Steel skew Student-t)
- [`JFSTUvSpec()`](https://madsyair.github.io/nimix/reference/JFSTUvSpec-class.md)
  : JFST mixture components (Jones-Faddy skew-t)
- [`SkewNormalMvSpec()`](https://madsyair.github.io/nimix/reference/SkewNormalMvSpec-class.md)
  : Skew multivariate Normal mixture components (Ferreira-Steel)
- [`SkewIStudentMvSpec()`](https://madsyair.github.io/nimix/reference/SkewIStudentMvSpec-class.md)
  : Skew multivariate independent-Student mixture components
  (Ferreira-Steel)
- [`SkewNormalMvOSpec()`](https://madsyair.github.io/nimix/reference/SkewNormalMvOSpec-class.md)
  : Skew multivariate Normal components with estimated orthogonal factor
- [`SkewIStudentMvOSpec()`](https://madsyair.github.io/nimix/reference/SkewIStudentMvOSpec-class.md)
  : Skew mv independent-Student components with estimated orthogonal
  factor
- [`SkewNormalMvOGenSpec()`](https://madsyair.github.io/nimix/reference/SkewNormalMvOGenSpec-class.md)
  : Skew multivariate Normal components with estimated O, general
  dimension
- [`SkewIStudentMvOGenSpec()`](https://madsyair.github.io/nimix/reference/SkewIStudentMvOGenSpec-class.md)
  : Skew mv independent-Student components with estimated O, general
  dimension
- [`orthogonalFactor()`](https://madsyair.github.io/nimix/reference/orthogonalFactor.md)
  : Orthogonal factor from Householder angles
- [`canonicaliseO()`](https://madsyair.github.io/nimix/reference/canonicaliseO.md)
  : Canonical representative of an FS orthogonal factor

## Neo-normal distribution functions

Numerically stable density, distribution, quantile, and RNG.

- [`dmsnburr()`](https://madsyair.github.io/nimix/reference/msnburr-distribution.md)
  [`pmsnburr()`](https://madsyair.github.io/nimix/reference/msnburr-distribution.md)
  [`qmsnburr()`](https://madsyair.github.io/nimix/reference/msnburr-distribution.md)
  [`rmsnburr()`](https://madsyair.github.io/nimix/reference/msnburr-distribution.md)
  : MSNBurr Distribution
- [`dmsnburr2a()`](https://madsyair.github.io/nimix/reference/msnburr2a-distribution.md)
  [`pmsnburr2a()`](https://madsyair.github.io/nimix/reference/msnburr2a-distribution.md)
  [`qmsnburr2a()`](https://madsyair.github.io/nimix/reference/msnburr2a-distribution.md)
  [`rmsnburr2a()`](https://madsyair.github.io/nimix/reference/msnburr2a-distribution.md)
  : MSNBurr-IIa Distribution
- [`dgmsnburr()`](https://madsyair.github.io/nimix/reference/gmsnburr-distribution.md)
  [`pgmsnburr()`](https://madsyair.github.io/nimix/reference/gmsnburr-distribution.md)
  [`qgmsnburr()`](https://madsyair.github.io/nimix/reference/gmsnburr-distribution.md)
  [`rgmsnburr()`](https://madsyair.github.io/nimix/reference/gmsnburr-distribution.md)
  : GMSNBurr Distribution
- [`dsep()`](https://madsyair.github.io/nimix/reference/sep-distribution.md)
  [`psep()`](https://madsyair.github.io/nimix/reference/sep-distribution.md)
  [`qsep()`](https://madsyair.github.io/nimix/reference/sep-distribution.md)
  [`rsep()`](https://madsyair.github.io/nimix/reference/sep-distribution.md)
  : Subbotin Exponential Power (SEP) Distribution
- [`dlep()`](https://madsyair.github.io/nimix/reference/lep-distribution.md)
  [`plep()`](https://madsyair.github.io/nimix/reference/lep-distribution.md)
  [`qlep()`](https://madsyair.github.io/nimix/reference/lep-distribution.md)
  [`rlep()`](https://madsyair.github.io/nimix/reference/lep-distribution.md)
  : Lunetta Exponential Power Distribution
- [`dfssn()`](https://madsyair.github.io/nimix/reference/fssn-distribution.md)
  [`pfssn()`](https://madsyair.github.io/nimix/reference/fssn-distribution.md)
  [`qfssn()`](https://madsyair.github.io/nimix/reference/fssn-distribution.md)
  [`rfssn()`](https://madsyair.github.io/nimix/reference/fssn-distribution.md)
  : Fernandez-Steel Skew Normal Distribution
- [`dfossep()`](https://madsyair.github.io/nimix/reference/fossep-distribution.md)
  [`pfossep()`](https://madsyair.github.io/nimix/reference/fossep-distribution.md)
  [`qfossep()`](https://madsyair.github.io/nimix/reference/fossep-distribution.md)
  [`rfossep()`](https://madsyair.github.io/nimix/reference/fossep-distribution.md)
  : Fernandez-Osiewalski-Steel Skew Exponential Power Distribution
- [`dfsst()`](https://madsyair.github.io/nimix/reference/fsst-distribution.md)
  [`pfsst()`](https://madsyair.github.io/nimix/reference/fsst-distribution.md)
  [`qfsst()`](https://madsyair.github.io/nimix/reference/fsst-distribution.md)
  [`rfsst()`](https://madsyair.github.io/nimix/reference/fsst-distribution.md)
  : Fernandez-Steel Skew t Distribution
- [`djfst()`](https://madsyair.github.io/nimix/reference/jfst-distribution.md)
  [`pjfst()`](https://madsyair.github.io/nimix/reference/jfst-distribution.md)
  [`qjfst()`](https://madsyair.github.io/nimix/reference/jfst-distribution.md)
  [`rjfst()`](https://madsyair.github.io/nimix/reference/jfst-distribution.md)
  : Jones-Faddy Skew-t Distribution
- [`dskewmvn()`](https://madsyair.github.io/nimix/reference/skewnormal-mv-distribution.md)
  [`rskewmvn()`](https://madsyair.github.io/nimix/reference/skewnormal-mv-distribution.md)
  : Ferreira-Steel skew multivariate Normal
- [`dskewmvit()`](https://madsyair.github.io/nimix/reference/skewistudent-mv-distribution.md)
  [`rskewmvit()`](https://madsyair.github.io/nimix/reference/skewistudent-mv-distribution.md)
  : Ferreira-Steel skew multivariate independent-Student
- [`dskewmvno()`](https://madsyair.github.io/nimix/reference/skewnormal-mv-o-distribution.md)
  [`rskewmvno()`](https://madsyair.github.io/nimix/reference/skewnormal-mv-o-distribution.md)
  : Ferreira-Steel skew multivariate Normal with estimated orthogonal
  factor
- [`dskewmvito()`](https://madsyair.github.io/nimix/reference/skewistudent-mv-o-distribution.md)
  [`rskewmvito()`](https://madsyair.github.io/nimix/reference/skewistudent-mv-o-distribution.md)
  : Ferreira-Steel skew multivariate independent-Student with estimated
  O

## Spatial mixtures (MRF)

Spatially constrained mixtures whose labels follow a Potts / Markov
random field over a neighbourhood graph.

- [`spatialWeights()`](https://madsyair.github.io/nimix/reference/spatialWeights.md)
  : Construct a SpatialWeightSpec from an adjacency/weight matrix
- [`gridAdjacency()`](https://madsyair.github.io/nimix/reference/gridAdjacency.md)
  : Rook/queen contiguity on a regular grid
- [`spacetimeAdjacency()`](https://madsyair.github.io/nimix/reference/spacetimeAdjacency.md)
  : Build a space-time adjacency for spatio-temporal mixtures
- [`nRegions()`](https://madsyair.github.io/nimix/reference/nRegions.md)
  : Number of regions in a spatial weight structure
- [`getAdjacency()`](https://madsyair.github.io/nimix/reference/getAdjacency.md)
  : Adjacency matrix of a spatial weight structure
- [`neighborsOf()`](https://madsyair.github.io/nimix/reference/neighborsOf.md)
  : Neighbours of one region
- [`MRFEngine()`](https://madsyair.github.io/nimix/reference/MRFEngine-class.md)
  : Markov random field engine (spatially constrained finite mixture)
- [`show(`*`<SpatialWeightSpec>`*`)`](https://madsyair.github.io/nimix/reference/SpatialWeightSpec-class.md)
  : Spatial neighbourhood structure for spatially constrained mixtures

## Model containers

- [`ClusterModel()`](https://madsyair.github.io/nimix/reference/ClusterModel.md)
  : Construct a ClusterModel
- [`RegressionMixModel()`](https://madsyair.github.io/nimix/reference/RegressionMixModel.md)
  : Construct a RegressionMixModel

## Bayesian workflow

Posterior predictive checks; convergence diagnostics print via
summary().

- [`ppCheck()`](https://madsyair.github.io/nimix/reference/ppCheck.md) :
  Posterior predictive check
- [`priorPredictive()`](https://madsyair.github.io/nimix/reference/priorPredictive.md)
  : Prior predictive check for a mixture model
- [`posteriorPredict()`](https://madsyair.github.io/nimix/reference/posteriorPredict.md)
  : Posterior predictive replicates
- [`drawsArray()`](https://madsyair.github.io/nimix/reference/drawsArray.md)
  : Posterior draws as an iterations x chains x parameters array
- [`ppcData()`](https://madsyair.github.io/nimix/reference/ppcData.md) :
  Observed data and posterior predictive replicates for graphical PPC

## Extension API (engine-facing generics)

The remaining S4 contract a component spec implements so a new family
plugs into every engine.

- [`runEngine()`](https://madsyair.github.io/nimix/reference/runEngine.md)
  : Run a mixture engine on a model (internal generic)

- [`componentDensity()`](https://madsyair.github.io/nimix/reference/componentDensity.md)
  : Component density evaluator (R-level, for posterior predictive
  checks)

- [`simulateParams()`](https://madsyair.github.io/nimix/reference/simulateParams.md)
  : Simulate component parameters from a prior (for inits / recovery
  tests)

- [`componentInits()`](https://madsyair.github.io/nimix/reference/componentInits.md)
  : Build dispersed initial values (engine-agnostic)

- [`buildConstants()`](https://madsyair.github.io/nimix/reference/buildConstants.md)
  : Assemble distribution-specific NIMBLE constants for the DPM engine

- [`buildDataList()`](https://madsyair.github.io/nimix/reference/buildDataList.md)
  :

  Shape the observed data into the NIMBLE `data` list

- [`extractParamTraces()`](https://madsyair.github.io/nimix/reference/extractParamTraces.md)
  : Parse raw cluster-parameter traces from the MCMC sample matrix

- [`relabelComponents()`](https://madsyair.github.io/nimix/reference/relabelComponents.md)
  : Permute cluster parameters and build the relabelled component
  summary

## Utilities

- [`nimixClearCache()`](https://madsyair.github.io/nimix/reference/nimixClearCache.md)
  : Clear the compiled-model cache
- [`dmvt_nimix()`](https://madsyair.github.io/nimix/reference/dmvt_nimix.md)
  [`rmvt_nimix()`](https://madsyair.github.io/nimix/reference/dmvt_nimix.md)
  : Multivariate Student-t log density (nimbleFunction)
- [`getEdges()`](https://madsyair.github.io/nimix/reference/getEdges.md)
  : Edge list of a spatial weight structure

## Datasets

- [`usStates2023`](https://madsyair.github.io/nimix/reference/usStates2023.md)
  : US state poverty and income, 2023 (SAIPE official estimates)
- [`usStateAdj`](https://madsyair.github.io/nimix/reference/usStateAdj.md)
  : Contiguity of the 48 contiguous US states + DC (official derivation)
- [`wdi2022`](https://madsyair.github.io/nimix/reference/wdi2022.md) :
  World Development Indicators, 2022 (country-level official statistics)

## Cluster profiling

- [`clusterProfile()`](https://madsyair.github.io/nimix/reference/clusterProfile.md)
  : Profile the clusters of a fitted mixture

## Model selection and ensembling

Predictive comparison (WAIC / PSIS-LOO) and weighted ensembles.

- [`nimixWAIC()`](https://madsyair.github.io/nimix/reference/nimixWAIC.md)
  : WAIC for a fitted mixture
- [`nimixLOO()`](https://madsyair.github.io/nimix/reference/nimixLOO.md)
  : PSIS-LOO for a fitted mixture
- [`modelSelect()`](https://madsyair.github.io/nimix/reference/modelSelect.md)
  : Compare mixture models by predictive fit
- [`ensembleFit()`](https://madsyair.github.io/nimix/reference/ensembleFit.md)
  : Ensemble several fitted mixtures
- [`predict(`*`<nimixEnsemble>`*`)`](https://madsyair.github.io/nimix/reference/predict-nimixEnsemble.md)
  : Weighted predictive density from a mixture ensemble
