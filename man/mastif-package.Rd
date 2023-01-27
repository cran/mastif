\name{mastif-package}
\alias{mastif-package}
\docType{package}
\concept{ecology}
\concept{Bayesian analysis}
\concept{seed production}
\concept{source detection}
\concept{hierarchical Bayes}
\concept{Gibbs sampling}
\concept{Hamiltonian Markov chain}
\concept{seed dispersal}
\concept{tree fecundity}
\concept{masting}
\concept{autoregressive model}
\concept{state-space model}
\concept{prediction}
\title{Mast Inference and Forecasting}
\description{
Seed production is estimated from censuses of trees and seed collections from traps.  From locations of known trees and seed traps, infers source strength, coefficients for predictor variables, and parameters for a dispersal kernel.  Fecundity is a state-space model allowing for random individual (tree) effects, random year effects and random AR(p) lag effects.  Estimates unknown redistribution of seed types to known species identities of trees.  Functions begin with 'mast' to avoid conflicts with other packages.
}
\details{
  \tabular{ll}{
  Package: \tab mastif\cr
  Type: \tab Package\cr
  Version: \tab 2.0\cr
  Date: \tab 2023-1-20\cr
  License: \tab GPL (>= 2)\cr
  URL: \tab http://sites.nicholas.duke.edu/clarklab/code/\cr
  }
  
The package mastif estimates fecundity of trees and dispersion of seed observed at seed traps, using information on locations of sources and detectors, and covariates that could explain source strength.  Data sets of this type are common and used to understand a range of processes related to seed dispersal, masting, environmental controls on reproduction, sex ratio, and allocation.

Posterior simulation is done by Gibbs sampling.  Analysis is done by these functions:

\code{\link{mastif}} fits model with Gibbs sampling.

\code{\link{mastSim}} simulates data for analysis by \code{mastif}.

\code{\link{mastFillCensus}} aligns sample years in tree census data with seed trap data.

\code{\link{mastClimate}} annotates tree data with covariates for fecundity modeling.

\code{\link{mastPlot}} generates plots of the output from \code{\link{mastif}}.
}


\author{
  Author: James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, JS, C Nunes, and B Tomasek. 2019. Masting as an unreliable resource: spatio-temporal host diversity merged with consumer movement, storage, and diet. Ecological Monographs, e01381.
}
\keyword{ package }
\seealso{
  \code{\link{mastif}}, 
  \code{\link{mastSim}}
  
A more detailed vignette is can be obtained with:

\code{browseVignettes('mastif')}

}
