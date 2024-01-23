\name{mastSpectralDensity}
\alias{mastSpectralDensity}
\title{Volatility and period for mast data}
\description{
Extracts time series attributes for tree or population fecundity.
}
\usage{
  mastSpectralDensity( x, maxPeriod = length(x)/2, PLOT = FALSE, ylim = NULL )
}
\arguments{
  \item{x}{\code{numeric vector} of sequential fecundity values.}
  
  \item{maxPeriod}{the number of frequencies/periods to include.}
  
  \item{PLOT}{\code{logical} to plot spectrum.}
  
  \item{ylim}{if \code{PLOT}, values for lower and upper bounds.}
}

\details{
Returns attributes of volatility and period for a sequence of fecundity values for a single tree (or population) in \code{x}, which may often be on a log scale.  

More detailed examples can be obtained with: \code{browseVignettes('mastif')}
}

\value{
Returns a \code{list} that includes \code{spect}, a \code{matrix} of power values ordered by frequency (1/period). To permit comparisons between series that differ in length, \code{totVar} (total variance) and \code{volatility} (period-weighted variance) are divided by the length of the series. Mean and standard deviation for the weighted period are \code{periodMu} and \code{periodSd}.
}

\examples{
\donttest{
d <- "https://github.com/jimclarkatduke/mast/blob/master/outputAbies.rdata?raw=True"
repmis::source_data( d )

# single Abies tree from fitted output$prediction$fecPred:

wt <- which( fecPred$treeID == "BAMT1-1" )
s  <- mastSpectralDensity( log(fecPred$fecEstMu[ wt ]), PLOT = TRUE )

# population year effects (log scale) for an ecoRegion_species in output$parameters$betaYrRand:

x <- betaYrRand['3_abiesAmabilis', ]    
x <- x[ x != 0 ]                        # ecoRegion_species vary in observation years
s <- mastSpectralDensity( x, PLOT = TRUE )
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Qiu, T, ..., and J.S. Clark. 2023. Mutualist dispersers and the global distribution of masting: mediation by climate and fertility. in review.
}
\seealso{

\code{\link{mastif}} for analysis

A more detailed vignette is can be obtained with:

\code{browseVignettes('mastif')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

