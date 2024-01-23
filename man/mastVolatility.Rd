\name{mastVolatility}
\alias{mastVolatility}
\title{Volatility and period for mast data, combining trees in populations}
\description{
Synthesis of volatility and period at the population scale.
}
\usage{
  mastVolatility( treeID, year, fec, minLength = 6, minFrequency = 1/20 )
}
\arguments{
  \item{treeID}{\code{vector} of labels for individual trees}
  
  \item{year}{\code{numeric vector} for year}
  
  \item{fec}{\code{numeric vector} for fecundity values}
  
  \item{minLength}{determines the minimum number of years to a tree to be included in population estimates}
  
  \item{minFrequency}{lowest frequency to include in volatility, period evaluation}
}

\details{
The three vectors \code{treeID, year, fec} are aligned by tree and year and, thus, of the same length. Tree fecundity values in the \code{numeric vector fec} can differ in number of years due to maturation times, deaths, and observation years.  Trees having fewer than \code{minLength} observations are omitted from the analysis. \code{minFrequency} is high enough to omit low frequencies that are missing in the shortest series to be compared.

More detailed examples can be obtained with: \code{browseVignettes('mastif')}
}

\value{
Returns a \code{list} that includes \code{stats}, which holds the period- and fecundity-weighted estimates of volatility and period at the population scale. The \code{matrix statsDensity} holds the means and standard deviations by period (1/frequency). The \code{matrix mastMatrix} holds for each tree the number of years, mean log fecundity, variance, volatility, and period mean and standard deviation. Returned as tree by frequency are \code{density} and \code{frequency}.
}

\examples{
\donttest{
d <- "https://github.com/jimclarkatduke/mast/blob/master/outputAbies.rdata?raw=True"
repmis::source_data( d )

# all trees in a plot:

wi  <- which( fecPred$plotSpec == 'BERK28 abiesGrandis' ) # tree-years in group
tmp <- mastVolatility( treeID = fecPred$treeID[wi], year = fecPred$year[wi], 
                       fec = fecPred$fecEstMu[wi], minLength = 10 )
period  <- 1/tmp$frequency
density <- tmp$density

plot( NA, xlim = range( period, na.rm = TRUE ), ylim = range( density, na.rm = TRUE ),
      xlab = 'Period (yr)', ylab = 'Density', log = 'xy' )
      
for( i in 1:nrow(density) )lines( period[i,], density[i, ], col = 'grey' )
  lines( tmp$statsDensity['Period', ], tmp$statsDensity['Mean', ], lwd = 2 )
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

