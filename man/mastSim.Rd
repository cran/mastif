\name{mastSim}
\alias{mastSim}
\title{Data simulation for mast model}
\description{
Simulates data for analysis by \code{\link{mastif}} in package \code{mastif}.
}
\usage{
  mastSim(sim)
}
\arguments{
  \item{sim}{\code{list} of inputs for simulation, see Details.}
}

\details{

The \code{list} \code{sim} contains the following:

\code{specNames}: \code{character vector} of species names.

\code{seedNames}: \code{character vector} of seed names.

\code{nyr = 5}:  average number of years for a plot

\code{ntree = 10}:  average number of trees in \code{specNames} on a plot

\code{plotWide = 100}:  diameter of plot

\code{nplot = 3}:  number of plots

\code{ntrap = 20}:  average number of seed traps on a plot

\code{meanDist = 25}:  mean dispersal (meters)

}
       
\value{

Returns an object of \code{\link{class} "mastif"}, a list containing the following components:

\item{distall}{\code{matrix} of distances between all seed traps (rows) and trees (columns).  Block diagonal structure by plot, padded with the value \code{10000} for tree/trap combinations on different plots.}

\item{distall}{seed trap by tree \code{matrix} of distances.}

\item{formulaFec}{\code{formula} for the fecundity model.}

\item{formulaRep}{\code{formula} for the maturation model.}

\item{plots}{\code{character vector} of plot names.}

\item{R}{species to seed type matrix.}

\item{seedData}{\code{data.frame} of seed data needed for \code{mastif}.}

\item{seedNames}{\code{character vector} of seed types.}

\item{sim}{inputs to \code{mastSim}.}

\item{specNames}{\code{character vector} of species.}

\item{treeData}{\code{data.frame} of tree data needed for \code{mastif}.}

\item{trueValues}{\code{list} holding true parameter values includes:

\code{fec}: \code{vector} of log fecundity states. 

\code{repr}: vector of maturation states.

\code{betaFec}: fecundity parameter values.

\code{betaRep}: maturation parameter values.

\code{upar}: dispersal parameter, m^2.
}

\item{xytrap}{\code{data.frame} of seed locations needed for \code{mastif}.}

\item{xytree}{\code{data.frame} of tree locations needed for \code{mastif}.}

\item{years}{\code{numeric vector} years.}

}


\examples{
\donttest{
# simulate data 
seedNames  <- specNames  <- 'acerRubr'
sim <- list(nyr = 10, ntree = 30, nplot = 5, ntrap = 40, 
            specNames = specNames, seedNames = seedNames)
inputs <- mastSim(sim)
output <- mastif( inputs = inputs, ng = 500, burnin = 200 ) 

# increase iterations, then plot:

output   <- mastif( inputs = output, ng = 2000, burnin = 1000 ) 
plotPars <- list(trueValues = inputs$trueValues)
mastPlot(output, plotPars)
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., C. Nunes, and B. Tomasek. 2019. Foodwebs based on unreliable foundations: spatio-temporal masting merged with consumer movement, storage, and diet. Ecological Monographs, in press.
}
\seealso{
\code{\link{mastSim}} simulates data

A more detailed vignette is can be obtained with:

\code{browseVignettes('mastif')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

