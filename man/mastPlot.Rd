\name{mastPlot}
\alias{mastPlot}
\title{Plot mast model}
\description{
Plots data fitted with \code{\link{mastif}} in package \code{mastif}.
}
\usage{
  mastPlot(output, plotPars = NULL)
  
}
\arguments{
  \item{output}{\code{list} of class \code{mast}.}
  
  \item{plotPars}{\code{list} of inputs for plotting can include:
  
  \code{outFolder}:  if \code{SAVEPLOTS = TRUE} path to saved files.
  
  \code{SAVEPLOTS = FALSE}: if \code{TRUE} writes plots to \code{.pdf} files.
  
  \code{trueValues}:  if passed from \code{mastSim} output then \code{trueValues} will be included in plots.
  
  \code{RMD}:  if set to \code{"pdf"} or \code{"html"} writes output to R markdown file in working directory.
  }
}

\details{

If \code{SAVEPLOTS = TRUE} plots are saved to files in \code{outFolder}. If \code{RMD = "pdf"}, output is written to a R markdown file that can be edited and knitted.  Maps are not included in this option.  Otherwise, plots are rendered to the screen.

More detailed vignettes can be obtained with:

\code{browseVignettes('mastif')}
}

\value{
Currently, there are graphical outputs.
}
\examples{
\donttest{
# simulate data 
seedNames  <- specNames  <- 'acerRubr'
sim <- list(nyr=10, ntree=30, nplot=5, 
            specNames = specNames, seedNames = seedNames)
inputs <- mastSim(sim)
output <- mastif( inputs = inputs, ng = 2000, burnin = 1000 )

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

\code{browseVignettes('mast')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

