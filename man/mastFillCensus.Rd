\name{mastFillCensus}
\alias{mastFillCensus}
\title{Interpolate census data for seed trap years}
\description{
Provides interpolated census data to include years when seed data are available.  This is used when tree sampling is at a lower frequency than seed-trap collections.
}
\usage{

  mastFillCensus(inputs, beforeFirst = 15, afterLast = 15, p = 0, verbose = FALSE)
  
}
\arguments{
  \item{inputs}{\code{list} used as input to \code{\link{mastif}} that must include two \code{character vector}s, (\code{specNames}, \code{seedNames}) and four \code{data.frame}s (\code{treeData}, \code{seedData}, \code{xytree}, \code{xytrap}).  See details.}
  
  \item{beforeFirst}{number of years before a tree is first observed in a census that it should be considered as potentially present.}
  
  \item{afterLast}{number of years after a tree is last observed in a census that it should be considered as potentially present.}
  
    \item{p}{if AR(p) model is used (in \code{yearEffect}), the number of lag terms.}
   
   \item{verbose}{if verbose = TRUE information is provided on filling progress.} 
}

\details{

Masting data sets contain tree census data, held in \code{treeData}, and seed trap data, help in \code{seedData}. Most studies monitor seed rain frequently (e.g., annual), while tree censuses occur at intervals of 2 to 5 years.  \code{mastFillCensus} 'fills in' the tree census so that the annual seed data can be used.  It is made available to the user so covariates can be added, e.g., with \code{\link{mastClimate}}.

\code{mastFillCensus} accepts the \code{list} of \code{inputs} used in \code{\link{mastif}}.  The missing years are inserted for each tree with interpolated diameters.  \code{inputs} is returned with objects updated to include the missing census years and modified slightly for analysis by \code{\link{mastif}}.

The function \code{mastFillCensus} is made accessible to the user, because covariates may be needed for the missing census years.  For example, models often include climate variables that change annually.  The version of \code{treeData} returned by \code{mastFillCensus} can be annotated with additional columns that can then be included in the model, as specified in \code{formulaFec}, \code{formulaRep}, and/or \code{randomEffect$formulaRan}.

\code{beforeFirst} and \code{afterLast} allow the user to control the assumptions about treatment of trees between (and before and after) tree census years. Seed trap data may begin before the first tree census or after the last tree census.  Trees may appear in the middle of the study due to ingrowth.  They may be lost to mortality.  In other words, census data can be left-, right-, and interval-censored.  

For the AR(p) model, values are imputed for \code{p} years before a tree is first observed and \code{p} years after the tree is last observed (\code{\link{mastif}}). These years are inserted by \code{mastFillCensus}, such that they too can then be annotated with covariate data.

More detailed vignettes can be obtained with:

\code{browseVignettes('mastif')}
}

\value{

\item{\code{inputs}}{\code{list} includes all inputs to the model, cleaned for \code{\link{mastif}}, with years inserted to \code{treeData}.}

}


\examples{
\donttest{
d <- "https://github.com/jimclarkatduke/mast/blob/master/liriodendronExample.rData?raw=True"
repmis::source_data(d)

inputs   <- list( specNames = specNames, seedNames = seedNames, 
                  treeData = treeData, seedData = seedData, 
                  xytree = xytree, xytrap = xytrap)
inputs <- mastFillCensus(inputs)

formulaFec <- formulaRep <- as.formula(~ diam)
                  
output <- mastif(inputs = inputs, formulaFec, formulaRep, ng = 1000, 
                 burnin = 400)
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., C. Nunes, and B. Tomasek. 2019. Foodwebs based on unreliable foundations: spatio-temporal masting merged with consumer movement, storage, and diet. Ecological Monographs, e01381.
}
\seealso{
\code{\link{mastSim}} simulates data

A more detailed vignette can be obtained with:

\code{browseVignettes('mastif')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

