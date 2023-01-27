\name{mastPriors}
\alias{mastPriors}
\title{Obtain prior parameter values for mastif from file}
\description{
Prior parameter values may be saved in a file by species or by genus. \code{mastPriors} looks for a species-level prior first.  If not found, it can substutitute a genus-level prior.
}
\usage{
  mastPriors(file, specNames, code, genus = 'NULL')
}
\arguments{
  \item{file}{\code{character} giving the path and file name where parameters are stored (see Details).}
  
  \item{specNames}{\code{character vector} if names in \code{treeData$species}.}
  
  \item{code}{\code{character} name of column holding \code{specNames} in parameter \code{file}.}
  
  \item{genus}{\code{character} name of genus held in "genus" column of parameter \code{file}.}
}

\details{
The \code{file} includes rows with genera, given in column "genus", or "species".  Species rows also have an entry for genus, with the species code given in the column named \code{code}.  Additional columns are names of prior parameters, including:

\code{priorDist}: mean parameter for dispersal kernel (m), related to kernel parameter \code{u} as \code{d <- pi*sqrt(u)/2}.  The estimated values for these parameters are found in \code{output$parameters$upars} and \code{output$parameters$dpars}, where \code{output} is an object fitted by \code{mastif}.

\code{minDist}: the lower bound for the mean parameter \code{d} of the dispersal kernel (m).  

\code{maxDist}: the upper bound for the mean parameter \code{d} of the dispersal kernel (m).  

\code{priorVDist}: variance on the mean parameter for dispersal kernel (m^2).  For large values, the prior distribution of \code{d} (and by variable change, \code{u}) becomes \code{dunif(d, minDist, maxDist)}.

\code{minDiam}: below this diameter trees of unknown status are assumed immature (cm). 

\code{maxDiam}: above this diameter trees of unknown status are assumed mature (cm).  

\code{maxFec}: maximum seeds per tree per year

More detailed vignettes can be obtained with: \code{browseVignettes('mastif')}
}

\value{
A \code{data.frame} with a row for each \code{specNames} and columns for prior parameter values.  Where \code{file} contains species-level parameter values, they will be used.  If a separate row in \code{file} holds genus-level parameters, with the entry for \code{code == 'NA'}, then genus-level parameters will be substituted.  In other words, these genus rows are default values.
}

\examples{
\donttest{
d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rdata?raw=True"
repmis::source_data(d)

# prior parameter values
pfile <- tempfile(fileext = '.txt')
d <- "https://github.com/jimclarkatduke/mast/blob/master/priorParameters.txt?raw=True"
download.file(d, destfile = pfile)

specNames <- c("pinuEchi","pinuRigi","pinuStro","pinuTaed","pinuVirg")
seedNames <- c(specNames, "pinuUNKN")

priorTable <- mastPriors(file = pfile, specNames, 
                         code = 'code4', genus = 'pinus')

inputs   <- list( specNames = specNames, seedNames = seedNames, 
                  treeData = treeData, seedData = seedData, 
                  xytree = xytree, xytrap = xytrap, 
                  priorTable = priorTable, seedTraits = seedTraits)

formulaRep <- as.formula( ~ diam )
formulaFec <- as.formula( ~ diam )        
  
output <- mastif(inputs = inputs, formulaFec, formulaRep,
                 ng = 1000, burnin = 400)
}
}
\author{
James S Clark, \email{jimclark@duke.edu}
}
\references{
Clark, J.S., C. Nunes, and B. Tomasek. 2019. Foodwebs based on unreliable foundations: spatio-temporal masting merged with consumer movement, storage, and diet. Ecological Monographs, e01381.
}
\seealso{

\code{\link{mastFillCensus}} to fill tree census

\code{\link{mastif}} for analysis

A more detailed vignette is can be obtained with:

\code{browseVignettes('mastif')}

website 'http://sites.nicholas.duke.edu/clarklab/code/'.
}

