\name{mastClimate}
\alias{mastClimate}
\title{Covariates for mast data}
\description{
Annotates \code{treeData} for \code{mastif} to include covariates.
}
\usage{
  mastClimate( file, plots, years, months = 1:12, FUN = 'mean', 
                         vname = '', normYr = c( 1990:2020 ), lastYear = 2021 )
}
\arguments{
  \item{file}{\code{character} giving the path and file name where covariate is stored in \code{mastif} covariate format (see Details).}
  
  \item{plots}{\code{character vector} from \code{treeData$plot}.}
  
  \item{years}{\code{numeric vector} from \code{treeData$year}.}
  
  \item{months}{\code{numeric vector} of months to include in annual values.}
  
  \item{FUN}{\code{character} specifies how monthly covariate values are combined into annual values.}
  
  \item{vname}{name to use for a variable in the model that comes from \code{file}.}
  
  \item{normYr}{years for climate norm for calculating anomalies.}
  
  \item{lastYear}{last data year to include.}
}

\details{
The version of \code{treeData} used in \code{mastif} can have additional tree years included when there are seed trap years that were not censused or when AR(p) effects extend observations to impute the p years before and after a tree was observed.  The \code{function mastFillCensus} makes this version of \code{treeData} available to the user.  The \code{function mastClimate} provides a quick way to add plot-year covariates to \code{treeData}. 

A covariate like minimum monthly temperature is stored in a \code{plot} by \code{year_month} format, where \code{rownames} of \code{file} are plot names matching \code{treeData$plot}, and \code{colnames} of \code{file} could be \code{2012_1, 2012_2, ...} for the 12 months in the year.  The \code{numeric vector months} holds the months to be included in the annual values, e.g., \code{c(3, 4)} for minimum winter temperatures during the period from March through April.  To find the minimum for this period, set \code{FUN} to \code{'min'}.

More detailed vignettes can be obtained with: \code{browseVignettes('mastif')}
}

\value{
A \code{numeric vector} equal in length to the number of rows in \code{treeData} that can be added as a \code{column} and included in \code{formulaFec}.
}

\examples{
\donttest{
d <- "https://github.com/jimclarkatduke/mast/blob/master/liriodendronExample.rData?raw=True"
repmis::source_data(d)

inputs   <- list( specNames = specNames, seedNames = seedNames, 
                  treeData = treeData, seedData = seedData, 
                  xytree = xytree, xytrap = xytrap)
                  
# interpolate census, add years for AR(p) model
inputs   <- mastFillCensus(inputs, p = 3)
treeData <- inputs$treeData               #now includes additional years

# include minimum spring temperature of previous year
cfile <- tempfile(fileext = '.csv')
d <- "https://github.com/jimclarkatduke/mast/blob/master/tmin.csv?raw=True"
download.file(d, destfile=cfile)

tyears <- treeData$year - 1
tplots <- treeData$plot
  
tmp <- mastClimate( file = cfile, plots = tplots, 
                    years = tyears, months = 1:4, FUN = 'min')
treeData$tminSprAnomaly <- tmp$x[,3]                             
inputs$treeData <- treeData
  
formulaRep <- as.formula( ~ diam )
formulaFec <- as.formula( ~ diam + tminSprAnomaly )        
  
inputs$yearEffect <- list(groups ='species', p = 3)  # AR(3) model, species are lag groups
  
output <- mastif(inputs = inputs, formulaFec, formulaRep, ng = 1000, burnin = 400)
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

