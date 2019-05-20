\name{mastMap}
\alias{mastMap}
\title{Map data and predictions for mastif model}
\description{
Maps dispersal data (trees and seed traps) with predictions.
}
\usage{
                        
  mastMap(mapList)
}
\arguments{
   \item{mapList}{\code{list} that includes \code{data.frame}s described in the \code{mastif} help page: \code{treeData}, \code{seedData}, \code{xytree}, \code{xytrap}. Alternatively, \code{mapList} can be an object of \code{class mastif} simulated by \code{mastSim} or fitted by \code{mast}.  If a \code{class mastif} object is used, it must have appended \code{mapPlot} (name of plot to map) and \code{mapYears} (years to map).  These and additional objects in \code{mapList} are summarized here:

   
      \code{mapPlot}: \code{character} for the plot name, which appears in \code{treeData$plot} and \code{seedData$plot}.

   \code{mapYears}: \code{numeric vector} containing years to plot, which appear in \code{treeData$year} and \code{seedData$year}.
  
   \code{treeSymbol}: \code{character} indicating a column name in \code{treeData} holding numeric values that scale the sizes of trees in the map, e.g., \code{diam} for diameter.
   
   \code{PREDICT}: if \code{TRUE}, then predicted seed will be mapped.  Must be predicted in \code{mastif}.
  
  \code{xlim}: \code{vector} with upper and lower plot limits.
  
  \code{ylim}: \code{vector} with upper and lower plot limits.
  
  \code{MAPTRAPS}: include seed traps in map.
  
  \code{MAPTREES}: include trees in map.
    
  \code{seedMax}: \code{numeric} value can be used to set maximum seed symbol/prediction contour.
  
    \code{fecMax}: \code{numeric} value can be used to set maximum fecundity symbol size.
  
    \code{mfrow}: \code{numeric vector} of length-2 to set arrangement of maps on plot, as in \code{pars()}.
    
  \code{LEGEND}: add \code{legend} with species names to map.
  
  \code{plotScale}: some level of control on size of plot is possible with values > 1.
  
  \code{treeScale}: relative scale for tree symbols, default value = 1.
    
  \code{trapScale}: relative scale for trap symbols, default value = 1.
  
  \code{SCALEBAR}: add a scale bar with distance to bottom of map.
  
  \code{scaleValue}: \code{numeric} value in meters used if \code{SCALEBAR = T}.
  
  \code{mfrow}: \code{numeric vector} of length-2, giving numbers of rows and columns for plots, see \code{help(pars)}. 
  
  \code{COLORSCALE}: if \code{TRUE}, then a color scale is added for predicted seed.
}
}

\details{

Generates of map of seed traps and trees, with symbols scaled to the sizes relative to seed counts in \code{sdata$seedNames} and \code{treeSymbol}.  Sizes are adjusted with \code{scaleTree} and \code{scaleTrap}.

If \code{PREDICT = T}, then predictions come in the object fitted in \code{mastif} with \code{predictList} used to specify prediction plots and years.  See the help page for \code{mastif}. 

More detailed vignettes can be obtained with:

\code{browseVignettes('mastif')}
}

\value{

Only graphical outputs.
}
\examples{
\donttest{
# simulate data (see \link{\code{mastSim}})
seedNames  <- specNames  <- 'acerRubr'
sim    <- list(nyr=10, ntree=30, nplot=5,
          specNames = specNames, seedNames = seedNames)
inputs <- mastSim(sim)

inputs$mapPlot <- 'p1'
inputs$ mapYears = inputs$years[1]

mastMap( inputs )

# for Pinus
d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rdata?raw=True"
repmis::source_data(d)

specNames <- c("pinuEchi","pinuRigi","pinuStro","pinuTaed","pinuVirg")
seedNames <- c(specNames, "pinuUNKN")

mapList <- list( treeData = treeData, seedData = seedData, 
                 specNames = specNames, seedNames = seedNames, 
                 xytree = xytree, xytrap = xytrap, mapPlot = 'DUKE_BW', 
                 mapYears = c(2004:2007), treeScale = .5, trapScale=1.2, 
                 plotScale = 1.2, LEGEND=T)
mastMap(mapList) 
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

