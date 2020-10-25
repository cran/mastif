## ----outform, echo=F----------------------------------------------------------
insertPlot <- function(file, caption){
#    outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
#  if(outputFormat == 'latex')
#    paste("![ ", caption, " ](", file, ")",sep="")
}
bigskip <- function(){
#  outputFormat = knitr::opts_knit$get("rmarkdown.pandoc.to")
#  if(outputFormat == 'latex')
#    "\\bigskip"
#  else
    "<br>"
}

## ----getFiles, eval = F, echo=F-----------------------------------------------
#  #library(Rcpp)
#  #library(RcppArmadillo)
#  Rcpp::sourceCpp('../RcppFunctions/cppFns.cpp')
#  source('../RFunctions/mastifFunctions.r')
#  library(RANN)

## ----simSetup0, eval = F------------------------------------------------------
#  seedNames  <- specNames  <- 'acerRubr'
#  sim <- list(nplot=5, nyr=10, ntree=30,  ntrap=40,
#                 specNames = specNames, seedNames = seedNames)

## ----sim0, eval = F-----------------------------------------------------------
#  inputs     <- mastSim(sim)        # simulate dispersal data
#  seedData   <- inputs$seedData     # year, plot, trap, seed counts
#  treeData   <- inputs$treeData     # year, plot, tree data
#  xytree     <- inputs$xytree       # tree locations
#  xytrap     <- inputs$xytrap       # trap locations
#  formulaFec <- inputs$formulaFec   # fecundity model
#  formulaRep <- inputs$formulaRep   # maturation model
#  trueValues <- inputs$trueValues   # true states and parameter values

## ----formulaFec, eval = F-----------------------------------------------------
#  formulaFec

## ----treeData0, eval = F------------------------------------------------------
#  head(treeData)

## ----xytree, eval = F---------------------------------------------------------
#  head(xytree, 5)

## ----treeData1, eval = F------------------------------------------------------
#  head(seedData)

## ----map1a, eval = F----------------------------------------------------------
#  dataTab <- table(treeData$plot,treeData$year)
#  
#  w <- which(dataTab > 0,arr.ind=T) # a plot-year with observations
#  w <- w[sample(nrow(w),1),]
#  
#  mapYears <- as.numeric( colnames(dataTab)[w[2]] )
#  mapPlot  <- rownames(dataTab)[w[1]]
#  inputs$mapPlot    <- mapPlot
#  inputs$mapYears   <- mapYears
#  inputs$treeSymbol <- treeData$diam
#  inputs$SCALEBAR   <- T
#  
#  mastMap(inputs)

## ----map3, eval=F-------------------------------------------------------------
#  inputs$treeSymbol <- trueValues$fec
#  inputs$treeScale  <- 2
#  inputs$trapScale  <- 1
#  
#  mastMap(inputs)

## ----map4, eval=F-------------------------------------------------------------
#  inputs$treeSymbol <- trueValues$repr
#  inputs$treeScale  <- .5
#  mastMap(inputs)

## ----hist0, eval = F----------------------------------------------------------
#  par( mfrow=c(2,1),bty='n', mar=c(4,4,1,1) )
#  seedData <- inputs$seedData
#  seedNames <- inputs$seedNames
#  
#  hist( as.matrix(seedData[,seedNames]) ,nclass=100,
#        xlab = 'seed count', ylab='per trap', main='' )
#  hist( trueValues$fec,nclass=100, xlab = 'seeds produced', ylab = 'per tree',
#        main = '')

## ----mast2, eval = F----------------------------------------------------------
#  output   <- mastif( inputs = inputs, ng = 3000, burnin = 500 )

## ----tabPars0, eval = F-------------------------------------------------------
#  summary( output )

## ----pars, eval = F-----------------------------------------------------------
#  plotPars <- list(trueValues = trueValues)
#  mastPlot(output, plotPars)

## ---- eval=F------------------------------------------------------------------
#  output   <- mastif( inputs = output, ng = 5000, burnin = 3000 )

## ----restart, eval=F----------------------------------------------------------
#  predList <- list( mapMeters = 3, plots = mapPlot, years = mapYears )
#  output   <- mastif( inputs = output, ng = 2000, burnin = 1000,
#                      predList = predList)

## ----mapout, eval = F---------------------------------------------------------
#  output$mapPlot    <- mapPlot
#  output$mapYears   <- mapYears
#  output$treeScale  <- 1.5
#  output$trapScale  <- .8
#  output$PREDICT <- T
#  output$scaleValue <- 20
#  output$plotScale  <- 1
#  output$COLORSCALE <- T
#  output$LEGEND     <- T
#  
#  mastMap(output)

## ----to rmd, eval=F-----------------------------------------------------------
#  plotPars$RMD <- 'pdf'
#  mastPlot(output, plotPars)

## ----simSetup, eval = F-------------------------------------------------------
#  specNames <- c('pinuTaeda','pinuEchi','pinuVirg')
#  seedNames <- c('pinuTaeda','pinuEchi','pinuVirg','pinuUNKN')
#  sim    <- list(nyr=4, ntree=25, nplot=10, ntrap=50, specNames = specNames,
#                    seedNames = seedNames)

## ----sim, eval = F------------------------------------------------------------
#  inputs <- mastSim(sim)        # simulate dispersal data
#  R      <- inputs$trueValues$R # species to seedNames probability matrix
#  round(R, 2)

## ----mast3, eval = F----------------------------------------------------------
#  output <- mastif( inputs = inputs, ng = 2000, burnin = 1000)

## ----tabPars, eval = F--------------------------------------------------------
#  summary( output )

## ----pars0, eval = F----------------------------------------------------------
#  plotPars <- list(trueValues = inputs$trueValues)
#  mastPlot(output, plotPars)

## ----again, eval=F------------------------------------------------------------
#  tab   <- with( inputs$seedData, table(plot, year) )
#  years <- as.numeric( colnames(tab)[tab[1,] > 0] ) # years for 1st plot
#  predList <- list( plots = 'p1', years = years )
#  output <- mastif( inputs = output, ng = 3000, burnin = 1500,
#                    predList = predList)
#  mastPlot(output)

## ---- eval=F------------------------------------------------------------------
#  specNames <- c('pinuTaed', 'pinuEchi')
#  
#  #seeds never differentiated:
#  seedNames <- c('pinuUNKN')
#  
#  #one species sometimes differentiated:
#  seedNames <- c('pinuTaed', 'pinuUNKN')
#  
#  #both species sometimes differentiated:
#  seedNames <- c('pinuTaed', 'pinuEchi', 'pinuUNKN')

## ----map21, eval=F------------------------------------------------------------
#  library(repmis)
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/liriodendronExample.rData?raw=True"
#  repmis::source_data(d)
#  mapList <- list( treeData = treeData, seedData = seedData,
#                   specNames = specNames, seedNames = seedNames,
#                   xytree = xytree, xytrap = xytrap, mapPlot = 'DUKE_EW',
#                   mapYears = 2011:2014, treeSymbol = treeData$diam,
#                   treeScale = .7, trapScale = 1.5, plotScale = 2,
#                   SCALEBAR=T, scaleValue=50)
#  mastMap(mapList)

## ----litu1, eval=F------------------------------------------------------------
#  head(treeData, 3)

## ----litu2, eval=F------------------------------------------------------------
#  head(seedData, 3)

## ----fit, eval=F--------------------------------------------------------------
#  formulaFec <- as.formula( ~ log(diam))    # fecundity model
#  formulaRep <- as.formula( ~ log(diam))    # maturation model
#  
#  inputs   <- list(specNames = specNames, seedNames = seedNames,
#                   treeData = treeData, seedData = seedData, xytree = xytree,
#                   xytrap = xytrap)
#  output <- mastif( formulaFec, formulaRep, inputs = inputs,  ng = 3000,
#                    burnin = 1000 )

## ----more, eval=F-------------------------------------------------------------
#  predList <- list( mapMeters = 10, plots = 'DUKE_EW', years = 2010:2015 )
#  output <- mastif( inputs = output, ng = 3000, burnin = 1000,
#                    predList = predList )

## ----plotmydata1, eval=F------------------------------------------------------
#  mastPlot(output)

## ----outpars, eval=F----------------------------------------------------------
#  summary( output )

## ----fitSum, eval=F-----------------------------------------------------------
#  output$fit

## ----nopred, eval=F-----------------------------------------------------------
#  #group plots in regions for year effects
#  region <- rep('sApps',nrow(treeData))
#  region[ as.character(treeData$plot) == 'DUKE_EW' ] <- 'piedmont'
#  
#  treeData$region <- region
#  
#  formulaFec   <- as.formula(~ diam)
#  formulaRep   <- as.formula( ~ diam )
#  yearEffect   <- list(groups = 'region')
#  randomEffect <- list(randGroups = 'treeID', formulaRan = as.formula( ~ 1 ) )
#  inputs <- list(specNames = specNames, seedNames = seedNames,
#                 treeData = treeData, seedData = seedData,
#                 xytree = xytree, xytrap = xytrap,
#                 priorDist = 28, priorVDist = 15, maxDist = 50, minDist = 15,
#                 minDiam = 25, maxF = 1e+6)
#  output <- mastif(inputs, formulaFec, formulaRep, ng = 2000, burnin = 1000,
#                 randomEffect = randomEffect, yearEffect = yearEffect )
#  mastPlot(output)

## ---- yr, eval=F--------------------------------------------------------------
#  yearEffect   <- list(groups = c('species','region'))   # year effects

## ---- eval=F------------------------------------------------------------------
#  inputs$priorList <- list(minDiam = 15, maxDiam = 60)

## ---- eval=F------------------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/priorParameters.txt?raw=True"
#  download.file(d, destfile="priorParameters.txt")
#  
#  specNames <- c("pinuEchi","pinuRigi","pinuStro","pinuTaed","pinuVirg")
#  seedNames <- c("pinuEchi","pinuRigi","pinuStro","pinuTaed","pinuVirg","pinuUNKN")
#  priorTable <- mastPriors("priorParameters.txt", specNames,
#                           code = 'code4', genus = 'pinus')

## ---- eval=F, echo=F----------------------------------------------------------
#  
#  # THIS BLOCK IS OMITTED
#  
#  load('../combinedSites/pinus.Rdata')
#  
#  pl <- c('CWT_118','CWT_218','DUKE_BW','DUKE_EW','MARS_F','MARS_P')
#  treeData <- treeData[treeData$plot %in% pl,]
#  seedData <- seedData[seedData$plot %in% pl,]
#  xytree <- xytree[xytree$plot %in% pl,]
#  xytrap <- xytrap[xytrap$plot %in% pl,]
#  
#  count <- as.matrix(seedData[,-c(1:5)])
#  count <- count[,colSums(count, na.rm=T) > 0]
#  count <- count[,-1]
#  seedData <- cbind(seedData[,1:5],count)
#  
#  treeData <- treeData[,1:6]
#  
#  seedNames <- colnames(count)
#  specNames <- sort(unique(treeData$species))
#  
#  save(treeData, seedData, xytree, xytrap,
#       specNames, seedNames, file='pinusExample.rData')
#  
#  #load('pinusExample.rData')

## ----priorBeta, eval=F--------------------------------------------------------
#  
#  ###############################
#  
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rdata?raw=True"
#  repmis::source_data(d)

## ---- eval=F------------------------------------------------------------------
#  formulaRep <- as.formula( ~ diam )
#  formulaFec <- as.formula( ~ diam + I(diam^2) )
#  
#  betaPrior <- list(pos = 'diam', neg = 'I(diam^2)')
#  
#  inputs <- list( treeData = treeData, seedData = seedData, xytree = xytree,
#                  xytrap = xytrap, specNames = specNames, seedNames = seedNames,
#                  betaPrior = betaPrior, priorTable = priorTable,
#                  seedTraits = seedTraits)
#  
#  output <- mastif(inputs = inputs, formulaFec, formulaRep,
#                   ng = 500, burnin = 200)
#  
#  # restart
#  output <- mastif(output, formulaFec, formulaRep,
#                   ng = 3000, burnin = 1000)
#  mastPlot(output)

## ----fill, eval=F-------------------------------------------------------------
#  
#  #d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rData?raw=True"
#  #repmis::source_data(d)
#  
#  # randomly remove years for this example:
#  years <- sort(unique(treeData$year))
#  sy    <- sample(years,5)
#  treeData <- treeData[treeData$year %in% sy,]
#  treeData[1:10,]
#  

## ---- removed, eval=F---------------------------------------------------------
#  inputs   <- list( specNames = specNames, seedNames = seedNames,
#                    treeData = treeData, seedData = seedData,
#                    xytree = xytree, xytrap = xytrap, priorTable = priorTable,
#                    seedTraits = seedTraits)
#  inputs <- mastFillCensus(inputs, beforeFirst=10, afterLast=10)
#  inputs$treeData[1:10,]

## ---- treeYr table, eval=F----------------------------------------------------
#  # original data
#  table(treeData$year)
#  
#  # filled census
#  table(inputs$treeData$year)
#  table(seedData$year)

## ---- treeOnly, eval=F--------------------------------------------------------
#  p <- 3
#  inputs   <- mastFillCensus(inputs, p = p)
#  treeData <- inputs$treeData

## ----arr, eval=F--------------------------------------------------------------
#  region = list(sApps = c("CWT_118", "MARS_F", "MARS_P"),
#                piedmont = c("DUKE_BW", "DUKE_EW","DUKE_FACE1","DUKE_FACE5","DUKE_FACE6"),
#                NE = c("HARV_BW","HARV_S"))
#  treeData$region <- 'sApps'
#  for(j in 1:length(region)){
#    wj <- which( as.character(treeData$plot) %in% region[[j]])
#    treeData$region[wj] <- names(region)[j]
#  }
#  inputs$treeData <- treeData
#  yearEffect <- list(groups = c('species','region'), p = p)

## ----def, eval=F--------------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/def.csv?raw=True"
#  download.file(d, destfile="def.csv")

## ----types, eval=F------------------------------------------------------------
#  treeData <- inputs$treeData
#  deficit  <- mastClimate( file = 'def.csv', plots = treeData$plot,
#                           years = treeData$year - 1, months = 6:8,
#                           FUN = 'sum', vname='def')
#  treeData <- cbind(treeData, deficit$x)
#  summary(deficit)

## ----covars, eval=F-----------------------------------------------------------
#  # include min winter temperature
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/tmin.csv?raw=True"
#  download.file(d, destfile="tmin.csv")
#  
#  # minimum winter temperature December through March of previous winter
#  
#  t1 <- mastClimate( file = 'tmin.csv', plots = treeData$plot,
#                     years = treeData$year - 1, months = 12, FUN = 'min',
#                     vname = 'tmin')
#  t2 <- mastClimate( file = 'tmin.csv', plots = treeData$plot,
#                     years = treeData$year, months = 1:3, FUN = 'min',
#                     vname = 'tmin')
#  tmin <- apply( cbind(t1$x[,1], t2$x[,1]), 1, min)
#  treeData$tminDecJanFebMar <- tmin
#  inputs$treeData <- treeData

## ---- runClim, eval=F---------------------------------------------------------
#  
#  formulaRep <- as.formula( ~ diam )
#  formulaFec <- as.formula( ~ diam + defJunJulAugAnom + tminDecJanFebMar )
#  output <- mastif(inputs = inputs, formulaFec, formulaRep,
#                   yearEffect = yearEffect,
#                   ng = 2500, burnin = 1000)

## ----ranEff0, eval=F----------------------------------------------------------
#  randomEffect <- list(randGroups = 'tree',
#                       formulaRan = as.formula( ~ diam ) )

## ----ranEff, eval=F-----------------------------------------------------------
#  formulaFec <- as.formula( ~ diam)    # fecundity model
#  formulaRep <- as.formula( ~ diam)    # maturation model
#  randomEffect <- list(randGroups = 'tree',
#                       formulaRan = as.formula( ~ diam) )
#  output <- mastif( formulaFec, formulaRep, inputs = inputs,
#                    ng = 2000, burnin = 1000,
#                    randomEffect = randomEffect )

## ----ranEff2, eval=F----------------------------------------------------------
#  output <- mastif( inputs = output, ng = 2000, burnin = 1000)
#  mastPlot(output)

## ----fitSum2, eval=F----------------------------------------------------------
#  output$fit

## ----regtab, eval=F-----------------------------------------------------------
#  with(treeData, colSums( table(plot, region)) )

## ----yrpl, eval=F-------------------------------------------------------------
#  yearEffect <- list(groups = 'region')

## ----fit3, eval=F-------------------------------------------------------------
#  inputs$treeData  <- treeData
#  output <- mastif(formulaFec, formulaRep, inputs = inputs,
#                   ng = 2500, burnin = 500,
#                   randomEffect = randomEffect, yearEffect = yearEffect)

## ----moreYR, eval=F-----------------------------------------------------------
#  predList <- list( mapMeters = 10, plots = 'DUKE_BW', years = 1998:2014 )
#  output <- mastif(inputs = output, predList = predList, ng = 3000, burnin = 1000)
#  mastPlot(output)

## ----map2, eval=F-------------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rdata?raw=True"
#  repmis::source_data(d)
#  
#  mapList <- list( treeData = treeData, seedData = seedData,
#                   specNames = specNames, seedNames = seedNames,
#                   xytree = xytree, xytrap = xytrap, mapPlot = 'DUKE_EW',
#                   mapYears = c(2007:2010), treeSymbol = treeData$diam,
#                   treeScale = .6, trapScale=1.4,
#                   plotScale = 1.2, LEGEND=T)
#  mastMap(mapList)

## ----fit0, eval=F-------------------------------------------------------------
#  formulaFec <- as.formula( ~ diam )   # fecundity model
#  formulaRep <- as.formula( ~ diam )            # maturation model
#  
#  yearEffect   <- list(groups = 'species', p = 4)   # AR(4)
#  randomEffect <- list(randGroups = 'tree',
#                       formulaRan = as.formula( ~ diam ) )
#  
#  inputs   <- list( specNames = specNames, seedNames = seedNames,
#                    treeData = treeData, seedData = seedData,
#                    xytree = xytree, xytrap = xytrap, priorDist = 20,
#                    priorVDist = 5, minDist = 15, maxDist = 30,
#                    minDiam = 12, maxDiam = 40,
#                    maxF = 1e+6, seedTraits = seedTraits)
#  output <- mastif(formulaFec, formulaRep, inputs = inputs, ng = 500,
#                 burnin = 100, yearEffect = yearEffect,
#                 randomEffect = randomEffect)

## ----moreAR, eval=F-----------------------------------------------------------
#  output <- mastif(inputs = output, ng = 2000, burnin = 1000)
#  plotPars <- list(MAPS = F)
#  mastPlot(output, plotPars = plotPars)

## ----moreYR2, eval=F----------------------------------------------------------
#  plots <- c('DUKE_EW','CWT_118')
#  years <- 1980:2025
#  predList <- list( mapMeters = 10, plots = plots, years = years )
#  output <- mastif(inputs = output, predList = predList, ng = 3000,
#                   burnin = 1000)

## ----yrPlot, eval=F-----------------------------------------------------------
#  mastPlot( output, plotPars = list(MAPS=F) )

## ----onemap, eval=F-----------------------------------------------------------
#  mapList <- output
#  mapList$mapPlot <- 'DUKE_EW'
#  mapList$mapYears <- c(2011:2012)
#  mapList$PREDICT <- T
#  mapList$treeScale <- 1.2
#  mapList$trapScale <- .8
#  mapList$LEGEND <- T
#  mapList$scaleValue <- 50
#  mapList$plotScale <- 2
#  mapList$COLORSCALE <- T
#  mapList$mfrow <- c(2,1)
#  
#  mastMap( mapList )

## ----onemap1, eval=F----------------------------------------------------------
#  mapList$mapPlot <- 'CWT_118'
#  mapList$mapYears <- 2015
#  mapList$PREDICT <- T
#  mapList$treeScale <- 1.5
#  mapList$trapScale <- .8
#  mapList$LEGEND <- T
#  mapList$scaleValue <- 50
#  mapList$plotScale <- 2
#  mapList$COLORSCALE <- T
#  mapList$mfrow <- c(1,1)
#  mastMap( mapList )

## ----outpars0, eval=F---------------------------------------------------------
#  summary( output )

