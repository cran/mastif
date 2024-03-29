## ----outform, echo=F----------------------------------------------------------
insertPlot <- function( file, caption ){
#    outputFormat = knitr::opts_knit$get( "rmarkdown.pandoc.to" )
#  if( outputFormat == 'latex' )
#    paste( "![ ", caption, " ]( ", file, " )", sep="" )
}
bigskip <- function( ){
#  outputFormat = knitr::opts_knit$get( "rmarkdown.pandoc.to" )
#  if( outputFormat == 'latex' )
#    "\\bigskip"
#  else
    "<br>"
}

## ----getFiles, eval = F, echo=F-----------------------------------------------
#  Rcpp::sourceCpp( '../RcppFunctions/cppFns.cpp' )
#  source( '../RFunctions/mastifFunctions.r' )
#  library( RANN )
#  library( robustbase )

## ----eval = F-----------------------------------------------------------------
#  library( mastif )

## ----simSetup0, eval = F------------------------------------------------------
#  seedNames  <- specNames  <- 'acerRubr'
#  sim <- list( nplot=5, nyr=10, ntree=30,  ntrap=40, specNames = specNames,
#               seedNames = seedNames )

## ----sim0, eval = F-----------------------------------------------------------
#  inputs     <- mastSim( sim )        # simulate dispersal data
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
#  head( treeData )

## ----xytree, eval = F---------------------------------------------------------
#  head( xytree, 5 )

## ----treeData1, eval = F------------------------------------------------------
#  head( seedData )

## ----map1a, eval = F----------------------------------------------------------
#  dataTab <- table( treeData$plot, treeData$year )
#  
#  w <- which( dataTab > 0, arr.ind=T ) # a plot-year with observations
#  w <- w[sample( nrow( w ), 1 ), ]
#  
#  mapYears <- as.numeric( colnames( dataTab )[w[2]] )
#  mapPlot  <- rownames( dataTab )[w[1]]
#  inputs$mapPlot    <- mapPlot
#  inputs$mapYears   <- mapYears
#  inputs$treeSymbol <- treeData$diam
#  inputs$SCALEBAR   <- T
#  
#  mastMap( inputs )

## ----map3, eval=F-------------------------------------------------------------
#  inputs$treeSymbol <- trueValues$fec
#  inputs$treeScale  <- 1.5
#  inputs$trapScale  <- 1
#  mastMap( inputs )

## ----map4, eval=F-------------------------------------------------------------
#  inputs$treeSymbol <- trueValues$repr
#  inputs$treeScale  <- .5
#  mastMap( inputs )

## ----hist0, eval = F----------------------------------------------------------
#  par( mfrow=c( 1, 2 ), bty='n', mar=c( 4, 4, 1, 1 ) )
#  seedData  <- inputs$seedData
#  seedNames <- inputs$seedNames
#  
#  hist( as.matrix( seedData[, seedNames] ) , nclass=100,
#        xlab = 'seed count', ylab='per trap', main='' )
#  hist( trueValues$fec, nclass=100, xlab = 'seeds produced', ylab = 'per tree', main = '' )

## ----mast2, eval = F----------------------------------------------------------
#  output   <- mastif( inputs = inputs, ng = 4000, burnin = 500 )

## ----tabPars0, eval = F-------------------------------------------------------
#  summary( output )

## ----pars, eval = F-----------------------------------------------------------
#  plotPars <- list( trueValues = trueValues )
#  mastPlot( output, plotPars )

## ----eval=F-------------------------------------------------------------------
#  output   <- mastif( inputs = output, ng = 5000, burnin = 3000 )

## ----restart, eval=F----------------------------------------------------------
#  output$predList <- list( mapMeters = 3, plots = mapPlot, years = mapYears )
#  output   <- mastif( inputs = output, ng = 2000, burnin = 1000 )

## ----mapout, eval = F---------------------------------------------------------
#  output$mapPlot    <- mapPlot
#  output$mapYears   <- mapYears
#  output$treeScale  <- 1.2
#  output$trapScale  <- .7
#  output$PREDICT    <- T
#  output$scaleValue <- 10
#  output$plotScale  <- 1
#  output$COLORSCALE <- T
#  output$LEGEND     <- T
#  
#  mastMap( output )

## ----to rmd, eval=F-----------------------------------------------------------
#  plotPars$RMD <- 'pdf'
#  mastPlot( output, plotPars )

## ----simSetup, eval = F-------------------------------------------------------
#  specNames <- c( 'pinuTaeda', 'pinuEchi', 'pinuVirg' )
#  seedNames <- c( 'pinuTaeda', 'pinuEchi', 'pinuVirg', 'pinuUNKN' )
#  sim    <- list( nyr=4, ntree=25, nplot=10, ntrap=50, specNames = specNames,
#                    seedNames = seedNames )

## ----sim, eval = F------------------------------------------------------------
#  inputs <- mastSim( sim )        # simulate dispersal data
#  R      <- inputs$trueValues$R   # species to seedNames probability matrix
#  round( R, 2 )

## ----mast3, eval = F----------------------------------------------------------
#  output <- mastif( inputs = inputs, ng = 2000, burnin = 1000 )

## ----tabPars, eval = F--------------------------------------------------------
#  summary( output )

## ----pars0, eval = F----------------------------------------------------------
#  plotPars <- list( trueValues = inputs$trueValues, RMAT = TRUE )
#  mastPlot( output, plotPars )

## ----again, eval=F------------------------------------------------------------
#  tab   <- with( inputs$seedData, table( plot, year ) )
#  mapPlot <- 'p1'
#  mapYears <- as.numeric( colnames( tab )[tab[1, ] > 0] )   # years for 1st plot
#  output$predList <- list( plots = mapPlot, years = mapYears )
#  output <- mastif( inputs = output, ng = 3000, burnin = 1500 )
#  mastPlot( output, plotPars = list( trueValues = inputs$trueValues, MAPS = TRUE )  )

## ----eval = F-----------------------------------------------------------------
#  output$PREDICT  <- T
#  output$LEGEND   <- T
#  output$mapPlot  <- mapPlot
#  output$mapYears <- mapYears
#  mastMap( output )

## ----eval=F-------------------------------------------------------------------
#  specNames <- c( 'pinuTaed', 'pinuEchi' )
#  
#  #seeds not differentiated:
#  seedNames <- c( 'pinuUNKN' )
#  
#  #one species sometimes differentiated:
#  seedNames <- c( 'pinuTaed', 'pinuUNKN' )
#  
#  #both species sometimes differentiated:
#  seedNames <- c( 'pinuTaed', 'pinuEchi', 'pinuUNKN' )

## ----map21, eval=F------------------------------------------------------------
#  library( repmis )
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/liriodendronExample.rData?raw=True"
#  repmis::source_data( d )
#  mapList <- list( treeData = treeData, seedData = seedData,
#                   specNames = specNames, seedNames = seedNames,
#                   xytree = xytree, xytrap = xytrap, mapPlot = 'DUKE_BW',
#                   mapYears = 2011:2014, treeSymbol = treeData$diam,
#                   treeScale = .7, trapScale = 1.5, plotScale = 2,
#                   SCALEBAR=T, scaleValue=50 )
#  mastMap( mapList )

## ----litu1, eval=F------------------------------------------------------------
#  head( treeData, 3 )

## ----litu2, eval=F------------------------------------------------------------
#  head( seedData, 3 )

## ----fit, eval=F--------------------------------------------------------------
#  formulaFec <- as.formula( ~ diam )    # fecundity model
#  formulaRep <- as.formula( ~ diam )    # maturation model
#  inputs   <- list( specNames = specNames, seedNames = seedNames,
#                   treeData = treeData, seedData = seedData, xytree = xytree,
#                   xytrap = xytrap )
#  output <- mastif( inputs, formulaFec, formulaRep, ng = 3000, burnin = 1000 )

## ----more, eval=F-------------------------------------------------------------
#  output$predList <- list( mapMeters = 10, plots = 'DUKE_EW', years = 2010:2015 )
#  output <- mastif( inputs = output, ng = 4000, burnin = 1000 )

## ----plotmydata1, eval=F------------------------------------------------------
#  mastPlot( output )

## ----outpars, eval=F----------------------------------------------------------
#  summary( output )

## ----fitSum, eval=F-----------------------------------------------------------
#  output$fit

## ----nopred, eval=F-----------------------------------------------------------
#  #group plots in regions for year effects
#  region <- rep( 'sApps', nrow( treeData ) )
#  region[ as.character( treeData$plot ) == 'DUKE_EW' ] <- 'piedmont'
#  
#  treeData$region <- region
#  
#  formulaFec   <- as.formula( ~ diam )
#  formulaRep   <- as.formula( ~ diam )
#  yearEffect   <- list( groups = 'region' )
#  randomEffect <- list( randGroups = 'treeID', formulaRan = as.formula( ~ 1 ) )
#  inputs <- list( specNames = specNames, seedNames = seedNames,
#                 treeData = treeData, seedData = seedData,
#                 xytree = xytree, xytrap = xytrap,
#                 priorDist = 28, priorVDist = 15, maxDist = 50, minDist = 15,
#                 minDiam = 25, maxF = 1e+6,
#                 randomEffect = randomEffect, yearEffect = yearEffect )
#  output <- mastif( inputs, formulaFec, formulaRep, ng = 2000, burnin = 1000 )
#  mastPlot( output )

## ----yr, eval=F---------------------------------------------------------------
#  yearEffect   <- list( groups = c( 'species', 'region' ) )   # year effects

## ----eval=F-------------------------------------------------------------------
#  inputs$priorList <- list( minDiam = 15, maxDiam = 60 )

## ----eval=F-------------------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/priorParametersPinus.csv?raw=True"
#  download.file( d, destfile="priorParametersPinus.csv" )
#  
#  specNames <- c( "pinuEchi", "pinuRigi", "pinuStro", "pinuTaed", "pinuVirg" )
#  seedNames <- c( "pinuEchi", "pinuRigi", "pinuStro", "pinuTaed", "pinuVirg", "pinuUNKN" )
#  priorTable <- mastPriors( "priorParametersPinus.csv", specNames,
#                           code = 'code4', genus = 'pinus' )

## ----eval=F, echo=F-----------------------------------------------------------
#  
#  # THIS BLOCK IS OMITTED
#  
#  load( '../combinedSites/pinus.Rdata' )
#  
#  pl <- c( 'CWT_118', 'CWT_218', 'DUKE_BW', 'DUKE_EW', 'MARS_F', 'MARS_P' )
#  treeData <- treeData[treeData$plot %in% pl, ]
#  seedData <- seedData[seedData$plot %in% pl, ]
#  xytree <- xytree[xytree$plot %in% pl, ]
#  xytrap <- xytrap[xytrap$plot %in% pl, ]
#  
#  count <- as.matrix( seedData[, -c( 1:5 )] )
#  count <- count[, colSums( count, na.rm=T ) > 0]
#  count <- count[, -1]
#  seedData <- cbind( seedData[, 1:5], count )
#  
#  treeData <- treeData[, 1:6]
#  
#  seedNames <- colnames( count )
#  specNames <- sort( unique( treeData$species ) )
#  
#  save( treeData, seedData, xytree, xytrap,
#       specNames, seedNames, file='pinusExample.rData' )
#  
#  #load( 'pinusExample.rData' )

## ----priorBeta, eval=F--------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rdata?raw=True"
#  repmis::source_data( d )

## ----eval=F-------------------------------------------------------------------
#  formulaRep <- as.formula( ~ diam )
#  formulaFec <- as.formula( ~ diam + I( diam^2 ) )
#  betaPrior  <- list( pos = 'diam', neg = 'I( diam^2 )' )
#  
#  inputs <- list( treeData = treeData, seedData = seedData, xytree = xytree,
#                  xytrap = xytrap, specNames = specNames, seedNames = seedNames,
#                  betaPrior = betaPrior, priorTable = priorTable,
#                  seedTraits = seedTraits )
#  output <- mastif( inputs = inputs, formulaFec, formulaRep,
#                   ng = 500, burnin = 200 )
#  
#  # restart
#  output <- mastif( output, formulaFec, formulaRep, ng = 5000, burnin = 1000 )
#  mastPlot( output )

## ----fill, eval=F-------------------------------------------------------------
#  # randomly remove years for this example:
#  years <- sort( unique( treeData$year ) )
#  sy    <- sample( years, 5 )
#  treeData <- treeData[treeData$year %in% sy, ]
#  treeData[1:10, ]

## ----removed, eval=F----------------------------------------------------------
#  inputs   <- list( specNames = specNames, seedNames = seedNames,
#                    treeData = treeData, seedData = seedData,
#                    xytree = xytree, xytrap = xytrap, priorTable = priorTable,
#                    seedTraits = seedTraits )
#  inputs <- mastFillCensus( inputs, beforeFirst=10, afterLast=10 )
#  inputs$treeData[1:10, ]

## ----treeYr table, eval=F-----------------------------------------------------
#  # original data
#  table( treeData$year )
#  
#  # filled census
#  table( inputs$treeData$year )
#  table( seedData$year )

## ----treeOnly, eval=F---------------------------------------------------------
#  p <- 3
#  inputs   <- mastFillCensus( inputs, p = p )
#  treeData <- inputs$treeData

## ----arr, eval=F--------------------------------------------------------------
#  region <- c( 'CWT', 'DUKE', 'HARV' )
#  treeData$region <- 'SCBI'
#  for( j in 1:length( region ) ){
#    wj <- which( startsWith( treeData$plot, region[j] ) )
#    treeData$region[wj] <- region[j]
#  }
#  inputs$treeData <- treeData
#  yearEffect <- list( groups = c( 'species', 'region' ), p = p )

## ----def, eval=F--------------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/def.csv?raw=True"
#  download.file( d, destfile="def.csv" )

## ----types, eval=F------------------------------------------------------------
#  treeData <- inputs$treeData
#  deficit  <- mastClimate( file = 'def.csv', plots = treeData$plot,
#                           years = treeData$year - 1, months = 6:8,
#                           FUN = 'sum', vname='def' )
#  treeData <- cbind( treeData, deficit$x )
#  summary( deficit )

## ----covars, eval=F-----------------------------------------------------------
#  # include min winter temperature
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/tmin.csv?raw=True"
#  download.file( d, destfile="tmin.csv" )
#  
#  # minimum winter temperature December through March of previous winter
#  
#  t1 <- mastClimate( file = 'tmin.csv', plots = treeData$plot,
#                     years = treeData$year - 1, months = 12, FUN = 'min',
#                     vname = 'tmin' )
#  t2 <- mastClimate( file = 'tmin.csv', plots = treeData$plot,
#                     years = treeData$year, months = 1:3, FUN = 'min',
#                     vname = 'tmin' )
#  tmin <- apply( cbind( t1$x[, 1], t2$x[, 1] ), 1, min )
#  treeData$tminDecJanFebMar <- tmin
#  inputs$treeData <- treeData

## ----runClim, eval=F----------------------------------------------------------
#  formulaRep <- as.formula( ~ diam )
#  formulaFec <- as.formula( ~ diam + defJunJulAugAnom + tminDecJanFebMar )
#  inputs$yearEffect <- yearEffect
#  output <- mastif( inputs = inputs, formulaFec, formulaRep, ng = 2500, burnin = 1000 )

## ----ranEff, eval=F-----------------------------------------------------------
#  formulaFec <- as.formula( ~ diam )    # fecundity model
#  formulaRep <- as.formula( ~ diam )    # maturation model
#  inputs$randomEffect <- list( randGroups = 'tree', formulaRan = as.formula( ~ 1 ) )
#  output <- mastif( inputs = inputs, formulaFec, formulaRep, ng = 2000, burnin = 1000 )

## ----ranEff2, eval=F----------------------------------------------------------
#  output <- mastif( inputs = output, ng = 4000, burnin = 1000 )
#  mastPlot( output )

## ----fitSum2, eval=F----------------------------------------------------------
#  output$fit

## ----regtab, eval=F-----------------------------------------------------------
#  with( treeData, colSums( table( plot, region ) ) )

## ----yrpl, eval=F-------------------------------------------------------------
#  yearEffect <- list( groups = c('species', 'region' ) )

## ----fit3, eval=F-------------------------------------------------------------
#  inputs$treeData      <- treeData
#  inputs$randomEffect  <- randomEffect
#  inputs$yearEffect    <- yearEffect
#  output <- mastif( inputs = inputs, formulaFec, formulaRep, ng = 2500, burnin = 500 )

## ----moreYR, eval=F-----------------------------------------------------------
#  predList <- list( mapMeters = 10, plots = 'DUKE_BW', years = 1998:2014 )
#  output$predList <- predList
#  output <- mastif( inputs = output, ng = 3000, burnin = 1000 )
#  mastPlot( output )

## ----map2, eval=F-------------------------------------------------------------
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/pinusExample.rdata?raw=True"
#  repmis::source_data( d )
#  
#  mapList <- list( treeData = treeData, seedData = seedData,
#                   specNames = specNames, seedNames = seedNames,
#                   xytree = xytree, xytrap = xytrap, mapPlot = 'DUKE_EW',
#                   mapYears = c( 2007:2010 ), treeSymbol = treeData$diam,
#                   treeScale = .6, trapScale=1.4,
#                   plotScale = 1.2, LEGEND=T )
#  mastMap( mapList )

## ----fit0, eval=F-------------------------------------------------------------
#  formulaFec <- as.formula( ~ diam )   # fecundity model
#  formulaRep <- as.formula( ~ diam )   # maturation model
#  
#  yearEffect   <- list( groups = 'species', p = 4 )   # AR( 4 )
#  randomEffect <- list( randGroups = 'tree',
#                       formulaRan = as.formula( ~ 1 ) )
#  
#  inputs   <- list( specNames = specNames, seedNames = seedNames,
#                    treeData = treeData, seedData = seedData,
#                    yearEffect = yearEffect,
#                    randomEffect = randomEffect,
#                    xytree = xytree, xytrap = xytrap, priorDist = 20,
#                    priorVDist = 5, minDist = 15, maxDist = 30,
#                    minDiam = 12, maxDiam = 40,
#                    maxF = 1e+6, seedTraits = seedTraits )
#  output <- mastif( inputs = inputs, formulaFec, formulaRep, ng = 500, burnin = 100 )

## ----moreAR, eval=F-----------------------------------------------------------
#  output <- mastif( inputs = output, ng = 3000, burnin = 1000 )
#  mastPlot( output, plotPars = plotPars )

## ----moreYR2, eval=F----------------------------------------------------------
#  plots <- c( 'DUKE_EW', 'CWT_118' )
#  years <- 1980:2025
#  output$predList <- list( mapMeters = 10, plots = plots, years = years )
#  output <- mastif( inputs = output, ng = 3000, burnin = 1000 )

## ----yrPlot, eval=F-----------------------------------------------------------
#  mastPlot( output, )

## ----onemap, eval=F-----------------------------------------------------------
#  mapList <- output
#  mapList$mapPlot <- 'DUKE_EW'
#  mapList$mapYears <- c( 2011:2012 )
#  mapList$PREDICT <- T
#  mapList$treeScale <- .5
#  mapList$trapScale <- .8
#  mapList$LEGEND <- T
#  mapList$scaleValue <- 50
#  mapList$plotScale  <- 2
#  mapList$COLORSCALE <- T
#  mapList$mfrow <- c( 2, 1 )
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
#  mapList$mfrow <- c( 1, 1 )
#  mastMap( mapList )

## ----outpars0, eval=F---------------------------------------------------------
#  summary( output )

## ----eval=F, echo = F---------------------------------------------------------
#  
#  # abies output west
#  load('prediction.rdata')
#  
#  fecPred <- prediction$fecPred
#  fecPred <- fecPred[ fecPred$matrEst > .5, ]
#  
#  ptab  <- table( fecPred$plot, fecPred$species )
#  ww    <- which( ptab > 1000, arr.ind = T )
#  keep  <- paste( rownames(ptab)[ww[,1]], colnames(ptab)[ww[,2]] )
#  
#  fecPred$plotSpec  <- paste(fecPred$plot, fecPred$species)
#  fecPred <- fecPred[ fecPred$plotSpec %in% keep, ]
#  specs   <- sort( unique( fecPred$species ) )
#  
#  yseq    <- seq( 0, 10, length = 100 )
#  intVal  <- matrix( 0, length(specs), 100 )
#  rownames( intVal ) <- specs
#  
#  par( mfrow = c(1, 3), bty = 'n', mar = c(2,4,1,1), omi = c(.5,.1,.1,.1) )
#  plot( NA, xlim = c(2, 10), ylim = c(.01, 1), xlab = 'Year', ylab = 'Density/yr',
#        log = 'xy')
#  
#  for( i in 1:length(keep) ){
#  
#    wi  <- which( fecPred$plotSpec == keep[i]  )
#    ci  <- fecPred$species[wi[1]]
#    col <- match( ci, specs )
#  
#    tmp <- mastVolatility( treeID = fecPred$treeID[wi], year = fecPred$year[wi],
#                           fec = fecPred$fecEstMu[wi] )
#    if( is.null(tmp) )next
#  
#    intVal[ col, ] <- intVal[ col, ] + dnorm( yseq, tmp$stats['Period', 1], tmp$stats['Period', 2] )
#  
#    dens <- tmp$statsDensity
#    lines( dens[ 'Period', ], dens[ 'Mean', ], lwd = 2, col = getColor( col, .4) )
#    lines( dens[ 'Period', ], dens[ 'Mean', ] - dens[ 'SD', ], lty = 2, col = getColor( col, .4) )
#    lines( dens[ 'Period', ], dens[ 'Mean', ] + dens[ 'SD', ], lty = 2, col = getColor( col, .4) )
#  }
#  title( 'a) Plot-species groups' )
#  
#  plot( NA, xlim = c(2, 6), ylim = c(0, .12), xlab = '', ylab = 'Density' )
#  for( i in 1:length(specs) ){
#    polygon( yseq, intVal[i,]/sum(intVal[i,]), border = i, col = getColor(i, .4) )
#  }
#  legend( 'topright', specs, text.col = c(1:length(specs)), bty = 'n', cex = .8 )
#  title( 'b) Period estimates' )
#  
#  # using year effects
#  
#  load('parameters.rdata')
#  
#  # if there are year effects in the model:
#  
#  betaYrRand <- parameters$betaYrRand
#  betaYrSE <- parameters$betaYrRandSE
#  
#  spec <- strsplit( rownames( betaYrRand ), '_' )
#  spec <- sapply( spec, function(x) x[2] )
#  specs <- sort( unique(spec) )
#  
#  mastMatrix <- matrix( 0, nrow(betaYrRand), 5 )
#  rownames(mastMatrix) <- rownames(betaYrRand)
#  colnames(mastMatrix) <- c( 'nyr', 'Variance', 'Volatility', 'Period Est', 'Period SD' )
#  
#  
#  #par( mfrow = c(1, 2), bty = 'n' )
#  plot( NA, xlim = c(2, 10), ylim = c(.002, .4), xlab = '', ylab = 'Density/yr',
#        log = 'xy')
#  
#  for(i in 1:nrow(betaYrRand)){
#  
#    wc <- which( betaYrRand[i,] != 0 )
#    if( length(wc) < 6 )next
#  
#    s <- spectralDensity( betaYrRand[i,wc] )
#    if( !is.matrix( s$spect ) )next
#  
#    mastMatrix[i, ] <- c( length(wc), s$totVar, s$volatility, s$periodMu, s$periodSd )
#  
#    period <- 1/s$spec[, 'frequency' ]
#    dens   <- s$spec[, 'spectralDensity' ]/length(wc)  # series vary in length
#  
#    col <- match( spec[i], specs )
#  
#    lines( period, dens, lwd = 2, col = getColor( col, .4) )
#  }
#  title( 'c) Year effects' )
#  mtext( 'Period (yr)', 1, line = 1, outer = T )
#  
#  keepRows   <- which(  is.finite(mastMatrix[,'Variance']) & mastMatrix[,'Variance'] != 0 )
#  keepCols   <- which( colSums( frequency, na.rm=T ) > 0 )
#  mastMatrix <- mastMatrix[ keepRows, ]
#  
#  save( fecPred, betaYrRand, file = 'outputAbies.rdata' )

## ----eval = F-----------------------------------------------------------------
#  getColor <- function( col, trans ){                  # transparent colors
#    tmp <- col2rgb( col )
#    rgb( tmp[ 1, ], tmp[ 2, ], tmp[ 3, ], maxColorValue = 255,
#         alpha = 255*trans, names = paste( col, trans, sep = '_' ) )
#  }
#  
#  getPlotLayout <- function( np, WIDE = TRUE ){
#  
#    # np - no. plots
#  
#    if( np == 1 )return( list( mfrow = c( 1, 1 ), left = 1, bottom = c( 1, 2 ) ) )
#    if( np == 2 ){
#      if( WIDE )return( list( mfrow = c( 1, 2 ), left = 1, bottom = c( 1, 2 ) ) )
#      return( list( mfrow = c( 2, 1 ), left = c( 1, 2 ), bottom = 2 ) )
#    }
#  
#    if( np == 3 ){
#      if( WIDE )return( list( mfrow = c( 1, 3 ), left = 1, bottom = c( 1:3 ) ) )
#      return( list( mfrow = c( 3, 1 ), left = 1:3, bottom = 3 ) )
#    }
#    if( np <= 4 )return( list( mfrow = c( 2, 2 ), left = c( 1, 3 ), bottom = c( 3:4 ) ) )
#    if( np <= 6 ){
#      if( WIDE )return( list( mfrow = c( 2, 3 ), left = c( 1, 4 ), bottom = c( 4:6 ) ) )
#      return( list( mfrow = c( 3, 2 ), left = c( 1, 3, 5 ), bottom = 5:6 ) )
#    }
#    if( np <= 9 )return( list( mfrow = c( 3, 3 ), left = c( 1, 4, 7 ), bottom = c( 7:9 ) ) )
#    if( np <= 12 ){
#      if( WIDE )return( list( mfrow = c( 3, 4 ), left = c( 1, 5, 9 ), bottom = c( 9:12 ) ) )
#      return( list( mfrow = c( 4, 3 ), left = c( 1, 4, 7, 10 ), bottom = 10:12 ) )
#    }
#    if( np <= 16 )return( list( mfrow = c( 4, 4 ), left = c( 1, 5, 9, 13 ),
#                              bottom = c( 13:16 ) ) )
#    if( np <= 20 ){
#      if( WIDE )return( list( mfrow = c( 4, 5 ), left = c( 1, 6, 11, 15 ),
#                              bottom = c( 15:20 ) ) )
#      return( list( mfrow = c( 5, 4 ), left = c( 1, 5, 9, 13 ), bottom = 17:20 ) )
#    }
#    if( np <= 25 )return( list( mfrow = c( 5, 5 ), left = c( 1, 6, 11, 15, 20 ),
#                              bottom = c( 20:25 ) ) )
#    if( np <= 25 ){
#      if( WIDE )return( list( mfrow = c( 5, 6 ), left = c( 1, 6, 11, 15, 20, 25 ),
#                              bottom = c( 25:30 ) ) )
#      return( list( mfrow = c( 6, 5 ), left = c( 1, 6, 11, 16, 21, 26 ), bottom = 26:30 ) )
#    }
#    if( np <= 36 ){
#      return( list( mfrow = c( 6, 6 ), left = c( 1, 7, 13, 19, 25, 31 ), bottom = c( 31:36 ) ) )
#    }
#    return( list( mfrow = c( 7, 6 ), left = c( 1, 7, 13, 19, 25, 31, 37 ), bottom = c( 37:42 ) ) )
#  }
#  
#    plotFec <- function( fec, groups = NULL, LOG = F ){
#  
#      if( is.null(groups) ){
#        groupID <- rep(1, nrow(fec) )
#        ngroup  <- 1
#      }else{
#        group   <- fec[, groups]
#        groups  <- sort( unique( group ) )
#        groupID <- match( group, groups )
#        ngroup  <- length( groups )
#      }
#  
#      if( LOG )fec$fecEstMu <- log10(fec$fecEstMu)
#  
#      xlim  <- range( fec$year )
#      ylim  <- range( fec$fecEstMu, na.rm = T )
#      mfrow <- getPlotLayout(ngroup)
#  
#      par( mfrow = mfrow$mfrow, bty = 'n', mar = c(4,4,1,1), omi = c( .5, .5, .2, .2) )
#  
#      for( k in 1:ngroup ){
#  
#        fk    <- fec[ groupID == k, ]
#        tree  <- fk$treeID
#        trees <- sort( unique( tree ) )
#        ntree <- length( trees )
#  
#        plot( NA, xlim = xlim, ylim = ylim, xlab = '', ylab = '' )
#  
#        for(j in 1:ntree){
#          fj <- fk[ tree == trees[j], ]
#          lines( fj$year, fj$fecEstMu, lwd = 1 )
#        }
#        title( groups[k] )
#      }
#      mtext( 'Year', 1, outer = T, cex = 1.2 )
#      ytext <- 'Seeds per tree'
#      if( LOG ) ytext <- 'Seeds per tree (log_10)'
#      mtext( ytext, 2, outer = T, cex = 1.2 )
#    }

## ----eval = F-----------------------------------------------------------------
#  
#  d <- "https://github.com/jimclarkatduke/mast/blob/master/outputAbies.rdata?raw=True"
#  repmis::source_data( d )
#  
#  specs   <- sort( unique( fecPred$species ) )     # accumulate period estimates
#  yseq    <- seq( 0, 10, length = 100 )
#  intVal  <- matrix( 0, length(specs), 100 )
#  weight  <- intVal
#  rownames( intVal ) <- specs
#  plotSpecs <- sort( unique( fecPred$plotSpec ) )  # label trees in a plot-species group
#  
#  # time series for one species:
#  plotFec( fec = fecPred[ fecPred$species == 'abiesGrandis', ], groups = 'species', LOG = T )
#  
#  par( mfrow = c(1, 2), bty = 'n', mar = c(3,4,1,1), omi = c(.5,.1,.1,.1) )
#  
#  plot( NA, xlim = c(1, 20), ylim = c(.01, 1), xlab = '', ylab = 'Density/nyr', log = 'xy')
#  
#  for( i in 1:length(plotSpecs) ){
#  
#    wi  <- which( fecPred$plotSpec == plotSpecs[i] ) # tree-years in group
#    ci  <- fecPred$species[wi[1]]
#    col <- match( ci, specs )                        # color by species
#    tmp <- mastVolatility( treeID = fecPred$treeID[wi], year = fecPred$year[wi],
#                           fec = fecPred$fecEstMu[wi] )
#    if( is.null(tmp) )next
#  
#    intVal[ col, ] <- intVal[ col, ] + dnorm( yseq, tmp$stats['Period', 1], tmp$stats['Period', 2] )
#    weight[ col, ] <- weight[ col, ] + length( wi )
#  
#    # density +/- 1 SD
#    dens <- tmp$statsDensity
#    lines( dens[ 'Period', ], dens[ 'Mean', ], lwd = 2, col = getColor( col, .4) )
#    lines( dens[ 'Period', ], dens[ 'Mean', ] - dens[ 'SD', ], lty = 2, col = getColor( col, .4) )
#    lines( dens[ 'Period', ], dens[ 'Mean', ] + dens[ 'SD', ], lty = 2, col = getColor( col, .4) )
#  }
#  title( 'a) Plot-species groups' )
#  legend( 'topright', specs, text.col = c(1:length(specs)), bty = 'n', cex = .8 )
#  
#  intVal <- intVal*weight/weight
#  
#  plot( NA, xlim = c(2, 8), ylim = c(0, .12), xlab = '', ylab = 'Density' )  # density of mean intervals
#  for( i in 1:length(specs) ){
#    polygon( yseq, intVal[i,]/sum(intVal[i,]), border = i, col = getColor(i, .4) )
#  }
#  title( 'b) Period estimates' )
#  mtext( 'Year', 1, outer = T, cex = 1.3 )

## ----eval = F-----------------------------------------------------------------
#  spec  <- strsplit( rownames( betaYrRand ), '_' ) # extract species from ecoregion_species rownames
#  spec  <- sapply( spec, function(x) x[2] )
#  specs <- sort( unique(spec) )
#  
#  mastMatrix <- matrix( 0, nrow(betaYrRand), 5 )   # store stats by group
#  rownames(mastMatrix) <- rownames(betaYrRand)
#  colnames(mastMatrix) <- c( 'nyr', 'Variance', 'Volatility', 'Period Est', 'Period SD' )
#  
#  plot( NA, xlim = c(2, 10), ylim = c(.002, .4), xlab = '', ylab = 'Density/yr', log = 'xy')
#  
#  for(i in 1:nrow(betaYrRand)){
#  
#    wc <- which( betaYrRand[i,] != 0 )
#    if( length(wc) < 6 )next
#  
#    s <- mastSpectralDensity( betaYrRand[i,wc] )
#    if( !is.matrix( s$spect ) )next
#  
#    mastMatrix[i, ] <- c( length(wc), s$totVar, s$volatility, s$periodMu, s$periodSd )
#  
#    period <- 1/s$spec[, 'frequency' ]
#    dens   <- s$spec[, 'spectralDensity' ]/length(wc)  # series vary in length
#    col    <- match( spec[i], specs )
#  
#    lines( period, dens, lwd = 2, col = getColor( col, .4) )
#  }
#  title( 'c) Year effects' )
#  mtext( 'Period (yr)', 1, line = 1, outer = T )
#  
#  keepRows   <- which(  is.finite(mastMatrix[,'Variance']) & mastMatrix[,'Variance'] != 0 )
#  #keepCols   <- which( colSums( frequency, na.rm=T ) > 0 )
#  mastMatrix <- mastMatrix[ keepRows, ]

