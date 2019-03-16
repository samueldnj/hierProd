# --------------------------------------------------------------------------
# plots.R
# 
# Plots functions script for the coastwide (single stock) hierarchical
# multispecies assessment model simulation-estimation experiments.
# 
# Author: Samuel D N Johnson
# Date: 8 September, 2016
#
# --------------------------------------------------------------------------


plotCVbarsMultiTable <- function( tabRoots = c("DoverBase", "DoverLoPE"),
                                  titles = c("Strong","Weak"),
                                  pars = c( "DnT", "Bmsy", "BnT",  "Umsy"  ),
                                  models = c("noJointPriors","qPriorOnly","qUpriors","UmsyPriorOnly"),
                                  modLabs = expression("Single-Stock","None",q,q/U,U)  )
{
  par( mfrow = c( 1, length(tabRoots) ), oma = c(3,2,2,1), mar = c(1,2,1,2) )
  for( tabIdx in 1:length( tabRoots ) ) plotCVbars( tableRoot = tabRoots[tabIdx],
                                                    pars = pars,
                                                    models = models,
                                                    legend = (tabIdx == length(tabRoots) ),
                                                    modLabs = modLabs )
  mtext( outer = T, side = 1, text = "Coefficient of Variation (%)", line = 2)
  mtext( outer = T, side = 3, text = titles, cex = 1.5, at = seq(1/length(tabRoots)/2, by = 1/length(tabRoots)))
}


plotCVbars <- function( tableRoot = "DoverBase",
                        pars = c( "BnT", "Bmsy",  "Umsy", "DnT" ),
                        models = c("noJointPriors","qPriorOnly","qUpriors"),
                        legend = FALSE,
                        parLabels = expression( B[T]/B[0], B[MSY], B[T], U[MSY] ),
                        modLabs = expression("Single-Stock","None",q,q/U)  )
{

  # Read in table
  tabFile <- paste("./project/Statistics/",tableRoot,".csv", sep = "")
  tab <- read.csv(tabFile, header = T)

  estSS <- paste("ss", pars, sep = "" )
  cvSS  <- paste(estSS, "CV", sep = "" )

  estMS <- paste("ms", pars, sep = "" )
  cvMS  <- paste(estMS, "CV", sep = "" )

  checkPDhess <- function( x, hess )
  {
    if( !any( hess ) ) return(NA)
    else return(x)
  }

  tab <-  tab %>%
          filter( mp %in% models ) %>%
          dplyr::select(  cvSS, cvMS, mp, ssAICc, 
                          msAICc, msHessPD, ssHessPD ) %>%
          group_by( mp ) %>%
          mutate( ssAICcSum = sum(ssAICc) ) %>%
          ungroup()

  tab[tab == 0] <- NA

  # Create colour palette
  cols <- brewer.pal( n = length(models)+1, name = "Dark2")

  nModels <- length(models) + 1

  # Now let's plot CVs
  plot( x = c(0,110), y = c(0,length(pars) ), type = "n", axes = F,
        xlab = "",
        ylab = "" )
    axis( side = 1 )
    axis( side = 2, at = 0.5 + 0:(length(pars)-1), labels = parLabels,
          las = 1 )
    for( pIdx in 1:length(pars) )
    {
      parCVss <- cvSS[grepl(pattern = pars[pIdx], x = cvSS )]
      parCVms <- cvMS[grepl(pattern = pars[pIdx], x = cvMS )]
      plotTab <-  tab %>%
                  group_by( mp ) %>%
                  summarise_all( funs(mean,sd) )

      meanCVss <- as.numeric(unique( plotTab[,paste(parCVss,"mean",sep = "_")] ))
      sdCVss <- as.numeric(unique( plotTab[,paste(parCVss,"sd",sep = "_")] ))

      rect( xleft = 0, xright = c(meanCVss) * 100,
            ybottom = (pIdx - 1) + nModels/(nModels+1) - 1/(nModels + 1) / 2,
            ytop =  (pIdx - 1)+nModels/(nModels+1) + 1/(nModels + 1)/2,
            col = cols[1], border = NA )
      segments( x0 = (c(meanCVss) - c(sdCVss)) * 100,
                x1 = (c(meanCVss) + c(sdCVss)) * 100,
                y0 = (pIdx - 1) + (nModels)/(nModels+1),
                lwd = 3, col = "black" )
      for( mIdx in 1:length(models) )
      {
        plotTabMod <-  plotTab %>%
                    filter( mp == models[mIdx] )

        meanCVms <- unlist(plotTabMod[,paste(parCVms,"mean",sep = "_")])
        sdCVms <- unlist(plotTabMod[,paste(parCVms,"sd",sep = "_")])

        if( is.na(meanCVms) | is.na(sdCVms) | length(meanCVms) == 0 | length(sdCVms) == 0 ) next

        rect( xleft = 0, xright = meanCVms * 100,
              ybottom = (pIdx - 1) + (nModels - mIdx)/(nModels+1) - 1/(nModels + 1) / 2,
              ytop =  (pIdx - 1)+(nModels - mIdx)/(nModels+1) + 1/(nModels + 1)/2,
              col = cols[mIdx + 1], border = NA )
        segments( x0 = (c(meanCVms) - c(sdCVms)) * 100,
                  x1 = (c(meanCVms) + c(sdCVms)) * 100,
                  y0 = (pIdx - 1) + (nModels - mIdx)/(nModels+1),
                  lwd = 3, col = "black" )
      }
    }
    if(legend) panLegend( x = "bottomright", legTxt = modLabs,
                          fill = cols, border = NA, bty = "n" )

}

plotStatTableGraphs <- function(  tableRoot = "allSame_fixedProcRE",
                                  resp = c("BnT","Umsy","Bmsy","Dep","HessPD","q_1","q_2","tau2_1","tau2_2"),
                                  axes = c("kappaMult","corr"),
                                  MARE = TRUE,
                                  RE = TRUE,
                                  groupPars = FALSE )
{
  # plotStatTableGraphs()
  # Plots the statistics table information and saves into
  # the figs directory.
  # inputs:   tableRoot = RE tableName root in statistics folder
  #           resp = response variables for plotting
  #           axes = 2-character giving factors to plot over
  # outputs:  NULL
  # side-eff: plots to quartz or figs directory

  # Create save directory
  saveDir <- file.path(getwd(),"project","figs",tableRoot)
  dir.create(saveDir)

  # First, plot MARE values
  # a place to save them
  MAREdir <- file.path(saveDir,"MARE")
  dir.create(MAREdir)
  # Now load the RE table
  MAREtabRoot  <- paste(tableRoot,"_MARE", sep = "")
  MAREtabFile  <- paste(tableRoot,"_MARE.csv", sep = "")
  MAREtabPath <- file.path(getwd(),"project","Statistics",MAREtabFile)
  MAREtab     <- read.csv( MAREtabPath, header=TRUE, stringsAsFactors=FALSE )
  nSpp        <- unique(MAREtab$nS)
  MPs         <- unique(MAREtab$mp)



  # Loop over complex sizes
  if( MARE )
  {
    for( sIdx in 1:length(nSpp) )
    {
      # Make folder for this complex size
      nSp           <- nSpp[sIdx]
      nSpFolderName <- paste(nSp,"stocks",sep = "_")
      nSpFolder     <- file.path(MAREdir,nSpFolderName)
      dir.create(nSpFolder)
      nSppMAREtab   <- MAREtab %>% filter( nS == nSp )
      spp           <- unique(nSppMAREtab$species) 

      # loop over response variables
      for( rIdx in 1:length(resp) )
      {
        # loop over species/stocks
        for( spIdx in 1:nSp )
        {
          response      <- resp[rIdx]
          sppName       <- spp[spIdx]
          plotFile      <- paste(response,"_",sppName,"_",nSp,"stocks_",MAREtabRoot,"_byMP.pdf", sep = "") 
          plotPath      <- file.path(nSpFolder,plotFile)
          pdf( file = plotPath, width = 17, height = 11 )
          plotMAREs(  tableName = MAREtabRoot, table = nSppMAREtab,
                      axes = c(axes,"mp"), nSp = nSp, spec = sppName, resp = response )
          dev.off()  
        }
        plotFile      <- paste(response,"_allStocks_",nSp,"stocks_",MAREtabRoot,"_byMP.pdf", sep = "") 
        plotPath      <- file.path(nSpFolder,plotFile)
        pdf( file = plotPath, width = 17, height = 11 )
        plotMAREs(  tableName = MAREtabRoot, table = nSppMAREtab,
                    axes = c(axes,"mp"), nSp = nSp, spec = spp, resp = response )
        dev.off()     
      }
    }
    # Now plot for all stock sizes aggregated
    folderName <- paste("allComplexSizes")
    folderPath <- file.path(MAREdir,folderName)
    dir.create( folderPath )
    spp        <- unique(MAREtab$species) 
    for( rIdx in 1:length(resp) )
      {
        for( spIdx in 1:max(nSpp) )
        {
          response      <- resp[rIdx]
          sppName       <- spp[spIdx]
          plotFile      <- paste(response,"_",sppName,"_allComplexSizes_",MAREtabRoot,"_byMP.pdf", sep = "") 
          plotPath      <- file.path(folderPath,plotFile)
          pdf( file = plotPath, width = 17, height = 11 )
          plotMAREs(  tableName = MAREtabRoot, table = MAREtab,
                      axes = c(axes,"mp"), nSp = NULL, spec = sppName, resp = response )
          dev.off()  
        }
        plotFile      <- paste(response,"_allStocks_allComplexSizes_",MAREtabRoot,"_byMP.pdf", sep = "") 
        plotPath      <- file.path(folderPath,plotFile)
        pdf( file = plotPath, width = 17, height = 11 )
        plotMAREs(  tableName = MAREtabRoot, table = MAREtab,
                    axes = c(axes,"mp"), nSp = NULL, spec = spp, resp = response )
        dev.off()     
      }

  }
  

  # Now plot relative error distributions
  REdir <- file.path(saveDir,"RE")
  dir.create(REdir)
  # Now load the RE table
  REtabRoot   <- paste(tableRoot,"_RE", sep = "")
  REtabFile   <- paste(tableRoot,"_RE.csv", sep = "")
  REtabPath   <- file.path(getwd(),"project","Statistics",REtabFile)
  if( RE )
  {
    REtab       <- read.csv( REtabPath, header=TRUE, stringsAsFactors=FALSE )
    nSpp        <- unique(REtab$nS)

    # Loop over complex sizes
    for( sIdx in 1:length(nSpp) )
    {
      # Make folder for this complex size
      nSp           <- nSpp[sIdx]
      nSpFolderName <- paste(nSp,"stocks",sep = "_")
      nSpFolder     <- file.path(REdir,nSpFolderName)
      dir.create(nSpFolder)
      nSppREtab   <- REtab %>% filter( nS == nSp )
      spp           <- unique(nSppREtab$species) 

      for( rIdx in 1:length(resp) )
      {
        # plot RE dists for each stock
        for( spIdx in 1:length(spp) )
        { 
          sppName       <- spp[spIdx]
          response      <- resp[rIdx]
          plotFile      <- paste(response,"_",sppName,"_",nSp,"stocks_",REtabRoot,"_byMP.pdf", sep = "") 
          plotPath      <- file.path(nSpFolder,plotFile)
          pdf( file = plotPath, width = 11, height = 17 )
          plotREdists(  tableName = REtabRoot, 
                        table = nSppREtab,
                        axes = c(axes,"mp"), 
                        nSp = nSp, 
                        spec = sppName, 
                        resp = response )
          dev.off()
        }
        # Now plot RE dists over the whole complex - might not be useful but
        # should check anyhow
        plotFile      <- paste(response,"_allStocks_",nSp,"stocks_",REtabRoot,"_byMP.pdf", sep = "") 
        plotPath      <- file.path(nSpFolder,plotFile)
        pdf( file = plotPath, width = 11, height = 17 )
        plotREdists(  tableName = REtabRoot, 
                      table = nSppREtab,
                      axes = c(axes,"mp"), 
                      nSp = nSp, 
                      spec = NULL, 
                      resp = response )
        dev.off()   
      }
    }

    # Now do it for all complex sizes
    nSp           <- max(nSpp)
    folderName    <- paste("allComplexSizes",sep = "_")
    folderPath    <- file.path(REdir,folderName)
    dir.create(folderPath)
    spp           <- unique(REtab$species) 

    for( rIdx in 1:length(resp) )
    {
      # plot RE dists for each stock
      for( spIdx in 1:length(spp) )
      { 
        sppName       <- spp[spIdx]
        response      <- resp[rIdx]
        plotFile      <- paste(response,"_",sppName,"_allComplexSizes_",REtabRoot,"_byMP.pdf", sep = "") 
        plotPath      <- file.path(folderPath,plotFile)
        pdf( file = plotPath, width = 11, height = 17 )
        plotREdists(  tableName = REtabRoot, 
                      table = REtab,
                      axes = c(axes,"mp"), 
                      nSp = NULL, 
                      spec = sppName, 
                      resp = response )
        dev.off()
      }
      # Now plot RE dists over the whole complex - might not be useful but
      # should check anyhow
      plotFile      <- paste(response,"_allStocks_allComplexSizes_",REtabRoot,"_byMP.pdf", sep = "") 
      plotPath      <- file.path(folderPath,plotFile)
      pdf( file = plotPath, width = 11, height = 17 )
      plotREdists(  tableName = REtabRoot, 
                    table = REtab,
                    axes = c(axes,"mp"), 
                    nSp = NULL, 
                    spec = NULL, 
                    resp = response )
      dev.off()   
    }
  }
  
  if( groupPars )
  {
    # Now plot group parameters
    MREtabRoot <- paste(tableRoot,"_MRE", sep = "" )
    for( nSp in nSpp )
    {
      for( mIdx in 1:length(MPs))
      {
        MP <- MPs[mIdx]
        plotFile <- paste( nSp, "nS_", MP,"_groupPars_", tableRoot,".pdf", sep = "" ) 
        plotPath <- file.path( saveDir, plotFile )
        pdf( file = plotPath, width = 11, height = 17 )
        plotGroupPars(  tableName = MREtabRoot,
                        axes = axes,
                        nSp = nSp,
                        MP = MP,
                        spec = "Stock1"
                     )  
        dev.off()
      } 
    }  
  }
  cat("Table plotting complete for ", tableRoot, ".\n", sep = "")
}


plotLatinSquareRect <- function( sqDim = 4, rectY = 3 )
{
  # Set up plotting area
  par( mfrow = c(2,1), mar = c(2,2,2,2), oma = c(1,1,1,1) )

  # First plot a latin square
  plot(x = c(0,1), y = c(0,1), type = "n", axes = F )
    rect( xleft = 0, xright = 1, ybottom = 0, ytop = 1 )
    for( i in 1:(sqDim-1) )
    {
      segments( x0 = 0, y0 = i/(sqDim), x1 = 1, y1 = i/(sqDim) )
      segments( y0 = 0, x0 = i/(sqDim), y1 = 1, x1 = i/(sqDim) )
    }
    centres <- seq( from = 1 / sqDim / 2, by = 1 / sqDim, length = sqDim)
    panLab( x = centres, y = centres[sample(1:sqDim,size = sqDim)],
            txt = "X", cex = 2  )
    mtext( side = 1, text = "(a)" )

  # Now the rectangle
  plot(x = c(0,1), y = c(0,rectY/sqDim), type = "n", axes = F )
    rect( xleft = 0, xright = 1, ybottom = 0, ytop = rectY/sqDim )
    for( i in 1:(sqDim-1) )
    {
      segments( x0 = 0, y0 = i/(sqDim), x1 = 1, y1 = i/(sqDim) )
      segments( y0 = 0, x0 = i/(sqDim), y1 = 1, x1 = i/(sqDim) )
    }
    browser()
    centres <- seq( from = 1 / sqDim / 2, by = 1 / sqDim, length = sqDim)
    panLab( x = centres, y = centres[sample(1:rectY,size = rectY)],
            txt = "X", cex = 1  )
    mtext( side = 1, text = "(b)" )
}



plotGroupPars <- function(  tableName = "rKq_msInc__MRE",
                            axes = c("qOM","UmsyOM"),
                            nSp = 10,
                            MP = "allJointPriors",
                            spec = "Stock1" )
{
  # plotGroupPars()
  # Plots group prior parameter estimates 
  # as a function of axes[1], grouped by axes[2].
  # inputs:   tableName = RE table in statistics folder
  #           nSp = number of species to restrict to
  #           axes = 2-numeric giving factor levels to plot over
  #           par = character name of group level parameter to plot
  # outputs:  NULL
  # side-eff: plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 

  if(!is.null(spec)) table <- table %>% filter(species == spec)

  # restrict to correct number of species, summarise group pars since we're
  # only using 2 axes
  table  <-   table %>%
              filter( mp == MP, nS %in% nSp, species %in% spec ) %>%
              group_by( kappaMult, corr, Umax, tUpeak ) %>%
              summarise_if(.predicate = "is.numeric",.funs="median") %>%
              ungroup()


  # Get levels for axes
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }
  # order levels
  x   <- unique(table %>% pull(axesCols[1]))
  x   <- x[order(x)]
  y   <- unique(table %>% pull(axesCols[2]))
  y   <- y[order(y)]

  groupPars <- c("Umsybar","tauq2","sigU2")
  
  # Create an array to hold info
  z <- array( NA, dim = c( length(x), length(y), length(groupPars), 3 ),
              dimnames = list(  paste(axes[1],x,sep = ""),
                                paste(axes[2],y,sep = ""),
                                groupPars,
                                c(0.025, 0.5, 0.975) ) )

  # Fill the array, looping over species and axes
  for( xIdx in 1:length(x) )
  {
    for( yIdx in 1:length(y) )
    {
      # Get x and y values
      xVal <- x[ xIdx ]
      yVal <- y[ yIdx ]
      # Get rows
      xRows <- which ( table[ ,axesCols[ 1 ] ] == xVal )
      yRows <- which ( table[ ,axesCols[ 2 ] ] == yVal )
      rows <- intersect( xRows, yRows )
      # browser()
      if( length(rows) == 0 ) next
      filtered <- table[rows,]
      if(nrow(filtered) > 1) filtered <- filtered[1,]
      # Fill array
      z[xIdx,yIdx,,2] <- as.numeric(filtered[1,groupPars])
      z[xIdx,yIdx,,1] <- as.numeric(filtered[1,paste(groupPars,"025",sep = "")])
      z[xIdx,yIdx,,3] <- as.numeric(filtered[1,paste(groupPars,"975",sep = "")])
    }
  }
  # Set up plotting colours and a jitter for different axes
  cols      <- brewer.pal( n = length(y), name = "Dark2" )
  xShift    <- 1/(length(y)+2)
  leftShift <- xShift*length(y)/2
  # Set up plotting environment
  par(mfrow = c(2,2), mar = c(2,2,1,1), oma = c(3,3,2,2) )
  # Loop over the group parameters
  for( p in 1:length(groupPars) )
  {
    plot( x = c(1-leftShift,length(x) + leftShift), 
          y = range(z[,,groupPars[p],], na.rm=T), axes = F, type = "n",
          xlab = "", ylab = "", main = groupPars[p] )
    axis( side = 1, at = 1:length(x), labels = x )
    axis( side = 2, las = 1 )
    for( yIdx in 1:length(y) )
    {
      lines(  x = 1:length(x) - leftShift + xShift*(yIdx-1), 
              y = z[,yIdx,groupPars[p],2], col = cols[yIdx], lwd = 0.5 )
      points( x = 1:length(x) - leftShift + xShift*(yIdx-1), 
              y = z[,yIdx,groupPars[p],2], col = cols[yIdx], pch = 16 )
      segments( x0 = 1:length(x) - leftShift + xShift*(yIdx-1),
                x1 = 1:length(x) - leftShift + xShift*(yIdx-1),
                y0 = z[,yIdx,groupPars[p],1], 
                y1 = z[,yIdx,groupPars[p],3], col = cols[yIdx] )
    }
  }
  
  panLegend(  legTxt = paste( axes[2], " = ", y ),
              lty = 1, pch = 16, cex = 1,
              lwd = 0.8, col = cols,
              x = 0.1, y = 0.9, bty = "n" )
  mtext( side = 1, outer = TRUE, text = axes[1], line = 1.5 )  
}


plotRastersPub <- function( saveName = "LSrasters.png",
                            tabRoot = "info_qREs_RP_MARE",
                            axes = c("nS","nDiff"),
                            titles = c("Complex Size","Number of Low Information Stocks"),
                            breakRange = c(-1,1),
                            MP = NULL,
                            resp = c("DeltaBnT","DeltaUmsy","DeltaDep"),
                            respTitles = expression( B[T], U[MSY], B[T]/B[0]),
                            save = T )
{
  saveFile <- file.path(".","project","figs",saveName)

  if(save) png( filename = saveFile, width = 1100, height = 850,
                type = "quartz")

  par( mfrow = c(2,length(resp)), mar = c(2,1.5,2,1), oma = c(4,6,3,1) )
  plotRasters(  tableName = paste(tabRoot,"Delta",sep = "_"),
                axes = axes, breakRange = breakRange, MP = MP, resp = resp,
                setPar = F, titles = respTitles )

  plotRasters(  tableName = paste(tabRoot,"Cplx",sep = "_"),
                axes = axes, breakRange = breakRange, MP = MP, resp = resp,
                setPar = F, titles = "" )  

  mtext( outer = T, text  = c("(a)", "(b)"), side = 2, line = 1,
          at = c(0.75, 0.25), las = 1)
  mtext( outer = T, side = 1, line = 2, text = titles[1], cex = 2 )
  mtext( outer = T, side = 2, line = 3, text = titles[2], las = 0, cex = 2 )
  if(save) dev.off()
  cat("Done!\n")
}





plotRasters <- function(  tableName = "allSame_fixedProcRE_MARE_Cplx",
                          axes = c("corr","kappaMult"),
                          resp = c("DeltaBnT","DeltaUmsy","Deltaq_1","Deltaq_2","DeltaDep","DeltaBmsy"),
                          spec = c("Stock1"),
                          MP = NULL,
                          nSp = c(4,7,10),
                          table = NULL,
                          breakRange = c(-1,1),
                          save = F,
                          titles = NULL,
                          setPar = T
                      )
{
  # plotRasterss()
  # Plots MARE (response variable) values as a function
  # of given axes (depends on experimental factors). 
  # Will group points using lines corresponding to the first axis (w) and 
  # plot across values of the second axes (x), with multi-panel rows 
  # corresponding to the third (y). Panel columns are 
  # single and Joint AMs.
  # inputs:   tableName = RE table in statistics folder
  #           axes = 3-char vector naming RE table columns
  #           resp = response variable of interest for RE values 
  #                  (usually a leading or derived AM parameter)
  #           table = optional df to replace loading file, useful for 
  #                   automated plotting to directory
  # outputs:  NULL

  # Read in table if none supplied
  if( is.null(table) )
  {
    # Load table
    fileName <- paste( tableName, ".csv", sep = "" )
    tablePath <- file.path ( getwd(), "project/Statistics", fileName )
    table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 
  }
  tab <- table
  # Restrict to subtables if requested
  if(!is.null(MP))    tab  <-   tab %>% filter( mp %in% MP ) 
  if(!is.null(spec))  tab  <-   tab %>% filter( species %in% spec ) 
  if(!is.null(nSp))   tab  <-   tab %>% filter( nS %in% nSp )

  # Now make rasters
  brickObj <- makeRasterFromTable( table = tab, axes = axes, colnames = resp )

  if( save )
  {
    outFile <- paste( tableRoot, "complexPubRasters.pdf", sep = "")
    outFile <- file.path("./project/figs/",outFile)
    pdf ( file = outFile, width = wdh, height = hgt )
  }

  if( setPar ) par(mfrow = c(length(resp),1), mar = c(2,2,3,1), oma = c(2,2,1,1) ) 
  colBreaks   <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  redScale    <- brewer.pal(8, "Reds")
  greenScale  <- brewer.pal(8, "Greens" )
  rgScale     <- c(redScale[8:1],greenScale)
  rrScale     <- c(redScale[8:1],redScale)

  # Now plot the raster brick layers
  for(k in 1:length(resp) )
  {
    plot( brickObj$brick[[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rgScale, breaks = colBreaks,
        axes = FALSE, cex.main = 1.5 )
      if(is.null(titles)) 
        mtext( text = names(brickObj$brick)[k], side = 3, line = 1, cex = 1.5)
        else mtext( text = titles[k], side = 3, line = 1, cex = 1.5)
      # if( k == 1 ) mtext ( text = "(a)", side = 2, las = 1, line = 2 )
      text( brickObj$brick[[k]], digits = 2, halo = T, cex = 2, axes = FALSE )
      axis( side = 2, at = 1:length(brickObj$y) - 0.5, labels = brickObj$y, las = 1, cex.axis = 1.5 )
      axis( side = 1, at = 1:length(brickObj$x) - 0.5, labels = brickObj$x, cex.axis = 1.5 ) 
      box(lwd = 2)
  }  

  if(save) dev.off()
}



plotMAREs <- function(  tableName = "allSame_RE_msIncr__MARE",
                        axes = c("corr","kappaMult","mp"),
                        resp = "Umsy",
                        spec = "Stock1",
                        MP = NULL,
                        nSp = 6,
                        table = NULL,
                        ssModel = NULL,
                        msModel = NULL
                      )
{
  # plotMAREs()
  # Plots MARE/MRE (response variable) values as a function
  # of given axes (depends on experimental factors). 
  # Will group points using lines corresponding to the first axis (w) and 
  # plot across values of the second axes (x), with multi-panel rows 
  # corresponding to the third (y). Panel columns are 
  # single and Joint AMs.
  # inputs:   tableName = RE table in statistics folder
  #           axes = 3-char vector naming RE table columns
  #           resp = response variable of interest for RE values 
  #                  (usually a leading or derived AM parameter)
  #           table = optional df to replace loading file, useful for 
  #                   automated plotting to directory
  # outputs:  NULL

  # Read in table if none supplied
  if( is.null(table) )
  {
    # Load table
    fileName <- paste( tableName, ".csv", sep = "" )
    tablePath <- file.path ( getwd(), "project/Statistics", fileName )
    table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 
  }
  tab <- table
  # Restrict to subtables if requested
  if(!is.null(MP))    tab  <-   tab %>% filter( mp == MP ) 
  if(!is.null(spec))  tab  <-   tab %>% filter( species %in% spec ) 
  if(!is.null(nSp))   tab  <-   tab %>% filter( nS %in% nSp )

  # # Run predictions if explanatory model supplied
  # predTab <- tab %>% mutate( ssPred = NA, msPred = NA )
  # if( !is.null(ssModel) )
  # {
  #   predTab <-  predTab %>% 
  #               mutate( ssPred = predict( ssModel, newdata = predTab ) )
  # }

  # if( !is.null(msModel) )
  # {
  #   predTab <-  predTab %>% 
  #               mutate( msPred = predict( msModel, newdata = predTab ) )
  # }

  # response colnames
  respCols <- paste(c("ss","ms"),resp,sep = "")


  # Get levels for axes
  axesCols <- integer(length=length(axes))
  # predAxesCols <- integer(length=length(axes))
  for ( i in 1:3)
  {
    axesCols[i]     <- which(names(tab) == axes[i])
    # predAxesCols[i] <- which(names(predTab) == axes[i])
  }
  # order levels
  w   <- unique(tab %>% pull(axesCols[1]))
  w   <- w[order(w)]
  x   <- unique(tab %>% pull(axesCols[2]))
  x   <- x[order(x)]
  y   <- unique(tab %>% pull(axesCols[3]))
  y   <- y[order(y)]

  # Create an array to hold info
  z <- array( NA, dim = c( length(w), length(x), length(y), 2 ),
              dimnames = list(  paste(axes[1],w,sep = ""),
                                paste(axes[2],x,sep = ""),
                                paste(axes[3],y,sep = ""),
                                respCols ) )

  # Create an array to hold info
  zPred <- array( NA, dim = c( length(w), length(x), length(y), 2 ),
                  dimnames = list(  paste(axes[1],w,sep = ""),
                                    paste(axes[2],x,sep = ""),
                                    paste(axes[3],y,sep = ""),
                                    respCols ) )


  # Fill the array, looping over species and axes
  for( wIdx in 1:length(w) )
  {
    for( xIdx in 1:length(x) )
    {
      for( yIdx in 1:length(y) )
      {
        # Get x and y values
        wVal <- w[ wIdx ]
        xVal <- x[ xIdx ]
        yVal <- y[ yIdx ]
        # Get rows
        wRows <- which ( tab[,axesCols[1] ] == wVal)
        xRows <- which ( tab[,axesCols[2] ] == xVal)
        yRows <- which ( tab[,axesCols[3] ] == yVal)
        # wRowsPred <- which ( predTab[,predAxesCols[1] ] == wVal)
        # xRowsPred <- which ( predTab[,predAxesCols[2] ] == xVal)
        # yRowsPred <- which ( predTab[,predAxesCols[3] ] == yVal)
        rows      <- intersect(intersect(wRows,xRows),yRows)
        # predRows  <- intersect(intersect(wRowsPred,xRowsPred),yRowsPred)
        if( length(rows) == 0 ) next
        filtered <- tab[rows,respCols]
        if( nrow(filtered) > 1 ) 
          filtered <-  apply(X = filtered, FUN = mean, MARGIN = 2, na.rm = T)
        # Fill array
        z[wIdx,xIdx,yIdx,respCols] <- as.numeric( filtered )
        # Now fill zPred array  
        # predTab <- predTab[predRows,c("ssPred","msPred")]
        # if(any(!is.na(predTab)))
        # {
        #   if( nrow(predTab > 1) ) predTab <- apply(X = predTab, FUN = mean, MARGIN = 2, na.rm = T)
        #   zPred[wIdx,xIdx,yIdx,respCols] <- as.numeric( predTab[ , c( "ssPred", "msPred" ) ] )
        # }  
      }
    }
  }
  if(any(grepl(pattern="msHessPD", x = respCols)) ) z[,,,"msHessPD"] <- z[,,,"msHessPD"] + 100

  zDelta <- (z[,,,1] / z[,,,2]) - 1

  # Because of missing factor levels for some stocks, the dimension
  # of zDelta might reduce too much for the plotting code below. This
  # will fix it.
  if( length( dim( zDelta ) ) < length( dim(z) ) - 1 )
  {
    missingDim            <- which( !(dim(z)[1:3] %in% dim(zDelta)) )
    newDims               <- numeric( length = 3 )
    newDims[missingDim]   <- 1
    newDims[-missingDim]  <- dim( zDelta )
    newDimNames           <- dimnames( z )[1:3]
    zDelta <- array(zDelta, dim = newDims, dimnames = newDimNames )
  }

  cols <- brewer.pal( n = length(w), name = "Dark2" )
  xShift <- 1 / (length(w) + 2)
  leftShift <- length(w)*xShift / 2

  
  # Plot contents. Use w values as line groupings, x values as
  # horizontal axis points, and y values as panel row groups
  par( mfcol = c(3, length(y)), mar = c(2,2,2,1), oma = c(3,3,2,2) )
  for( yIdx in 1:length(y) )
  {
    # Plot SS model
    plot( x = c(1 - leftShift,length(x)+leftShift), 
          y = c(0,max(z, na.rm = T )), type = "n",
          axes = FALSE )
      mtext(side = 3,, text =  y[yIdx])
      if( yIdx == length(y) ) mtext( side = 4, text = "Single Species Model", cex = 0.8)
      if( yIdx == 1)  mtext( side = 2, outer = F, text = "MARE (%)", cex = 1.2, line = 2 )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      for( wIdx in 1:length(w) )
      {
        rect( xleft = 1:length(x) - leftShift + xShift*(wIdx-1), ybottom = 0, 
              xright = 1:length(x) - leftShift + xShift*(wIdx), ytop = z[ wIdx, , yIdx, respCols[1] ],
              col = cols[wIdx] )
      }
    # plot MS model
    plot( x = c(1 - leftShift,length(x)+leftShift), 
          y = c(0,max(z,na.rm = T)), type = "n",
          axes = FALSE )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      if( yIdx == length(y) ) mtext( side = 4, text = "Joint Model", cex = 0.8)
      if( yIdx == 1)  mtext( side = 2, outer = F, text = "MARE (%)", cex = 1.2, line = 2 )
      for( wIdx in 1:length(w) )
      {
        rect( xleft = 1:length(x) + xShift*(wIdx-1) - leftShift, ybottom = 0, 
              xright = 1:length(x) + xShift*(wIdx) - leftShift , ytop = z[ wIdx, , yIdx, respCols[2] ],
              col = cols[wIdx] )
      }

    # plot Delta statistic
    plot( x = c(1 - leftShift,length(x)+leftShift), 
          y = range(zDelta,na.rm = T), type = "n",
          axes = FALSE )
      if( yIdx == length(y) ) panLegend(  legTxt = paste( axes[1], " = ", w ),
                                  lty = 1, pch = 16, cex = 1,
                                  lwd = 0.8, col = cols,
                                  x = 0.1, y = 0.7, bty = "n" )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      if( yIdx == length(y) ) mtext( side = 4, text = "Delta Values", cex = 0.8)
      if( yIdx == 1)  mtext( side = 2, outer = F, text = "Delta", cex = 1.2, line = 2 )
      for( wIdx in 1:length(w) )
      {
        xShift <- 1/(length(w) + 2)
        rect( xleft = 1:length(x) + xShift*(wIdx-1) - leftShift, ybottom = 0, 
              xright = 1:length(x) + xShift*(wIdx) - leftShift , ytop = zDelta[ wIdx, , yIdx ],
              col = cols[wIdx] )
      }
  }
  mtext( side = 3, outer = T, text = tableName, cex = 1.2 )
  mtext( side = 1, outer =T, text = axes[2], cex = 1.2 )
  mtext( side = 1, outer = T, text = resp, cex = 1, col = "grey50",
          adj = 0.9, line = 1.3 )
  mtext( side = 1, outer = T, text = spec, cex = 1, col = "grey50",
          adj = 0.1, line = 1.3 )
}

plotREdists <- function(  tableName = "allSame_RE_msIncr__RE",
                          axes = c("corr","kappaMult","mp"),
                          resp = "Umsy",
                          spec = "Stock1",
                          MP = NULL,
                          nSp = 6,
                          table = NULL,
                          ssModel = NULL,
                          msModel = NULL,
                          yLim = c(-2,2)
                        )
{
  # plotREdists()
  # Plots relative error distributions as  a function
  # of given axes (depends on experimental factors). 
  # Will group points using lines corresponding to the first axis (w) and 
  # plot across values of the second axes (x), with multi-panel rows 
  # corresponding to the third (y). Panel columns are 
  # single and Joint AMs.
  # inputs:   tableName = RE table in statistics folder
  #           axes = 3-char vector naming RE table columns
  #           resp = response variable of interest for RE values 
  #                  (usually a leading or derived AM parameter)
  #           table = optional df to replace loading file, useful for 
  #                   automated plotting to directory
  # outputs:  NULL

  # Read in table if none supplied
  if( is.null(table) )
  {
    # Load table
    fileName <- paste( tableName, ".csv", sep = "" )
    tablePath <- file.path ( getwd(), "project/Statistics", fileName )
    table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 
  }
  tab <- table
  # Restrict to subtables if requested
  if(!is.null(MP))    tab  <-   tab %>% filter( mp == MP ) 
  if(!is.null(spec))  tab  <-   tab %>% filter( species %in% spec ) 
  if(!is.null(nSp))   tab  <-   tab %>% filter( nS %in% nSp )

  # Run predictions if explanatory model supplied
  predTab <- tab %>% mutate( ssPred = NA, msPred = NA )
  if( !is.null(ssModel) )
  {
    predTab <-  predTab %>% 
                mutate( ssPred = predict( ssModel, newdata = predTab ) )
  }

  if( !is.null(msModel) )
  {
    predTab <-  predTab %>% 
                mutate( msPred = predict( msModel, newdata = predTab ) )
  }

  # response colnames
  respCols <- paste(c("ss","ms"),resp,sep = "")


  # Get levels for axes
  axesCols <- integer(length=length(axes))
  predAxesCols <- integer(length=length(axes))
  for ( i in 1:3)
  {
    axesCols[i]     <- which(names(tab) == axes[i])
    predAxesCols[i] <- which(names(predTab) == axes[i])
  }
  # order levels
  w   <- unique(tab %>% pull(axesCols[1]))
  w   <- w[order(w)]
  x   <- unique(tab %>% pull(axesCols[2]))
  x   <- x[order(x)]
  y   <- unique(tab %>% pull(axesCols[3]))
  y   <- y[order(y)]

  # Create an array to hold info
  z <- array( NA, dim = c( length(w), length(x), length(y), 2, 3 ),
              dimnames = list(  paste(axes[1],w,sep = ""),
                                paste(axes[2],x,sep = ""),
                                paste(axes[3],y,sep = ""),
                                respCols, c("0.025","0.5","0.975") ) )

  # Create an array to hold info
  zPred <- array( NA, dim = c( length(w), length(x), length(y), 2 ),
                  dimnames = list(  paste(axes[1],w,sep = ""),
                                    paste(axes[2],x,sep = ""),
                                    paste(axes[3],y,sep = ""),
                                    respCols ) )

  # Fill the array, looping over species and axes
  for( wIdx in 1:length(w) )
  {
    for( xIdx in 1:length(x) )
    {
      for( yIdx in 1:length(y) )
      {
        # Get x and y values
        wVal <- w[ wIdx ]
        xVal <- x[ xIdx ]
        yVal <- y[ yIdx ]
        # Get rows
        wRows <- which ( tab[,axesCols[1] ] == wVal)
        xRows <- which ( tab[,axesCols[2] ] == xVal)
        yRows <- which ( tab[,axesCols[3] ] == yVal)
        wRowsPred <- which ( predTab[ , predAxesCols[ 1 ] ] == wVal)
        xRowsPred <- which ( predTab[ , predAxesCols[ 2 ] ] == xVal)
        yRowsPred <- which ( predTab[ , predAxesCols[ 3 ] ] == yVal)
        rows      <- intersect(intersect(wRows,xRows),yRows)
        predRows  <- intersect(intersect(wRowsPred,xRowsPred),yRowsPred)
        if( length(rows) == 0 ) next
        filtered <- tab[rows,respCols]
        quantiles <- apply( X = filtered, FUN = quantile, 
                            probs = c(0.025,0.5, 0.975), MARGIN = 2,
                            na.rm = T)
        medPred <- apply( X = predTab[predRows,c("ssPred","msPred")],
                          FUN = median, na.rm = T, MARGIN = 2 )
        # Fill array
        # browser()
        z[wIdx,xIdx,yIdx,,] <- t(quantiles)
        zPred[wIdx,xIdx,yIdx,] <- as.numeric(medPred)
      }
    }
  }

  cols <- brewer.pal( n = length(w), name = "Dark2" )
  xShift <- 1 / (length(w) + 2)
  leftShift <- length(w)*xShift / 2
  
  if( is.null(yLim) ) yLim <- range(z,na.rm = T)

  # Plot contents. Use w values as line groupings, x values as
  # horizontal axis points, and y values as panel row groups
  par( mfrow = c(length(y), 2), mar = c(2,2,2,1), oma = c(3,3,2,2) )
  for( yIdx in 1:length(y) )
  {
    # Plot SS model
    plot( x = c(1-leftShift,length(x)+leftShift), 
          y = yLim, type = "n",
          axes = FALSE )
      if( yIdx == 1 ) panLegend(  legTxt = paste( axes[1], " = ", w ),
                                  lty = 1, pch = 16, cex = 1,
                                  lwd = 0.8, col = cols,
                                  x = 0.1, y = 0.7, bty = "n" )
      if( yIdx == 1 ) mtext( side = 3, text = "Single Species Model", cex = 0.8)
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      # axis( side = 3, at = 1:length(x), labels = FALSE )
      abline(h = 0, lty = 2, lwd = 0.8)
      for( wIdx in 1:length(w) )
      {
        lines(  x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[1], 2], col = cols[wIdx],
                lwd = 0.6, lty = 3 )
        points( x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[1], 2 ], col = cols[wIdx],
                pch = 16, cex = 2 )
        segments( x0 = 1:length(x) - leftShift + xShift*(wIdx-1), y0 = z[ wIdx, , yIdx, respCols[1], 1],
                  x1 = 1:length(x) - leftShift + xShift*(wIdx-1), y1 = z[ wIdx, , yIdx, respCols[1], 3],
                  col = cols[wIdx], lwd = 3 )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )
    # plot MS model
    plot( x = c(1-leftShift,length(x)+leftShift), 
          y = yLim, type = "n",
          axes = FALSE )
      axis( side = 1, at = 1:length(x),labels = x )
      axis( side = 2, las = 1 )
      # axis( side = 3, at = 1:length(x), labels = FALSE)
      abline(h = 0, lty = 2, lwd = 0.8)
      if( yIdx == 1 ) mtext( side = 3, text = "Joint Model", cex = 0.8)
      for( wIdx in 1:length(w) )
      {
        lines(  x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[2], 2], col = cols[wIdx],
                lwd = 1, lty = 3 )
        points( x = 1:length(x) - leftShift + xShift*(wIdx-1), y = z[ wIdx, , yIdx, respCols[2], 2 ], col = cols[wIdx],
                pch = 16, cex = 2 )
        segments( x0 = 1:length(x) - leftShift + xShift*(wIdx-1), y0 = z[ wIdx, , yIdx, respCols[2], 1],
                  x1 = 1:length(x) - leftShift + xShift*(wIdx-1), y1 = z[ wIdx, , yIdx, respCols[2], 3],
                  col = cols[wIdx], lwd = 3 )
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[3], " = ", y[yIdx], sep = "" ) )
  }
  mtext( side = 3, outer = T, text = tableName, cex = 1.2 )
  mtext( side = 1, outer =T, text = axes[2], cex = 1.2 )
  mtext( side = 2, outer = T, text = "Relative Error (%)", cex = 1.2, line = 1 )
  mtext( side = 1, outer = T, text = resp, cex = 1, col = "grey50",
          adj = 0.9, line = 1.3 )
  mtext( side = 1, outer = T, text = MP, cex = 1, col = "grey50",
          adj = 0.1, line = 1.3 )
}

plotREtableSlice <- function( tableName = "Fhist_pub_MARE",
                              nSp = 5,
                              axes = c("tUpeak","Umax"),
                              specVals = c( Dover = 0.30,
                                            English = 0.40,
                                            Rock = 0.35,
                                            Petrale = 0.25,
                                            Arrowtooth = 0.15 ),
                              par = "Umsy"
                            )
{
  # plotREtableSlice()
  # Plots MARE/MRE (response variable) values as a function
  # of given specVals (OM values of q, Umsy, etc). 
  # Will group lines corresponding to the first axis and create panel rows 
  # corresponding to the second. Panel columns are single and Joint
  # AMs.
  # inputs:   tableName = RE table in statistics folder
  #           nSp = number of species to restrict to
  #           axes = 2-char vector naming RE table columns
  #           specVals = named numeric of OM quantities to use as plot x axis
  #           par = estimated parameter of interest for MRE/MARE values
  # outputs:  NULL
  # side-eff: plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE ) 

  # restrict to correct number of species
  table  <-   table %>%
              filter( nS == nSp )

  specVals <- specVals[1:nSp]
  specVals <- specVals[order(specVals)]


  # Get levels for axes
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }
  # order levels
  x   <- unique(table %>% pull(axesCols[1]))
  x   <- x[order(x)]
  y   <- unique(table %>% pull(axesCols[2]))
  y   <- y[order(y)]

  mps <- unique( table %>% pull(mp) )

  # Create an array to hold info
  z <- array( NA, dim = c( length(x), length(y), length(specVals), 2, length(mps) ),
              dimnames = list( x, y, names(specVals), c("ss","ms"), mps ) )

  # Fill the array, looping over species and axes
  for(sIdx in 1:length(specVals) )
  {
    spec <- names( specVals )[ sIdx ]
    specTab <- table %>% filter( species == spec )
    for( xIdx in 1:length(x) )
    {
      for( yIdx in 1:length(y) )
      {
        # Get x and y values
        xVal <- x[ xIdx ]
        yVal <- y[ yIdx ]
        # Get rows
        xRows <- which ( specTab[,axesCols[1] ] == xVal)
        yRows <- which ( specTab[,axesCols[2] ] == yVal)
        rows <- intersect(xRows,yRows)
        if( length(rows) == 0 ) next
        # Fill array
        parColNames <- paste(c("ss","ms"),par, sep = "")
        z[xIdx,yIdx,sIdx,,] <- as.numeric(specTab[rows,parColNames])
      }
    }
  }

  cols <- brewer.pal( n = length(x), name = "Dark2" )
  
  # Plot contents. Use x values as line groupings, y values as
  # panel groupings, and species productivities as x values
  par( mfrow = c(length(y), 2), mar = c(2,2,2,1), oma = c(3,3,2,2) )
  for( yIdx in 1:length(y) )
  {
    # Plot SS model
    plot( x = c(1,length(specVals)), 
          y = c( min( z ), max( z ) ), type = "n",
          axes = FALSE )
      if( yIdx == 1 ) panLegend(  legTxt = paste( axes[1], " = ", x ),
                                  lty = 1, pch = 16, cex = 1,
                                  lwd = 0.8, col = cols,
                                  x = 0.1, y = 0.7, bty = "n" )
      if( yIdx == 1 ) mtext( side = 3, text = "Single Species Model", cex = 0.8)
      axis( side = 1, at = 1:length(specVals),labels = specVals )
      axis( side = 2, las = 1 )
      axis( side = 3, at = 1:length(specVals), labels = names(specVals))
      abline(h = 0, lty = 2, lwd = 0.8)
      for( m in 1:length(mps) )
      {
        for( xIdx in 1:length(x) )
              {
                lines(  x = 1:length(specVals), y = z[ xIdx, yIdx, , "ss", m ], col = cols[xIdx],
                        lwd = 0.8, lty = m )
                points( x = 1:length(specVals), y = z[ xIdx, yIdx, , "ss", m ], col = cols[xIdx],
                        pch = 16, cex = 0.5 )
              }
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[2], " = ", y[yIdx], sep = "" ) )

    plot( x = c(1,length(specVals)), 
          y = c( min( z ), max( z ) ), type = "n",
          axes = FALSE )
      axis( side = 1, at = 1:length(specVals),labels = specVals )
      axis( side = 2, las = 1 )
      axis( side = 3, at = 1:length(specVals), labels = names(specVals))
      abline(h = 0, lty = 2, lwd = 0.8)
      if( yIdx == 1 ) mtext( side = 3, text = "Joint Model", cex = 0.8)
      for( m in 1:length(mps) )
      {
        for( xIdx in 1:length(x) )
        {
          lines(  x = 1:length(specVals), y = z[ xIdx, yIdx, , "ms", m ], col = cols[xIdx],
                  lwd = 0.8, lty = m )
          points( x = 1:length(specVals), y = z[ xIdx, yIdx, , "ms", m ], col = cols[xIdx],
                  pch = 16, cex = 0.5 )
        }  
      }
      panLab( x = 0.8, y = 0.8, 
              txt = paste( axes[2], " = ", y[yIdx], sep = "" ) )
  }
  mtext( side = 3, outer = T, text = tableName, cex = 1.2 )
  mtext( side = 1, outer = T, text = par, cex = 1, col = "grey80",
          adj = 0.9, line = 1.5 )
}

plotFhist <- function ( sims = 1:2,
                        nTraj = 60,
                        seed = 2,
                        save = FALSE,
                        fileName = "depPlot.pdf",
                        reverse = FALSE,
                        nS = 3,
                        ... )
{
  # plotFhist()
  # This function plots traces of biomass trajectories
  # for different fishing histories. 
  # This one's a little path specific, and only needs
  # to be pointed to the correct sims
  # inputs      sims=numeric vector indicating sims in project dir
  
  .loadSim(sims[1])

  # Pull parameters
  Bmsy  <- blob$opMod$Bmsy
  nT    <- max( blob$opMod$lYear - blob$opMod$fYear + 1 )
  nReps <- blob$ctrl$signifReps
  if( is.null(nS) )
    nS    <- blob$opMod$nS

  if( reverse ) sIndices <- nS:1
  else sIndices <- 1:nS

  # create an array for holding biomass
  depBlob <- array (  NA, dim = c(length(sims),nTraj,nS,nT),
                      dimnames = list(sims,1:nTraj,blob$ctrl$speciesNames[1:nS],1:nT) )

  Fblob   <- array (  NA, dim = c(length(sims),nS,nT),
                      dimnames = list(sims,blob$ctrl$speciesNames[1:nS],1:nT) ) 

  # Also collect Fhist info
  initDep   <- blob$opMod$initDep
  Umax      <- numeric( length(sims) )

  # set random seed, select traces
  set.seed(seed)

  traces <- sample( 1:nReps, nTraj )

  # Save first sim's biomass and Fhist info
  for( sIdx in 1:nS )
  {
    depBlob[ 1, , sIdx , ]   <- blob$om$Bt[ traces, sIdx, ] / Bmsy[sIdx] / 2
    Fblob[ 1, sIdx,  ]       <- blob$om$Ut[ 1, sIdx, ]
  }
  Umax[ 1 ]         <- blob$opMod$Umult[2]


  

  # Now loop over the other three and collect the info
  for( simIdx in 2:length(sims) )
  {
    .loadSim( sims[simIdx] )
    for( sIdx in 1:nS )
    {
      depBlob[ simIdx, , sIdx , ]   <- blob$om$Bt[ traces, sIdx, ] / Bmsy[sIdx] / 2
      Fblob[ simIdx, sIdx,  ]       <- blob$om$Ut[ 1, sIdx, ]
    }
    Umax[ simIdx ]       <- blob$opMod$Umult[2]
  }

  ordUmax <- unique( Umax ) [ order( unique( Umax  ) ) ]


  nrows <- length( ordUmax )
  ncols <- nS
  # calculate median depletion trajectory
  medDep <- apply ( X = depBlob, FUN = median, MARGIN = c(1,3,4), na.rm = T )

  # If saving to figs folder
  if( save )
  {
    outFile <- file.path(getwd(),"project","figs", fileName )
    pdf( file = outFile, ... )
  } 
  # Make plot labels
  labs <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)" )
  counter <- 1

  # now make x axis year ticks
  years <- seq(1984,2017,by=5)
  years <- c(years,2017)


  par( mfrow = c( nrows, ncols ), mar = c(0,0,0,0), oma = c(5,6,3,5) )
  # Loop over U and t values
  for( U in ordUmax )
  {
    for ( sIdxNum in 1:length(sIndices) )
    { 
      sIdx <- sIndices[sIdxNum]
      sim <- which ( Umax == U )
      plot( x = c(1984,2017), y = c(0,1.2), type = "n", axes = FALSE,
            xlab = "", ylab = "" )
      mfg <- par("mfg")
      if( mfg[1] == mfg[3] )
        axis( side = 1, at = years )
      if( mfg[2] == mfg[4] )
        axis( side = 4, las = 1 )
      if( mfg[2] == 1 )
        axis( side = 2, las = 1 )
      box()
      if( length(sim) == 0 ) next
        for( traj in 1:nTraj ) 
          lines( x = 1984:2017, 
                 y = depBlob[ sim, traj, sIdx, ],
                 col = "grey60", lwd = 0.5 )
        lines( x = 1984:2017, y = medDep[ sim, sIdx, ], lwd = 3, lty = 2 )
        lines( x = 1984:2017, y = Fblob[ sim, sIdx, ], lty = 1, lwd = 3)

        panLab( txt = labs[counter], x = 0.8, y = 0.9, cex = 2 )
      counter <- counter + 1
    }
  }
  mtext( side = 2, outer = TRUE, text = expression(Relative~Biomass~B[t]/B[0]), cex = 1.5, line = 3 )
  mtext( side = 1, outer = TRUE, text = "Year", cex = 1.5, line = 3 )
  mtext( side = 4, outer = TRUE, text = "Harvest Rate", cex = 1.5, line = 3 )
  if( save ) dev.off()
}


plotUPriorSens <- function (  tableName = "priorSens_pub_MRE",
                              nSp = 2,
                              betaU = c(.01,.3),
                              showLegend = FALSE )
{
  # plotUPriorSens()
  # Reads in a stats table for prior sensitivity analyses
  # and plots the performance
  # inputs:     tableName = character name of stats table
  #             par=paremeter that is affected by sens analysis
  # outputs:    NULL
  # side-effs:  plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )
  
  # restrict to the table with the 
  table     <-  table %>%
                # filter( fixqbar == TRUE ) %>%
                filter( nS == nSp ) %>%
                group_by( mp, species) %>%
                summarise(  msUmsy = mean(msUmsy),
                            msUmsy025 = mean(msUmsy025),
                            msUmsy975 = mean(msUmsy975),
                            UmsyOM = mean(UmsyOM),
                            sigU2 = mean(sigU2),
                            sigU2P2 = mean(sigU2P2),
                            Umsybar = mean(Umsybar),
                            Umsybar025 = mean(Ubar025),
                            Umsybar975 = mean(Ubar975),
                            sigU2025 = mean(sigU2025),
                            sigU2975 = mean(sigU2975),
                            ssUmsy = mean(ssUmsy),
                            ssUmsy025 = mean(ssUmsy025),
                            ssUmsy975 = mean(ssUmsy975),
                            s2U       = mean(s2Umsy) ) %>%
                mutate( msUmsy      = UmsyOM * msUmsy + UmsyOM,
                        msUmsy025   = UmsyOM * msUmsy025 + UmsyOM,
                        msUmsy975   = UmsyOM * msUmsy975 + UmsyOM,
                        ssUmsy      = UmsyOM * ssUmsy + UmsyOM,
                        ssUmsy025   = UmsyOM * ssUmsy025 + UmsyOM,
                        ssUmsy975   = UmsyOM * ssUmsy975 + UmsyOM ) %>%
                dplyr::select(  msUmsy, msUmsy025, msUmsy975, 
                                ssUmsy, ssUmsy025, ssUmsy975, 
                                Umsybar, Umsybar025, Umsybar975,
                                mp, species, UmsyOM, 
                                sigU2, sigU2025, sigU2975, s2U,
                                sigU2P2 )

  # plotting quantities
  s2U <- unique(table$s2U)
  s2U <- s2U[ order(s2U) ]

  species     <- unique( table$species )
  nS          <- length( species )

  # Species plotting locations
  specX       <- seq (-0.3,0.3, length=nS )
  fullX       <- c()
  for( k in 1:length(s2U) ) fullX <- c( fullX, k + specX ) 

  # pch for each species
  pchSpec <- c(0,1,2,5,6)
  
  textLab <- c("(a)","(b)","(c)","(d)")

  par( mfrow = c(length(betaU),1), mar = c(1,1,2,1), oma = c(3,3,0,0) )
  for( bIdx in 1:length(betaU) )
  {
    beta <- betaU[bIdx]
    # plotMat <- matrix(NA, nrow=nS*length(s2U), ncol = 7 )
    plotTable <-  table %>% filter( sigU2P2 == beta ) 
    # browser()
    Umsybar <- plotTable %>% group_by(s2U) %>% 
            dplyr::select(s2U,Umsybar,Umsybar025,Umsybar975)

    sigU2 <- plotTable %>% group_by(s2U) %>% 
            dplyr::select(s2U,sigU2,sigU2025,sigU2975)

    Uspec <- plotTable %>% group_by(s2U) %>% 
            dplyr::select(s2U,msUmsy,msUmsy025,msUmsy975,UmsyOM,ssUmsy,ssUmsy025,ssUmsy975 )
    
    plot( x = c( 0, length(s2U) + 1 ), y = c( 0.9*min(plotTable$msUmsy025, plotTable$UmsyOM), 1.1*max(plotTable$msUmsy975, plotTable$UmsyOM) ), 
          type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      polygon(  x = c( fullX, rev(fullX) ),
                y = c(Umsybar$Umsybar025, rev(Umsybar$Umsybar975) ),
                col = "grey90", border = NA )
      lines( x = fullX, y = Umsybar$Umsybar, lty = 2, lwd = 2 )
      abline( h = exp(-1.23), lty = 3, lwd = 2 )
      points( x = fullX, y = Uspec$msUmsy, pch = pchSpec[1:nS], cex = 1.5 )
      points( x = fullX+0.05, y = Uspec$ssUmsy, pch = pchSpec[1:nS], cex = 1.5 )
      points( x = fullX, y = Uspec$UmsyOM, pch = 20, cex = 1.5 )
      segments( x0 = fullX, y0 = as.numeric( Uspec$msUmsy025 ),
                y1 = as.numeric( Uspec$msUmsy975 ), 
                col = "grey40", lwd = 2 )
      segments( x0 = fullX+0.05, y0 = as.numeric( Uspec$ssUmsy025 ),
                y1 = as.numeric( Uspec$ssUmsy975 ), 
                col = "grey40", lwd = 2, lty = 4 )
      panLab( x = 0.05, y = 0.8, txt = textLab[bIdx] )
      if( showLegend & bIdx == 1)
        panLegend(  x = 0.05, y = 0.1, legTxt = c( "q OM", "q Est", "qbar"),
                    bty = "n", pch = c( 1, 16, NA ), 
                    lty = c( NA, NA, 2 ),
                    lwd = c( NA, NA, 2 ),
                    cex = c( 1.5, 1.5, NA ) )
      axis( side = 1, at = 1:length(s2U), labels = s2U )
      axis( side = 2, las = 1 )
      if( bIdx > 1 ) axis( side = 3, at = 1:length(s2U), labels = FALSE )
  }

  mtext ( side = 2, text = "Optimal exploitation rate (Umsy)", las = 0 , outer = TRUE,
          line = 2 )

  mtext ( side = 1, text = "Mean Hyperprior Variance (v_U)", 
          las = 0 , outer = TRUE,
          line = 2 )
}


plotqPriorSens <- function (  tableName = "priorSens_pub_MRE",
                              nSp = 2,
                              betaq = c(.01,.3),
                              showLegend = FALSE )
{
  # plotqPriorSens()
  # Reads in a stats table for prior sensitivity analyses
  # and plots the performance
  # inputs:     tableName = character name of stats table
  #             par=paremeter that is affected by sens analysis
  # outputs:    NULL
  # side-effs:  plots to quartz

  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )

  DERPAcols <- brewer.pal(n = 5, name = "Set1" )

  # restrict to the table with the 
  table <-      table %>%
                # filter( fixqbar == TRUE ) %>%
                filter( nS == nSp, tauq2P2 %in% betaq ) %>%
                group_by( mp, species) %>%
                summarise(  msq = mean(msq),
                            qOM = mean(qOM),
                            tauq2 = mean(tauq2),
                            qbar = mean(qbar),
                            msq025 = mean(msq025),
                            msq975 = mean(msq975),
                            qbar025 = mean(qbar025),
                            qbar975 = mean(qbar975),
                            tauq2025 = mean(tauq2025),
                            tauq2975 = mean(tauq2975),
                            tauq2P2 = mean(tauq2P2),
                            s2q     = mean(s2q),
                            ssq     = mean(ssq),
                            ssq025  = mean(ssq025),
                            ssq975  = mean(ssq975) ) %>%
                mutate( msq     = qOM * msq + qOM,
                        msq025  = qOM * msq025 + qOM,
                        msq975  = qOM * msq975 + qOM,
                        ssq     = qOM * ssq + qOM,
                        ssq025  = qOM * ssq025 + qOM,
                        ssq975  = qOM * ssq975 + qOM ) %>%
                dplyr::select(  msq, msq025, msq975, 
                                qbar, qbar025, qbar975,
                                mp, species, qOM, s2q, tauq2P2,
                                tauq2, tauq2025, tauq2975, species,
                                ssq025, ssq975, ssq )
  # plotting quantities
  s2q <- unique(table$s2q)
  s2q <- s2q[ order(s2q) ]

  species     <- unique( table$species )
  nS          <- length( species )

  # Species plotting locations
  specX       <- seq (-0.3,0.3, length=nS )
  fullX       <- c()
  for( k in 1:length(s2q) ) fullX <- c( fullX, k + specX ) 

  # pch for each species
  pchSpec <- 1:5
  
  textLab <- c("(a)","(b)","(c)","(d)")

  par( mfrow = c(length(betaq),1), mar = c(1,1,2,1), oma = c(3,3,0,0) )
  for( bIdx in 1:length(betaq) )
  {
    beta <- betaq[bIdx]
    # plotMat <- matrix(NA, nrow=nS*length(s2q), ncol = 7 )
    plotTable <-  table %>% filter( tauq2P2 == beta ) 
    # browser()
    qbar <- plotTable %>% group_by(s2q) %>% 
            dplyr::select(s2q,qbar,qbar025,qbar975)

    tauq2 <- plotTable %>% group_by(s2q) %>% 
            dplyr::select(s2q,tauq2,tauq2025,tauq2975)

    qSpec <- plotTable %>% group_by(s2q) %>% 
            dplyr::select(s2q,msq,msq025,msq975,qOM,ssq,ssq025,ssq975 )
    
    plot( x = c( 0, length(s2q) + 1 ), y = c( 0.9*min(plotTable$msq025, plotTable$qOM), 1.1*max(plotTable$msq975, plotTable$qOM) ), 
          type = "n",
          axes = FALSE, xlab = "", ylab = "" )
      polygon(  x = c( fullX, rev(fullX) ),
                y = c(qbar$qbar025, rev(qbar$qbar975) ),
                col = "grey90", border = NA )
      lines( x = fullX, y = qbar$qbar, lty = 2, lwd = 2 )
      abline( h = exp(-0.51), lty = 3, lwd = 2 )
      points( x = fullX, y = qSpec$msq, pch = pchSpec[1:nS], cex = 1.5,
              col = DERPAcols )
      points( x = fullX+0.05, y = qSpec$ssq, pch = pchSpec[1:nS], cex = 1.5,
              col = DERPAcols )
      points( x = fullX, y = qSpec$qOM, pch = 20, cex = 1.5, col = DERPAcols )
      segments( x0 = fullX, y0 = as.numeric( qSpec$msq025 ),
                y1 = as.numeric( qSpec$msq975 ), 
                col = DERPAcols, lwd = 2 )
      segments( x0 = fullX+0.05, y0 = as.numeric( qSpec$ssq025 ),
                y1 = as.numeric( qSpec$ssq975 ), 
                col = DERPAcols, lwd = 2, lty = 4 )
      panLab( x = 0.05, y = 0.8, txt = textLab[bIdx] )
      if( showLegend & bIdx == 1)
        panLegend(  x = 0.05, y = 0.1, legTxt = c( "q OM", "q Est", "qbar"),
                    bty = "n", pch = c( 1, 16, NA ), 
                    lty = c( NA, NA, 2 ),
                    lwd = c( NA, NA, 2 ),
                    cex = c( 1.5, 1.5, NA ) )
      axis( side = 1, at = 1:length(s2q), labels = s2q )
      axis( side = 2, las = 1 )
      if( bIdx > 1 ) axis( side = 3, at = 1:length(s2q), labels = FALSE )
  }

  mtext ( side = 2, text = "Catchability", las = 0 , outer = TRUE,
          line = 2 )

  mtext ( side = 1, text = "Mean Hyperprior Variance (v_q)", 
          las = 0 , outer = TRUE,
          line = 2 )
}


plotObsRasters <- function  ( tableName = "obsErr_pub_MARE",
                              pars = c("BnT","Umsy","q","Dep"),
                              wdh = 18, hgt = 18,
                              breakRange = c(-1,1),
                              compMean = TRUE,
                              hess = FALSE,
                              ... 
                            )
{
  # plotObsRasters()
  # Reads in a statistics table produced by .statTableXXX()
  # and produces plots of performance contours.
  # inputs:     tableName=charactre vector of file name root
  #             pars=character vector of parameters to plot contours for
  #             wdh=numeric width (inches) of plot
  #             hgt=numeric height (inches) of plot
  #             breakRange=range of colour breaks for raster plots
  #             compMean=Switch for computing the mean (TRUE) or 
  #                       sum (FALSE) for the complex contour
  #             ...=args to pass to plotObsCompContour()
  # outputs:    NA
  # side eff:   Prints pdf plots of each MP and nS combination
  #             to getwd()/project/figs/<tableName>/
  
  # Load table
  fileName <- paste( tableName, ".csv", sep = "" )
  tablePath <- file.path ( getwd(), "project/Statistics", fileName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )

  tabFolder <- file.path(getwd(),"project","figs",tableName )
  dir.create(tabFolder)
  # Now get the list of MPs
  MPs <- unique( table$mp )
  # browser()
  # Loop over MPs, and plot one set of raster fields for
  # each of the other levels
  for (mp in MPs)
  {
    mpFolder  <- file.path(tabFolder,mp)
    dir.create(mpFolder)
    mpTab     <- table %>% filter( mp == mp )
    nS        <- unique(mpTab$nS)
    levels    <- list (nS = nS )
    grid      <- expand.grid(levels)
    # browser()
    for( g in 1:nrow(grid) )
    {
      # Set up output file
      outFile <- paste( grid$nS[g], 
                        mp,
                        "obsErrDelta.pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsCompContour( tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          tcex = 1.5,
                          breakRange = breakRange,
                          compMean = compMean,
                          ... )
      dev.off()
      outFile <- paste( grid$nS[g], 
                        mp,
                        "obsErr_MS.pdf", sep = "")
      outFile <- file.path(mpFolder,outFile)
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsModContour(  tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          model = "ms", tcex = 1.5,
                          breakRange = breakRange,
                          compMean = compMean,
                          hess = hess )
      dev.off()
      outFile <- paste( grid$nS[g], mp,
                        "obsErr_SS.pdf", sep = "")
      outFile <- file.path( mpFolder, outFile )
      pdf ( file = outFile, width = wdh, height = hgt )
      plotObsModContour(  tableName = fileName,
                          mpLabel = mp,
                          pars = pars,
                          nSp = grid[g,"nS"],
                          model = "ss", tcex =1.5,
                          breakRange = breakRange,
                          compMean = compMean,
                          hess = hess )
      dev.off()
    }
  }
  # No return
}


plotObsModContour <- function (   tableName   = "obsErr_MARE.csv", 
                                  pars        = c("BnT", "Umsy", "q", "Dep"),
                                  nSp         = 2,
                                  method      = "deviance",
                                  compMean    = FALSE,                                
                                  tcex        = 1,
                                  breakRange  = c(-1,1),
                                  mpLabel     = "corrTRUE_SigPriorIG",
                                  model       = "ss",
                                  devLabels   = FALSE,
                                  hess        = FALSE 
                              )
{
  # plotObsModContour()
  # Plots the observation model experiment MSE contours for
  # estimation performance of given pars by species
  # and for the complex.
  # inputs:     tableName=charactre vector of file name root
  #             pars=character vector of parameters to plot contours for
  #             nSp=number of species (scenario selector)
  #             method=method for comparing ms and ns models,
  #                     deviance => log2(ss/ms),
  #                     difference => abs(ss) - abs(ms)
  #             compMean=compute mean for complex (TRUE), or sum (FALSE)?
  #             tcex=cex for text calls (prints numbers over raster cells)
  #             breakRange=range of colour breaks for raster plots
  #             mpLabel=character vector selecting MP
  #             model=single species "ss" or multi-species "ms"
  # outputs:    NA
  # side eff:   Prints plots of model estimation performance contours for
  #             selected pars and MP

  # Load stat table
  tablePath <- file.path ( getwd(), "project/Statistics", tableName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )
  # reduce to the correct MP
  table <-  table  %>% 
            filter( mp == mpLabel, nS == nSp ) %>%
            mutate( CV = round( sqrt( exp( tau2OM ) -1 ), digits = 2 ) )

  # browser()
  
  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=nS)

  # Axis ticks:
  CVs   <- unique( table[ , "CV" ] )
  CVs   <- CVs[ order( CVs ) ]
  # scenario name values
  tau2  <- unique( table[ , "tau2OM" ] ) 
  tau2  <- tau2[ order( tau2 ) ]

  if( hess ) pars <- c( pars, "HessPD" )
 
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim=c(length(CVs),length(CVs),length(pars)),
                dimnames=list(CVs[length(CVs):1],CVs,pars) )

    for (xIdx in 1:length(CVs))
    {
      for (yIdx in 1:length(CVs))
      {
        # Get x and y values
        xVal <- tau2[xIdx]
        yVal <- tau2[yIdx]
        # Get rows
        scenName <- paste( "nS", nS, "_CVc(", xVal, ",", yVal, ")", sep = "" )
        # Get rows
        rows <- which ( specTab[,"scenario"] == scenName )
        # browser()
        if( length(rows) == 0 ) next
        colNames <- paste(model,pars,sep="")
        z[length(CVs)-yIdx+1,xIdx,] <- as.numeric(specTab[rows,colNames])
      }
    } 
    # browser()
    specList[[s]] <- brick( z, xmn = 0, xmx = length(CVs), ymn = 0, ymx=length(CVs) )
  }

  names(specList) <- species

  par(mfrow = c((nS+1),length(pars) ),  mar = c(1,1,1,1), oma = c(5,5,2,2) )
  for (s in 1:nS)
  {
    # Create a raster from the stat table info
    plotBrick <- specList[[s]]
    # Make titles for the plots
    specName <- species[s]

    titles <- names(plotBrick)
    l <- length(titles)
    if( !hess ) l <- l + 1

    colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
    colBreaksHess <- seq(0,200,length = 8)
    colScale  <- brewer.pal(8, "Reds")
    colScale      <- c(colScale[8:1],colScale)
    
    # plot the rasters, if first species plot main titles
    if ( s == 1 )
    {
      for(k in 1:(l-1))
      {
        plotBrick[[ k ]][plotBrick[[k]] > 10 ] <- 99
        plot( plotBrick[[ k ]],main = titles[ k ], legend = FALSE, 
              asp=NA, ylab = "", xlab = "", axes = FALSE, 
              col = colScale, breaks = colBreaks )
          if( k == 1 ) mtext ( side = 2, text = specName, line = 2, las = 0)
          text( plotBrick[[ k ]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      }
      if( hess )
      {
        plot( plotBrick[[ l ]],main = titles[ l ], legend = FALSE, 
              asp=NA, ylab = "", xlab = "", axes = FALSE, 
              col = colScale[9:16], breaks = colBreaksHess )
          text( plotBrick[[ l ]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )  
      }
    }
    else {
      # browser()
      for(k in 1:( l - 1 ) )
      {
        plotBrick[[ k ]][plotBrick[[k]] > 10 ] <- 99
        plot( plotBrick[[k]],main = "", legend = FALSE, 
              asp=NA, ylab = "", xlab = "", axes = FALSE, 
              col = colScale, breaks = colBreaks )
          text( plotBrick[[k]], digits = 2, halo = T, cex = tcex ) 
          if( k == 1 ) mtext ( side = 2, text = specName, line = 2, las = 0)
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      }
      if( hess )
      {
        plot( plotBrick[[l]],main = "", legend = FALSE, 
                    asp=NA, ylab = "", xlab = "", axes = FALSE, 
                    col = colScale[9:16], breaks = colBreaksHess )
          text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
          axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      }
    }   
  }

  # browser()
  # Now create an average/sum raster brick for the complex
  # fudge factor for hess logical flag
  m <- l
  if( !hess ) m <- l-1
  for(k in 1:m)
  {
    # use leftover plotBrick in memory to compute complex
    plotBrick[[k]] <- 0
    for( s in 1:length(species) )
    {
      plotBrick[[k]] <- plotBrick[[k]] + specList[[s]][[k]]
    }
    if( compMean ) plotBrick[[k]] <- plotBrick[[k]]/length(species)
  }

  titles <- names(plotBrick)
  
  # browser()
  
  for(k in 1:m )
  {
    plotBrick[[ k ]][plotBrick[[k]] > 10 ] <- 99
    plot( plotBrick[[k]],main = "", legend = FALSE, 
          asp=NA, ylab="", xlab = "", axes = FALSE, 
          col = colScale, breaks = colBreaks )
      if( k == 1 ) mtext ( side = 2, text = "Complex", line = 2, las = 0)
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )
      axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
      axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )   
  }
  if( hess )
  {
    plot( plotBrick[[l]],main = "", legend = FALSE, 
          asp=NA, ylab="", xlab = "", axes = FALSE,
          col = colScale[9:16], breaks = nSp*colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex )   
      axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
      axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )   
  }
  

  mtext( side = 1, text = "Dover Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  mtext( side = 2, text = "English Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  if( devLabels )
  {
    grid.text(  paste("Raw error ", model, " model.", sep = ""), 
                x=unit(0.5, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  tableName,
                x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  mpLabel,
                x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)  
  }
  
  # grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
              # x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)
}


# plotObsCompContour()
# Plots the correlation comparative experiment performance contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex. The plots are panel plots of Rasters, so need
# a little love as far as sizing goes.
# inputs: table=character name of table generated by .statTableXXX
#         mpLabel=character name of MP label
#         axes=character 2-vector for axes of raster plots
#         nSp=integer number of species in scenario
#         method="deviance" for log2(SS/MS), difference for SS-MS
#         compMean=compute the mean for the complex?
#         tUtr=integer value for tUtrough in Fhist scenarios
#         tcex=cex of text in raster cells
plotObsCompContour <- function (  tableName   = "obsErr_MARE.csv", 
                                  pars        = c("BnT", "Umsy", "q", "Dep"),
                                  nSp         = 2,
                                  method      = "deviance",
                                  compMean    = FALSE,                                
                                  tcex        = 2,
                                  breakRange  = c(-2,2),
                                  mpLabel     = "corrTRUE_SigPriorIG",
                                  legend      = FALSE,
                                  devLabels   = FALSE 
                                  )
{
  # Load stat table
  tablePath <- file.path ( getwd(),"project/Statistics",tableName)
  table <- read.csv (tablePath, header=TRUE, stringsAsFactors=FALSE)
  # reduce to the correct MP
  table <- table %>% filter(mp == mpLabel )

  # Filter by other 
  table <- table %>% 
           filter( nS == nSp ) %>%
           mutate ( CV = round( sqrt( exp(tau2OM) - 1 ), digits = 2 ) )


  # browser()
  # Now make a pallette for the colours
  colScaleHi    <- brewer.pal( 8, "Greens" )
  colScaleLo    <- brewer.pal( 8, "Reds" )
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colScale      <- c(colScaleLo[8:1],colScaleHi)

  # calculate response
  # Deviance is log2(ss/ms)
  if (method == "deviance" )
  {
    table <- table %>% mutate(  BnT     = log2(abs(ssBnT/msBnT)),
                                Umsy    = log2(abs(ssUmsy/msUmsy)),
                                q       = log2(abs(ssq/msq)),
                                Dep     = log2(abs(ssDep/msDep)),
                                Bmsy    = log2(abs(ssBmsy/msBmsy)) )
  }   
  # Difference is just ss - ms
  if ( method == "difference")
  {
    table <- table %>% mutate(  BnT     = abs(ssBnT) - abs(msBnT),
                                Umsy    = abs(ssUmsy) - abs(msUmsy),
                                q       = abs(ssq) - abs(msq),
                                Dep     = abs(ssDep) - abs(msDep),
                                Bmsy    = log2(abs(ssBmsy/msBmsy)) )
  }

  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=length(species))

  CVs   <- unique( table[ , "CV" ] )
  CVs   <- CVs[ order( CVs ) ]

  tau2  <- unique( table[ , "tau2OM" ] ) 
  tau2  <- tau2[ order( tau2 ) ]
  # browser()
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim = c( length( CVs ), length( CVs ), length(pars) ),
                dimnames=list(  CVs[ length( CVs ):1 ], CVs,
                                pars ) )

    for (xIdx in 1:length(CVs))
    {
      for (yIdx in 1:length(CVs))
      {
        # Get x and y values
        xVal <- tau2[xIdx]
        yVal <- tau2[yIdx]
        scenName <- paste( "nS", nS, "_CVc(", xVal, ",", yVal, ")", sep = "" )
        # Get rows
        rows <- which ( specTab[,"scenario"] == scenName )
        # browser()
        if (length(rows) == 0 ) next
        # Recover MSE comparison
        z[ length( CVs ) - yIdx + 1, xIdx, ] <- as.numeric(specTab[rows,pars])
      }
    } 
    specList[[s]] <- brick( z, xmn = 0, xmx = length(CVs), ymn = 0, ymx=length(CVs) )
  }
  names(specList) <- species

  # Start plotting

  par(mfrow = c(nS+1,length(pars)), mar = c(1.5,1.5,1.5,1.5), oma = c(5,5,2,2) )
  for (s in 1:nS)
  {
    # Create a raster from the stat table info 
    plotBrick <- specList[[s]]
    # Make titles for the plots
    titles <- names(plotBrick)
    # browser()
    # plot the rasters, without legends
    for ( k in 1:length(titles) )
    {
      plotBrick[[ k ]][plotBrick[[k]] > 100 ] <- 99
      plotBrick[[ k ]][plotBrick[[k]] < -100 ] <- -99 
      if(s == 1) plot( plotBrick[[k]],main = titles[k], legend = FALSE, 
                        breaks = colBreaks, col=colScale, asp=NA, axes = FALSE )
      else plot( plotBrick[[k]],main = "", legend = FALSE, axes = FALSE, 
                  breaks = colBreaks, col=colScale, asp=NA )
      if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
      axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
      axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
    }
  }

  # browser()
  # Now create the complex plot brick - idea: group by scenario
  # and MP, there should be a unique set of species for each, then
  # apply the complex rule
  compTable <- table %>%  group_by( scenario, mp ) %>%
                          summarise(  ssBnT = sum( ssBnT ),
                                      msBnT = sum( msBnT ),
                                      ssUmsy = sum( ssUmsy ),
                                      msUmsy = sum( msUmsy ),
                                      ssBmsy = sum( ssBmsy ),
                                      msBmsy = sum( msBmsy ),
                                      ssq = sum( ssq ),
                                      msq = sum( msq ),
                                      ssHessPD = sum( ssHessPD ),
                                      msHessPD = sum( msHessPD ),
                                      # ssMSY = sum( ssMSY ),
                                      # msMSY = sum( msMSY ),
                                      ssDep = sum( ssDep ),
                                      msDep = sum( msDep ),
                                      corr = mean(corr),
                                      kappaMult = mean(kappaMult),
                                      tUpeak = mean(tUpeak),
                                      tUtrough = mean(tUtrough),
                                      Umax = mean(Umax),
                                      tau2OM = mean (tau2OM),
                                      CV = mean ( CV ) )

  # Now make the comparison based on the compMean logical flag
  if( method == "deviance" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = log2( abs(ssBnT) / abs(msBnT) ),
                          Umsy = log2( abs(ssUmsy) / abs(msUmsy) ),
                          q = log2( abs(ssq) / abs(msq) ),
                          Bmsy = log2( abs(ssBmsy) / abs(msBmsy) ),
                          Dep = log2 ( abs(ssDep) / abs(msDep) ) 
                        )
  }

  if( method == "difference" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = abs(ssBnT) - abs(msBnT) ,
                          Umsy = abs(ssUmsy) - abs(msUmsy) ,
                          q = abs(ssq) - abs(msq) ,
                          Bmsy = abs(ssBmsy) - abs(msBmsy) ,
                          Dep = abs(ssDep) - abs(msDep) 
                        )
    if( compMean )
      compTable <-  compTable %>% 
                    mutate( BnT = BnT/nS,
                            Umsy = Umsy/nS,
                            q = q/nS,
                            # MSY = MSY/nS,
                            Dep = Dep/nS )
  }

  # Make a raster brick with the selected pars
  z <- array( NA, dim = c( length( CVs ), length( CVs ), length(pars) ),
              dimnames=list(  CVs[ length(CVs):1], CVs,
                              pars ) )

  # browser()
  for (xIdx in 1:length(CVs))
  {
    for (yIdx in 1:length(CVs))
    {
      # Get x and y values
      xVal <- tau2[xIdx]
      yVal <- tau2[yIdx]
      # Get rows
      scenName <- paste( "nS", nS, "_CVc(", xVal, ",", yVal, ")", sep = "" )
      # Get rows
      rows <- which ( compTable[,"scenario"] == scenName )
      if( length(rows) == 0 ) next
      # Recover MSE comparison
      z[length(CVs) - yIdx + 1,xIdx,] <- as.numeric(compTable[rows,pars])
    }
  } 
  compBrick <- brick( z, xmn = 0, xmx = length(CVs), ymn = 0, ymx=length(CVs) )

  for(k in 1:length(titles))
  {
    compBrick[[ k ]][compBrick[[k]] > 10 ] <- 99
    compBrick[[ k ]][compBrick[[k]] < -10 ] <- -99
    plot( compBrick[[k]],main = "", legend = FALSE, 
          asp=NA, col = colScale, breaks = colBreaks,
          axes = FALSE )
      if( k == 1 ) mtext ( text = "Complex", side = 2, las = 0, line = 2 )
    text( compBrick[[k]], digits = 2, halo = T, cex = tcex )
    axis( side = 2, at = 1:length(CVs) - 0.5, labels = CVs, las = 1 )
    axis( side = 1, at = 1:length(CVs) - 0.5, labels = CVs )
  }

  if( legend )
  {
    # Plot the legends
    plot( plotBrick[[1]], legend.only=TRUE, col=colScale, fill=colScale,
          legend.width=1, legend.shrink=0.75,
          axis.args = list( at = colBreaks,
                            labels = round(colBreaks,2),
                            cex.axis = 0.6),
          zlim = c(-3,3),
          legend.args = list( text="log2(ssMSE/msMSE)", 
                              side=4, font=2, 
                              line=3, cex=0.75, las=0) )
  
  }
  
  # Labels
  mtext( side = 1, text = "Dover Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  mtext( side = 2, text = "English Sole Survey CV", outer = TRUE, las = 0,
          line = 3 )
  if( devLabels )
  {
    grid.text(  "Relative errors ss/ms", 
                x=unit(0.5, "npc"), y=unit(0.95, "npc"), rot=0)
    grid.text(  tableName,
                x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  mpLabel,
                x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
  }
  # grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
  # grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
  #             x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)

}

plotComplexPubRasters <- function(  tableRoot = "info_qREs_RP",
                                    axes = c("nDiff","nS"),
                                    pars = c("Umsy","BnT","Dep"),
                                    wdh = 14, hgt = 18,
                                    breakRange = c(-1,1),
                                    tcex = 1.5,
                                    legend = FALSE,
                                    mpLabel   = "allJointPriors",
                                    hess = FALSE,
                                    nSp       = 5,
                                    tUtr      = 15,
                                    negCorr   = FALSE,
                                    devLabels = FALSE,
                                    save      = FALSE
                                  )
{
  # Load in MARE and MRE table
  browser()
  MAREfile  <- paste( tableRoot, "_MARE.csv", sep = "" )
  MREfile   <- paste( tableRoot, "_MRE.csv", sep = "" )

  MAREpath  <- file.path ( getwd(), "project/Statistics", MAREfile )
  MREpath   <- file.path ( getwd(), "project/Statistics", MREfile )

  MREtable  <- read.csv( MREpath, header = TRUE, stringsAsFactors = FALSE )
  MAREtable <- read.csv( MAREpath, header = TRUE, stringsAsFactors = FALSE )

  # reduce to the correct subset of sims
  MREtable <- MREtable  %>% 
              filter( mp == mpLabel,
                      nS == nSp,
                      tUtrough == tUtr,
                      lastNegCorr == negCorr ) %>%
              group_by( scenario, mp ) %>%
              summarise( ssBnT = mean( ssBnT ),
                         msBnT = mean( msBnT ),
                         ssUmsy = mean( ssUmsy ),
                         msUmsy = mean( msUmsy ),
                         ssq = mean( ssq ),
                         msq = mean( msq ),
                         ssDep = mean( ssDep ),
                         msDep = mean( msDep ),
                         corr = mean(corr),
                         kappaMult = mean(kappaMult) )
              
  
  MAREtable <-  MAREtable  %>% 
                filter( mp == mpLabel,
                        nS == nSp,
                        tUtrough == tUtr,
                        lastNegCorr == negCorr ) %>%
                group_by( scenario, mp ) %>%
                summarise(  ssBnT = mean( ssBnT ),
                            msBnT = mean( msBnT ),
                            ssUmsy = mean( ssUmsy ),
                            msUmsy = mean( msUmsy ),
                            ssq = mean( ssq ),
                            msq = mean( msq ),
                            ssDep = mean( ssDep ),
                            msDep = mean( msDep ),
                            corr = mean(corr),
                            kappaMult = mean(kappaMult) )
  # Compute Delta values              
  MAREtable <-  MAREtable %>%
                mutate( BnT     = log2(abs(ssBnT/msBnT)),
                        Umsy    = log2(abs(ssUmsy/msUmsy)),
                        q       = log2(abs(ssq/msq)),
                        Dep     = log2(abs(ssDep/msDep)) )

  # Now start pulling rasters
  msColNames      <- paste( "ms", pars, sep = "" )
  ssColNames      <- paste( "ss", pars, sep = "" )
  DeltaColNames   <- pars

  # use makeRasterFromTable() to create raster bricks
  ssMRE <- makeRasterFromTable( MREtable, axes, ssColNames )
  msMRE <- makeRasterFromTable( MREtable, axes, msColNames )
  ssMARE <- makeRasterFromTable( MAREtable, axes, ssColNames )
  msMARE <- makeRasterFromTable( MAREtable, axes, msColNames )
  Delta <- makeRasterFromTable( MAREtable, axes, DeltaColNames )

  if( save )
  {
    outFile <- paste( tableRoot, "complexPubRasters.pdf", sep = "")
    outFile <- file.path("./project/figs/",outFile)
    pdf ( file = outFile, width = wdh, height = hgt )
  }
  
  par(mfrow = c(5,length(pars)), mar = c(1,1,1,1), oma = c(5,5,3,2) ) 
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colBreaksHess <- seq(0,200,length = 8)
  redScale    <- brewer.pal(8, "Reds")
  greenScale  <- brewer.pal(8, "Greens" )
  rgScale     <- c(redScale[8:1],greenScale)
  rrScale     <- c(redScale[8:1],redScale)

  titles <- c("B2017","Umsy","q","D2017")

  # First plot Delta
  for(k in 1:length(pars) )
  {
    plot( Delta[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rgScale, breaks = colBreaks,
        axes = FALSE, cex.main = 1.5 )
      mtext( text = titles[k], side = 3, line = 2, cex = 1.5)
      if( k == 1 ) mtext ( text = "(a)", side = 2, las = 1, line = 2 )
      text( Delta[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(Delta$y) - 0.5, labels = ssMARE$y, las = 1 )
      axis( side = 1, at = 1:length(Delta$x) - 0.5, labels = Delta$x ) 
  }  
  # Then MARE plots
  for(k in 1:length(pars) )
  {
    plot( ssMARE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(b)", side = 2, las = 1, line = 2 )
      text( ssMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(ssMARE$y) - 0.5, labels = ssMARE$y, las = 1 )
      axis( side = 1, at = 1:length(ssMARE$x) - 0.5, labels = ssMARE$x ) 
  }  
  for(k in 1:length(pars) )
  {
    plot( msMARE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(c)", side = 2, las = 1, line = 2 )
      text( msMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(msMARE$y) - 0.5, labels = msMARE$y, las = 1 )
      axis( side = 1, at = 1:length(msMARE$x) - 0.5, labels = msMARE$x ) 
  }  
  # Then MRE plots
  for(k in 1:length(pars) )
  {
    plot( ssMRE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(d)", side = 2, las = 1, line = 2 )
      text( ssMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(ssMRE$y) - 0.5, labels = ssMRE$y, las = 1 )
      axis( side = 1, at = 1:length(ssMRE$x) - 0.5, labels = ssMRE$x ) 
  }  
  for(k in 1:length(pars) )
  {
    plot( msMRE[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
        axes = FALSE )
      if( k == 1 ) mtext ( text = "(e)", side = 2, las = 1, line = 2 )
      text( msMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(msMRE$y) - 0.5, labels = msMRE$y, las = 1 )
      axis( side = 1, at = 1:length(msMRE$x) - 0.5, labels = msMRE$x ) 
  }  
  mtext(  side = 1, outer = TRUE, text = "Correlation in species effects",
          line = 2 )
  mtext(  side = 2, outer = TRUE, text = "Relative magnitude of shared effects",
          line = 3 )
  if( save ) dev.off()
}


plotComplexPubRasters2 <- function(  tableRoot = "RE_coarse_pub",
                                    axes = c("corr","kappaMult"),
                                    pars = c("BnT","Bmsy","Umsy","Dep"),
                                    wdh = 14, hgt = 18,
                                    breakRange = c(-1,1),
                                    tcex = 1.5,
                                    legend = FALSE,
                                    mpLabel   = "pubAM",
                                    hess = FALSE,
                                    nSp       = 5,
                                    tUtr      = 15,
                                    negCorr   = FALSE,
                                    devLabels = FALSE,
                                    save      = FALSE
                                  )
{
  # Load in MARE and MRE tables
  # browser()
  MAREfile  <- paste( tableRoot, "_MARE.csv", sep = "" )
  MREfile   <- paste( tableRoot, "_MRE.csv", sep = "" )

  MAREpath  <- file.path ( getwd(), "project/Statistics", MAREfile )
  MREpath   <- file.path ( getwd(), "project/Statistics", MREfile )

  MREtable  <- read.csv( MREpath, header = TRUE, stringsAsFactors = FALSE )
  MAREtable <- read.csv( MAREpath, header = TRUE, stringsAsFactors = FALSE )

  # reduce to the correct subset of sims
  MREtable <- MREtable  %>% 
              filter( mp == mpLabel,
                      nS == nSp,
                      tUtrough == tUtr,
                      lastNegCorr == negCorr ) %>%
              group_by( scenario, mp ) %>%
              summarise( ssBnT = mean( ssBnT ),
                         msBnT = mean( msBnT ),
                         ssBmsy = mean( ssBmsy ),
                         msBmsy = mean( msBmsy ),
                         ssUmsy = mean( ssUmsy ),
                         msUmsy = mean( msUmsy ),
                         ssq = mean( ssq ),
                         msq = mean( msq ),
                         ssDep = mean( ssDep ),
                         msDep = mean( msDep ),
                         corr = mean(corr),
                         kappaMult = mean(kappaMult) )
              
  
  MAREtable <-  MAREtable  %>% 
                filter( mp == mpLabel,
                        nS == nSp,
                        tUtrough == tUtr,
                        lastNegCorr == negCorr ) %>%
                group_by( scenario, mp ) %>%
                summarise(  ssBnT = mean( ssBnT ),
                            msBnT = mean( msBnT ),
                            ssBmsy = mean( ssBmsy ),
                            msBmsy = mean( msBmsy ),
                            ssUmsy = mean( ssUmsy ),
                            msUmsy = mean( msUmsy ),
                            ssq = mean( ssq ),
                            msq = mean( msq ),
                            ssDep = mean( ssDep ),
                            msDep = mean( msDep ),
                            corr = mean(corr),
                            kappaMult = mean(kappaMult) )
  # Compute Delta values              
  MAREtable <-  MAREtable %>%
                mutate( BnT     = log2(abs(ssBnT/msBnT)),
                        Umsy    = log2(abs(ssUmsy/msUmsy)),
                        Bmsy    = log2(abs(ssBmsy/msBmsy)),
                        q       = log2(abs(ssq/msq)),
                        Dep     = log2(abs(ssDep/msDep)) )

  # Now start pulling rasters
  msColNames      <- paste( "ms", pars, sep = "" )
  ssColNames      <- paste( "ss", pars, sep = "" )
  DeltaColNames   <- pars

  # use makeRasterFromTable() to create raster bricks
  ssMRE <- makeRasterFromTable( MREtable, axes, ssColNames )
  msMRE <- makeRasterFromTable( MREtable, axes, msColNames )
  ssMARE <- makeRasterFromTable( MAREtable, axes, ssColNames )
  msMARE <- makeRasterFromTable( MAREtable, axes, msColNames )
  Delta <- makeRasterFromTable( MAREtable, axes, DeltaColNames )

  if( save )
  {
    outFile <- paste( tableRoot, "complexPubRasters.pdf", sep = "")
    outFile <- file.path("./project/figs/",outFile)
    pdf ( file = outFile, width = wdh, height = hgt )
  }
  
  par(mfrow = c(1,length(pars)), mar = c(1,1,1,1), oma = c(5,5,3,2) ) 
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colBreaksHess <- seq(0,200,length = 8)
  redScale    <- brewer.pal(8, "Reds")
  greenScale  <- brewer.pal(8, "Greens" )
  rgScale     <- c(redScale[8:1],greenScale)
  rrScale     <- c(redScale[8:1],redScale)

  titles <- c("B2017","B0","Umsy","D2017")

  # First plot Delta
  for(k in 1:length(pars) )
  {
    plot( Delta[[1]][[k]], legend = FALSE, 
        asp=NA, ylab="", xlab = "", col = rgScale, breaks = colBreaks,
        axes = FALSE, cex.main = 1.5 )
      mtext( text = titles[k], side = 3, line = 2, cex = 1.5)
      # if( k == 1 ) mtext ( text = "(a)", side = 2, las = 1, line = 2 )
      text( Delta[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
      axis( side = 2, at = 1:length(Delta$y) - 0.5, labels = ssMARE$y, las = 1 )
      axis( side = 1, at = 1:length(Delta$x) - 0.5, labels = Delta$x ) 
  }  
  # # Then MARE plots
  # for(k in 1:length(pars) )
  # {
  #   plot( ssMARE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(b)", side = 2, las = 1, line = 2 )
  #     text( ssMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(ssMARE$y) - 0.5, labels = ssMARE$y, las = 1 )
  #     axis( side = 1, at = 1:length(ssMARE$x) - 0.5, labels = ssMARE$x ) 
  # }  
  # for(k in 1:length(pars) )
  # {
  #   plot( msMARE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(c)", side = 2, las = 1, line = 2 )
  #     text( msMARE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(msMARE$y) - 0.5, labels = msMARE$y, las = 1 )
  #     axis( side = 1, at = 1:length(msMARE$x) - 0.5, labels = msMARE$x ) 
  # }  
  # # Then MRE plots
  # for(k in 1:length(pars) )
  # {
  #   plot( ssMRE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(d)", side = 2, las = 1, line = 2 )
  #     text( ssMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(ssMRE$y) - 0.5, labels = ssMRE$y, las = 1 )
  #     axis( side = 1, at = 1:length(ssMRE$x) - 0.5, labels = ssMRE$x ) 
  # }  
  # for(k in 1:length(pars) )
  # {
  #   plot( msMRE[[1]][[k]], legend = FALSE, 
  #       asp=NA, ylab="", xlab = "", col = rrScale, breaks = colBreaks,
  #       axes = FALSE )
  #     if( k == 1 ) mtext ( text = "(e)", side = 2, las = 1, line = 2 )
  #     text( msMRE[[1]][[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
  #     axis( side = 2, at = 1:length(msMRE$y) - 0.5, labels = msMRE$y, las = 1 )
  #     axis( side = 1, at = 1:length(msMRE$x) - 0.5, labels = msMRE$x ) 
  # }  
  mtext(  side = 1, outer = TRUE, text = "Inter-species Correlation",
          line = 2 )
  mtext(  side = 2, outer = TRUE, text = "Shared Environmental Effects",
          line = 3 )
  if( save ) dev.off()
}


pull <- function(x,y) 
{
  # pull()
  # Pulls a single column from a tbl_df object and returns
  # it as a vector. Usefule for working with dplyr tables
  # inputs:   x = tbl_df object
  #           y = column name/number
  x[,if(is.name(substitute(y))) deparse(substitute(y)) 
      else y, drop = FALSE][[1]]
}

# makeRasterFromTable()
# Creates a raster brick from a table, which is either a stat
# table created in .statTableXXX() or a pared down version 
# inside a function
# inputs:   table=data.frame of stats
#           axes=variables that make up the x and y axes of the raster
#           colnames=names of brick layers
# ouputs:   raster brick made up of table entries
makeRasterFromTable <- function ( table, 
                                  axes, 
                                  colnames )
{
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }

  x   <- unique(table %>% pull(axesCols[1]))
  x   <- x[order(x)]
  y   <- unique(table %>% pull(axesCols[2]))
  y   <- y[order(y)]


  # First, create an array to turn into the raster
  z <- array( NA, dim=c(length(y),length(x),length(colnames)),
                dimnames = list( y[ length( y ):1 ], x, colnames ) )
  for (xIdx in 1:length(x))
  {
    for (yIdx in 1:length(y))
    {
      # Get x and y values
      xVal <- x[xIdx]
      yVal <- y[yIdx]
      # Get rows
      xRows <- which ( table[,axesCols[1] ] == xVal)
      yRows <- which ( table[,axesCols[2] ] == yVal)
      rows <- intersect(xRows,yRows)
      if( length(rows) == 0 ) next
      # browser()
      z[length(y)-yIdx+1,xIdx,] <- apply(X = as.matrix(table[rows,colnames]), FUN = mean, MARGIN = 2, na.rm = T)
    }
  }
  # create raster brick
  raster <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )  
  # output raster brick and axis vectors
  out <- list()
  out$brick <- raster
  out$x <- x
  out$y <- y

  out
}

# plotModContour()
# Plots the correlation  experiment MSE contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex.
plotModContour <- function (  tableName = "RE_coarse_Bmsy_MARE.csv", 
                              mpLabel   = "corrTRUE_SigPriorIG",
                              axes      = c("corr","kappaMult"),
                              pars      = c("BnT","Umsy","q"),
                              nSp       = 2,
                              model     = "ss",
                              compMean  = FALSE,
                              tUtr      = 15,
                              negCorr   = FALSE,
                              tcex      = 2,
                              breakRange= c(-1,1),
                              devLabels = FALSE,
                              hess      = FALSE )
{
  # Load stat table
  tablePath <- file.path ( getwd(), "project/Statistics", tableName )
  table <- read.csv( tablePath, header=TRUE, stringsAsFactors=FALSE )
  # reduce to the correct MP
  table <- table  %>% filter( mp == mpLabel )
  # Filter by other 
  if( "nS" %in% names(table) ) 
    table <- table %>% 
             filter( nS == nSp )
  if( "tUtrough" %in% names(table) ) 
    table <- table %>% 
             filter( tUtrough == tUtr )
  if( "lastNegCorr" %in% names(table) ) 
    table <- table %>% 
             filter( lastNegCorr == negCorr )

  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=nS)

  if ( length(axes) != 2) return ("Error! length(axes) != 2.\n")
  # Pick out columns for axes
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }
  # Add hessianPD column if requested
  if ( hess ) pars <- c( pars, "HessPD" )

  x   <- unique(table[,axesCols[1] ])
  x   <- x[order(x)]
  y   <- unique(table[,axesCols[2] ])
  y   <- y[order(y)]

  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    z <- array( NA, dim=c(length(y),length(x),length(pars)),
                dimnames = list( y[ length( y ):1 ], x, pars ) )

    for (xIdx in 1:length(x))
    {
      for (yIdx in 1:length(y))
      {
        # Get x and y values
        xVal <- x[xIdx]
        yVal <- y[yIdx]
        # Get rows
        xRows <- which ( specTab[,axesCols[1] ] == xVal)
        yRows <- which ( specTab[,axesCols[2] ] == yVal)
        rows <- intersect(xRows,yRows)
        if( length(rows) == 0 ) next
        colNames <- paste(model,pars,sep="")
        z[length(y)-yIdx+1,xIdx,] <- as.numeric(specTab[rows,colNames])
      }
    } 
    # browser()
    specList[[s]] <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )
  }

  names(specList) <- species
  
  par(mfrow = c((nS+1),length(pars)),  mar = c(1,1,1,1), oma = c(5,5,2,2) )
  for (s in 1:nS)
  {
    # Create a raster from the stat table info
    plotBrick <- specList[[s]]
    # Make titles for the plots
    specName <- species[s]

    titles <- names(plotBrick)
    if( hess ) l <- length(titles) else l <- length(titles) + 1

    colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
    colBreaksHess <- seq(0,200,length = 8)
    colScale  <- brewer.pal(8, "Reds")
    colScale      <- c(colScale[8:1],colScale)
    # plot the rasters, if first species plot main titles
    if ( s == 1 )
    {
      for(k in 1:(l-1))
      {
        plot( plotBrick[[k]], main = titles[k], legend = FALSE, 
              asp=NA, ylab="", xlab = "", col = colScale, breaks = colBreaks,
              axes = FALSE )
          if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
          text( plotBrick[[k]], digits = 2, halo = T, cex = tcex, axes = FALSE )
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
      }
      if( hess )
      {
        plot( plotBrick[[l]],main = titles[l], legend = FALSE, 
            asp=NA, ylab = "", xlab = "", col = colScale[9:16], 
            breaks = colBreaksHess,
            axes = FALSE )
          text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )    
      }
      
    }
    else {
      for(k in 1:(l-1))
      {
        plot( plotBrick[[k]],main = "", legend = FALSE, 
              asp=NA, ylab = "", xlab = "", col = colScale, breaks = colBreaks,
              axes = FALSE )
          if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
          text( plotBrick[[k]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
      }
      if( hess )
      {
        plot( plotBrick[[l]],main = "", legend = FALSE, 
              asp=NA, ylab = "", xlab = "", col = colScale[9:16], 
              breaks = colBreaksHess, axes = FALSE )
          text( plotBrick[[l]], digits = 2, halo = T, cex = tcex ) 
          axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
          axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
      }
      
    }   
  }

  # Now create an average/sum raster brick for the complex
  # Some fudging due to the hessian logical flag
  if( hess ) m <- l else m <- l-1
  for( k in 1:m )
  {
    # use leftover plotBrick in memory to compute complex
    plotBrick[[k]] <- 0
    for( s in 1:length(species) )
    {
      plotBrick[[k]] <- plotBrick[[k]] + specList[[s]][[k]]
    }
    if( compMean ) plotBrick[[k]] <- plotBrick[[k]]/length(species)
  }

  
  titles <- names(plotBrick)

  # Get colour breaks, based on range of data
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colBreaksHess <- seq(0,200,length = 8)
  colScale  <- brewer.pal(8, "Reds")
  colScale      <- c(colScale[8:1],colScale)

  for( k in 1:( l - 1 ) )
  {
    plot( plotBrick[[k]],main = "", legend = FALSE, 
          asp=NA, ylab = "", xlab = "", col = colScale, breaks = colBreaks,
          axes = FALSE )
      if( k == 1 ) mtext ( text = "Complex", side = 2, las = 0, line = 2 )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
      axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
      axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
  }
  if ( hess )
  {
    plot( plotBrick[[l]],main = "", legend = FALSE, 
          asp=NA, ylab = "", xlab = "", col = colScale, breaks = nSp * colBreaksHess )
      text( plotBrick[[l]], digits = 2, halo = T, cex = tcex )  
      axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
      axis( side = 1, at = 1:length(x) - 0.5, labels = x )   
  }
  
  # Labels
  mtext( side = 1, text = axes[1], outer = TRUE, cex = 2, line = 3, las = 0 )
  mtext( side = 2, text = axes[2], outer = TRUE, cex = 2, line = 3, las = 0)
  if( devLabels )
  {
    grid.text(  paste("Raw error ", model, " model.", sep = ""), 
                x=unit(0.5, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  tableName,
                x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  mpLabel,
                x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
    grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
                x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
    grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
                x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)  
  }
  
}



# plotCorrCompContour()
# Plots the correlation comparative experiment performance contours for
# estimation of BnT and Umsy (abundance and productivity) by species
# and for the complex. The plots are panel plots of Rasters, so need
# a little love as far as sizing goes.
# inputs: table=character name of table generated by .statTableXXX
#         mpLable=character name of MP label
#         axes=character 2-vector for axes of raster plots
#         nSp=integer number of species in scenario
#         method="deviance" for log2(SS/MS), difference for SS-MS
#         compMean=compute the mean for the complex?
#         tUtr=integer value for tUtrough in Fhist scenarios
#         negCorr=logical value of simCtl par lastNegCorr in corr scenarios
#         tcex=cex of text in taster cells
plotCompContour <- function ( tableName   = "RE_coarse_pub_MARE.csv", 
                              mpLabel     = "pubAM",
                              axes        = c("corr","kappaMult"),
                              pars        = c("BnT","Umsy","q","Dep"),
                              nSp         = 2,
                              method      = "deviance",
                              compMean    = FALSE,
                              tUtr        = 15,
                              negCorr     = FALSE,
                              tcex        = 2,
                              breakRange  = c(-2,2),
                              devLabels   = FALSE,
                              legend      = FALSE  
                              )
{
  # browser()
  # Load stat table
  tablePath <- file.path ( getwd(),"project/Statistics",tableName)
  table <- read.csv (tablePath, header=TRUE, stringsAsFactors=FALSE)
  # reduce to the correct MP
  table <- table %>% filter(mp == mpLabel )

  # Filter by other 
  if( "nS" %in% names(table) ) 
    table <- table %>% 
             filter( nS == nSp )
  if( "tUtrough" %in% names(table) ) 
    table <- table %>% 
             filter( tUtrough == tUtr )
  if( "lastNegCorr" %in% names(table) ) 
    table <- table %>% 
             filter( lastNegCorr == negCorr )

  # browser()
  # Now make a pallette for the colours
  colScaleHi    <- brewer.pal( 8, "Greens" )
  colScaleLo    <- brewer.pal( 8, "Reds" )
  colBreaks     <- c(-100,seq(breakRange[1],0,length=8),seq(0,breakRange[2],length=8)[-1],100)
  colScale      <- c(colScaleLo[8:1],colScaleHi)

  # calculate response
  # Deviance is log2(ss/ms)
  if (method == "deviance" )
  {
    table <- table %>% mutate(  BnT     = log2(abs(ssBnT/msBnT)),
                                Umsy    = log2(abs(ssUmsy/msUmsy)),
                                q       = log2(abs(ssq/msq)),
                                Dep     = log2(abs(ssDep/msDep)) )
  }   
  # Difference is just ss - ms
  if ( method == "difference")
  {
    table <- table %>% mutate(  BnT     = abs(ssBnT) - abs(msBnT),
                                Umsy    = abs(ssUmsy) - abs(msUmsy),
                                q       = abs(ssq) - abs(msq),
                                Dep     = abs(ssDep) - abs(msDep) )
  }

  species   <- as.character(unique( table$species))
  nS <- length(species)
  specList  <- vector( mode="list", length=length(species))

  # Pick out columns for axes
  # browser()
  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(table) == axes[i])
  }

  x   <- unique(table[,axesCols[1] ])
  x   <- x[order(x)]
  y   <- unique(table[,axesCols[2] ])
  y   <- y[order(y)]
  for (s in 1:nS)
  {
    spec <- species[s]
    specTab <- table %>% filter (species == spec )
    # browser()
    z <- array( NA, dim = c( length( y ), length( x ), length(pars) ),
                dimnames=list(  y[ length( y ):1 ], x,
                                pars ) )

    for (xIdx in 1:length(x))
    {
      for (yIdx in 1:length(y))
      {
        # Get x and y values
        xVal <- x[xIdx]
        yVal <- y[yIdx]
        # Get rows
        xRows <- which ( specTab[ ,axesCols[1] ] == xVal)
        yRows <- which ( specTab[ ,axesCols[2] ] == yVal)
        rows <- intersect(xRows,yRows)
        # Recover MSE comparison
        # browser()
        z[ length( y ) - yIdx + 1, xIdx, ] <- as.numeric(specTab[rows,pars])
      }
    } 
    specList[[s]] <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )
  }
  names(specList) <- species

  # Start plotting

  par(mfrow = c(nS+1,length(pars)), mar = c(1,1,1,1), oma = c(5,5,2,2) )
  for (s in 1:nS)
  {
    # Create a raster from the stat table info 
    plotBrick <- specList[[s]]
    # Make titles for the plots
    titles <- names(plotBrick)
    # browser()
    # plot the rasters, without legends
    for ( k in 1:length(titles) )
    {
      if(s == 1) plot( plotBrick[[k]], main = titles[k], legend = FALSE, 
                        breaks = colBreaks, col=colScale, asp=NA,
                        xlab = "", ylab = "", axes = FALSE )
      else plot( plotBrick[[k]],main = "", legend = FALSE, 
                  breaks = colBreaks, col=colScale, asp=NA,
                  xlab = "", ylab = "", axes = FALSE )
      if( k == 1 ) mtext ( text = species[s], side = 2, las = 0, line = 2 )
      text( plotBrick[[k]], digits = 2, halo = T, cex = tcex )  
      axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
      axis( side = 1, at = 1:length(x) - 0.5, labels = x )  
    }
  }

  # browser()
  # Now create the complex plot brick - idea: group by scenario
  # and MP, there should be a unique set of species for each, then
  # apply the complex rule
  compTable <- table %>%  group_by( scenario, mp ) %>%
                          summarise(  ssBnT = sum( ssBnT ),
                                      msBnT = sum( msBnT ),
                                      ssUmsy = sum( ssUmsy ),
                                      msUmsy = sum( msUmsy ),
                                      ssq = sum( ssq ),
                                      msq = sum( msq ),
                                      ssHessPD = sum( ssHessPD ),
                                      msHessPD = sum( msHessPD ),
                                      # ssMSY = sum( ssMSY ),
                                      # msMSY = sum( msMSY ),
                                      ssDep = sum( ssDep ),
                                      msDep = sum( msDep ),
                                      corr = mean(corr),
                                      kappaMult = mean(kappaMult),
                                      tUpeak = mean(tUpeak),
                                      tUtrough = mean(tUtrough),
                                      Umax = mean(Umax) )

  # Now make the comparison based on the compMean logical flag
  if( method == "deviance" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = log2( abs(ssBnT) / abs(msBnT) ),
                          Umsy = log2( abs(ssUmsy) / abs(msUmsy) ),
                          q = log2( abs(ssq) / abs(msq) ),
                          # MSY = log2( ssMSY / msMSY ),
                          Dep = log2 ( abs(ssDep) / abs(msDep) ) 
                        )
  }

  if( method == "difference" )
  {
    compTable <-  compTable %>% 
                  mutate( BnT = abs(ssBnT) - abs(msBnT) ,
                          Umsy = abs(ssUmsy) - abs(msUmsy) ,
                          q = abs(ssq) - abs(msq) ,
                          # MSY = abs(ssMSY) - abs(msMSY) ,
                          Dep = abs(ssDep) - abs(msDep) 
                        )
    if( compMean )
      compTable <-  compTable %>% 
                    mutate( BnT = BnT/nS,
                            Umsy = Umsy/nS,
                            q = q/nS,
                            # MSY = MSY/nS,
                            Dep = Dep/nS )
  }

  axesCols <- integer(length=2)
  for ( i in 1:2)
  {
    axesCols[i] <- which(names(compTable) == axes[i])
  }

  # Make a raster brick with the selected pars
  z <- array( NA, dim = c( length( y ), length( x ), length(pars) ),
              dimnames=list(  y[ length(y):1], x,
                              pars ) )
  # browser()
  for (xIdx in 1:length(x))
  {
    for (yIdx in 1:length(y))
    {
      # Get x and y values
      xVal <- x[xIdx]
      yVal <- y[yIdx]
      # Get rows
      xRows <- which ( compTable[,axesCols[1] ] == xVal)
      yRows <- which ( compTable[,axesCols[2] ] == yVal)
      rows <- intersect(xRows,yRows)
      # Recover MSE comparison
      z[length(y) - yIdx + 1,xIdx,] <- as.numeric(compTable[rows,pars])
    }
  } 
  compBrick <- brick( z, xmn = 0, xmx = length(x), ymn = 0, ymx = length(y) )

  for(k in 1:length(titles))
  {
    plot( compBrick[[k]],main = "", legend = FALSE, 
          asp=NA, col = colScale, breaks = colBreaks,
          xlab = "", ylab = "", axes = FALSE )
      if( k == 1 ) mtext ( text = "Complex", side = 2, las = 0, line = 2 )
    text( compBrick[[k]], digits = 2, halo = T, cex = tcex ) 
    axis( side = 2, at = 1:length(y) - 0.5, labels = y, las = 1 )
    axis( side = 1, at = 1:length(x) - 0.5, labels = x )   
  }

  # Plot the legends
  if ( legend )
  {
    plot( plotBrick[[1]], legend.only=TRUE, col=colScale, fill=colScale,
          legend.width=1, legend.shrink=0.75,
          axis.args = list( at = colBreaks,
                            labels = round(colBreaks,2),
                            cex.axis = 0.6),
          zlim = c(-3,3),
          legend.args = list( text="log2(ssMSE/msMSE)", 
                              side=4, font=2, 
                              line=3, cex=0.75, las=0) )
  
  }
  
  # Labels
  mtext( side = 1, text = axes[1], outer = TRUE, cex = 2, line = 3, las = 0 )
  mtext( side = 2, text = axes[2], outer = TRUE, cex = 2, line = 3, las = 0)
  if ( devLabels )
  {
    grid.text(  "Relative errors ss/ms", 
              x=unit(0.5, "npc"), y=unit(0.95, "npc"), rot=0)
    grid.text(  tableName,
                x=unit(0.15, "npc"), y=unit(0.98, "npc"), rot=0)
    grid.text(  mpLabel,
                x=unit(0.75, "npc"), y=unit(0.02, "npc"), rot=0)
    grid.text(  paste( "lastNegCorr = ", negCorr, sep = "" ),
                x=unit(0.15, "npc"), y=unit(0.02, "npc"), rot=0)
    grid.text(  paste( "tUtrough = ", tUtr, sep = "" ),
                x=unit(0.15, "npc"), y=unit(0.05, "npc"), rot=0)  
  }
}

# plotRepScan()
# Scans through BCU plots for every replicate in a simulation
# blob. Good for checking the single reps 
# inputs:   sim=number of the simulation in the project folder (alphabetical)
#           rep=optional number of replicate to start scanning from
#           est=character indicating which esitimates for plotting (MLE or MCMC)
plotRepScan <- function ( sim =1, rep = 1)
{
  # First load the simulation
  .loadSim(sim[1])

  # Recover number of reps
  nReps <- blob$ctrl$nReps

  if (rep > nReps ) 
  {
    cat ( "(Error MSG) Starting rep too high, this sim has ",
          nReps, " replicates.", sep = "")
    return()
  }

  for (r in rep:nReps)
  {
    plotBC(rep = r, sims = sim )
    invisible(readline(prompt="Press [enter] for the next plot."))
  }
  return(invisible())
}

# plotMCMCpar()
# Function that will plot MCMC output from ADMB models for a nominated
# parameter, simulation and replicate. Uses the coda package.
# inputs:   rep=replicate number; 
#           par=character name of parameter
#           sim=number indicating blob to load from project dir
# output:   NULL
# usage:    post-sim, showing MCMC performance
plotMCMCpar <- function ( rep=1, par="Bmsy", sim=1 )
{
  # Check if blob is loaded, if not, load the
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctrl$nReps,sep="")

  # Recover blob elements for plotting
  nS    <- blob$opMod$nS
  nT    <- blob$opMod$nT

  # Species names
  specNames <- blob$ctl$specNames

  # Recover MCMC output for each model and species
  parMCoutSS <- t(blob$am$ss$mcPar[rep,,,as.character(par)])
  parMCoutMS <- t(blob$am$ms$mcPar[rep,,,as.character(par)])
  # There's an issue with having fewer MC replicates than expected
  # so reduce the length of mc output 
  nonNAss <- numeric(nS)
  for (s in 1:nS)
  { 
    nonNAss[s] <- length ( parMCoutSS[,s][!is.na(parMCoutSS[,s])])
  }
  nonNAms <- length(parMCoutMS[!is.na(parMCoutMS)])/nS
  # Set up mc output thinning
  parMCoutSS <- as.mcmc(parMCoutSS[1:min(nonNAss),])
  parMCoutMS <- as.mcmc(parMCoutMS[1:nonNAms,])

  # Now to plot each model's posterior distributions
  titleSS <- paste ( "SS MCMC for ", par, sep = "")
  titleMS <- paste ( "MS MCMC for ", par, sep = "")

  # plot SS trace and posteriors
  plot(parMCoutSS,las=1)
  mtext(side=3,text=titleSS,outer=TRUE,padj=2)

  # create a new device and plot MS trace and posteriors
  dev.new()
  plot(parMCoutMS,las=1)
  mtext(side=3,text=titleMS,outer=TRUE,padj=2)  

  return()
}

# plotMCMCspecies()
# Function that will plot MCMC output from ADMB models for a nominated
# species and replicate from the loaded blob. Uses the coda package.
# inputs:   rep=replicate number to start from; 
#           spec=number indicating the species
#           sim=number indicating simulation to load if no blob present
# output:   NULL
# usage:    post-sim, showing MCMC performance
plotMCMCspecies <- function ( rep=1, spec=1, sim = 1 )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)
 
  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctl$nReps,sep="")

  # Recover blob elements for plotting
  nS    <- blob$opMod$nS
  nT    <- blob$opMod$nT

  # Species names
  if (!is.null(blob$ctrl$specNames))
  {
    specNames <- blob$ctrl$specNames
    spec <- specNames[spec]
  }
  
  # Recover MCMC output for each model and species
  parMCoutSS <- blob$am$ss$mcPar[rep,spec,,]
  parMCoutMS <- blob$am$ms$mcPar[rep,spec,,]
  # There's an issue with having fewer MC replicates than expected
  # so reduce the length of mc output 
  nonNAss <- length(parMCoutSS[!is.na(parMCoutSS)])/ncol(parMCoutSS)
  nonNAms <- length(parMCoutMS[!is.na(parMCoutMS)])/ncol(parMCoutMS)
  # Set up mc output thinning
  parMCoutSS <- as.mcmc(parMCoutSS[1:nonNAss,])
  parMCoutMS <- as.mcmc(parMCoutMS[1:nonNAms,])

  # Now to plot each model's posterior distributions
  titleSS <- paste ( "SS MCMC for ", spec, sep = "")
  titleMS <- paste ( "MS MCMC for ", spec, sep = "")

  # plot SS trace and posteriors
  plot(parMCoutSS,las=1)
  mtext(side=3,text=titleSS,outer=TRUE,padj=2)

  # create a new device and plot MS trace and posteriors
  dev.new()
  plot(parMCoutMS,las=1)
  mtext(side=3,text=titleMS,outer=TRUE,padj=2)  

  return()
}

# plotCIbioScan()
# First load the simulation
plotCIbioScan <- function ( sim =1, rep = 1 )
{
  .loadSim(sim)

  # Recover number of reps
  nReps <- blob$ctrl$nReps

  if (rep > nReps ) 
  {
    cat ( "(Error MSG) Starting rep too high, this sim has ",
          nReps, " replicates.", sep = "")
    return()
  }

  for (r in rep:nReps)
  {
    plotCIbio( r )
    invisible( readline( prompt = "Press [enter] for the next plot." ) )
  }
  return(invisible())
}

# plotCICbio()
# Function that plots true and posterior biomass for each species for a
# given simulation and replicate.
# inputs:   rep=replicate number
#           sim=number indicating blob to load from project dir
#           quant=numeric of percentiles to be calculated
plotCIbio <- function ( rep = 1, sim=1)
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Recover control pars
  nS <- blob$opMod$nS
  nT <- blob$opMod$nT
  species <- blob$ctrl$speciesNames[1:nS]

  browser()

  # load biomass trajectories
  omBio <- blob$om$Bt[rep,,]
  msBio <- array( NA, dim = c( nS, nT, 3 ), 
                  dimnames = list( species, 1:nT, c("lCI", "val", "uCI") ) )
  ssBio <- array( NA, dim = c( nS, nT, 3 ), 
                  dimnames = list( species, 1:nT, c("lCI", "val", "uCI") ) )

  # Fill MS model biomass
  if( !is.null( blob$am$ms$CIs[[ rep ]] ) )
  {
    msCIs <- blob$am$ms$CIs[[ rep ]]
  } 
  else {
    msStdErr <- blob$am$ms$sdrep[[ rep ]]
    if( is.null( msStdErr ) ) 
      msCIs <- NA
    else {  
      msStdErr <- summary( msStdErr )
      colnames( msStdErr ) <- c( "val", "se" )
      msCIs <-  msStdErr %>%
                as.data.frame() %>%
                mutate( par = rownames(msStdErr),
                        lCI = val - 1.96*se,
                        uCI = val + 1.96*se ) %>%
                dplyr::select( par, val, se, lCI, uCI )
    }
  }
  # browser()
  if( !is.na(msCIs) )
  {
    Btrows <- which( msCIs$par == "lnBt" )
    msBio[ , , 1 ] <- matrix( exp(msCIs[ Btrows , "uCI" ]), nrow = nS, ncol = nT, byrow = FALSE )
    msBio[ , , 2 ] <- matrix( exp(msCIs[ Btrows, "val" ]), nrow = nS, ncol = nT, byrow = FALSE )
    msBio[ , , 3 ] <- matrix( exp(msCIs[ Btrows, "lCI" ]), nrow = nS, ncol = nT, byrow = FALSE ) 
  }

  # Now fill SS model biomas
  if( !is.null(blob$am$ss$CIs) ) 
  {
    ssCIs <- blob$am$ss$CIs[[rep]]
  }
  else {
    # Create CIs from sd report objects
    ssCIs <- vector( mode = "list", length = 3 )
    ssSDreps <- blob$am$ss$sdrep[[rep]]
    # The following makes up for some coding errors in earlier sims
    if( length( ssSDreps ) == 1 ) k <- nS else k <- 1
    # loop over valid sd objects
    for( s in k:nS )
    {
      ssStdErr <- ssSDreps[[ s ]]
      # Save single species sd rep obj
      if( is.null(ssStdErr) ) 
        ssCIs[[ s ]] <- NA
      else {
      ssStdErr <- summary( ssStdErr )
      colnames( ssStdErr ) <- c("val","se")
      ssCIs[[ s ]] <-   ssStdErr %>%
                        as.data.frame() %>%
                        mutate( par = rownames(ssStdErr),
                                lCI = val - 1.96*se,
                                uCI = val + 1.96*se ) %>%
                                dplyr::select( par, val, se, lCI, uCI )    
      }
    }
  }
  for( s in 1:nS )
  {
    if( is.null( ssCIs[[ s ]] ) ) next
    if( is.na( ssCIs[[ s ]] ) ) next
    Btrows <- which( ssCIs[[ s ]]$par == "lnBt" )
    ssBio[ s, , 1 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "uCI" ]), nrow = 1, ncol = nT )
    ssBio[ s, , 2 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "val" ]), nrow = 1, ncol = nT )
    ssBio[ s, , 3 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "lCI" ]), nrow = 1, ncol = nT )  
  }
  

  # Fit diagnostics (hessian)
  hpdSS <- blob$am$ss$hesspd[rep,]
  hpdMS <- blob$am$ms$hesspd[rep]

  # Set colours
  ssCol <- "steelblue"
  msCol <- "salmon"
  # Create transparent fill colours for polys
  ssFill <- col2rgb(ssCol)/255
  ssFill <- rgb(ssFill[1],ssFill[2],ssFill[3],alpha=0.5)
  msFill <- col2rgb(msCol)/255
  msFill <- rgb(msFill[1],msFill[2],msFill[3],alpha=0.5)

  # Now plot!
  par(mfrow = c(nS,1), mar = c(2,4,1,1), oma=c(4,4,1,1))
  for (s in 1:nS)
  {
    # Compute max B value
    yLim <- 2*max (omBio[s,])
    # Blank plot
    plot (  x=c(1,nT), y =c(1,yLim), type = "n", xlab = "", ylab = "Biomass (t)",
            las=1 )
    # Create polygon vertices
    xPoly <- c(1:nT,nT:1)
    yPolySS <- c(ssBio[s,1:nT,1],ssBio[s,nT:1,3])
    yPolyMS <- c(msBio[s,1:nT,1],msBio[s,nT:1,3])
    if (s == 1) panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                            fill=c(ssFill,msFill), border=c(NA,NA), 
                            cex=c(0.7), bty="n" )
    if ( hpdSS[s] ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=0.7)
    if ( hpdMS ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=0.7)
    # plot CI bands
    # SS
    polygon (x = xPoly, y = yPolySS,border=NA,col = ssFill,density=NA)
    # MS
    polygon (x = xPoly, y = yPolyMS,border=NA,col = msFill,density=NA)
    # Plot Medians
    lines ( x = 1:nT, y = ssBio[s,,2], col = ssCol,lty=2,lwd=2 )
    lines ( x = 1:nT, y = msBio[s,,2], col = msCol,lty=2,lwd=2 )
    # Plot true biomass
    lines ( x=1:nT, y = omBio[s,],col="black", lwd = 2)
  }
}

plotBenvelopes <- function( sim = 1,
                              fYear = 1988,
                              est = FALSE,
                              labSize = 2, 
                              tickSize = 1.2,
                              devLabels = TRUE,
                              alpha = .8 )
{
  # First, load blob
  .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS <- blob$opMod$nS
  nT <- blob$opMod$nT

  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt
  Bmsy  <- blob$opMod$Bmsy

  # Single species model
  ssBt  <- blob$am$ss$Bt

  # Multispecies model
  msBt  <- blob$am$ms$Bt
  msBt[ msBt < 0 ] <- NA

  years <- fYear:(fYear + nT - 1)

  for( s in 1:nS )
  {
    omBt[,s,]  <- omBt[,s,]/Bmsy[s]/2
    ssBt[,s,]  <- ssBt[,s,]/Bmsy[s]/2
    msBt[,s,]  <- msBt[,s,]/Bmsy[s]/2
  }  
  
  # Set colours
  ssCol <- "steelblue"
  msCol <- "salmon"
  # Create transparent fill colours for polys
  ssFill <- col2rgb(ssCol)/255
  ssFill <- rgb(ssFill[1],ssFill[2],ssFill[3],alpha=0.4)
  msFill <- col2rgb(msCol)/255
  msFill <- rgb(msFill[1],msFill[2],msFill[3],alpha=0.4)

  ## CHECK THAT THESE ARRAY CALCS ARE WORKING ##
  # Now, reduce to envelopes
  # OM values
  omBt  <- apply( X = omBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  # SS AM 
  ssBt  <- apply( X = ssBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  
  # MS AM 
  msBt  <- apply( X = msBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  
  # Set up plotting environment
  par ( mfrow = c(nS,1), mar = c(1,2,2,0), oma = c(4,4,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # loop over species and plot biomass
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,1.1), type = "n",
              ylim = c(0,1.1), ylab = "", axes=FALSE, las = 1, xlab = ""  )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      mtext( side = 2, text = specNames[s], line = 2.5, las = 0 )
      # OM
      polygon(  x = c(years,rev(years)),
                y = c(omBt[1,s,],rev(omBt[3,s,])),
                border = NA, col = "grey80" )
      if( est )
      {
        # SS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( ssBt[ 1, s, ], rev( ssBt[ 3, s, ] ) ),
                 border = NA, col = ssFill )
        # MS envelop
        polygon( x = c( years, rev( years ) ),
                 y = c( msBt[ 1, s, ], rev( msBt[ 3, s, ] ) ),
                 border = NA, col = msFill )
        # median
        lines( x = years, y = ssBt[ 2, s, ], col = ssCol, lwd = 2, lty = 2 )
        lines( x = years, y = msBt[ 2, s, ], col = msCol, lwd = 2, lty = 2 )

      }
      lines( x = years, y = omBt[ 2, s, ], lwd = 2, lty = 2 )
      abline( h = (1 - alpha / 2), lty = 3, lwd = .8 )
    
  }
  mtext( text = "Depletion", side = 2, line = 2.5, las = 0, outer = T, cex = 1.5 )
  mtext( text = "Year", outer = TRUE, side = 1, padj = 1.5, cex = 1.5, line = 2)
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.9),padj=2,col="grey50",cex=0.8 )
}


plotUshrinkage <- function( sim = 1, rep = 1, SS = T, MS = T, MSdist = T, OM = T,
                            devLabels = T )
{
  # load the simulation
  .loadSim(sim)

  # Now grab OM U values
  nStocks <- blob$opMod$nS
  U_om    <- blob$opMod$Umsy[1:nStocks]
  
  # grab both SS and MS estimated Umsy
  estUmsy <- array( 0, dim = c( nStocks, 2, 3 ),
                        dimnames = list( 1:nStocks, c("ss","ms"), c("LB", "MLE", "UB" ) )
                  )
  # Save estimate tables
  SS_CIs <- blob$am$ss$CIs[[rep]]
  MS_CIs <- blob$am$ms$CIs[[rep]]
  # First grab single stock estimates

  for( sIdx in 1:nStocks )
  {
    stockCIs <- SS_CIs[[sIdx]]
    estUmsy[sIdx, "ss", ] <- as.numeric(stockCIs[stockCIs$par == "Umsy",c(4,2,5)])
  }
  estUmsy[,"ms",] <- as.matrix(MS_CIs[MS_CIs$par == "Umsy",c(4,2,5)])

  Umsybar <- blob$am$ms$Umsybar[rep,]
  sigU2   <- blob$am$ms$sigU2[rep,]

  mpLabel <- blob$ctrl$mpLabel

  # Set up plotting window
  estUmsy[!is.finite(estUmsy)] <- NA
  URange  <- range( estUmsy, U_om )
  xLim    <- c( .8*URange[1],1.2*URange[2] )

  UestNorm <- dnorm(seq(xLim[1],xLim[2], length = 100), mean = Umsybar, sd = sqrt(sigU2) )

  cols <- brewer.pal(n = nStocks+1, name= "PuOr")

  yPoints <- seq(0,max(UestNorm),length = nStocks+1)
  yPoints <- yPoints[2:4]

  yJitter <- seq(-.1, .1, length = nStocks)  
  


  plot( x = xLim, y = c( 0, 1.2*max(UestNorm) ), type = "n", axes = F, xlab = "", ylab = "" )
    x       <- seq( xLim[1], xLim[2], length = 100 )
    axis( side = 1 )
    if( SS )
    {
      points( x = estUmsy[,"ss","MLE"], y = yPoints[1] + yJitter, pch = 17, cex = 2, col = cols )
      segments( x0 = estUmsy[,"ss","LB"], x1 = estUmsy[,"ss","UB"], y0 = yPoints[1] + yJitter,
                col = cols, lwd = 2 )
    }
    if( MS )
    {
      points( x = estUmsy[,"ms","MLE"], y = yPoints[3] + yJitter, pch =16, cex = 2, col = cols )
      segments( x0 = estUmsy[,"ms","LB"], x1 = estUmsy[,"ms","UB"], y0 = yPoints[3] + yJitter,
                col = cols, lwd = 2 )

    }
    if( MSdist )
    {
      lines( x = x, y = UestNorm, lty = 2, lwd = 3, col = "grey50" )
      abline( v = Umsybar, col ="grey50" )
    }
    if( OM ) points( x = U_om, y = yPoints[2] + yJitter, col = cols, pch = 18, cex = 2 )
    if( (SS & MS) & (!OM) )
      segments( x0 = estUmsy[,"ss","MLE"], x1 = estUmsy[,"ms","MLE"],
                y0 = yPoints[1] + yJitter, y1 = yPoints[3] + yJitter,
                col = cols)
    if( SS & MS & OM )
    {
      segments( x0 = estUmsy[,"ss","MLE"], x1 = U_om,
                y0 = yPoints[1] + yJitter, y1 = yPoints[2] + yJitter,
                col = cols )
      segments( x0 = estUmsy[,"ms","MLE"], x1 = U_om,
                y0 = yPoints[3] + yJitter, y1 = yPoints[2] + yJitter,
                col = cols)
    }


  mtext( side = 1, text = "Productivity", outer = F, line = 2 )
  if( devLabels ) mtext( side = 1, line = 2, adj = 1, text = mpLabel )

  
}


# plotqDists()
# Does what it says, plots to quartz the estimated and simulated q distributions
# for the complex in a given simulation and replicate.
# inputs:   sim = simulation number
#           rep = replicate number
# ouputs:   NULL
# side-eff: plots to quartz
# source: SDNJ
plotqDists <- function( sim = 1, rep = 1, surveys = c(1,2), devLabels = TRUE,
                        OM = TRUE, MS = TRUE, SS = TRUE )
{
  # load the simulation
  .loadSim(sim)

  # Now grab opMod q distribution pars
  qSurvM  <- blob$opMod$qSurvM
  qSurvSD <- blob$opMod$qSurvSD
  nSurv   <- blob$opMod$nSurv
  nStocks <- blob$opMod$nS
  # grab simulated q 
  simq    <- blob$om$q_os[rep,,]
  # grab both SS and MS estimated q
  estqSS  <- blob$am$ss$q_os[rep,,]
  estqMS  <- blob$am$ms$q_os[rep,,]
  # And MS estimated distribution
  estqSurvM   <- blob$am$ms$qbar[rep,]
  estqSurvSD  <- sqrt(blob$am$ms$tauq2[rep,])

  SS_CIs <- blob$am$ss$CIs[[rep]]
  MS_CIs <- blob$am$ms$CIs[[rep]]


  mpLabel <- blob$ctrl$mpLabel

  # Split plotting window
  par( mfrow = c(length(surveys),1), mar = c(1,1,1,1), oma = c(3,3,1,1) )
  # Calculate range of q values
  qRange  <- range( simq, estqSS, estqMS, na.rm = T )
  xLim    <- c( .5*qRange[1],1.5*qRange[2] )
  for( survIdx in surveys )
  {
    # Now make simulated and estimated distributions
    x       <- seq( xLim[1], xLim[2], length = 100 )
    qlNorm    <- dlnorm( x, meanlog = log(qSurvM[survIdx]), sdlog = qSurvSD[survIdx] )
    qlNormEst <- dlnorm( x, meanlog = log(estqSurvM[survIdx]), sdlog = estqSurvSD[survIdx] )
    # Initialise plots
    plot( x = xLim, y = c( 0, 1.2*max(qlNorm) ), type = "n", axes = F, xlab = "", ylab = "" )
      axis( side = 1 )
      panLab( x = .8, y = 0.8, txt = paste("Survey ", survIdx, sep = "") )
      # Plot distributions
      if( OM ) 
        lines( x = x, y = qlNorm, lty = survIdx + 1, lwd = 3 )
      if( MS )
        lines( x = x, y = qlNormEst, lty = survIdx + 1, lwd = 3, col = "grey60" )
      # Plot simulated survey mean q
      if( OM )
      {
        abline( v = qSurvM[survIdx], lty = survIdx + 1 )
        # Now plot the simulated values for each stock
        text( x = simq[survIdx,], y = (1:nStocks)*max(qlNorm)/(nStocks + 1),
              labels = 1:nStocks )
      }
      
      # Joint SS and MS estimates with a thin line
      if( SS & MS )
        segments( x0 = estqSS[survIdx,], y0 = (1:nStocks)*max(qlNorm)/(nStocks + 1),
                  x1 = estqMS[survIdx,], y1 = (1:nStocks)*max(qlNorm)/(nStocks + 1),
                  lwd = .8, lty = 4 )
      # Plot SS and MS estimates
      if( SS )
        points( x = estqSS[survIdx,], y = (1:nStocks)*max(qlNorm)/(nStocks + 1), pch = 16 )
      if( MS )
      {
        points( x = estqMS[survIdx,], y = (1:nStocks)*max(qlNorm)/(nStocks + 1), pch = 17 )
        abline( v = estqSurvM[survIdx], lty = survIdx + 1 + nSurv )
      }
  }
  mtext( side = 1, text = "Catchability", outer = T, line = 1 )
  if( devLabels ) mtext( side = 1, line = 2, adj = 1, text = mpLabel )

  
} 

# Plots simulation tulip plots of biomass, catch and harvest
# rate
plotBCUenvelopes <- function( sim = 1,
                              fYear = 1988,
                              est = FALSE,
                              labSize = 2, 
                              tickSize = 1.2,
                              devLabels = TRUE )
{
  # First, load blob
  .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS <- blob$opMod$nS

  # Year indexing
  if(!is.null(blob$opMod$fYear)) 
  {
    nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
    fYear <- min(blob$opMod$fYear)
  } else nT <- blob$opMod$nT

  years <- fYear:(fYear + nT - 1)


  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt
  Ct    <- blob$om$Ct
  Ut    <- blob$om$Ut

  # Single species model
  ssBt  <- blob$am$ss$Bt

  # Multispecies model
  msBt  <- blob$am$ms$Bt
  msBt[msBt < 0] <- NA

  years <- fYear:(fYear + nT - 1)

  # Estimated Ut
  ssUt <- Ct / ssBt
  msUt <- Ct / msBt

  ## Now scale quantities to eqbm values
  # Pull Bmsy and Umsy from the blob
  Umsy <- blob$opMod$Umsy
  Bmsy <- blob$opMod$Bmsy
  MSY  <- Umsy*Bmsy

  
  for( s in 1:nS )
  {
    omBt[,s,]  <- omBt[,s,]/Bmsy[s]/2
    Ct[,s,]    <- Ct[,s,]/MSY[s]
    Ut[,s,]    <- Ut[,s,]/Umsy[s]
    ssBt[,s,]  <- ssBt[,s,]/Bmsy[s]/2
    msBt[,s,]  <- msBt[,s,]/Bmsy[s]/2
    msUt[,s,]    <- msUt[,s,]/Umsy[s]
    ssUt[,s,]    <- ssUt[,s,]/Umsy[s]
  }  
  
  # Set colours
  ssCol <- "steelblue"
  msCol <- "salmon"
  # Create transparent fill colours for polys
  ssFill <- col2rgb(ssCol)/255
  ssFill <- rgb(ssFill[1],ssFill[2],ssFill[3],alpha=0.4)
  msFill <- col2rgb(msCol)/255
  msFill <- rgb(msFill[1],msFill[2],msFill[3],alpha=0.4)

  ## CHECK THAT THESE ARRAY CALCS ARE WORKING ##
  # Now, reduce to envelopes
  # OM values
  omBt  <- apply( X = omBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  Ct    <- apply( X = Ct, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  omUt  <- apply( X = Ut, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.5), na.rm = TRUE )
  # SS AM 
  ssBt  <- apply( X = ssBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  ssUt  <- apply( X = ssUt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )

  # MS AM 
  msBt  <- apply( X = msBt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )
  msUt  <- apply( X = msUt, FUN = quantile, MARGIN = c(2,3), 
                  probs = c(0.025,0.5,0.975), na.rm = TRUE )

  # Set up plotting environment
  par ( mfrow = c(3,nS), mar = c(1,3,2,0), oma = c(3,3,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # loop over species and plot biomass
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,1.1), type = "n",
              ylim = c(0,1.1), ylab = "", axes=FALSE, las = 1, xlab = "" ,
              main = specNames[s] )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      # OM
      polygon(  x = c(years,rev(years)),
                y = c(omBt[1,s,],rev(omBt[3,s,])),
                border = NA, col = "grey80" )
      if( est )
      {
        # SS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( ssBt[ 1, s, ], rev( ssBt[ 3, s, ] ) ),
                 border = NA, col = ssFill )
        # MS envelop
        polygon( x = c( years, rev( years ) ),
                 y = c( msBt[ 1, s, ], rev( msBt[ 3, s, ] ) ),
                 border = NA, col = msFill )
        # median
        lines( x = years, y = ssBt[ 2, s, ], col = ssCol, lwd = 2, lty = 2 )
        lines( x = years, y = msBt[ 2, s, ], col = msCol, lwd = 2, lty = 2 )

      }
      lines( x = years, y = omBt[ 2, s, ], lwd = 2, lty = 2 )
    if( s == 1 ) mtext( text = "Depletion", side = 2, line = 2.5, las = 0 )
  }

  # Now catch
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,2.1), type = "n",
              ylim = c(0,2.1), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      # OM
      polygon(  x = c(years,rev(years)),
                y = c(Ct[1,s,],rev(Ct[3,s,])),
                border = NA, col = "grey80" )
      lines( x = years, y = Ct[ 2, s, ], lwd = 2, lty = 2 )
    if( s == 1 ) mtext( text = "Ct/MSY", side = 2, line = 2.5, las = 0 )
  }

  # Lastly, exploitation rate
  for( s in 1:nS )
  {
    plot    ( x = c(fYear,max(years)), y = c(0,2.2), type = "n",
              ylim = c(0,2.2), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
      if( est )
      {
        # SS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( ssUt[ 1, s, ], rev( ssUt[ 3, s, ] ) ),
                 border = NA, col = ssFill )
        # MS envelope
        polygon( x = c( years, rev( years ) ),
                 y = c( msUt[ 1, s, ], rev( msUt[ 3, s, ] ) ),
                 border = NA, col = msFill )
        # medians
        lines( x = years, y = ssUt[ 2, s, ], lwd = 2, lty = 2, col = ssCol )
        lines( x = years, y = msUt[ 2, s, ], lwd = 2, lty = 2, col = msCol )

      }
      # OM
      lines( x = years, y = omUt[ s, ], lwd = 2, lty = 2 )
    if( s == 1 ) mtext( text = "Ut/Umsy", side = 2, line = 2.5, las = 0 )
  }
  mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5)
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.9),padj=2,col="grey50",cex=0.8 )
}

# plotBC()
# Plots a given replicate's true and estimated time series of biomass and 
# catch for all stocks, with sims as columns
# inputs:   sims=numeric indicating blobs to load (alpha order in project folder)
#           legend=logical indicating plotting of legend (SS/MS line cols)
#           data=logical indicating whether to plot survey data
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBenv <- function(  sims=1, legend=TRUE,
                        data = FALSE, labSize = 2, tickSize = 1.2,
                        fYear = 1988, devLabels = TRUE,
                        scale = FALSE, est = TRUE,
                        titles = NULL, density = NULL )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  .loadSim(sims[1])

  # Recover blob elements for plotting
  nS <- blob$opMod$nS
  nO <- blob$opMod$nSurv

  # Set up plot window
  par ( mfcol = c(nS,length(sims)), mar = c(1,3,2,0), oma = c(4,4,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )

  for( simIdx in 1:length(sims) )
  {
    .loadSim( sims[simIdx] )
  

    # Create a stamp from scenario and mp name
    scenario  <- blob$ctrl$scenario
    mp        <- blob$ctrl$mp
    stamp     <- paste(scenario,":",mp,sep="")

    
    

    # Species names
    specNames <- blob$ctrl$speciesNames

    # True OM quantities
    omBt  <- apply( X = blob$om$Bt, MARGIN = c(2,3), FUN = quantile, probs = c(0.05, 0.5, 0.95), na.rm = T )
    # Single species model
    ssBt <- blob$am$ss$Bt
    ssBt[ssBt<0] <- NA
    ssBt  <- apply( X = ssBt, MARGIN = c(2,3), FUN = quantile, probs = c(0.05, 0.5, 0.95), na.rm = T )
    # Multispecies model
    msBt <- blob$am$ms$Bt
    msBt[msBt<0] <- NA
    msBt  <- apply( X = msBt, MARGIN = c(2,3), FUN = quantile, probs = c(0.05, 0.5, 0.95), na.rm = T )



    # Year indexing
    if(!is.null(blob$opMod$fYear)) 
    {
      nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
      fYear <- min(blob$opMod$fYear)
    } else nT <- blob$opMod$nT
    # browser()
    years <- fYear:(fYear + nT - 1)

    # Pull Bmsy and Umsy from the blob
    Bmsy <- blob$opMod$Bmsy

    if( scale )
    {
      # Now loop over species and scale
      for( s in 1:nS )
      {
        omBt[s,]  <- omBt[s,]/Bmsy[s]/2
        ssBt[s,]  <- ssBt[s,]/Bmsy[s]/2
        msBt[s,]  <- msBt[s,]/Bmsy[s]/2
      }  
    }
    
    # Set colours for each model
    ssCol <- "steelblue"
    msCol <- "salmon"

    # Plot biomass, actual and estimated, including 2 index series,
    # scaled by model estimated q
    maxBt <- 1.1*max(omBt,na.rm=T)
    for ( s in 1:nS )
    {
      plot    ( x = c(fYear,max(years)), y = c(0,maxBt), type = "n",
                ylim = c(0,maxBt), ylab = "", axes=FALSE, xlab = "" ,
                main = "" )
        if( !is.null(titles) & s == 1 ) title( main = titles[simIdx] )
        if(s == nS)axis( side = 1, las = 0, cex.axis = labSize )
        axis( side = 2, las = 1, cex.axis = labSize )
      if (s == 1 & est & legend) 
        panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                    col=c(ssCol,msCol), lty = c(2,3), 
                    lwd = c(2,2), cex=c(1), bty="n" )
        panLab( x = .5, y = .9, txt = paste("Stock ", s, sep = "" ), cex = labSize )
 
      polygon   ( x = c(years,rev(years)), y = c(omBt[1,s,],rev(omBt[3,s,])), border = NA, 
                  col = "grey70" )
      lines( x = years, y = omBt[2,s,], col = "black", lwd = 2 )
      if(est)
      {
        ssColPoly <- adjustcolor( ssCol, alpha.f = 0.4)
        msColPoly <- adjustcolor( msCol, alpha.f = 0.4)
        polygon(  x = c(years, rev(years)), y = c(ssBt[1,s,],rev(ssBt[3,s,])),
                  border = NA, col = ssColPoly, angle = 45, density = density,
                  lty = 1, lwd = 2 )
        polygon(  x = c(years, rev(years)), y = c(msBt[1,s,],rev(msBt[3,s,])),
                  border = NA, col = msColPoly, angle = -45, density = density,
                  lty = 1, lwd = 2 )
        lines   ( x = years, y = ssBt[2,s,], col = ssCol, lwd = 2, lty = 2 )
        lines   ( x = years, y = msBt[2,s,], col = msCol, lwd = 2, lty = 3 )
      }
    }
    mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5, cex = labSize )
    mtext ( text = "Biomass (kt)", outer = TRUE, side = 2, line = 2, las = 0, cex = labSize )
    if( devLabels )
      mtext ( text = c(stamp),side=1, outer = TRUE, 
              at = c(0.8),padj=2,col="grey50",cex=0.8 )
  }
}

makeSimNumTable <- function()
{
  simList <- list.files("./project/", full.names = TRUE)
  simList <- simList[grepl(pattern = "sim", x = simList)]
  # Count simulations, make a table of scenarios/MPs
  nSims <- length(simList)
  simNumTable <- matrix(NA, nrow = nSims, ncol = 4 )
  colnames(simNumTable) <- c("simNum","scenario","mp","nS")

  simNumTable <- as.data.frame(simNumTable)
  simNumTable$simNum <- 1:nSims

  for( i in 1:nSims )
  {
    .loadSim(i)
    simNumTable[i,"scenario"] <- blob$ctrl$scenarioName
    simNumTable[i,"mp"] <- blob$ctrl$mpLabel
    simNumTable[i,"nS"] <- blob$opMod$nS
  }
  simNumTable
}

dumpBCsim <- function(  simPath = "./project",
                        prefix = "pubBase",
                        MPs = c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors" ),
                        MPlabels = expression("Single Stock","None", q, r, q/r ),
                        simNumTable = NULL, devLabels = TRUE )
{
  # Need to read in all the sims in the folder, and then
  # tabulate so we can plot out the sets of MPs
  # List directories
  simList <- list.files(simPath, full.names = TRUE)
  simList <- simList[grepl(pattern = "sim", x = simList)]
  if(simPath != "./project")
  {
    file.copy( from = simList, to = "./project/", recursive = TRUE)
  }

  plotPath <- file.path("./project/figs",prefix)
  if(!dir.exists(plotPath))
    dir.create(plotPath)

  plotPath <- file.path(plotPath,"BCsim")
  if(!dir.exists(plotPath))
    dir.create(plotPath)  

  if(is.null(simNumTable))
    simNumTable <- makeSimNumTable()

  scenarios <- unique(simNumTable$scenario)
  if(is.null(MPs))
    MPs     <- unique(simNumTable$mp)

  for( sIdx in 1:length(scenarios) )
  {
    scenLabel <- scenarios[sIdx]
    scenPath <- file.path(plotPath, scenLabel)
    if(!dir.exists(scenPath))
      dir.create(scenPath)

    subTable <- simNumTable %>%
                filter( scenario == scenLabel )

    mpOrder <- numeric(length(MPs))
    for( mIdx in 1:length(MPs))
      mpOrder[mIdx] <- subTable[which(subTable[,"mp"] == MPs[mIdx] ),"simNum"]

    nS <- unique(subTable$nS)

    # Now, loop over rep numbers and plot
    plotBCsimReps(  sims = mpOrder,
                    saveFileRoot = file.path(scenPath,"BCsim"),
                    legend=FALSE,
                    data = FALSE,
                    CIs = TRUE,
                    scale = FALSE,
                    titles = MPlabels, MPtitles = FALSE,
                    labSize = 2, tickSize = 2, 
                    devLabels = devLabels, maxBt = rep(45, nS) )

  }

  if(simPath != "./project")
  {
    simList <- list.files("./project", full.names = FALSE)
    simList <- simList[grepl(pattern = "sim", x = simList)]
    for(k in 1:length(simList))
      system( paste("rm -r ", simList[k], sep = "") )
  }
    cat("BCsim() plots for ", prefix," complete\n", sep = "")
}

dumpPerfMetrics <- function(  tabNameRoot = "pubBase", stockLabel = "Stock1",
                              vars = c("Umsy","BnT","Bmsy","Dep","q_1","q_2"),
                              varLabels = expression(U[MSY], B[T], B[MSY], D[T], q[11], q[21]),
                              MPs = c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors" ),
                              MPlabels = expression("None", q, r, q/r ),
                              bw = FALSE )                              
{
  graphics.off()
  # First, create a directory to dump out all the plots
  batchPlotPath <- file.path("./project/figs",tabNameRoot)
  if(!dir.exists(batchPlotPath) )
    dir.create(batchPlotPath)

  # Read in tables
  tablePath <- "./project/statistics" 

  # Pred Interval
  PItab     <- paste(tabNameRoot,"_PI.csv", sep = "" )
  PItab     <- read.csv(file.path(tablePath,PItab), header = T, stringsAsFactors = FALSE)

  # Interval Coverage
  ICtab     <- paste(tabNameRoot,"_IC.csv", sep = "" )
  ICtab     <- read.csv(file.path(tablePath,ICtab), header = T, stringsAsFactors = FALSE)

  # MREs
  MREtab     <- paste(tabNameRoot,"_MRE.csv", sep = "" )
  MREtab     <- read.csv(file.path(tablePath,MREtab), header = T, stringsAsFactors = FALSE)

  # MAREs
  MAREtab     <- paste(tabNameRoot,"_MARE.csv", sep = "" )
  MAREtab     <- read.csv(file.path(tablePath,MAREtab), header = T, stringsAsFactors = FALSE)

  scenarios <- unique(PItab$scenario)
  if(is.null(MPs))
    MPs       <- unique(PItab$mp)

  perfPlotPath <- file.path(batchPlotPath, "perfMetrics" )
  if(!dir.exists(perfPlotPath) )
    dir.create(perfPlotPath)

  # Loop over scenarios, plot 1 for each scenario
  for( scenIdx in 1:length(scenarios) )
  {
    scenLabel <- scenarios[scenIdx]
    scenPlotPath <- file.path(perfPlotPath,paste(scenLabel,".pdf",sep = ""))
    pdf( file = scenPlotPath, width = 8.5, height = round(8.5*1.3/2) )
    par(mfcol = c(length(vars), length(MPs) ), oma =c(6,6,6,7), mar = c(0,0,0,0),
        lab.cex = 1.5, cex.axis = 1.5 )
    # MPs are columns, and vars are rows
    # we are going row by row, so loop over vars first
    for( mpIdx in 1:length(MPs) )
    { 
      for(varIdx in 1:length(vars))
      {
        varName <- vars[varIdx]

        # Check if run finished
        subTab <- PItab %>%
                  filter( mp == MPs[mpIdx], scenario == scenLabel)
        if( nrow(subTab) < 2 )
        {
          plot( x = c(0,1), y = c(0,1),
                type = "n", axes = F )
          next
        }

        plotPerfMetrics(  scenLabel = scenLabel,
                          tabNameRoot = tabNameRoot,
                          mpLabel = MPs[mpIdx],
                          variable = varName,
                          stockLabel = stockLabel,
                          PItab = PItab,
                          ICtab = ICtab,
                          MREtab = MREtab,
                          MAREtab = MAREtab,
                          axisLabs = FALSE,
                          bw = bw )
        mfg <- par("mfg")
        # Plot column header if in top row
        if(mfg[1] == 1)
          mtext(side = 3, line = 2, text = MPlabels[mpIdx],cex = 1.5)
        # Plot variable label in right hand margin
        if(mfg[2] == mfg[4])
          mtext( side = 4, line = 1, text = varLabels[varIdx], cex = 1.5,
                 las = 1 )
        
      }
    }
    mtext( side = 1, text = expression(Q(theta)), outer = T, line = 4, font = 2, cex = 1.5)
    mtext( side = 2, text = "Density", outer = T, line = 4, font = 2, cex = 1.5 )
    # mtext( side = 3, text = "AM Configuration", outer = T, line = 1, font = 2, cex = 1.2 )
    # mtext( side = 4, text = "Variable", outer = T, line = 1.5, font = 2, cex = 1.5 )
    dev.off()
  }


  cat("Perf metrics for ", tabNameRoot," complete\n", sep = "")
}

plotPerfMetrics <- function(  scenLabel = "base_1way", tabNameRoot = "pubBase",
                              mpLabel = "UmsyPriorOnly", variable = "Bmsy", 
                              stockLabel = "Stock1", axisLabs = TRUE,
                              PItab = NULL,
                              ICtab = NULL,
                              MREtab = NULL,
                              MAREtab = NULL,
                              bw = FALSE )
{
  tablePath <- "./project/statistics"

  # Read in various perf tables if required, filter to given scenario
  if(is.null(PItab))
  {
    PItab     <- paste(tabNameRoot,"_PI.csv", sep = "" )
    PItab     <- read.csv(file.path(tablePath,PItab), header = T, stringsAsFactors = FALSE)
  }
  # Filter to scenario/stock/MP
  PItab <-  PItab %>%
            dplyr::filter(  scenario == scenLabel,
                            mp == mpLabel,
                            species == stockLabel ) 

  if(is.null(ICtab))
  {
    ICtab     <- paste(tabNameRoot,"_IC.csv", sep = "" )
    ICtab     <- read.csv(file.path(tablePath,ICtab), header = T, stringsAsFactors = FALSE) 
  }
  # Filter to scenario/stock/MP
  ICtab <-  ICtab %>%
            dplyr::filter(  scenario == scenLabel,
                            mp == mpLabel,
                            species == stockLabel ) 

  if(is.null(MREtab))
  {
    MREtab    <- paste(tabNameRoot,"_MRE.csv", sep = "" )
    MREtab    <- read.csv(file.path(tablePath,MREtab), header = T, stringsAsFactors = FALSE )
  }
  # Filter to scenario/stock/MP
  MREtab <- MREtab %>%
            dplyr::filter(  scenario == scenLabel,
                            mp == mpLabel,
                            species == stockLabel ) 
  if(is.null(MAREtab))
  {
    MAREtab   <- paste(tabNameRoot,"_MARE.csv", sep = "" )
    MAREtab   <- read.csv(file.path(tablePath,MAREtab), header = T, stringsAsFactors = FALSE)
  }
  MAREtab <-  MAREtab %>%
              dplyr::filter(  scenario == scenLabel,
                              mp == mpLabel,
                              species == stockLabel ) 

  
  ssVar <- paste("ss",variable, sep = "")
  msVar <- paste("ms",variable, sep = "")

  if(!bw)
    cols <- brewer.pal(3, "Dark2")
  if(bw)
    cols <- c("grey25","grey25")


  # get prediction integral values
  ssPIhist <- hist(PItab[,ssVar], plot = FALSE, breaks = seq(0,1,length = 21), na.rm = T )
  msPIhist <- hist(PItab[,msVar], plot = FALSE, breaks = seq(0,1,length = 21), na.rm = T )
  ssPIdens <- density(PItab[,ssVar], na.rm = T )
  msPIdens <- density(PItab[,msVar], na.rm = T )

  # Get IC
  ssIC <- format(round(ICtab[1,ssVar],3),nsmall = 3)
  msIC <- format(round(ICtab[1,msVar],3),nsmall = 3)

  # Get MARE values
  ssMARE <- format(round(MAREtab[,ssVar],3),nsmall = 3)
  msMARE <- format(round(MAREtab[,msVar],3),nsmall = 3)

  # Get MRE values
  ssMRE <- format(round(MREtab[,ssVar], 3),nsmall = 3)
  msMRE <- format(round(MREtab[,msVar], 3),nsmall = 3)

  # get total reps
  totReps <- ICtab[,"totReps"]
  msTries <- format(signif(ICtab[,"msTries"],2),nsmall = 3)
  ssTries <- format(signif(ICtab[,"ssTries"],2),nsmall = 3)


  plot( x = c(0,1), y = c(0,3),
        type = "n", axes = F, xlab = "",
        ylab = "" )
    if(axisLabs)
    {
      mtext( side = 1, line = 2, text = "P(X|hat(x),se(X))")
      mtext( side = 2, line = 2, text = "Density")
    }
    mfg <- par("mfg")
    if( mfg[2] == 1 )
      axis( side = 2, las = 1, at = seq(0,2.5,by=.5))
    # Plot axis if on the bottom row
    if( mfg[1] == mfg[3] )
      axis( side = 1 )

    box()
    plot(ssPIhist, add =T, col = alpha(cols[1],alpha = .5), freq = FALSE,
          angle = 45, density = 15, lwd = 2 )
    plot(msPIhist, add =T, col = alpha(cols[2],alpha = .5), freq =FALSE,
          lwd = 2 )
    lines(ssPIdens, col = cols[1], lwd = 3, lty = 2 )
    lines(msPIdens, col = cols[2], lwd = 3, lty = 1 )
    # Plot performance metrics
    rect( xleft = 0.68, xright = 1.2, ybottom = 1.85, ytop = 3.2,
          col = "white", border = "black" )
    text( x = c(0.855,0.985), y = 2.9, label = c("SS","MS"), font = 2, cex = 1 )
    text( x = c(0.74,0.855,0.985), y = 2.6, font = 2, label = c("IC", ssIC, msIC), cex = 1)
    text( x = c(0.74,0.855,0.985), y = 2.3, font = 2, label = c("MARE", ssMARE, msMARE), cex = 1 )
    text( x = c(0.74,0.855,0.985), y = 2.0, font = 2, label = c("MRE", ssMRE, msMRE), cex = 1 )
    box()

    # Plot number of tries
    # text(x = 0.8, y = 2.8, label = paste( "totReps = ", totReps, sep = ""))
    # text(x = 0.8, y = 2.6, label = paste( "msTries = ", msTries, sep = ""))
    # text(x = 0.8, y = 2.4, label = paste( "ssTries = ", ssTries, sep = ""))


}

corrTitles <- c("corr(Bt) = 0.1", "corr(Bt) = 0.9", "corr(Ct) = 0.1", "corr(Ct) = 0.9" )

mpTitles <- expression("Single-stock", q, U, q/U)





# plotBCfit()
# Plot fit, CIs and data for a set of simulations. This is made for sim objects 
# that were run on realData, so some of the code is a little path specific.
# Single stock fits are plotted as the first column, taken from the first
# sim object. Takes the first replicate.
# inputs:   sims=numeric indicating blobs to load (alpha order in project folder)
#           data=logical indicating plotting data scaled by estimated q
#           CIs=logical indicating plotting CIs (1se) as polygons around estimate
#           scale= logical indicating whether catch is scaled by estimated MSY
#                   and biomass is scaled by estimated B0
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBCfit <- function(  sims=1, legend=FALSE,
                        data = FALSE,
                        CIs = FALSE,
                        scale = FALSE,
                        titles = expression("","None",q,U[MSY],q/U[MSY]), 
                        MPtitles = FALSE,
                        labSize = 1.5, tickSize = 1.5, 
                        devLabels = TRUE, maxBt = c(80,40,40),
                        saveName = "BCfit_largePE.pdf",
                        savePlot = FALSE )
{

  blobs <- vector( mode = "list", length = length(sims))
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  for( simIdx in 1:length(sims) )
  {
    .loadSim(sims[simIdx])
    blobs[[simIdx]] <- blob
  }

  # Recover blob elements for plotting
  nStocks <- blob$opMod$nS
  nS <- nStocks
  nO <- blob$opMod$nSurv

  if( MPtitles ) titles <- numeric( length = length(sims)+1 )

  titles[1] <- "Single-Stock"
  
  # Set up plot window
  if(savePlot)
    pdf( file.path("./project/figs",saveName), width = 11, height = 8.5 )

  par ( mfcol = c(nStocks,length(sims)+1), 
        mar = c(0,0,0,0), oma = c(5,5,4,2),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )


  for( simIdx in 1:length(sims) )
  {
    blob <- blobs[[simIdx]]

    # Create a stamp from scenario and mp name
    scenario  <- blob$ctrl$scenario
    mp        <- blob$ctrl$mp
    stamp     <- paste(scenario,":",mp,sep="")

    # Species names
    specNames <- blob$ctrl$speciesNames

    # True OM quantities
    Ct    <- blob$om$Ct[1,,]
    I_ost <- blob$om$I_ost[1,,,]

    # Year indexing
    if(!is.null(blob$opMod$fYear)) 
    {
      nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
      fYear <- min(blob$opMod$fYear)
    } else nT <- blob$opMod$nT

    sT <- blob$opMod$fYear - fYear + 1
    years <- fYear:(fYear + nT - 1)

    # Groom I_ost to remove unused data
    I_ost[I_ost < 0] <- NA
    for( s in 1:nS ) 
    { 
      fYear_s <- blob$opMod$fYear[s]
      I_ost[,s,1:(fYear_s-fYear)] <- NA
    }

    # Leading par estimates
    # Single Stock
    ssBmsy  <- blob$am$ss$Bmsy
    ssUmsy  <- blob$am$ss$Umsy
    ssq     <- blob$am$ss$q_os[1,,]
    ssMSY   <- blob$am$ss$msy

    # Multi-stock
    msBmsy  <- blob$am$ms$Bmsy
    msUmsy  <- blob$am$ms$Umsy
    msq     <- blob$am$ms$q_os[1,,]  
    msMSY   <- blob$am$ms$msy

    ssfitRep <- blob$am$ss$fitrep[[1]]
    msfitRep <- blob$am$ms$fitrep[[1]]

    # arrays to hold CIs
    ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
    msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

    msCIs <- blob$am$ms$CIs[[ 1 ]]
    ssCIs <- blob$am$ss$CIs[[1]]
    # Populate CIs
    if( !any(is.na(msCIs)) )
    {
      Btrows <- which( msCIs$par == "lnBt" )
      msBio[ , , 1 ] <- matrix( exp(msCIs[ Btrows, "uCI" ]), nrow = nS, ncol = nT, byrow = FALSE )
      msBio[ , , 2 ] <- matrix( exp(msCIs[ Btrows, "val" ]), nrow = nS, ncol = nT, byrow = FALSE )
      msBio[ , , 3 ] <- matrix( exp(msCIs[ Btrows, "lCI" ]), nrow = nS, ncol = nT, byrow = FALSE ) 
    }

    msBio[ msBio == -1 ] <- NA
    
    for( s in 1:nS )
    {

      if( is.null( ssCIs[[ s ]] ) ) next
      if( any(is.na( ssCIs[[ s ]] )) ) next
      Btrows <- which( ssCIs[[ s ]]$par == "lnBt" )
      ssBio[ s, sT[s]:nT, 1 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "uCI" ]), nrow = 1, ncol = nT-sT[s]+1 )
      ssBio[ s, sT[s]:nT, 2 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "val" ]), nrow = 1, ncol = nT-sT[s]+1 )
      ssBio[ s, sT[s]:nT, 3 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "lCI" ]), nrow = 1, ncol = nT-sT[s]+1 )  
    }

    if( MPtitles ) titles[simIdx+1] <- mp
    
    
    # Recover diagnostics for the fits
    hpdSS <- blob$am$ss$hesspd[1,]
    grdSS <- blob$am$ss$maxGrad[1,]
    hpdMS <- blob$am$ms$hesspd[1]
    grdMS <- blob$am$ms$maxGrad[1]

    # Set colours for each model
    ssCol <- "steelblue"
    msCol <- "salmon"
    legText <- c()
    legPch  <- c()
    legCol <- c()
    legLty  <- c()
    legLwd  <- c()

    require(scales)
    polyCol <- alpha("grey70", alpha = .5)
    rectCol <- "grey40"
    lineCol <- "black"
    pointCol <- "grey30"
  
    if( data ) 
    {
      legText <- c(legText,"Survey 1", "Survey 2")
      legPch  <- c(legPch,1,2)
      legCol <- c(legCol,"grey50","grey50")
      legLty  <- c(legLty,NA,NA)
      legLwd  <- c(legLwd,NA,NA)
    }


    vertLines <- seq(fYear,max(years),by = 10)

    # Plot SS model if first sim
    if( simIdx == 1 )
    { 
      for ( s in 1:nStocks )
      {
        plot( x = c(fYear,max(years) + 2), y = c(0,maxBt[s]), type = "n",
              ylim = c(0,maxBt[s]), ylab = "", axes=FALSE, xlab = "" ,
              main = "" )
          mfg <- par("mfg")
          if( !is.null(titles) & mfg[1] == 1 ) 
            mtext( side = 3, text = titles[simIdx], line = 1, cex = 1.2 )
          
          # plot catch
          rect( xleft = years-.3, xright = years+.3,
                ybottom = 0, ytop = Ct[s,], col = rectCol, border = NA )
          if( mfg[1] == mfg[3])
            axis( side = 1, las = 0, cex.axis = 1.5, at = vertLines )
          if( mfg[2] == 1 )
            axis( side = 2, las = 1, cex.axis = 1.5 )
          abline(v = vertLines, lty = 2, lwd = .8)
          box()

        if (s == 2 & legend) 
          panLegend ( x=0.2,y=1,legTxt=legText,
                      col=legCol, lty = legLty, 
                      lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
          
       
        # Add developer labels
        if( devLabels )
        {
          if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",cex=1.1)
        }
        # Add confidence intervals
        if( CIs )
        {
          yearsPoly <- c(years[sT[s]:nT],rev(years[sT[s]:nT]))
          polygon(  x = yearsPoly, y = c(ssBio[s,sT[s]:nT,1],rev(ssBio[s,sT[s]:nT,3])),
                    col = polyCol, border = NA )
        }
        panLab( x = .2, y = .9, txt = specNames[s],
                  cex = labSize )
        # Plot data
        if( data )
          for( o in 1:nO ) 
            points( x = years, y = I_ost[o,s,]/ssq[o,s], pch = o + 15, 
                    cex = 3, col = pointCol )
        
        # Add point estimates
        lines( x = years[sT[s]:nT], y = ssBio[s,sT[s]:nT,2], lwd = 4 )
      }
    }

    # Plot biomass, actual and estimated, including 2 index series,
    # scaled by model estimated q
    for ( s in 1:nStocks )
    {
      plot( x = c(fYear,max(years) + 2), y = c(0,maxBt[s]), type = "n",
            ylim = c(0,maxBt[s]), ylab = "", axes=FALSE, xlab = "" ,
            main = "" )
        mfg <- par("mfg")
        if( !is.null(titles) & mfg[1] == 1 ) 
          mtext( side = 3, text = titles[simIdx+1], line = 1, cex = 1.2 )
        # plot catch
        rect( xleft = years-.4, xright = years+.4,
              ybottom = 0, ytop = Ct[s,], col = rectCol, border = NA )
        if( mfg[1] == mfg[3])
            axis( side = 1, las = 0, cex.axis = 1.5, at = vertLines )
          if( mfg[2] == 1 )
            axis( side = 2, las = 1, cex.axis = 1.5 )
        abline(v = vertLines, lty = 2, lwd = .8)
        box()
      if (s == 2 & legend) 
        panLegend ( x=0.2,y=1,legTxt=legText,
                    col=legCol, lty = legLty, 
                    lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
        # panLab( x = .8, y = .9, txt = specNames[s],
        #         cex = labSize )
      
      # Add developer labels
      if( devLabels )
      {
        if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",cex=1.1)  
      }
      # Add confidence intervals
      if( CIs )
      {
        yearsPoly <- c(years[sT[s]:nT],rev(years[sT[s]:nT]))
        polygon(  x = yearsPoly, y = c(msBio[s,sT[s]:nT,1],rev(msBio[s,sT[s]:nT,3])),
                  col = polyCol, border = NA )
      }
      # Plot data
      if( data )
        for( o in 1:nO ) 
          points( x = years, y = I_ost[o,s,]/msq[o,s], pch = o+15, 
                  cex = 3, col = pointCol )
      

      # Add point estimates
      lines( x = years[sT[s]:nT], y = msBio[s,sT[s]:nT,2], lwd = 4 )

      
    }
    mtext ( text = "Year", outer = TRUE, side = 1, cex = labSize, line = 3)
    mtext ( text = "Catch and Biomass (kt)", outer = TRUE, side = 2, line = 3, 
            las = 0, cex = labSize )
    
  }
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.9),padj=2,col="grey50",cex=0.8 )

  if(savePlot)
    dev.off()
}

# plotBCfit()
# Plot fit, CIs and data for a set of simulations. This is made for sim objects 
# that were run on realData, so some of the code is a little path specific.
# Single stock fits are plotted as the first column, taken from the first
# sim object. Takes the first replicate.
# inputs:   sims=numeric indicating blobs to load (alpha order in project folder)
#           data=logical indicating plotting data scaled by estimated q
#           CIs=logical indicating plotting CIs (1se) as polygons around estimate
#           scale= logical indicating whether catch is scaled by estimated MSY
#                   and biomass is scaled by estimated B0
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBCsim <- function(  sims=1, rep = 1, legend=FALSE,
                        data = FALSE,
                        CIs = FALSE,
                        scale = FALSE,
                        titles = expression("Single Stock","None", q, r, q/r ), 
                        MPtitles = FALSE,
                        labSize = 2, tickSize = 2, 
                        devLabels = TRUE, maxBt = NULL )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  .loadSim(sims[1])

  # Recover blob elements for plotting
  nStocks <- blob$opMod$nS
  nS <- nStocks
  nO <- blob$opMod$nSurv

  # Set up plot window

  par ( mfcol = c(nStocks,length(sims)+1), mar = c(1,2,2,2), oma = c(4.5,4.5,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )

  if( MPtitles ) titles <- character( length = length(sims)+1 )

  titles[1] <- "Single-Stock"

  for( simIdx in 1:length(sims) )
  {

    .loadSim( sims[simIdx] )
  

    # Create a stamp from scenario and mp name
    scenario  <- blob$ctrl$scenario
    mp        <- blob$ctrl$mp
    stamp     <- paste(scenario,":",mp,sep="")

    # Species names
    specNames <- paste("Stock ", 1:nS, sep = "" )

    # True OM quantities
    Ct    <- blob$om$Ct[rep,,]
    I_ost <- blob$om$I_ost[rep,,,]

    # Year indexing
    if(!is.null(blob$opMod$fYear)) 
    {
      nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
      fYear <- min(blob$opMod$fYear)
    } else nT <- blob$opMod$nT

    sT <- blob$opMod$fYear - fYear + 1
    years <- fYear:(fYear + nT - 1)

    # Groom I_ost to remove unused data
    I_ost[I_ost < 0] <- NA
    for( s in 1:nS ) 
    { 
      fYear_s <- blob$opMod$fYear[s]
      I_ost[,s,1:(fYear_s-fYear)] <- NA
    }

    omBt    <- blob$om$Bt[rep,,] 

    # Leading par estimates
    # Single Stock
    ssBmsy  <- blob$am$ss$Bmsy[rep,]
    ssUmsy  <- blob$am$ss$Umsy[rep,]
    ssq     <- blob$am$ss$q_os[rep,,]
    ssMSY   <- blob$am$ss$msy[rep,]

    # Multi-stock
    msBmsy  <- blob$am$ms$Bmsy[rep,]
    msUmsy  <- blob$am$ms$Umsy[rep,]
    msq     <- blob$am$ms$q_os[rep,,]  
    msMSY   <- blob$am$ms$msy[rep,]


    # arrays to hold CIs
    ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
    msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

    msCIs <- blob$am$ms$CIs[[ rep ]]
    ssCIs <- blob$am$ss$CIs[[ rep ]]
    # Populate CIs
    if( !is.null(msCIs) )
    {
      if( !any(is.na( msCIs )) )
      {
        Btrows <- which( msCIs$par == "lnBt" )    
        
        msBio[ , , 1 ] <- matrix( exp(msCIs[ Btrows, "uCI" ]), nrow = nS, ncol = nT, byrow = FALSE )
        msBio[ , , 2 ] <- matrix( exp(msCIs[ Btrows, "val" ]), nrow = nS, ncol = nT, byrow = FALSE )
        msBio[ , , 3 ] <- matrix( exp(msCIs[ Btrows, "lCI" ]), nrow = nS, ncol = nT, byrow = FALSE ) 
      }
    }

    msBio[ msBio == -1 ] <- NA
    
    for( s in 1:nS )
    {

      if( is.null( ssCIs[[ s ]] ) ) next
      if( any(is.na( ssCIs[[ s ]] )) ) next
      Btrows <- which( ssCIs[[ s ]]$par == "lnBt" )
      ssBio[ s, sT[s]:nT, 1 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "uCI" ]), nrow = 1, ncol = nT-sT[s]+1 )
      ssBio[ s, sT[s]:nT, 2 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "val" ]), nrow = 1, ncol = nT-sT[s]+1 )
      ssBio[ s, sT[s]:nT, 3 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "lCI" ]), nrow = 1, ncol = nT-sT[s]+1 )  
    }

    if( MPtitles ) titles[simIdx+1] <- mp
    
    
    # Recover diagnostics for the fits
    hpdSS <- blob$am$ss$hesspd[rep,]
    grdSS <- blob$am$ss$maxGrad[rep,]
    hpdMS <- blob$am$ms$hesspd[rep]
    grdMS <- blob$am$ms$maxGrad[rep]

    # Set colours for each model
    ssCol <- "grey50"
    msCol <- "grey50"
    legText <- c()
    legPch  <- c()
    legCol <- c()
    legLty  <- c()
    legLwd  <- c()

    require(scales)
    polyCol <- alpha("grey70", alpha = .5)
    msPolyCol <- "grey70"
    ssPolyCol <- "grey70"


  
    if( data ) 
    {
      legText <- c(legText,"Survey 1", "Survey 2")
      legPch  <- c(legPch,1,2)
      legCol <- c(legCol,"grey50","grey50")
      legLty  <- c(legLty,NA,NA)
      legLwd  <- c(legLwd,NA,NA)
    }
    if(is.null(maxBt))
      maxBt <- rep(NA, nS)

    # Plot SS model if first sim
    if( simIdx == 1 )
    { 
      for ( s in 1:nStocks )
      {
        if(is.na(maxBt[s]))
          maxBt[s] <- max(ssBio[s,sT[s]:nT,], msBio[s,sT[s]:nT,], omBt[s,sT[s]:nT], na.rm = T)
        plot( x = c(fYear,max(years)), y = c(0,maxBt[s]), type = "n",
              ylim = c(0,maxBt[s]), ylab = "", axes=FALSE, xlab = "" ,
              main = "" )
          if( !is.null(titles) & s == 1 ) title( main = titles[simIdx] )
          # Plot confidence intervals
          if( CIs )
          {
            yearsPoly <- c(years[sT[s]:nT],rev(years[sT[s]:nT]))
            polygon(  x = yearsPoly, y = c(ssBio[s,sT[s]:nT,1],rev(ssBio[s,sT[s]:nT,3])),
                      col = ssPolyCol, border = NA )
          }
          # plot catch
          rect( xleft = years-.4, xright = years+.4,
                ybottom = 0, ytop = Ct[s,], col = "grey10", border = NA )
          if(s == nS) axis( side = 1, las = 0, cex.axis = tickSize )
          axis( side = 2, las = 1, cex.axis = tickSize )
        if (s == 2 & legend) 
          panLegend ( x=0.2,y=1,legTxt=legText,
                      col=legCol, lty = legLty, 
                      lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
          panLab( x = .2, y = .9, txt = specNames[s],
                  cex = labSize )
       
        # Add developer labels
        if( devLabels )
        {
          if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",cex=1.1)
        }
        # Add confidence intervals
        # Plot OM Bt
        lines( x = years[sT[s]:nT], y = omBt[s,sT[s]:nT], lwd = 3 )
        # Add point estimates
        lines( x = years[sT[s]:nT], y = ssBio[s,sT[s]:nT,2], lwd = 3, col = ssCol, lty = 2 )

         # Plot data
        if( data )
        {
          for( o in 1:nO ) 
            points( x = years, y = I_ost[o,s,]/ssq[o,s], pch = o, cex = 1.5, col = "grey10" )
        }
      }
    }

    # Plot biomass, actual and estimated, including 2 index series,
    # scaled by model estimated q
    for ( s in 1:nStocks )
    {
      plot( x = c(fYear,max(years)), y = c(0,maxBt[s]), type = "n",
            ylim = c(0,maxBt[s]), ylab = "", axes=FALSE, xlab = "" ,
            main = "" )
        if( !is.null(titles) & s == 1 ) title( main = titles[simIdx+1] )
        # plot catch
        rect( xleft = years-.4, xright = years+.4,
              ybottom = 0, ytop = Ct[s,], col = "grey10", border = NA )
        if(s == nS) axis( side = 1, las = 0, cex.axis = tickSize )
        axis( side = 2, las = 1, cex.axis = tickSize )
      if (s == 2 & legend) 
        panLegend ( x=0.2,y=1,legTxt=legText,
                    col=legCol, lty = legLty, 
                    lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
        # panLab( x = .8, y = .9, txt = specNames[s],
        #         cex = labSize )
      
      # Add developer labels
      if( devLabels )
      {
        if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",cex=1.1)  
      }
      # Add confidence intervals
      if( CIs )
      {
        yearsPoly <- c(years[sT[s]:nT],rev(years[sT[s]:nT]))
        polygon(  x = yearsPoly, y = c(msBio[s,sT[s]:nT,1],rev(msBio[s,sT[s]:nT,3])),
                  col = msPolyCol, border = NA )
      }
      # Plot data
      if( data )
      {
        for( o in 1:nO ) 
          points( x = years, y = I_ost[o,s,]/msq[o,s], pch = o, cex = 1.5, col = "grey10" )
      }
      # Plot OM Bt
      lines( x = years[sT[s]:nT], y = omBt[s,sT[s]:nT], lwd = 3 )
      # Add point estimates
      lines( x = years[sT[s]:nT], y = msBio[s,sT[s]:nT,2], lwd = 3, col = msCol, lty = 2 )

      
    }
    mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5, cex = labSize)
    mtext ( text = "Biomass (kt)", outer = TRUE, side = 2, line = 2, las = 0, cex = labSize )
    
  }
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.9),padj=2,col="grey50",cex=0.8 )
}


# plotBCfit()
# Plot fit, CIs and data for a set of simulations. This is made for sim objects 
# that were run on realData, so some of the code is a little path specific.
# Single stock fits are plotted as the first column, taken from the first
# sim object. Takes the first replicate.
# inputs:   sims=numeric indicating blobs to load (alpha order in project folder)
#           data=logical indicating plotting data scaled by estimated q
#           CIs=logical indicating plotting CIs (1se) as polygons around estimate
#           scale= logical indicating whether catch is scaled by estimated MSY
#                   and biomass is scaled by estimated B0
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBCsimReps <- function(  sims=1, legend=FALSE,
                            saveFileRoot = "BCsim",
                            data = FALSE,
                            CIs = TRUE,
                            scale = FALSE,
                            titles = expression("Single Stock","None", q, r, q/r ), 
                            MPtitles = FALSE,
                            labSize = 2, tickSize = 2, 
                            devLabels = TRUE, maxBt = NULL )
{
  blobList <- vector(mode = "list", length = length(sims))

  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  for( simIdx in 1:length(sims))
  {
    .loadSim(sims[simIdx])
    blobList[[simIdx]] <- blob
  }
  
  # Recover blob elements for plotting
  nStocks <- blobList[[1]]$opMod$nS
  nS <- nStocks
  nO <- blobList[[1]]$opMod$nSurv


  # Set up plot window
  for(rep in 1:100)
  {
    saveFile <- paste(saveFileRoot,"_rep", rep,".pdf",sep = "")
    pdf( saveFile, width = 16, height = 9 )
    par ( mfcol = c(nStocks,length(sims)+1), mar = c(0,0,0,0), oma = c(6.5,6.5,4,2.5),
          las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )

    if( MPtitles ) titles <- character( length = length(sims)+1 )

    titles[1] <- "Single-Stock"

    for( simIdx in 1:length(sims) )
    {
      blob <- blobList[[simIdx]]
      # Create a stamp from scenario and mp name
      scenario  <- blob$ctrl$scenario
      mp        <- blob$ctrl$mp
      stamp     <- paste(scenario,":",mp,sep="")

      # Species names
      specNames <- paste("Stock ", 1:nS, sep = "" )

      # True OM quantities
      Ct    <- blob$om$Ct[rep,,]
      I_ost <- blob$om$I_ost[rep,,,]

      # Year indexing
      if(!is.null(blob$opMod$fYear)) 
      {
        nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
        fYear <- min(blob$opMod$fYear)
      } else nT <- blob$opMod$nT

      sT <- blob$opMod$fYear - fYear + 1
      years <- fYear:(fYear + nT - 1)

      # Groom I_ost to remove unused data
      I_ost[I_ost < 0] <- NA
      for( s in 1:nS ) 
      { 
        fYear_s <- blob$opMod$fYear[s]
        I_ost[,s,1:(fYear_s-fYear)] <- NA
      }

      omBt    <- blob$om$Bt[rep,,] 

      # Leading par estimates
      # Single Stock
      ssBmsy  <- blob$am$ss$Bmsy[rep,]
      ssUmsy  <- blob$am$ss$Umsy[rep,]
      ssq     <- blob$am$ss$q_os[rep,,]
      ssMSY   <- blob$am$ss$msy[rep,]

      # Multi-stock
      msBmsy  <- blob$am$ms$Bmsy[rep,]
      msUmsy  <- blob$am$ms$Umsy[rep,]
      msq     <- blob$am$ms$q_os[rep,,]  
      msMSY   <- blob$am$ms$msy[rep,]


      # arrays to hold CIs
      ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
      msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

      msCIs <- blob$am$ms$CIs[[ rep ]]
      ssCIs <- blob$am$ss$CIs[[ rep ]]
      # Populate CIs
      if( !is.null(msCIs) )
      {
        if( !any(is.na( msCIs )) )
        {
          Btrows <- which( msCIs$par == "lnBt" )    
          
          msBio[ , , 1 ] <- matrix( exp(msCIs[ Btrows, "uCI" ]), nrow = nS, ncol = nT, byrow = FALSE )
          msBio[ , , 2 ] <- matrix( exp(msCIs[ Btrows, "val" ]), nrow = nS, ncol = nT, byrow = FALSE )
          msBio[ , , 3 ] <- matrix( exp(msCIs[ Btrows, "lCI" ]), nrow = nS, ncol = nT, byrow = FALSE ) 
        }
      }

      msBio[ msBio == -1 ] <- NA
      
      for( s in 1:nS )
      {

        if( is.null( ssCIs[[ s ]] ) ) next
        if( any(is.na( ssCIs[[ s ]] )) ) next
        Btrows <- which( ssCIs[[ s ]]$par == "lnBt" )
        ssBio[ s, sT[s]:nT, 1 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "uCI" ]), nrow = 1, ncol = nT-sT[s]+1 )
        ssBio[ s, sT[s]:nT, 2 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "val" ]), nrow = 1, ncol = nT-sT[s]+1 )
        ssBio[ s, sT[s]:nT, 3 ] <- matrix( exp(ssCIs[[ s ]][ Btrows, "lCI" ]), nrow = 1, ncol = nT-sT[s]+1 )  
      }

      if( MPtitles ) titles[simIdx+1] <- mp
      
      
      # Recover diagnostics for the fits
      hpdSS <- blob$am$ss$hesspd[rep,]
      grdSS <- blob$am$ss$maxGrad[rep,]
      hpdMS <- blob$am$ms$hesspd[rep]
      grdMS <- blob$am$ms$maxGrad[rep]

      # Set colours for each model
      ssCol <- "grey50"
      msCol <- "grey50"
      legText <- c()
      legPch  <- c()
      legCol <- c()
      legLty  <- c()
      legLwd  <- c()

      require(scales)
      polyCol <- alpha("grey70", alpha = .5)
      msPolyCol <- "grey70"
      ssPolyCol <- "grey70"


    
      if( data ) 
      {
        legText <- c(legText,"Survey 1", "Survey 2")
        legPch  <- c(legPch,1,2)
        legCol <- c(legCol,"grey50","grey50")
        legLty  <- c(legLty,NA,NA)
        legLwd  <- c(legLwd,NA,NA)
      }
      if(is.null(maxBt))
        maxBt <- rep(NA, nS)

      # Plot SS model if first sim
      if( simIdx == 1 )
      { 
        for ( s in 1:nStocks )
        {
          if(is.na(maxBt[s]))
            maxBt[s] <- max(ssBio[s,sT[s]:nT,], msBio[s,sT[s]:nT,], omBt[s,sT[s]:nT], na.rm = T)
          if(is.na(maxBt[s]))
            browser()
          plot( x = c(fYear,max(years)), y = c(0,maxBt[s]), type = "n",
                ylim = c(0,maxBt[s]), ylab = "", axes=FALSE, xlab = "" ,
                main = "" )
            mfg <- par("mfg")
            # Plot axes if asked for
            if( mfg[1] == mfg[3] )
              axis( side = 1, cex.axis = tickSize )
            if( mfg[2] == 1 )
              axis( side = 2, las = 1, cex.axis = tickSize )
            box()
            if( !is.null(titles) & s == 1 ) 
              mtext( side = 3, text = titles[simIdx], font = 2, cex = 1.5 )
            # Plot confidence intervals
            if( CIs )
            {
              yearsPoly <- c(years[sT[s]:nT],rev(years[sT[s]:nT]))
              polygon(  x = yearsPoly, y = c(ssBio[s,sT[s]:nT,1],rev(ssBio[s,sT[s]:nT,3])),
                        col = ssPolyCol, border = NA )
            }
            # plot catch
            rect( xleft = years-.4, xright = years+.4,
                  ybottom = 0, ytop = Ct[s,], col = "grey10", border = NA )
            if(s == nS) axis( side = 1, las = 0, cex.axis = tickSize )
            
          if (s == 2 & legend) 
            panLegend ( x=0.2,y=1,legTxt=legText,
                        col=legCol, lty = legLty, 
                        lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
            panLab( x = .25, y = .9, txt = specNames[s],
                    cex = labSize )
         
          # Add developer labels
          if( devLabels )
          {
            if( !is.na(hpdSS[s]) ) 
              if( hpdSS[s] )
                panLab (x=0.9,y=0.9,txt="h",cex=1.1)
          }
          # Add confidence intervals
          # Plot OM Bt
          lines( x = years[sT[s]:nT], y = omBt[s,sT[s]:nT], lwd = 3 )
          # Add point estimates
          lines( x = years[sT[s]:nT], y = ssBio[s,sT[s]:nT,2], lwd = 3, col = ssCol, lty = 2 )

           # Plot data
          if( data )
          {
            for( o in 1:nO ) 
              points( x = years, y = I_ost[o,s,]/ssq[o,s], pch = o, cex = 1.5, col = "grey10" )
          }
        }
      }

      # Plot biomass, actual and estimated, including 2 index series,
      # scaled by model estimated q
      for ( s in 1:nStocks )
      {
        if(is.na(maxBt[s]))
            browser()
        plot( x = c(fYear,max(years)), y = c(0,maxBt[s]), type = "n",
              ylim = c(0,maxBt[s]), ylab = "", axes=FALSE, xlab = "" ,
              main = "" )
          mfg <- par("mfg")
          # Plot axes if asked for
          if( mfg[1] == mfg[3] )
            axis( side = 1, cex.axis = tickSize )
          if( mfg[2] == 1 )
            axis( side = 2, las = 1, cex.axis = tickSize )
          box()
          if( !is.null(titles) & s == 1 ) 
            mtext( side = 3, text = titles[simIdx+1], font = 2, cex = 1.5 )
          # plot catch
          rect( xleft = years-.4, xright = years+.4,
                ybottom = 0, ytop = Ct[s,], col = "grey10", border = NA )
        if (s == 2 & legend) 
          panLegend ( x=0.2,y=1,legTxt=legText,
                      col=legCol, lty = legLty, 
                      lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
          # panLab( x = .8, y = .9, txt = specNames[s],
          #         cex = labSize )
        
        # Add developer labels
        if( devLabels )
        {
          if ( hpdMS & !is.na(hpdMS) ) 
            panLab (x=0.9,y=0.85,txt="h",cex=1.1)  
        }
        # Add confidence intervals
        if( CIs )
        {
          yearsPoly <- c(years[sT[s]:nT],rev(years[sT[s]:nT]))
          polygon(  x = yearsPoly, y = c(msBio[s,sT[s]:nT,1],rev(msBio[s,sT[s]:nT,3])),
                    col = msPolyCol, border = NA )
        }
        # Plot data
        if( data )
        {
          for( o in 1:nO ) 
            points( x = years, y = I_ost[o,s,]/msq[o,s], pch = o, cex = 1.5, col = "grey10" )
        }
        # Plot OM Bt
        lines( x = years[sT[s]:nT], y = omBt[s,sT[s]:nT], lwd = 3 )
        # Add point estimates
        lines( x = years[sT[s]:nT], y = msBio[s,sT[s]:nT,2], lwd = 3, col = msCol, lty = 2 )

        
      }
      mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5, cex = labSize,
               line = 2.5)
      mtext ( text = "Biomass (kt)", outer = TRUE, side = 2, line = 4, las = 0, cex = labSize )
      
    }
    if( devLabels )
      mtext ( text = c(stamp),side=1, outer = TRUE, 
              at = c(0.9),padj=2,col="grey50",cex=0.8,
              line = 3 )
    dev.off()
  }
}



# plotBC()
# Plots a given replicate's true and estimated time series of biomass and 
# catch for all stocks, with sims as columns
# inputs:   rep=replicate number for plotting; 
#           sims=numeric indicating blobs to load (alpha order in project folder)
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBC <- function( rep = 1, sims=1, legend=TRUE,
                    data = FALSE, labSize = 2, tickSize = 2,
                    fYear = 1988, devLabels = TRUE, CIs = FALSE,
                    scale = FALSE, est = TRUE,
                    titles = NULL, MPtitles = TRUE,
                    nStocks = NULL, MS = TRUE, SS = TRUE   )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  .loadSim(sims[1])

  # Recover blob elements for plotting
  nS <- blob$opMod$nS
  nO <- blob$opMod$nSurv

  if( !is.null( nStocks) ) nStocks <- min( nS, nStocks ) else nStocks <- nS
  # Set up plot window

  par ( mfcol = c(nStocks,length(sims)), mar = c(1,2,2,0), oma = c(4.5,4.5,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )

  if( MPtitles ) titles <- numeric( length = length(sims) )

  for( simIdx in 1:length(sims) )
  {
    .loadSim( sims[simIdx] )
  

    # Create a stamp from scenario and mp name
    scenario  <- blob$ctrl$scenario
    mp        <- blob$ctrl$mp
    stamp     <- paste(scenario,":",mp,sep="")
    repCount  <- paste("Replicate ",rep,"/",blob$ctrl$nReps,sep="")

    # Species names
    specNames <- blob$ctrl$speciesNames

    # True OM quantities
    omBt  <- blob$om$Bt[rep,,]
    Ct    <- blob$om$Ct[rep,,]
    I_ost <- blob$om$I_ost[rep,,,]
    Ut    <- blob$om$Ut[rep,,]
    q_os  <- blob$om$q_os[rep,,]
    if(any(is.na(q_os))) q_os <- matrix(1,nrow=nO, nS)
    
    # Single species model
    ssBt  <- blob$am$ss$Bt[rep,,]
    ssq   <- blob$am$ss$q_os[rep,,]

    # Multispecies model
    msBt  <- blob$am$ms$Bt[rep,,]
    msBt[msBt < 0] <- NA
    msq   <- blob$am$ms$q_os[rep,,]  

    if( MPtitles ) titles[simIdx] <- mp

    # Year indexing
    if(!is.null(blob$opMod$fYear)) 
    {
      nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
      fYear <- min(blob$opMod$fYear)
    } else nT <- blob$opMod$nT
    # browser()
    sT <- blob$opMod$fYear - fYear + 1
    years <- fYear:(fYear + nT - 1)

    # Groom I_ost to remove unused data
    I_ost[I_ost < 0] <- NA
    for( s in 1:nS ) 
    { 
      fYear_s <- blob$opMod$fYear[s]
      I_ost[,s,1:(fYear_s-fYear)] <- NA
    }

    # Estimated Ut
    ssUt <- Ct / ssBt
    msUt <- Ct / msBt

    if( CIs )
    {
      # arrays to hold CIs
      ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
      msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

      msCIs <- blob$am$ms$CIs[[ rep ]]
      ssCIs <- blob$am$ss$CIs[[rep]]

      # browser()
      if( !any(is.na(msCIs)) )
      {
        Btrows <- which( msCIs$par == "Bt" )
        msBio[ , , 1 ] <- matrix( msCIs[ Btrows , "uCI" ], nrow = nS, ncol = nT, byrow = FALSE )
        msBio[ , , 2 ] <- matrix( msCIs[ Btrows, "val" ], nrow = nS, ncol = nT, byrow = FALSE )
        msBio[ , , 3 ] <- matrix( msCIs[ Btrows, "lCI" ], nrow = nS, ncol = nT, byrow = FALSE ) 
      }
      for( s in 1:nS )
      {

        if( is.null( ssCIs[[ s ]] ) ) next
        if( any(is.na( ssCIs[[ s ]] )) ) next
        Btrows <- which( ssCIs[[ s ]]$par == "Bt" )
        ssBio[ s, sT[s]:nT, 1 ] <- matrix( ssCIs[[ s ]][ Btrows, "uCI" ], nrow = 1, ncol = nT-sT[s]+1 )
        ssBio[ s, sT[s]:nT, 2 ] <- matrix( ssCIs[[ s ]][ Btrows, "val" ], nrow = 1, ncol = nT-sT[s]+1 )
        ssBio[ s, sT[s]:nT, 3 ] <- matrix( ssCIs[[ s ]][ Btrows, "lCI" ], nrow = 1, ncol = nT-sT[s]+1 )  
      }
    }

    # Pull Bmsy and Umsy from the blob
    Umsy <- blob$opMod$Umsy
    Bmsy <- blob$opMod$Bmsy
    MSY  <- Umsy*Bmsy

    if( scale )
    {
      # Now loop over species and scale
      for( s in 1:nS )
      {
        omBt[s,]  <- omBt[s,]/Bmsy[s]/2
        Ct[s,]    <- (Ct[s,]/MSY[s])
        Ut[s,]    <- Ut[s,]/Umsy[s]
        I_ost[,s,]<- I_ost[,s,]/Bmsy[s]/2
        ssBt[s,]  <- ssBt[s,]/Bmsy[s]/2
        msBt[s,]  <- msBt[s,]/Bmsy[s]/2
        msUt[s,]    <- msUt[s,]/Umsy[s]
        ssUt[s,]    <- ssUt[s,]/Umsy[s]
        if(CIs)
        {
          ssBio[s,,]<- ssBio[s,,]/Bmsy[s]/2
          msBio[s,,]<- msBio[s,,]/Bmsy[s]/2  
        }
        Bmsy[s]   <- 1
      }  
    }
    
    
    # Recover diagnostics for the fits
    hpdSS <- blob$am$ss$hesspd[rep,]
    grdSS <- blob$am$ss$maxGrad[rep,]
    hpdMS <- blob$am$ms$hesspd[rep]
    grdMS <- blob$am$ms$maxGrad[rep]

    # Set colours for each model
    ssCol <- "steelblue"
    msCol <- "salmon"
    legText <- c()
    legPch  <- c()
    legCol <- c()
    legLty  <- c()
    legLwd  <- c()
    
    if( est ) 
    {
      legText <- c(legText,"ss","ms")
      legPch  <- c(legPch,NA,NA)
      legCol <- c(legCol,ssCol,msCol)
      legLty  <- c(legLty,2,3)
      legLwd  <- c(legLwd,3,3)
    }
    if( data ) 
    {
      legText <- c(legText,"Survey 1", "Survey 2")
      legPch  <- c(legPch,1,2)
      legCol <- c(legCol,"grey50","grey50")
      legLty  <- c(legLty,NA,NA)
      legLwd  <- c(legLwd,NA,NA)
    }



    # Plot biomass, actual and estimated, including 2 index series,
    # scaled by model estimated q
    maxBt <- 1
    for ( s in 1:nStocks )
    {
      if( !scale ) maxBt <- 1.1*max ( maxBt, omBt[s,], 2*Bmsy[s] ,na.rm=TRUE)
      if( (!scale) & CIs ) maxBt <- max ( maxBt, msBio[s,,], ssBio[s,,], na.rm = T )
      if( scale | data ) maxBt <- max(2,maxBt)
      if(!MS & !SS & all(q_os == 1)) maxBt <- 1.1*max(I_ost, na.rm = T)
      plot    ( x = c(fYear,max(years)), y = c(0,maxBt), type = "n",
                ylim = c(0,maxBt), ylab = "", axes=FALSE, xlab = "" ,
                main = "" )
        if( !is.null(titles) & s == 1 ) title( main = titles[simIdx] )
        # plot catch
        rect( xleft = years-.4, xright = years+.4,
              ybottom = 0, ytop = Ct[s,], col = "grey80", border = NA )
        if(s == nS) axis( side = 1, las = 0, cex.axis = tickSize )
        axis( side = 2, las = 1, cex.axis = tickSize )
      if (s == 2 & legend) 
        panLegend ( x=0.2,y=1,legTxt=legText,
                    col=legCol, lty = legLty, 
                    lwd = legLwd, pch = legPch, cex=c(2), bty="n" )
        panLab( x = .8, y = .9, txt = specNames[s],
                cex = labSize )
      # Plot data
      if( data )
      {
        for( o in 1:nO ) 
        {
          points( x = years, y = I_ost[o,s,]/q_os[o,s], pch = o, cex = 1.5, col = "grey50" )
          if( SS ) 
            points( x = years, y = I_ost[o,s,]/ssq[o,s], pch = o, cex = 1.5, col = "grey50" )
          if( MS ) 
            points( x = years, y = I_ost[o,s,]/msq[o,s], pch = o, cex = 1.5, col = "grey50" )

        }
      }
      # Add developer labels
      if( devLabels )
      {
        if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=1.1)
        if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=1.1)  
      }
      # Add OM biomass
      lines   ( x = years, y = omBt[s,], col = "black", lwd = 3)
      if(est)
      {
        if( SS ) 
          lines( x = years, y = ssBt[s,], col = ssCol, lwd = 3, lty = 2 )
        if( MS ) 
          lines( x = years, y = msBt[s,], col = msCol, lwd = 3, lty = 3 )
        if( CIs )
        {
          if( SS )
          {
            lines( x = years[sT[s]:nT], y = ssBio[s,sT[s]:nT,1], col = ssCol, lwd = 1, lty = 2)
            lines( x = years[sT[s]:nT], y = ssBio[s,sT[s]:nT,3], col = ssCol, lwd = 1, lty = 2)  
          }
          if( MS )
          {
            lines( x = years[sT[s]:nT], y = msBio[s,sT[s]:nT,1], col = msCol, lwd = 1, lty = 3)
            lines( x = years[sT[s]:nT], y = msBio[s,sT[s]:nT,3], col = msCol, lwd = 1, lty = 3)  
          }
          
        }  
      }
      
    }
    mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5, cex = labSize)
    mtext ( text = "Biomass (kt)", outer = TRUE, side = 2, line = 2, las = 0, cex = labSize )
    if( devLabels )
      mtext ( text = c(stamp,repCount),side=1, outer = TRUE, 
              at = c(0.9,0.1),padj=2,col="grey50",cex=0.8 )
  }
}

# plotBCU()
# Function that will open a supplied saved blob object and plot a 
# given replicate's true and estimated time series of biomass, 
# catch, exploitation and IoA data. 
# inputs:   rep=replicate number for plotting; 
#           sim=number indicating blob to load (alpha order in project folder)
#           folder=name of folder/blob file (supercedes sim number)
# output:   NULL
# usage:    post-simulation run, plotting performance
plotBCU <- function ( rep = 1, sim=1, legend=TRUE,
                      data = FALSE, labSize = 2, tickSize = 1.2,
                      fYear = 1988, devLabels = TRUE, CIs = FALSE,
                      scale = FALSE )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")
  repCount  <- paste("Replicate ",rep,"/",blob$ctrl$nReps,sep="")

  # Recover blob elements for plotting
  nS <- blob$opMod$nS
  

  # Species names
  specNames <- blob$ctrl$speciesNames

  # True OM quantities
  omBt  <- blob$om$Bt[rep,,]
  Ct    <- blob$om$Ct[rep,,]
  It    <- blob$om$It[rep,,]
  Ut    <- blob$om$Ut[rep,,]
  
  # Single species model
  ssBt  <- blob$am$ss$Bt[rep,,]
  # ssq   <- blob$am$ss$q[rep,]

  # Multispecies model
  msBt  <- blob$am$ms$Bt[rep,,]
  msBt[msBt < 0] <- NA
  # msq   <- blob$am$ms$q[rep,]  

  # Year indexing
  if(!is.null(blob$opMod$fYear)) 
  {
    nT <- blob$opMod$lYear - min(blob$opMod$fYear) + 1
    fYear <- min(blob$opMod$fYear)
  } else nT <- blob$opMod$nT
  # browser()
  years <- fYear:(fYear + nT - 1)

  # Estimated Ut
  ssUt <- Ct / ssBt
  msUt <- Ct / msBt

  if( CIs )
  {
    # arrays to hold CIs
    ssBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )
    msBio <- array( NA, dim = c(nS,nT,3), dimnames = list(1:nS,1:nT,c("uCI","Bst","lCI")) )

    # Fill MS model biomass
    if( !is.null( blob$am$ms$CIs[[ rep ]] ) )
    {
      msCIs <- blob$am$ms$CIs[[ rep ]]
    } 
    else {
      msStdErr <- blob$am$ms$sdrep[[ rep ]]
      if( is.null( msStdErr ) ) 
        msCIs <- NA
      else {  
        msStdErr <- summary( msStdErr )
        colnames( msStdErr ) <- c( "val", "se" )
        msCIs <-  msStdErr %>%
                  as.data.frame() %>%
                  mutate( par = rownames(msStdErr),
                          lCI = val - 1.645*se,
                          uCI = val + 1.645*se ) %>%
                  dplyr::select( par, val, se, lCI, uCI )
      }
    }
    # browser()
    if( !is.na(msCIs) )
    {
      Btrows <- which( msCIs$par == "Bt" )
      msBio[ , , 1 ] <- matrix( msCIs[ Btrows , "uCI" ], nrow = nS, ncol = nT, byrow = FALSE )
      msBio[ , , 2 ] <- matrix( msCIs[ Btrows, "val" ], nrow = nS, ncol = nT, byrow = FALSE )
      msBio[ , , 3 ] <- matrix( msCIs[ Btrows, "lCI" ], nrow = nS, ncol = nT, byrow = FALSE ) 
    }

    # Now fill SS model biomas
    if( !is.null(blob$am$ss$CIs) ) 
    {
      ssCIs <- blob$am$ss$CIs[[rep]]
    }
    else {
      # Create CIs from sd report objects
      ssCIs <- vector( mode = "list", length = 3 )
      ssSDreps <- blob$am$ss$sdrep[[rep]]
      # The following makes up for some coding errors in earlier sims
      if( length( ssSDreps ) == 1 ) k <- nS else k <- 1
      # loop over valid sd objects
      for( s in k:nS )
      {
        ssStdErr <- ssSDreps[[ s ]]
        # Save single species sd rep obj
        if( is.null(ssStdErr) ) 
          ssCIs[[ s ]] <- NA
        else {
        ssStdErr <- summary( ssStdErr )
        colnames( ssStdErr ) <- c("val","se")
        ssCIs[[ s ]] <-   msStdErr %>%
                          as.data.frame() %>%
                          mutate( par = rownames(ssStdErr),
                                  lCI = val - 1.645*se,
                                  uCI = val + 1.645*se ) %>%
                                  dplyr::select( par, val, se, lCI, uCI )    
        }
      }
    }
    for( s in 1:nS )
    {
      if( is.null( ssCIs[[ s ]] ) ) next
      if( is.na( ssCIs[[ s ]] ) ) next
      Btrows <- which( ssCIs[[ s ]]$par == "Bt" )
      ssBio[ s, , 1 ] <- matrix( ssCIs[[ s ]][ Btrows, "uCI" ], nrow = 1, ncol = nT )
      ssBio[ s, , 2 ] <- matrix( ssCIs[[ s ]][ Btrows, "val" ], nrow = 1, ncol = nT )
      ssBio[ s, , 3 ] <- matrix( ssCIs[[ s ]][ Btrows, "lCI" ], nrow = 1, ncol = nT )  
    }
  
  }

  # Pull Bmsy and Umsy from the blob
  Umsy <- blob$opMod$Umsy
  Bmsy <- blob$opMod$Bmsy
  MSY  <- Umsy*Bmsy

  if( scale )
  {
    # Now loop over species and scale
    for( s in 1:nS )
    {
      omBt[s,]  <- omBt[s,]/Bmsy[s]/2
      Ct[s,]    <- Ct[s,]/MSY[s]
      Ut[s,]    <- Ut[s,]/Umsy[s]
      ssBt[s,]  <- ssBt[s,]/Bmsy[s]/2
      msBt[s,]  <- msBt[s,]/Bmsy[s]/2
      msUt[s,]    <- msUt[s,]/Umsy[s]
      ssUt[s,]    <- ssUt[s,]/Umsy[s]
      if(CIs)
      {
        ssBio[s,,]<- ssBio[s,,]/Bmsy[s]/2
        msBio[s,,]<- msBio[s,,]/Bmsy[s]/2  
      }
    }  
  }
  
  
  # Recover diagnostics for the fits
  hpdSS <- blob$am$ss$hesspd[rep,]
  grdSS <- blob$am$ss$maxGrad[rep,]
  hpdMS <- blob$am$ms$hesspd[rep]
  grdMS <- blob$am$ms$maxGrad[rep]

  # Set colours for each model
  ssCol <- "steelblue"
  msCol <- "salmon"

  # Set up plot window
  par ( mfrow = c(3,nS), mar = c(1,3,2,0), oma = c(3,3,2,0.5),
        las = 1, cex.lab = labSize, cex.axis=tickSize, cex.main=labSize )
  # Plot biomass, actual and estimated, including 2 index series,
  # scaled by model estimated q
  maxBt <- 1.1*max(omBt,na.rm=T)
  if(CIs) maxBt <- max(maxBt,2)
  for ( s in 1:nS )
  {
    if( !scale ) maxBt <- 1.1*max ( omBt[s,] ,na.rm=TRUE)
    if( (!scale) & CIs ) maxBt <- max ( maxBt, msBio[s,,], ssBio[s,,], na.rm = T )
    plot    ( x = c(fYear,max(years)), y = c(0,maxBt), type = "n",
              ylim = c(0,maxBt), ylab = "", axes=FALSE, las = 1, xlab = "" ,
              main = specNames[s] )
      axis( side = 1, las = 0 )
      axis( side = 2, las = 1 )
    if (s == 1) panLegend ( x=0.2,y=1,legTxt=c("ss","ms"),
                            col=c(ssCol,msCol), lty = c(2,3), 
                            lwd = c(2,2), cex=c(1), bty="n" )
    if( scale )
    {
      ylabBio <- "(a)"
      labLas <- 1
    } else {
      ylabBio <- "Biomass (kt)"
      labLas <- 0
    }
    if( s == 1 ) mtext( text = ylabBio, side = 2, line = 2.5, las = labLas )
    if( devLabels )
    {
      if ( hpdSS[s] & !is.na(hpdSS[s]) ) panLab (x=0.9,y=0.9,txt="h",col=ssCol,cex=1.1)
      if ( hpdMS & !is.na(hpdMS) ) panLab (x=0.9,y=0.85,txt="h",col=msCol,cex=1.1)  
    }
    if( data ) points  ( x = years, y = It[s,]/ssq[s], pch = 2, cex = 0.6, col="grey70" )
    if( data ) points  ( x = years, y = It[s,]/msq[s], pch = 5, cex = 0.6, col="grey70" )
    lines   ( x = years, y = omBt[s,], col = "black", lwd = 2)
    lines   ( x = years, y = ssBt[s,], col = ssCol, lwd = 2, lty = 2 )
    lines   ( x = years, y = msBt[s,], col = msCol, lwd = 2, lty = 3 )
    if( CIs )
    {
      lines ( x = years, y = ssBio[s,,1], col = ssCol, lwd = 1, lty = 2)
      lines ( x = years, y = ssBio[s,,3], col = ssCol, lwd = 1, lty = 2)
      lines ( x = years, y = msBio[s,,1], col = msCol, lwd = 1, lty = 3)
      lines ( x = years, y = msBio[s,,3], col = msCol, lwd = 1, lty = 3)
    }
  }
  # Now plot catch
  maxCt <- max(Ct, na.rm=T)
  for ( s in 1:nS )
  {
    if( !scale ) maxCt <- max ( Ct[s,])
    plot  ( x = years, y = Ct[s,], col = "blue", lwd = 2, type = "l",
            ylim = c(0,maxCt), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1 )
      axis( side = 2, las = 1 )
    if( scale )
    {
      ylabCt <- "(b)"
      labLas <- 1
    } else {
      ylabCt <- "Catch (kt)"
      labLas <- 0
    }
    if( s == 1 ) mtext( text = ylabCt, side = 2, line = 2.5, las = labLas )

  }
  # Now F
  maxUt <- max(Ut, na.rm=T)
  for ( s in 1:nS )
  {
    if(!scale) maxUt <- max ( Ut[s,])
    plot  ( x = years, y = Ut[s,], col = "black", lwd = 2, type = "l",
            ylim = c(0,maxUt), ylab = "", axes=FALSE, las = 1, xlab = "" )
      axis( side = 1 )
      axis( side = 2, las = 1 )
    if( scale )
    {
      ylabUt <- "(c)"
      labLas <- 1
    } else 
    {
      ylabUt <- "Exploitation Rate"
      labLas <- 0
    }
    if( s == 1 ) mtext( text = ylabUt, side = 2, line = 2.5, las = labLas )
    lines (x = years, y = ssUt[s,], col= ssCol, lwd = 2, lty = 2 )
    lines (x = years, y = msUt[s,], col= msCol, lwd = 2, lty = 3 )
  }
  mtext ( text = "Year", outer = TRUE, side = 1, padj = 1.5)
  if( devLabels )
    mtext ( text = c(stamp,repCount),side=1, outer = TRUE, 
            at = c(0.9,0.1),padj=2,col="grey50",cex=0.8 )
}


mpNamesPerfPlot <- c( noJointPriors ="None", 
                      YearEffOnly = expression(epsilon[t]), 
                      qPriorOnly = expression(q), 
                      UmsyPriorOnly = expression(r), 
                      qUpriors = expression( q / r ), 
                      qYEpriors = expression(q / epsilon[t]),
                      UmsyYEpriors = expression(r / epsilon[t]), 
                      allJointPriors = expression(q / r / epsilon[t]) )



dumpStockPerf <- function(  simPath = file.path("./project/pubBase_2018-09-10/project"),
                            prefix = "pubBase",
                            MPs = c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors" ),
                            MPlabels = c( noJointPriors = "None", 
                                          qPriorOnly = expression(q), 
                                          UmsyPriorOnly = expression(r), 
                                          qUpriors = expression(q/r) ),
                            simNumTable = NULL )
{
  # Need to read in all the sims in the folder, and then
  # tabulate so we can plot out the sets of MPs
  # List directories
  simList <- list.files(simPath, full.names = TRUE)
  simList <- simList[grepl(pattern = "sim", x = simList)]
  if(simPath != "./project")
  {
    file.copy( from = simList, to = "./project/", recursive = TRUE)
  }

  plotPath <- file.path("./project/figs",prefix)
  if(!dir.exists(plotPath))
    dir.create(plotPath)

  plotPath <- file.path(plotPath,"stockPerf")
  if(!dir.exists(plotPath))
    dir.create(plotPath)
  if(is.null(simNumTable))
    simNumTable <- makeSimNumTable()

  scenarios <- unique(simNumTable$scenario)
  if(is.null(MPs))
    MPs     <- unique(simNumTable$mp)

  for( sIdx in 1:length(scenarios) )
  {
    scenLabel <- scenarios[sIdx]
    plotFile <- file.path(plotPath, paste(scenLabel,".pdf",sep = "") )

    subTable <- simNumTable %>%
                filter( scenario == scenLabel )

    mpOrder <- numeric(length(MPs))
    for( mIdx in 1:length(MPs))
      mpOrder[mIdx] <- subTable[which(subTable[,"mp"] == MPs[mIdx] ),"simNum"]

    pdf(plotFile, width = 5, height = 8 )
    plotStockPerfMultiSim(  pars = c("Umsy","BnT","dep","Bmsy","q_os"), 
                            sims = mpOrder, spec = 1, nSurv = 2,
                            devLabels = TRUE,
                            title = TRUE, plotMARE = FALSE,
                            mpNames = MPlabels,
                            labSize = 1 )
    dev.off()
  }
  cat("Completed plotting stock performance plot for ", prefix, "\n", sep = "" )
}

# plotSimPerf()
# A function to plot simulation-estimation performance for a whole
# collection of replicates. Resulting plot is a dot and bar graph
# showing relative error distributions for nominated parameters
# inputs:   sim = numeric indicator of simulation in project folder
#           folder = optional character indicating sim folder name
#           pars = nominated estimated leading and derived parameters 
# outputs:  NULL
plotStockPerfMultiSim <- function ( pars = c("Umsy","BnT","Bmsy","dep","q_os"), 
                                    sims = 1, spec = 1, nSurv = 2,
                                    devLabels = TRUE,
                                    title = TRUE, plotMARE = FALSE,
                                    mpNames = mpNamesPerfPlot,
                                    labSize = 1 )
{
  # Create a wrapper function for generating quantiles
  quantWrap <- function ( entry = 1, x = ssRE, ... )
  {
    if (is.null(dim(x[[entry]]))) x[[entry]] <- matrix(x[[entry]],ncol=1)
    xDim <- dim(x[[entry]])
    margins <- 2:length(xDim)
    quants <- apply ( X = x[[entry]], FUN = quantile, MARGIN = margins, ... )
    if(class(quants) == "array")
      return( aperm(quants) )
    else return( t(quants) )
  }  

  blobList <- vector(mode = "list", length=length(sims))
  # Loop over sims and load their MLEs
  for( sIdx in 1:length(sims) )
  { 
    sim <- sims[sIdx]
    .loadSim(sim)

    blobList[[sIdx]]$mp       <- blob$ctrl$mp
    blobList[[sIdx]]$scenario <- blob$ctrl$scenario
    blobList[[sIdx]]$stamp    <- paste( blob$ctrl$scenario, blob$ctrl$mp, sep = ":")
    blobList[[sIdx]]$nO       <- blob$opMod$nSurv

    blobList[[sIdx]]$ssQuant  <- lapply ( X = seq_along(blob$am$ss$err.mle[pars]), 
                                          FUN = quantWrap, x=blob$am$ss$err.mle[pars], 
                                          probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    names( blobList[[sIdx]]$ssQuant ) <- names( blob$am$ss$err.mle[pars] )

    
    blobList[[sIdx]]$msQuant  <- lapply ( X = seq_along(blob$am$ms$err.mle[pars]), 
                                          FUN = quantWrap, x=blob$am$ms$err.mle[pars], 
                                          probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

    names( blobList[[sIdx]]$msQuant ) <- names( blob$am$ms$err.mle[pars] )


  }

  # Multisurvey parameters
  if( "q_os" %in% pars )
  {
    q_os <- paste( "q_", 1:nSurv, sep = "" )
    pars <- c( pars[pars != "q_os"], q_os ) 
  }

  if( "tau2_o" %in% pars )
  {
    tau2_o <- paste( "tau2_", 1:nSurv, sep = "" )
    pars <- c( pars[pars != "tau2_o"], tau2_o ) 
  }

  # Create an array of quantile values 
  # (dims = species,percentiles,parameters,models)
  quantiles <- array ( NA, dim = c(length(sims),3,length(pars),2))
  dimnames (quantiles) <- list (  character(length(sims)),
                                  c("2.5%","50%","97.5%"), 
                                  pars,
                                  c("ss","ms") )


  MPs   <- c("Single Stock")


  for( sIdx in 1:length(sims) )
  {
    # populate table
    for (p in 1:length(pars))
    {
      parName <- pars[p]
      # Change behaviour if survey specific stuff is involved
      # SDNJ: Check backwards compatibility for master branch
      if( grepl( "q_", parName ) | grepl("tau2_", parName ) )
      {
        if( grepl( "q_", parName ) ) parName <- "q_os"
        if( grepl("tau2_", parName ) ) parName <- "tau2_o"  
        surveyNo <- as.numeric(str_split(pars[p],"_")[[1]][2])
        # Now fill in the survey specific deets
        # REs
        quantiles[sIdx,,pars[p],1] <- blobList[[sIdx]]$ssQuant[ parName ][[1]][spec,surveyNo,]
        quantiles[sIdx,,pars[p],2] <- blobList[[sIdx]]$msQuant[ parName ][[1]][spec,surveyNo,]
        # # MARE
        # MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s,surveyNo]
        # MARE[s,pars[p],2] <- msAQuant[[ parName ]][s,surveyNo]
      } else {
        quantiles[sIdx,,pars[p],1] <- blobList[[sIdx]]$ssQuant[ parName ][[1]][spec,]
        quantiles[sIdx,,pars[p],2] <- blobList[[sIdx]]$msQuant[ parName ][[1]][spec,]
        # # MARE
        # MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s]
        # MARE[s,pars[p],2] <- msAQuant[[ parName ]][s]
      }
    }

    MPs <- c(MPs,blobList[[sIdx]]$mp)
  }

  if( !is.null(mpNames) ) 
  {
    # Replace each term in MPs with the corresponding mpNames expression,
    # however, we can't run comparisons in MPs after replacing one, as
    # comparisons of expressions are not allowed
    mpNum <- numeric(length(mpNames))
    for( nomIdx in 1:length(mpNames) )
    { 
      mpNum[nomIdx] <- which(MPs == names(mpNames)[nomIdx])
    }
    for( num in 1:length(mpNum))
      MPs[mpNum[num]] <- mpNames[num]
  } 

  # Pull out RE dists
  med <- quantiles[ ,"50%",, ]
  q975 <- quantiles[ ,"97.5%",, ]
  q025 <- quantiles[ ,"2.5%",, ]



  # Set up plotting environment
  par( mfrow = c(length(pars),1), mar = c(0,4,0,2), oma = c(4,7,1,8) )
  xlim  <- c(-1,1.5)  
  # Loop over pars
  for( pIdx in 1:length(pars) )
  {
    par <- pars[pIdx]

    plot( x = xlim, y = 3*c(0,length(MPs)+1), type = "n", axes = F, ylab = "" )
      if( pIdx == length(pars) ) axis( side = 1, cex.axis = labSize )
      abline( v = 0, lty = 3, lwd = 1 )
      axis( side = 2 , at = 3*1:length(MPs), labels = MPs, las = 1, cex.axis = 1.2 )
      segments( x0 = q025[1,par,"ss"], y0 = 3, x1 = q975[1,par,"ss"], y1 = 3, lty=2, col="grey60", lwd = 4 )
      points( x = med[1,par,"ss"], y = 3, pch = 17, cex = 2 )
      segments( x0 = q025[,par,"ms"], y0 = 3*(2:length(MPs)), x1 = q975[,par,"ms"], y1 = 3*(2:length(MPs)), lty=1, col="grey60", lwd = 4 )
      points( x = med[,par,"ms"], y = 3*(2:length(MPs)), pch = 16, cex = 2 )
      if( par == "BnT" )  parLab <- expression(B[T])
      if( par == "dep" )  parLab <- expression(B[T]/B[0])
      if( par == "Bmsy" ) parLab <- expression(B[MSY])
      if( par == "Umsy" ) parLab <- expression(U[MSY])
      if( grepl(x = par, pattern = "q_") )
      { 
        oIdx <- str_split(par,"_")
        oIdx <- unlist(oIdx)[2]
        if( oIdx == 1 )
          parLab <- expression(q[1,1])
        if( oIdx == 2 )
          parLab <- expression(q[2,1])
      }
      if( is.null(parLab) ) parLab <- par
      mtext( side = 4, text = parLab, las = 1, line = 2, cex = 1.5 )
      
  }

  if( title )
  {
    mtext ( text = "Relative Error", side = 1, outer = TRUE, line = 3, cex = 1.5)
    mtext ( text = "AM configuration", side = 2, outer = TRUE, line = 5, cex = 1.5)
    # mtext ( text = "Parameter", side = 4, outer = TRUE, line =  cex = 1.5, las = 0)
  }

}


# plotSimPerf()
# A function to plot simulation-estimation performance for a whole
# collection of replicates. Resulting plot is a dot and bar graph
# showing relative error distributions for nominated parameters
# inputs:   sim = numeric indicator of simulation in project folder
#           folder = optional character indicating sim folder name
#           pars = nominated estimated leading and derived parameters 
# outputs:  NULL
plotSimPerf <- function ( pars = c("Bmsy","Umsy","q_os","dep","BnT","tau2_o","totRE"), sim=1, 
                          est="MLE", devLabels = TRUE,
                          title = TRUE, plotMARE = FALSE )
{
  # Blob should be loaded in global environment automatically,
  # if not, load first one by default (or whatever is nominated)
  if (!exists(x="blob",where=1)) .loadSim(sim)

  # Create a stamp from scenario and mp name
  scenario  <- blob$ctrl$scenario
  mp        <- blob$ctrl$mp
  stamp     <- paste(scenario,":",mp,sep="")

  # Recover blob elements for plotting
  nS        <- blob$opMod$nS
  specNames <- blob$ctrl$speciesNames[1:nS]
  nO        <- ifelse(!is.null(blob$opMod$nSurv), blob$opMod$nSurv, 1)

  # Recover relative error distributions
  if (est == "MLE")
  {
    ssRE <- blob$am$ss$err.mle[pars]
    msRE <- blob$am$ms$err.mle[pars]
  }

  
  # Create a wrapper function for generating quantiles
  quantWrap <- function ( entry = 1, x = ssRE, ... )
  {
    if (is.null(dim(x[[entry]]))) x[[entry]] <- matrix(x[[entry]],ncol=1)
    xDim <- dim(x[[entry]])
    margins <- 2:length(xDim)
    quants <- apply ( X = x[[entry]], FUN = quantile, MARGIN = margins, ... )
    if(class(quants) == "array")
      return( aperm(quants) )
    else return( t(quants) )
  }

  # browser()
  ssARE <- lapply(X = ssRE, FUN = abs)
  msARE <- lapply(X = msRE, FUN = abs)

  # generate quantiles
  ssQuant <- lapply ( X = seq_along(ssRE), FUN = quantWrap, x=ssRE, 
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  msQuant <- lapply ( X = seq_along(msRE), FUN = quantWrap, x=msRE, 
                      probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

  ssAQuant <- lapply (  X = seq_along(ssARE), FUN = quantWrap, x=ssARE, 
                        probs = c(0.5), na.rm = TRUE)
  msAQuant <- lapply (  X = seq_along(msARE), FUN = quantWrap, x=msARE, 
                        probs = c(0.5), na.rm = TRUE)

  # get names of parameters
  ssPars <- names ( ssRE )
  msPars <- names ( msRE )
  names ( ssQuant) <- ssPars
  names ( msQuant) <- msPars

  names ( ssAQuant) <- ssPars
  names ( msAQuant) <- msPars

  # Multisurvey parameters
  if( "q_os" %in% pars )
  {
    q_os <- paste( "q_", 1:nO, sep = "" )
    pars <- c( pars[pars != "q_os"], q_os ) 
  }

  if( "tau2_o" %in% pars )
  {
    tau2_o <- paste( "tau2_", 1:nO, sep = "" )
    pars <- c( pars[pars != "tau2_o"], tau2_o ) 
  }

  # Create an array of quantile values 
  # (dims = species,percentiles,parameters,models)
  quantiles <- array ( NA, dim = c(nS,3,length(pars),2))
  dimnames (quantiles) <- list (  specNames,
                                  c("2.5%","50%","97.5%"), 
                                  pars,
                                  c("ss","ms") )

  MARE <- array ( NA, dim = c(nS,length(pars),2),
                  dimnames = list ( specNames, pars,
                                    c("ss","ms") ) )


  # populate table
  for ( s in 1:nS )
  {
    for (p in 1:length(pars))
    {
      parName <- pars[p]
      # Change behaviour if survey specific stuff is involved
      # SDNJ: Check backwards compatibility for master branch
      if( grepl( "q_", parName ) | grepl("tau2_", parName ) )
      {
        if( grepl( "q_", parName ) ) parName <- "q_os"
        if( grepl("tau2_", parName ) ) parName <- "tau2_o"  
        surveyNo <- as.numeric(str_split(pars[p],"_")[[1]][2])
        # Now fill in the survey specific deets
        # REs
        quantiles[s,,pars[p],1] <- ssQuant[ parName ][[1]][s,surveyNo,]
        quantiles[s,,pars[p],2] <- msQuant[ parName ][[1]][s,surveyNo,]
        # MARE
        MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s,surveyNo]
        MARE[s,pars[p],2] <- msAQuant[[ parName ]][s,surveyNo]
      } else {
        quantiles[s,,pars[p],1] <- ssQuant[ parName ][[1]][s,]
        quantiles[s,,pars[p],2] <- msQuant[ parName ][[1]][s,]
        # MARE
        MARE[s,pars[p],1] <- ssAQuant[[ parName ]][s]
        MARE[s,pars[p],2] <- msAQuant[[ parName ]][s]
      }
    }
  }
  # Set plotting window
  par (mfrow = c(1,nS), mar = c(1,0,2,1), oma = c(3,0,2,0) )

  for ( s in 1:nS )
  {
    med <- quantiles[s,"50%",,]
    mare <- MARE[s,,]
    q975 <- quantiles[s,"97.5%",,]
    q025 <- quantiles[s,"2.5%",,]

    # Plot main dotchart
    plotTitle <- paste(specNames[s])
    FitCounter <- paste( "nReps = ", sum(blob$goodReps[,s]), sep = "" )
    if ( s == 1) dotchart ( x = t(med), xlim = c(-1.5,1.5), main = "",
                            pch = 16, cex = 1.2 )
    else dotchart ( x = t(med), xlim = c(-1.5,1.5), main = "",
                    pch = 16, cex = 1.2)
      mtext( side = 3, text = plotTitle, cex = 1.5, las = 0, line = 2 )
      mtext( side = 3, text = FitCounter, cex = 1, las = 0, line = 1 )
      # axis ( side = 1, at = seq(-1.5,1.5,by=0.5) )
    
    # Now add segments
    for ( p in length(pars):1)
    {
      parIdx <- length(pars) - p + 1
      plotY <- 2 * p + 2 * (p-1)
      segments( q025[pars[parIdx],"ss"],plotY-1,q975[pars[parIdx],"ss"],plotY-1,lty=1,col="grey60", lwd=3)  
      segments( q025[pars[parIdx],"ms"],plotY,q975[pars[parIdx],"ms"],plotY,lty=1,col="grey60", lwd=3)  
      if(plotMARE) 
      {
        points(x = mare[parIdx,"ss"],y = plotY-1, pch = 2, cex = 1.2 )
        points(x = mare[parIdx,"ms"],y = plotY, pch = 2, cex = 1.2 )
      }
    }

    abline ( v = 0, lty = 3, lwd = 0.5 )

  } 
  if( title )
  {
    mtext ( text = "Relative Error", side = 1, outer = TRUE, line = 2, cex = 1.5)
    mtext ( text = "Parameter", side = 2, outer = TRUE, line = 1, cex = 1.5)
  }
  
  if( devLabels )
    mtext ( text = c(stamp),side=1, outer = TRUE, 
            at = c(0.75),padj=1,col="grey50" )

}

