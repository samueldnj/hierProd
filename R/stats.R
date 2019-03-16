# --------------------------------------------------------------------------
# stats.R
# 
# Statistics functions script for the coastwide (single stock) 
# simulation-estimation experiments comparing single species and 
# multispecies assessment models.
# 
# Author: Samuel D N Johnson
# Date: 7 September, 2016
#
# --------------------------------------------------------------------------


fitTablePub <- function(  fitTab = "DoverAssessPubFinal.csv",
                          pars = c(  "Umsy", "BnT", "Bmsy", "DnT", "U_Umsy", "AICc" ),
                          models = c("qPriorOnly","UmsyPriorOnly","qUpriors"),
                          saveFile = "DoverAssessPub_AICTable.csv" )
{
  # First, read in the table
  tabPath <- file.path(getwd(),"project/statistics",fitTab)
  table <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE )  

  # Reduce to only wanted MPs
  table <- table %>% filter( mp %in% models )

  # Now split by scenario
  scenarios <- unique(table$scenario)

  # Count stocks
  stocks <- unique( table$stock )
  nStocks <- length( stocks )

  baseFitMtx <- matrix( NA, ncol = 5, nrow = length(models) + 1)
  colnames(baseFitMtx) <- c("Parameter", "Model", stocks )
  baseFitTab <- as.data.frame( baseFitMtx )

  # Column names for models
  ssPars  <- paste( "ss", pars, sep = "" )
  ssParSE <- paste( ssPars, "se", sep = "")

  msPars  <- paste( "ms", pars, sep = "" )
  msParSE <- paste( msPars, "se", sep = "")

  # Create a list to hold a fit table for each scenario
  scenList <- vector(mode="list", length = 2)
  names(scenList) <- scenarios

  for( scenIdx in 1:length( scenarios ) )
  {
    scenTable <-  table %>%
                  filter( scenario == scenarios[scenIdx] )

    parList <- vector( mode = "list", length = length(pars) + 1 )
    names(parList) <- pars

    for( pIdx in 1:length(pars) )
    {
      # Now start populating a parTab
      parTab <- baseFitTab
      parTab$Parameter <- pars[pIdx]
      parTab$Model <- c( "Single-stock", models )

      # Loop over models
      for( mIdx in 1:length(models) )
      {
        modTab <- filter( scenTable, mp == models[mIdx] )
        if( mIdx == 1 )
        {
          # Populate single stock model
          for( sIdx in 1:nStocks )
          {
            # Now pull stock spec estimates
            stockTab <- filter( modTab, stock == stocks[sIdx] )
            parEst  <- round(stockTab[,ssPars[pIdx]], digits = 2)
            if( pars[pIdx] != "AICc" ) 
            {
              parSE   <- round(stockTab[,ssParSE[pIdx]], digits =2)
              parEst  <- paste( parEst, " (",parSE,")", sep = "")
            }
            parTab[1,stocks[sIdx]] <- parEst        
          }
        }
        # Now pull MS model estimates
        # Populate single stock model
        for( sIdx in 1:nStocks )
        {
          # Now pull stock spec estimates
          stockTab <- filter( modTab, stock == stocks[sIdx] )
          parEst  <- round(stockTab[,msPars[pIdx]], digits = 2)
          if( pars[pIdx] != "AICc" ) 
          {
            parSE   <- round(stockTab[,msParSE[pIdx]], digits = 2)
            parEst  <- paste( parEst, " (",parSE,")", sep = "")
          }
          # browser()
          parTab[mIdx + 1,stocks[sIdx]] <- parEst        
        }
      }
      parList[[pIdx]] <- parTab  
    }
    scenList[[scenIdx]] <- do.call("rbind",parList)
  }
  outTable <- do.call("cbind",scenList)
  savePath <- file.path(".","project","Statistics",saveFile)
  write.csv( outTable, savePath )
}


checkCorr <- function(  tabName = "allSame_randProcCorr_MRE.csv",
                        cols = c("ssq_1","ssq_2","ssBmsy","ssBnT","msq_1","msq_2","msBmsy","msBnT"),
                        spec = NULL, MP = NULL, nSp = NULL )
{
  # first read in the table
  tabPath <- file.path(getwd(),"project/statistics",tabName)
  table <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE ) 

  # restrict to correct number of species
  # Restrict to subtables if requested
  if(!is.null(MP))    table  <-   table %>% filter( mp %in% MP ) 
  if(!is.null(spec))  table  <-   table %>% filter( species %in% spec ) 
  if(!is.null(nSp))   table  <-   table %>% filter( nS %in% nSp )

  # Restrict to columns
  table <- table[,cols]

  # Now compute the correlation
  corrMtx <- cor( table, use = "pairwise.complete.obs" )

  # Plot a pairs plot for shits and giggles
  pairs(table)

  corrMtx
}


.makeDeltaCols <- function( tabName = "pubBase_MARE" )
{
  # read in table
  tabFile     <- paste( tabName, ".csv", sep  = "" )
  tabPath     <- file.path(getwd(),"project","Statistics",tabFile)
  tab         <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE )

  # Takes the mean of numeric columns, or the first value
  # of other columns (will pick stock1 in most cases)
  meanTop <- function(x)
  {
    if(is.numeric(x)) out <- mean(x)
    else out <- x[1]

    out
  }

  # summarise complex table to get mean MARE values
  # and stock1's other factor levels
  complexTab <- tab %>%
                group_by(scenario, mp) %>%
                summarise_all( .funs=meanTop) %>%
                ungroup() %>%
                mutate( DeltaBmsy = (ssBmsy / msBmsy) - 1,
                        DeltaBnT  = (ssBnT / msBnT) - 1,
                        DeltaUmsy = (ssUmsy / msUmsy) - 1,
                        Deltaq_1  = (ssq_1 / msq_1) - 1,
                        Deltaq_2  = (ssq_2 / msq_2) - 1,
                        DeltaDep  = (ssDep / msDep)  - 1)


  tab <-  tab %>%
          mutate( DeltaBmsy = (ssBmsy / msBmsy) - 1,
                  DeltaBnT  = (ssBnT / msBnT) - 1,
                  DeltaUmsy = (ssUmsy / msUmsy) - 1,
                  Deltaq_1  = (ssq_1 / msq_1) - 1,
                  Deltaq_2  = (ssq_2 / msq_2) - 1,
                  DeltaDep  = (ssDep / msDep)  - 1)

  # Write complex delta table
  cplxTabName   <- paste(tabName,"cplx",sep = "_")
  cplxFileName  <- paste(cplxTabName,".csv",sep = "")
  cplxSavePath  <- file.path(getwd(),"project","Statistics",cplxFileName)
  write.csv( x = complexTab, file = cplxSavePath )

  # write stock delta table
  tabFile       <- paste( tabName, "_Delta.csv", sep = "" )
  tabPath       <- file.path(getwd(),"project","Statistics",tabFile)
  write.csv( x = tab, file = tabPath )

  cat("Delta columns appended to ", tabPath, ",\n",
      "Complex tab file ", cplxTabName, " created at ", cplxSavePath, "\n", sep = "" )
}

# Fits a GLM to the supplied stat table, using the named columns as explanatory 
# and response variables. 
metaModels <- function( tabName = "allSame_infoScenarios_MARE_cplx",
                        multiResp = c("BnT", "Umsy","q_1","q_2","Dep","Bmsy"),
                        singleResp = c("DeltaBnT","DeltaUmsy","Deltaq_1","Deltaq_2","DeltaDep","DeltaBmsy"),
                        spec = c("Stock1"),
                        expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
                        sig = .1, intercept = TRUE,
                        scaled = TRUE, saveOut = TRUE, interactions = FALSE, dropTest = FALSE,
                        tabSavePath = "./project/Statistics/" )
{
  # Create an rDataFile output name
  rDataName <- tabName
  ssResp <- paste("ss", multiResp, sep = "" )
  msResp <- paste("ms", multiResp, sep = "" )
  
  # Now some lists to hold meta-models
  fitList <- vector( mode = "list", length = 3 )
  names(fitList) <- c("ss","ms","group")
  # responses that have single species and multispecies versions
  respList <- vector(mode = "list", length = length(multiResp))
  names( respList ) <- multiResp
  # Group parameters (MS prior mean/var, convergence stats )
  groupList <- vector(mode = "list", length = length(singleResp))
  if(!is.null(singleResp)) names(groupList) <- singleResp
  fitList$ss <- respList
  fitList$ms <- respList
  fitList$group <- groupList
  
  # Now loop over response variables and fit models
  for( rIdx in 1:length(multiResp) )
  {
    resp <- multiResp[rIdx]
    ssResp <- paste( "ss", resp, sep = "" )
    msResp <- paste( "ms", resp, sep = "" )
    fitList$ss[[rIdx]] <- .DASEmodelSelection(  tabName = tabName,
                                                resp = ssResp,
                                                expVars = expVars,
                                                sig = sig,
                                                intercept = intercept,
                                                scaled = scaled,
                                                abs = FALSE,
                                                spec = spec,
                                                interactions = interactions,
                                                dropTest = dropTest )
    fitList$ms[[rIdx]] <- .DASEmodelSelection(  tabName = tabName,
                                                resp = msResp,
                                                expVars = expVars,
                                                sig = sig,
                                                intercept = intercept,
                                                scaled = scaled,
                                                abs = FALSE,
                                                spec = spec,
                                                interactions = interactions,
                                                dropTest = dropTest )
  }
  if(!is.null(singleResp))
  {
    for( rIdx in 1:length(singleResp) )
    {
      resp <- singleResp[rIdx]
      fitList$group[[rIdx]] <- .DASEmodelSelection( tabName = tabName,
                                                    resp = resp,
                                                    expVars = expVars,
                                                    sig = sig,
                                                    intercept = intercept,
                                                    scaled = scaled,
                                                    abs = FALSE,
                                                    interactions = interactions,
                                                    dropTest = dropTest,
                                                    spec = spec )
    }  
  } 
  # Create directories if needed
  if(!dir.exists(tabSavePath))
    dir.create(tabSavePath)
  
  savePath <- file.path( tabSavePath, rDataName ) 
  if(!dir.exists(savePath))
    dir.create( savePath )
  if(saveOut) save( fitList, 
                    file = file.path( savePath,
                                      paste(rDataName,".RData", sep = "" ) ) )

  # Now output the aic tables
  # First, the ss/ms pars
  for( par in multiResp)
  {
    # recover ss and ms AIC tables
    AICss <- fitList$ss[[par]]$AICrank
    AICms <- fitList$ms[[par]]$AICrank
    # Main effects first
    AICssFile <- paste(rDataName,par,"AICss.csv",sep = "")
    AICmsFile <- paste(rDataName,par,"AICms.csv",sep = "")
    
    write.csv( x = AICss, file = file.path(savePath,AICssFile ) )
    write.csv( x = AICms, file = file.path(savePath,AICmsFile ) )
  }
  # then the group pars
  if( !is.null(singleResp))
  {
    for( pIdx in 1:length(singleResp) )
    {
      par <- singleResp[pIdx]
      # recover ss and ms AIC tables
      AIC <- fitList$group[[pIdx]]$AICrank
      # Main effects first
      AICfile <- paste(rDataName,par,"AIC.csv",sep = "")
      
      write.csv( x = AIC, file = file.path(savePath,AICfile ) )
    }
  }

  # Now organise the top ranked models into a single table
  pull1AIC <- function(par, list, prefix = "")
  {
    out <- list[[par]]$AICrank[1,]
    out$par <- paste(prefix,par,sep = "")

    out
  }

  # Now organise the top ranked models into a single table
  # pulling the SEs this time
  pull1AIC.se <- function(par, list, prefix = "")
  {
    out <- list[[par]]$AICrank.se[1,]
    out$par <- paste(prefix,par,sep = "")

    out
  }

  # Now organise the top ranked models into a single table
  # pulling the SEs this time
  pull1AIC.eff <- function(par, list, prefix = "")
  {
    out <- list[[par]]$AICrank.eff[1,]
    out$par <- paste(prefix,par,sep = "")

    out
  }

  ssTop <- lapply(X = multiResp, FUN = pull1AIC, list = fitList$ss, prefix = "ss" )
  ssTop <- do.call(what = "rbind", args = ssTop )

  ssTop.se <- lapply(X = multiResp, FUN = pull1AIC.se, list = fitList$ss, prefix = "ss" )
  ssTop.se <- do.call(what = "rbind", args = ssTop.se )

  ssTop.eff <- lapply(X = multiResp, FUN = pull1AIC.eff, list = fitList$ss, prefix = "ss" )
  ssTop.eff <- do.call(what = "rbind", args = ssTop.eff )

  msTop <- lapply(X = multiResp, FUN = pull1AIC, list = fitList$ms, prefix = "ms" )
  msTop <- do.call(what = "rbind", args = msTop )

  msTop.se <- lapply(X = multiResp, FUN = pull1AIC.se, list = fitList$ms, prefix = "ms" )
  msTop.se <- do.call(what = "rbind", args = msTop.se )

  msTop.eff <- lapply(X = multiResp, FUN = pull1AIC.eff, list = fitList$ms, prefix = "ms" )
  msTop.eff <- do.call(what = "rbind", args = msTop.eff )

  topRank       <- rbind( ssTop, msTop )
  topRank.se    <- rbind( ssTop.se, msTop.se )
  topRank.eff   <- rbind( ssTop.eff, msTop.eff )

  if( !is.null(singleResp) )
  {
    groupTop <- lapply(X = singleResp, FUN = pull1AIC, list = fitList$group, prefix = "" )
    groupTop <- do.call(what = "rbind", args = groupTop )

    groupTop.se <- lapply(X = singleResp, FUN = pull1AIC.se, list = fitList$group, prefix = "" )
    groupTop.se <- do.call(what = "rbind", args = groupTop.se )

    groupTop.eff <- lapply(X = singleResp, FUN = pull1AIC.eff, list = fitList$group, prefix = "" )
    groupTop.eff <- do.call(what = "rbind", args = groupTop.eff )

    topRank     <- rbind( topRank, groupTop )
    topRank.se  <- rbind( topRank.se, groupTop.se )
    topRank.eff <- rbind( topRank.eff, groupTop.eff )
  }

  # Add an se tag to each parameter name
  topRank.se$par <- paste(topRank.se$par,"se",sep = "")

  # Now interleave them
  idx <- rbind( 1:nrow(topRank.eff), (nrow(topRank.se)+1):(2*nrow(topRank.se)) )
  topRank.interleave <- rbind(topRank.eff,topRank.se)[idx,]

  # Save
  # Interleaved SEs version
  topRankFile <- paste( rDataName, "topRank_interleave.csv", sep = "" )
  write.csv( x = topRank.interleave, file = file.path(tabSavePath,topRankFile) )

  # Inline SEs version
  topRankFile <- paste( rDataName, "topRank.csv", sep = "" )
  write.csv( x = topRank, file = file.path(tabSavePath,topRankFile) )

  fitList
  
}

# .DASEmodelSelection()
# Uses the following procedure to select candidate models for
# analysing simulation model output.
# Basic idea, uses strong heredity of main effects: 
#   1. Fit all main effects, check for significance
#   2. Choose main effects model with lowest AIC value from among significant models, 
#       repeat 1 with all possible interactions and second order effects from that model
#   3. Choose among first order, fo+interactions, fo+second order and
#       fo+so+interactions models using AICc
# Models from out$ms$sel and out$ss$sel can be used in plotOMSpecResponse 
# to plot the glm fits next to the observed values.
.DASEmodelSelection <- function(  tabName = "obsErrNew_MRE",
                                  resp = "ssq",
                                  spec = "Stock1",
                                  expVars = c("nHi","CVlo","CVhi","nS","mp"),
                                  sig = 0.1,
                                  intercept = FALSE,
                                  abs = FALSE,
                                  scaled = TRUE,
                                  interactions = FALSE,
                                  dropTest = TRUE )
{
  # First, fit main effects
  mainEff <- .DASEexperiment( tabName = tabName,
                              resp = resp, spec = spec,
                              scaled = scaled,
                              baseRHS = NULL,
                              abs = abs,
                              newExp = expVars,
                              intercept = intercept )

  if( interactions )
  {
    # Then, look at the ANOVA table for the 
    # factors, since some may be qualitative and have 
    # insignificant levels (small effect sizes)
    mainDrop <- lapply( X = mainEff$GLM, FUN = drop1, test = "LRT" )
    mainDrop <- mainDrop[[length(mainDrop)]]

    insig <- which( mainDrop[-1,5] >= sig )
    insig <- rownames(mainDrop[-1,])[insig]
    
    dropExp <- names(insig)

    signifExp <- expVars[ !(expVars %in% dropExp) ]

    # Function to make interactions out of main effects
    makeIntEff <- function ( main = expVars )
    {
      interactions <- character(length = 0)
      for( i in 1:(length(main)-1) )
      {
        for( j in (i+1):length(main) )
        {
          interactions <- c( interactions, paste( main[i], main[j], sep = ":" ) )
        }
      }
      interactions
    }
    # Now paste main effects together
    signifRHS <- paste(signifExp, collapse = " + ")

    # make interactions and second order effects
    intEff <- makeIntEff ( signifExp )
    # remove qualitative factors for second order effects
    if( "mp" %in% signifExp ) signifExpSO <- signifExp[ signifExp != "mp" ]
    else signifExpSO <- signifExp
    soEff  <- paste( "I(", signifExpSO, "^2)", sep = "" )

    # And run the experiment again
    fullModels <- .DASEexperiment(  tabName = tabName, resp = resp,
                                    spec = spec, scaled = FALSE,
                                    baseRHS = signifRHS,
                                    newExp = c(intEff, soEff) )

    allModels <- c(mainEff$GLM[-length(mainEff$GLM)],fullModels$GLM)
  } else allModels <- mainEff$GLM
  if( grepl(pattern = "Delta", x = resp ) ) scale = 1 else scale = 100
  # Now rank AICc values in each model
  aicRank     <- AICrank( allModels, sig = sig, scale = scale, drop = dropTest )

  aicRank
}


# Automate ranking by AIC values, append effect sizes
# to the output.
AICrank <- function ( modelList, sig, scale, drop = TRUE )
{
  # Count models
  nModels   <- length( modelList )
  summList  <- lapply( X = modelList, FUN = summary )
  if( drop ) 
    dropList  <- lapply( X = modelList, FUN = drop1, test = "LRT" )

  # Pull out model with all coefficients for effect size recording
  fullModel <- modelList[[ nModels ]]
  fullSumm  <- summList[[ nModels ]]
  coefFull  <- coef( fullModel )
  coefFSumm <- coef( fullSumm )
  nCoef     <- nrow( coefFSumm )

  # Create a data frame to hold the info
  rank.df           <- matrix(NA, ncol = 4, nrow = nModels )
  colnames(rank.df) <- c( "Number", "AICc", "deltaAICc", "maxPr" )
  rank.df           <- as.data.frame(rank.df)

  # Make effect size df with nModels rows, and a matching df for standard errors
  effect.df             <- matrix ( NA, ncol = nrow( coefFSumm ) + 1, nrow = nModels )
  colnames(effect.df)   <- c("Number", dimnames(coefFSumm)[[1]])
  effect.df             <- as.data.frame(effect.df)
  se.df                 <- effect.df
  combined.df           <- effect.df
  
  # Add an idenitifying number (useful for looking through model objects)
  rank.df[,"Number"] <- 1:nModels
  effect.df[,"Number"] <- 1:nModels
  se.df[,"Number"] <- 1:nModels
  combined.df[,"Number"] <- 1:nModels
  for(k in 1:nModels)
  {
    nPar <- nrow(summList[[k]]$coefficients)
    nObs <- length(summList[[k]]$deviance.resid)
    rank.df[k,"AICc"] <- summList[[k]]$aic + nPar*(nPar + 1) / (nObs - nPar - 1)
    if(drop)
      rank.df[k,"maxPr"] <- max(dropList[[k]][-1,5])
    coefSumm <- coef(summList[[k]])
    coefSumm[,1:2] <- round(coefSumm[,1:2] * scale,2)
    effects <- coefSumm[,1]
    effects.se <- coefSumm[,2]
    for( eIdx in 1:nrow( coefSumm ) )
    {
      effName <- dimnames(coefSumm)[[1]][eIdx]
      effect.df[k,effName] <- effects[eIdx]
      se.df[k,effName] <- effects.se[eIdx]
      combined.df[k,effName] <- paste( format(effects[eIdx],nsmall = 2), " (", format(effects.se[eIdx],nsmall = 2),")", sep = "")
    }
  }

  # Remove insignificant models if desired
  if( drop )
    rank.df <- rank.df  %>%
               filter( maxPr < sig )

  # Join effect sizes to ranks, and remove insignificant models
  rank.comb.df <- rank.df %>% 
                  left_join( y = combined.df, by = "Number" ) 

  rank.eff.df <-  rank.df %>%
                  left_join( y = effect.df, by = "Number" )                  

  rank.se.df <- rank.df %>% 
                left_join( y = se.df, by = "Number" ) 


  
  rank.comb.df[,"deltaAICc"] <- rank.comb.df[,"AICc"] - min(rank.comb.df[,"AICc"])
  rank.comb.df <- rank.comb.df[ order(rank.comb.df[,"deltaAICc"]), ]

  rank.eff.df[,"deltaAICc"] <- rank.eff.df[,"AICc"] - min(rank.eff.df[,"AICc"])
  rank.eff.df <- rank.eff.df[ order(rank.eff.df[,"deltaAICc"]), ]

  rank.se.df[,"deltaAICc"]  <- rank.se.df[,"AICc"] - min(rank.se.df[,"AICc"])
  rank.se.df                <- rank.se.df[ order(rank.se.df[,"deltaAICc"]), ]
  
  return( list( AICrank = rank.comb.df,
                AICrank.se = rank.se.df,
                AICrank.eff = rank.eff.df,
                models = modelList[ rank.df[ , "Number" ] ] )
        )
}


# .DASEexperiment()
# Function to apply the DASE methodology to the DASEex.csv table, my
# first attempt at discovering interactions between model inputs. Following
# exposition in Ch2 of Kleijnen, 2008: Design and Analysis of Simulation Experiments
# inputs:   tabName = character giving the root of the statTable
#           par = character indicating response variable parameter
#           species = character name of species to look at for effect
#           scaled = logical indicating whether to scale factor levels identically
#           rhs = character containing rhs of formula object for linear meta-model
# outputs:  fit = lm output for metamodel
.DASEexperiment <- function(  tabName = "lowK_rKq_RE",
                              resp = "q",
                              spec = "Dover",
                              scaled = TRUE,
                              baseRHS = NULL,
                              abs = FALSE,
                              intercept = FALSE,
                              newExp = c( "qOM", "UmsyOM", "BmsyOM" ) )
{
  tabFile <- paste(tabName, ".csv", sep = "")
  tabPath <- file.path(getwd(),"project/Statistics",tabFile)
  table <- read.csv( tabPath, header=TRUE, stringsAsFactors=FALSE ) 

  if( abs )
  {
    table[,resp] <- abs( table[,resp])
  }

  # restrict to nominated species/stock
  if(!is.null(spec)) 
    table <-  table %>%
              filter( species %in% spec )

  # HACK: reducing to the four MPs we want (throw away YE)
  table <-  table %>%
            filter( mp %in% c("noJointPriors","qPriorOnly","UmsyPriorOnly","qUpriors"))

  

  # Functions to rearrange the table using dplyr
  medTop <- function(x)
  {
    if(is.numeric(x)) out <- median(x)
    if(length(unique(x)) == 1) out <- x[1]

    out
  }

  # Take median of entries (MRE,MARE values), will just replicate the table
  # if using MRE/MARE data
  medTab <- table %>% 
            group_by(scenario,mp,species) %>%
            summarise_all(.funs=medTop) %>%
            ungroup()
  
  scaleNum <- function(x)
  {
    if(length(unique(x)) == 1) return(x[1])
    ranX <- range(x)
    gradX <- 2 / (ranX[2] - ranX[1])
    out <- gradX * x + (1 - gradX*ranX[2]) 
    out
  }
  pull    <- function(x,n) {x[[n]]}
  pullMP  <- function(x)
  {
    if(length(x) == 1) out <- x
    if(length(x) == 2) out <- x[[1]]
    if(length(x) == 3) out <- paste(x[[1]],x[[2]],sep = "_")
    out
  }

  cutScenario <- function(x,n)
  {
    x <- str_split(x,"_")
    x <- sapply(X = x, FUN = pull, n)
    x
  }

  cutMP <- function(x)
  {
    x <- str_split(x,"_")
    x <- sapply(X = x, FUN = pullMP )
    x
  }


  # medTab <- file.path(getwd(),"project/statistics",medTab)
  # medTab <- read.csv( medTab, header=TRUE, stringsAsFactors=FALSE )   
  # medTab <- medTab %>% filter( species == spec )

  # Scale inputs
  if( scaled )
  {

    table <-  table %>%
              mutate( nS = scaleNum(nS),
                      Umax = scaleNum(Umax),
                      fYear = scaleNum(fYear),
                      initDep = scaleNum(initDep),
                      nDiff = scaleNum(nDiff) )
  }

  # Function to fit GLMs given...
  fitGLM <- function( rhs = expVars[1],
                      resp = ssPar,
                      dat = table )
  {
    form    <- as.formula( paste( resp, "~", rhs, sep = "" ) )
    glmObj  <- glm( formula = form, data = dat, na.action = "na.omit" )

    return(glmObj)
  }

  # Function to make interactions out of main effects
  makeIntEff <- function ( main = expVars )
  {
    interactions <- character(length = 0)
    for( i in 1:(length(main)-1) )
    {
      for( j in (i+1):length(main) )
      {
        interactions <- c( interactions, paste( main[i], main[j], sep = ":" ) )
      }
    }
    interactions
  }

  # browser()
  table$mp <- as.factor(table$mp)
  table$mp <- relevel(table$mp, ref = "noJointPriors")

  # Generate the sets of possible effects
  # First, take the power set of main effects and collapse into formulae
  newExpCombos  <- powerset( newExp )
  newExpCombos  <- lapply( X = newExpCombos, FUN = paste, collapse = " + " )
  rhsList       <- paste( baseRHS, newExpCombos, sep = " + " )
  if( !intercept ) rhsList <- paste( rhsList, "0", sep = " + ")
  # Fit to all possible main effects for both models
  GLM       <- lapply( X = rhsList, FUN = fitGLM, resp = resp, dat = table )

  # A function to lapply to the glmObjects and compute the AIC
  medResids   <- function( glmObj, data, resp )
  {
    pred <- predict( glmObj, newdata = data )
    resids <- data[,resp] - pred
    resids
  }

  # Now run predict the medTab from the models
  resids    <- lapply( X = GLM, FUN = medResids, data = medTab, resp = resp )


  return( list (  GLM = GLM, resids = resids,
                  medTab = medTab ) )
}

# .makeMetaModelTables <- function(tabNameRoot)
# {
#   makeNewTabCols( tableRoot = "allSame_infoScenarios" )
#   tabName <- paste(tabNameRoot,"_MARE",sep = "")
#   .makeDeltaCols( tabName = tabName)

#   # Complex level pars
#   metaModels( tabName = paste(tabName,"_cplx",sep = ""),
#               multiResp = c("BnT", "Umsy","q_1","q_2","Dep","Bmsy"),
#               singleResp = c("DeltaBnT","DeltaUmsy","Deltaq_1","Deltaq_2","DeltaDep","DeltaBmsy"),
#               spec = c("Stock1"),
#               expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
#               sig = .05, intercept = TRUE,
#               scaled = TRUE, saveOut = TRUE, interactions = FALSE )

#   # Then stock level pars for stock1, with Delta values
#   metaModels( tabName = paste(tabName,"_Delta",sep = ""),
#               multiResp = c("BnT", "Umsy","q_1","q_2","Dep","Bmsy"),
#               singleResp = c("DeltaBnT","DeltaUmsy","Deltaq_1","Deltaq_2","DeltaDep","DeltaBmsy"),
#               spec = c("Stock1"),
#               expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
#               sig = .05, intercept = TRUE,
#               scaled = TRUE, saveOut = TRUE, interactions = FALSE )
# }


#.statTables()
# Wrapper for .statTableXXX functions, produces stacked tables of stats 
# for a group of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
# returns:  NULL
.statTables <- function(  sims=1,tabNameRoot = "statTable", par = F,
                          nCores = detectCores()-1,
                          IC = TRUE,
                          PI = TRUE)
{ 

  if(!dir.exists("./project/Statistics"))
    dir.create("./project/Statistics")

  if( par ) cluster <- makeCluster(nCores)

  

  # raw RE distributions
  REtable <- .statTableRE(  sims,paste(tabNameRoot,"_RE.csv", sep = ""),
                            par = par, clust = cluster )

  # raw RE distributions
  if(IC)
    ICtable <- .statTableIC(  sims,paste(tabNameRoot,"_ICraw.csv", sep = ""),
                              par = par, clust = cluster )

  # raw RE distributions
  if(PI)
    PItable <- .statTablePI(  sims,paste(tabNameRoot,"_PI.csv", sep = ""),
                              par = par, clust = cluster )

  if( par ) stopCluster( cluster )

  # Calculate run stats ( number of failed reps, number of repeated tries
  # for convergence on good reps )
  
  runStats <- REtable %>%
              group_by(scenario, mp, species) %>%
              summarise(  totReps = max(rep) ) %>%
              ungroup()

  if(PI)
  {
    browser()
    runStats2 <-  PItable %>%
                  group_by(scenario, mp, species) %>%
                  summarise(  msTries = mean(msTries),
                              ssTries = mean(ssTries) ) %>%
                  ungroup()

    runStats <- runStats %>%
                left_join(runStats2)
  }

  # MRE - summarised from RE
  MREtable <- REtable %>%
              dplyr::select(  scenario, mp,
                              species, nS,
                              fYear, Umax, initDep,
                              ssBnT, msBnT,
                              ssBmsy, msBmsy,
                              ssUmsy, msUmsy,
                              ssDep, msDep,
                              ssq_1, ssq_2,
                              msq_1, msq_2 ) %>%
              group_by(scenario, mp, species ) %>%
              summarise_all( .funs = .statTableSummarise ) %>%
              left_join( runStats )

  MREname   <- paste(tabNameRoot,"_MRE.csv", sep = "" )
  MREpath   <- file.path(getwd(),"project","Statistics",MREname)
  write.csv( MREtable, file = MREpath )

  expVarsTable <- MREtable %>%
                  filter(species == MREtable[1,"species"] ) %>%
                  dplyr::select( initDep,fYear,Umax )


  # Summarise interval coverage table
  if(IC)
  {
    ICtable <-  ICtable %>%
                    dplyr::select(  scenario, mp,
                                    species, nS,
                                    fYear, Umax, initDep,
                                    ssBnT, msBnT,
                                    ssBmsy, msBmsy,
                                    ssUmsy, msUmsy,
                                    ssDep, msDep,
                                    ssq_1, ssq_2,
                                    msq_1, msq_2 ) %>%
                    group_by( scenario, mp, species ) %>%
                    summarise_if(is.logical, .funs = funs(mean) ) %>%
                    mutate(nS = n()) %>%
                    ungroup() %>%
                    left_join( runStats ) %>%
                    left_join( expVarsTable)

    
    ICtabName <- paste( tabNameRoot, "_IC.csv", sep = "" )
    ICpath    <- file.path(getwd(),"project","Statistics",ICtabName)
    write.csv( ICtable, file = ICpath )
  }

  
  # MARE - summarised from RE
  MAREtable <-  REtable %>%
                dplyr::select(  scenario, mp,
                                species, nS,
                                fYear, Umax, initDep,
                                ssBnT, msBnT,
                                ssBmsy, msBmsy,
                                ssUmsy, msUmsy,
                                ssDep, msDep,
                                ssq_1, ssq_2,
                                msq_1, msq_2 ) %>%
                group_by( scenario, mp, species ) %>%
                summarise_all( .funs = .statTableSummariseAbs )%>%
                left_join( runStats )

  MAREname   <- paste(tabNameRoot,"_MARE.csv", sep = "" )
  MAREpath   <- file.path(getwd(),"project","Statistics",MAREname)
  write.csv( MAREtable, file = MAREpath )

  return(NULL)
}

#.statTableRE()
# Wrapper for .simStatRE, produces stacked tables of relative errors 
# for a group of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTableRE <- function (sims=1,tabName = "REstatTable.csv", par = FALSE, clust )
{ 
  # call function
  if( par ) tableList <- parLapply ( cl = clust, X = sims, fun = .simStatRE )
  else tableList <- lapply ( X = sims, FUN = .simStatRE )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}

#.statTableIC()
# Wrapper for .simStatIC, produces stacked tables of relative errors 
# for a group of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTableIC <- function (sims=1,tabName = "ICstatTable.csv", par = FALSE, clust )
{ 
  # call function
  if( par ) tableList <- parLapply ( cl = clust, X = sims, fun = .simStatIC )
  else tableList <- lapply ( X = sims, FUN = .simStatIC )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}

#.statTablePI()
# Wrapper for .simStatPI, produces stacked tables of relative errors 
# for a group of simulations.
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.statTablePI <- function (sims=1,tabName = "PIstatTable.csv", par = FALSE, clust )
{ 
  # call function
  if( par ) tableList <- parLapply ( cl = clust, X = sims, fun = .simStatPI )
  else tableList <- lapply ( X = sims, FUN = .simStatPI )

  # now make the table and return
  statTable <-  do.call("rbind",tableList)

  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( statTable, file = savePath )
  statTable
}


#.fitTable()
# Wrapper for .simFits, produces stacked tables of estimates
# and SEs for multiple fit runs
# inputs:   sims=integer vector indicating simulations in ./project/
# outputs:  statTable=table of statistics for a project/group
# usage:    to produce output for a project and create .csv tables of 
#           performance statistics
# side-eff: creates tables of statistics in ./project/stats/
.fitTable <- function (sims=1,tabName = "fitTable.csv", par = FALSE, clust )
{ 
  # call function
  if( par ) tableList <- parLapply ( cl = clust, X = sims, fun = .simFits )
  else tableList <- lapply ( X = sims, FUN = .simFits )

  # now make the table and return
  fitTable <-  do.call("rbind",tableList)
  savePath <- file.path(getwd(),"project","Statistics",tabName)
  write.csv ( fitTable, file = savePath )
  fitTable
}

# .statTableSummarise() & .statTableSummariseAbs()
# A custom function that will take the columns of the RE stat table
# and return a summarised version, using summarise_all(). Numeric
# columns have the median returned, single valued columns have the top
# value returned, and logical columns with multiple values have 
# their sum taken (this is for hessPD). 
# .statTableSummariseAbs() takes the absolute value of numeric columns
.statTableSummarise <- function(x)
{
  if(is.numeric(x)) 
  {
    out <- median(as.numeric(x), na.rm = T)
    return(out)
  }
  if(length(unique(x)) == 1) 
  {

    out <- x[1]
    return(out)
  }
  if(is.logical(x))
  {
    out <- sum(x,na.rm = T)
    return(out)
  }
  browser()
}

.statTableSummariseAbs <- function(x)
{
  if(is.numeric(x)) 
  {
    out <- median(abs(as.numeric(x)), na.rm = T)
    return(out)
  }
  if(length(unique(x)) == 1) 
  {

    out <- x[1]
    return(out)
  }
  if(is.logical(x))
  {
    out <- sum(x,na.rm = T)
    return(out)
  }
  out
}


# .simStatRE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simFits <- function ( sim=1 )
{
  # First, load blob
  source("tools.R")
  .loadSim(sim)

  om    <- blob$om
  opMod <- blob$opMod
  pars  <- blob$opMod$pars
  ss    <- blob$am$ss
  ms    <- blob$am$ms
  
  # Control info
  nS      <- blob$opMod$nS
  nSurv   <- blob$opMod$nSurv
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps
  fYear   <- blob$opMod$fYear
  lYear   <- blob$opMod$lYear

  goodReps <- blob$goodReps
  goodReps <- which(as.logical(apply(X =goodReps, FUN = prod, MARGIN = 1)))

  nT <- lYear - fYear + 1
  sT <- fYear - min(fYear) + 1
  maxT <- max(nT)

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","stock","ssAIC", "ssAICc","msAIC","msAICc",
                  "ssBnT","ssBnTse","msBnT","msBnTse",
                  "ssUmsy","ssUmsyse","msUmsy","msUmsyse",
                  "ssBmsy","ssBmsyse","msBmsy","msBmsyse",
                  "ssDnT","ssDnTse","msDnT","msDnTse",
                  "ssU_Umsy","ssU_Umsyse", "msU_Umsy", "msU_Umsyse",
                  "msHessPD", "ssHessPD","nReps", 
                  "nS", "nSurv")

  q_surv        <- paste( "q_", 1:nSurv, sep = "" )
  q_surv.se     <- paste( q_surv, "se", sep = "" )
  q_surv.seSS   <- paste( c("ss"), q_surv.se, sep = "" )
  q_surv.seMS   <- paste( c("ms"), q_surv.se, sep = "" )
  q_survSS      <- paste( c("ss"), q_surv, sep = "" )
  q_survMS      <- paste( c("ms"), q_surv, sep = "" )
  tau2_surv     <- paste( "tau2_", 1:nSurv, sep = "" )
  tau2_survSS   <- paste( c("ss"), tau2_surv, sep = "" )
  tau2_survMS   <- paste( c("ms"), tau2_surv, sep = "" )
  qbar_surv     <- paste( "qbar_", 1:nSurv, sep = "" )
  tauq2_surv    <- paste( "tauq2_", 1:nSurv, sep = "" )


  colLabels <- c( colLabels, 
                  q_survSS, q_survMS, q_surv.seSS, q_surv.seMS,
                  tau2_survSS, tau2_survMS, 
                  qbar_surv,
                  tauq2_surv )
  
  fitTable <- matrix( NA, nrow = nS*(length(goodReps)), ncol = length( colLabels ) )
  

  colnames(fitTable)    <- colLabels
  fitTable              <- as.data.frame(fitTable)

  # Start filling stat table
  # First, OM pars and labels
  fitTable$scenario        <- blob$ctrl$scenarioName
  fitTable$mp              <- blob$ctrl$mpLabel
  fitTable$stock           <- species[1:nS]
  fitTable$nS              <- opMod$nS
  fitTable$fYear           <- fYear[1:nS]
  fitTable$nSurv           <- nSurv
  fitTable$sigmaPriorCode  <- blob$assess$SigmaPriorCode
  fitTable$lnUPriorCode    <- blob$assess$lnUPriorCode
  fitTable$lnqPriorCode    <- blob$assess$lnqPriorCode
  fitTable$initBioCode     <- blob$assess$initBioCode[1:nS]


  # Now values that change with the replicate
  for( sIdx in 1:nS )
  {
    # Compute AICc for each model
    # survey q values, survey obs err var, 
    # 3 leading pars (B, r, PE)
    # initial biomass if estimated
    ssK <- nSurv*2 + 3 + blob$assess$initBioCode[sIdx]

    # q_os values, tau_o, B and r for each stock, initBioCodes, sigma for PE
    msK <- nS * nSurv + nSurv + 2*nS + sum(blob$assess$initBioCode[1:nS]) + 1
    # Add prior q if estiamted
    if( blob$assess$lnqPriorCode == 1 ) msK <- msK + 2 * nSurv 
    # Add prior U if estimated
    if( blob$assess$lnUPriorCode == 1 ) msK <- msK + 2
    # Add shared YE if estimated
    if( blob$assess$estYearEff == TRUE ) msK <- msK + 1

    # Compute AIC
    ssAIC <- 2*ssK + 2*ss$fitrep[[1]][[sIdx]]$nll
    msAIC <- 2*msK + 2*ms$fitrep[[1]]$nll

    nSampSS <- nT[sIdx]
    nSampMS <- sum(nT)
    # Compute AICc
    ssAICc <- ssAIC + 2*ssK*(ssK+1) / (nSampSS - ssK - 1)
    msAICc <- msAIC + 2*msK*(msK+1) / (nSampMS - msK - 1)

    if( length(ssAICc) == 0) ssAICc <- NA
    if( length(msAICc) == 0) msAICc <- NA

    fitTable[sIdx, c("ssAIC","ssAICc")] <- c(ssAIC, ssAICc)
    fitTable[sIdx, c("msAIC","msAICc")] <- c(msAIC, msAICc)

    # For those that exist, pull estimates and ses from CIs object
    ssCI <- ss$CIs[[1]][[sIdx]]
    msCI <- ms$CIs[[1]]

    if( is.null(ssCI) ) ssCI <- NA
    if( is.null(msCI) ) msCI <- NA

    # SS model
    # Check if converged
    if( !is.na(ssCI) )
    {
      # First, pull BnT
      ssBtrows <- which( ssCI$par == "lnBnT" )
      ssBnTrow <- max(ssBtrows)   
      fitTable[sIdx,c("ssBnT","ssBnTse")]        <- ssCI[ssBnTrow,c("val","se")]
      
      # Productivity
      ssUmsyRow <- which( ssCI$par == "lnUmsy" )
      fitTable[sIdx,c("ssUmsy","ssUmsyse")]        <- ssCI[ssUmsyRow,c("val","se")]
      
      # Optimal Biomass
      ssBmsyRow <- which( ssCI$par == "lnBmsy" )
      fitTable[sIdx,c("ssBmsy","ssBmsyse")]        <- ssCI[ssBmsyRow,c("val","se")]
      
      # Current depletion
      ssDnTRow <- which( ssCI$par == "lnDnT" )
      fitTable[sIdx,c("ssDnT","ssDnTse")]        <- ssCI[ssDnTRow,c("val","se")]

      # Current U_Umsy
      relUrows <- which( ssCI$par == "lnU_UmsyT" )
      ssrelUrow <- max( relUrows )
      fitTable[sIdx,c("ssU_Umsy","ssU_Umsyse")]        <- ssCI[ssrelUrow,c("val","se")]      
    }
      

    # Check if MS model converged
    if(!is.na(msCI))
    {
      # MS model
      # biomass
      msBtrows <- which( msCI$par == "lnBnT" )
      msBnTrow <- max(msBtrows) - nS + sIdx
      fitTable[sIdx,c("msBnT","msBnTse")]        <- msCI[msBnTrow,c("val","se")]
      # Prod
      msUmsyRow <- max(which( msCI$par == "lnUmsy" )) - nS + sIdx
      fitTable[sIdx,c("msUmsy","msUmsyse")]        <- msCI[msUmsyRow,c("val","se")]
      # opt biomass
      msBmsyRow <- max(which( msCI$par == "lnBmsy" )) - nS + sIdx    
      fitTable[sIdx,c("msBmsy","msBmsyse")]        <- msCI[msBmsyRow,c("val","se")]
      # current depletion
      msDnTRow <- max(which( msCI$par == "lnDnT" )) - nS + sIdx
      fitTable[sIdx,c("msDnT","msDnTse")]        <- msCI[msDnTRow,c("val","se")] 
      # relative fishing mortality (U/Umsy)
      relUrows <- which( msCI$par == "lnU_UmsyT" )
      msrelUrow <- max(relUrows) - nS + sIdx
      fitTable[sIdx,c("msU_Umsy","msU_Umsyse")]        <- msCI[msrelUrow,c("val","se")]

    }
  
    # Loop over surveys for catchability values
    for( oIdx in 1:nSurv )
    {
      sscolPar <- paste( "ssq_", oIdx, sep = "" )
      sscolSE  <- paste( sscolPar,"se", sep = "" )

      mscolPar <- paste( "msq_", oIdx, sep = "" )
      mscolSE  <- paste( mscolPar,"se", sep = "" )

      if( !is.na(ssCI) )
      {
        ssqRow  <- which(ssCI$par == "lnq_os" )[oIdx]
        fitTable[sIdx, c(sscolPar,sscolSE)]      <- ssCI[ssqRow, c("val","se")]  
      }
      if( !is.na(msCI) )
      {
        msqRow  <- max(which(msCI$par == "lnq_os" ) ) - nS * nSurv + (sIdx-1)*nSurv + oIdx
        fitTable[sIdx, c(mscolPar,mscolSE)]      <- msCI[msqRow, c("val","se")]  
      }
    } 

    # Avoid the uncertainty in these estimates - these aren't really meaningful
    fitTable[seq(sIdx,nS*length(goodReps),by = nS), "msHessPD"]    <- ms$hesspd[goodReps]
    fitTable[seq(sIdx,nS*length(goodReps),by = nS), "ssHessPD"]    <- ss$hesspd[goodReps,sIdx]
    fitTable[seq(sIdx,nS*length(goodReps),by = nS), tau2_survSS]   <- ss$tau2_o[goodReps,,sIdx]
    fitTable[seq(sIdx,nS*length(goodReps),by = nS), tau2_survMS]   <- ms$tau2_o[goodReps,]
    fitTable[seq(sIdx,nS*length(goodReps),by = nS), qbar_surv]     <- ms$qbar_o[goodReps,]
    fitTable[seq(sIdx,nS*length(goodReps),by = nS), tauq2_surv]    <- ms$tauq2_o[goodReps,]

  }

  fitTable <- fitTable %>%
              group_by( mp ) %>%
              mutate( ssAICcSum = sum(ssAICc),
                      ssBnTCV = sqrt(exp(ssBnTse^2) - 1) * ssHessPD,
                      msBnTCV = sqrt( exp(msBnTse^2) - 1 )* msHessPD,
                      ssUmsyCV = sqrt( exp(ssUmsyse^2) - 1 )* ssHessPD,
                      msUmsyCV = sqrt( exp(msUmsyse^2) - 1 )* msHessPD,
                      ssBmsyCV = sqrt( exp(ssBmsyse^2) - 1 )* ssHessPD,
                      msBmsyCV = sqrt( exp(msBmsyse^2) - 1 )* msHessPD,
                      ssDnTCV = sqrt( exp(ssDnTse^2) - 1 )* ssHessPD,
                      msDnTCV = sqrt( exp(msDnTse^2) - 1 )* msHessPD)
  
  # return
  fitTable
}


# .simStatRE()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of mean squared error BnT and Umsy
# usage:    in lapply to produce stats for a group of simulations
.simStatRE <- function ( sim=1 )
{
  # First, load blob
  source("tools.R")
  .loadSim(sim)

  om    <- blob$om
  opMod <- blob$opMod
  pars  <- blob$opMod$pars
  ss    <- blob$am$ss
  ms    <- blob$am$ms
  
  # Control info
  nS      <- blob$opMod$nS
  nT      <- blob$opMod$nT
  nSurv   <- blob$opMod$nSurv
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  goodReps  <- apply( X = success, FUN = prod, MARGIN = 1)
  goodReps  <- which(goodReps == 1)
  nGood     <- length(goodReps)

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "corr","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssDep","msDep",
                  "msHessPD", "ssHessPD","nReps",
                  "Umax", "tUpeak", 
                  "nS", "nSurv",
                  "UmsyOM", "BmsyOM", "rep",
                  "ssTries","msTries")

  q_surv        <- paste( "q_", 1:nSurv, sep = "" )
  q_survOM      <- paste( q_surv, "OM", sep = "" )
  q_survSS      <- paste( c("ss"), q_surv, sep = "" )
  q_survMS      <- paste( c("ms"), q_surv, sep = "" )
  tau2_surv     <- paste( "tau2_", 1:nSurv, sep = "" )
  tau2_survOM   <- paste( tau2_surv, "OM", sep = "" )
  tau2_survSS   <- paste( c("ss"), tau2_surv, sep = "" )
  tau2_survMS   <- paste( c("ms"), tau2_surv, sep = "" )
  qbar_surv     <- paste( "qbar_", 1:nSurv, sep = "" )
  qbar_survOM   <- paste( qbar_surv, "OM", sep = "" )
  tauq2_surv    <- paste( "tauq2_", 1:nSurv, sep = "" )
  tauq2_survOM  <- paste( tauq2_surv, "OM", sep = "" )


  colLabels <- c( colLabels, 
                  q_survSS, q_survMS, q_survOM, 
                  tau2_survSS, tau2_survMS, tau2_survOM, 
                  qbar_surv, qbar_survOM,
                  tauq2_surv, tauq2_survOM )
  
  statTable <- matrix( NA, nrow = nS*nGood, ncol = length( colLabels ) )
  

  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, OM pars and labels
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- apply(X = success, FUN = sum, MARGIN = 2)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUpeak          <- opMod$tUpeak
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$fYear           <- opMod$fYear[1:nS]
  statTable$initDep         <- opMod$initDep[1:nS]
  statTable$nSurv           <- nSurv
  statTable$fixProc         <- blob$ctrl$fixProc
  statTable$corrTargVar     <- blob$opMod$corrTargVar
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$initBioCode     <- blob$assess$initBioCode[1:nS]

  # Now values that change with the replicate
  for(gIdx in 1:nGood)
  {
    r <- goodReps[gIdx]
    for( survIdx in 1:nSurv )
    {
      statTable[ (gIdx-1)*nS+(1:nS), q_survOM[survIdx] ]      <- blob$om$q_os[r,survIdx,1:nS] 
      statTable[ (gIdx-1)*nS+(1:nS), qbar_survOM[survIdx] ]   <- blob$opMod$qSurvM[survIdx]
      statTable[ (gIdx-1)*nS+(1:nS), tauq2_survOM[survIdx] ]  <- blob$opMod$qSurvSD[survIdx]^2
      statTable[ (gIdx-1)*nS+(1:nS), tau2_survOM[survIdx] ]   <- blob$opMod$tauSurv[survIdx]^2
    }
    
    statTable[(gIdx-1)*nS+(1:nS),"rep"]           <- r
    statTable[(gIdx-1)*nS+(1:nS),"ssBnT"]         <- (ss$err.mle$BnT[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msBnT"]         <- (ms$err.mle$BnT[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssUmsy"]        <- (ss$err.mle$Umsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msUmsy"]        <- (ms$err.mle$Umsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssBmsy"]        <- (ss$err.mle$Bmsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msBmsy"]        <- (ms$err.mle$Bmsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssDep"]         <- (ss$err.mle$dep[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msDep"]         <- (ms$err.mle$dep[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssHessPD"]      <- (ss$hesspd[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msHessPD"]      <-  ms$hesspd[r]
    statTable[(gIdx-1)*nS+(1:nS),"sigU2"]         <- ( blob$am$ms$sigU2[r])
    statTable[(gIdx-1)*nS+(1:nS),"Umsybar"]       <- ( blob$am$ms$Umsybar[r])
    statTable[(gIdx-1)*nS+(1:nS), q_survSS]       <- ss$err.mle$q_os[r,,]
    statTable[(gIdx-1)*nS+(1:nS), q_survMS]       <- ms$err.mle$q_os[r,,]
    statTable[(gIdx-1)*nS+(1:nS), tau2_survSS]    <- ( ss$err.mle$tau2_o[r,,])
    statTable[(gIdx-1)*nS+(1:nS), tau2_survMS]    <- ( ms$err.mle$tau2_o[r,,])
    statTable[(gIdx-1)*nS+(1:nS), qbar_surv]      <- ( ms$err.mle$qbar_o[r,])
    statTable[(gIdx-1)*nS+(1:nS), tauq2_surv]     <- ( ms$err.mle$tauq2_o[r,])
    if(!is.null(ss$nTries) )
      statTable[(gIdx-1)*nS+(1:nS),"ssTries"]       <- (ss$nTries[r,] )
    if(!is.null(ms$nTries) ) 
      statTable[(gIdx-1)*nS+(1:nS),"msTries"]       <-  ms$nTries[r]

  }
  
  # return
  statTable
}

# .simStatIC()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of interval coverage
# usage:    in lapply to produce stats for a group of simulations
.simStatIC <- function ( sim=1 )
{
  # First, load blob
  source("tools.R")
  .loadSim(sim)

  om    <- blob$om
  opMod <- blob$opMod
  pars  <- blob$opMod$pars
  ss    <- blob$am$ss
  ms    <- blob$am$ms
  
  # Control info
  nS      <- blob$opMod$nS
  nT      <- blob$opMod$nT
  nSurv   <- blob$opMod$nSurv
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  goodReps  <- apply( X = success, FUN = prod, MARGIN = 1)
  goodReps  <- which(goodReps == 1)
  nGood     <- length(goodReps)

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "corr","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssDep","msDep",
                  "msHessPD", "ssHessPD","nReps",
                  "Umax", "tUpeak", 
                  "nS", "nSurv",
                  "UmsyOM", "BmsyOM", "rep",
                  "msTries","ssTries")

  q_surv        <- paste( "q_", 1:nSurv, sep = "" )
  q_survOM      <- paste( q_surv, "OM", sep = "" )
  q_survSS      <- paste( c("ss"), q_surv, sep = "" )
  q_survMS      <- paste( c("ms"), q_surv, sep = "" )


  colLabels <- c( colLabels, 
                  q_survSS, q_survMS, q_survOM )
  
  statTable <- matrix( NA, nrow = nS*nGood, ncol = length( colLabels ) )
  

  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, OM pars and labels
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- apply(X = success, FUN = sum, MARGIN = 2)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUpeak          <- opMod$tUpeak
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$fYear           <- opMod$fYear[1:nS]
  statTable$initDep         <- opMod$initDep[1:nS]
  statTable$nSurv           <- nSurv
  statTable$fixProc         <- blob$ctrl$fixProc
  statTable$corrTargVar     <- blob$opMod$corrTargVar
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$initBioCode     <- blob$assess$initBioCode[1:nS]

  # Now values that change with the replicate
  for(gIdx in 1:nGood)
  {
    r <- goodReps[gIdx]
    for( survIdx in 1:nSurv )
    {
      statTable[ (gIdx-1)*nS+(1:nS), q_survOM[survIdx] ]      <- blob$om$q_os[r,survIdx,1:nS] 
    }
    
    statTable[(gIdx-1)*nS+(1:nS),"rep"]           <- r
    statTable[(gIdx-1)*nS+(1:nS),"ssBnT"]         <- (ss$intCov$BnT[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msBnT"]         <- (ms$intCov$BnT[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssUmsy"]        <- (ss$intCov$Umsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msUmsy"]        <- (ms$intCov$Umsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssBmsy"]        <- (ss$intCov$Bmsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msBmsy"]        <- (ms$intCov$Bmsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssDep"]         <- (ss$intCov$dep[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msDep"]         <- (ms$intCov$dep[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssHessPD"]      <- (ss$hesspd[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msHessPD"]      <-  ms$hesspd[r]
    statTable[(gIdx-1)*nS+(1:nS),"ssTries"]       <- (ss$nTries[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msTries"]       <-  ms$nTries[r]
    statTable[(gIdx-1)*nS+(1:nS), q_survSS]       <- ss$intCov$q_os[r,,]
    statTable[(gIdx-1)*nS+(1:nS), q_survMS]       <- ms$intCov$q_os[r,,]

  }
  
  # return
  statTable
}


# .simStatIC()
# Produces a statistics table for leading pars in a simulation from
# the output produced by a runSimEst() call
# inputs:   sim=int indicating which simulation to compute stats for
# outputs:  statTable=data.frame of interval coverage
# usage:    in lapply to produce stats for a group of simulations
.simStatPI <- function ( sim=1 )
{
  # First, load blob
  source("tools.R")
  .loadSim(sim)

  om    <- blob$om
  opMod <- blob$opMod
  pars  <- blob$opMod$pars
  ss    <- blob$am$ss
  ms    <- blob$am$ms
  
  # Control info
  nS      <- blob$opMod$nS
  nT      <- blob$opMod$nT
  nSurv   <- blob$opMod$nSurv
  species <- blob$ctrl$speciesName[1:nS]
  nReps   <- blob$ctrl$nReps

  # get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  success <- blob$goodReps

  goodReps  <- apply( X = success, FUN = prod, MARGIN = 1)
  goodReps  <- which(goodReps == 1)
  nGood     <- length(goodReps)

  # First, create a data.frame of NAs with a row for each of MRE,MARE
  colLabels <- c( "scenario","mp","species","kappaTrue",
                  "SigmaTrue", "kappaMult", "corr","ssBnT","msBnT","ssUmsy","msUmsy",
                  "ssBmsy","msBmsy","ssDep","msDep",
                  "msHessPD", "ssHessPD","nReps",
                  "Umax", "tUpeak", 
                  "nS", "nSurv",
                  "UmsyOM", "BmsyOM", "rep","msTries","ssTries")

  q_surv        <- paste( "q_", 1:nSurv, sep = "" )
  q_survOM      <- paste( q_surv, "OM", sep = "" )
  q_survSS      <- paste( c("ss"), q_surv, sep = "" )
  q_survMS      <- paste( c("ms"), q_surv, sep = "" )


  colLabels <- c( colLabels, 
                  q_survSS, q_survMS, q_survOM )
  
  statTable <- matrix( NA, nrow = nS*nGood, ncol = length( colLabels ) )
  

  colnames(statTable)   <- colLabels
  statTable             <- as.data.frame(statTable)

  # get multiplier for shared effects
  kappaMult <- opMod$kappaMult
  if (is.null(opMod$kappaMult)) kappaMult <- 1

  # Start filling stat table
  # First, OM pars and labels
  statTable$scenario        <- blob$ctrl$scenarioName
  statTable$mp              <- blob$ctrl$mpLabel
  statTable$species         <- species[1:nS]
  statTable$kappaTrue       <- ifelse(is.null(opMod$kappa2),sqrt(opMod$pars$kappa2),sqrt(opMod$kappa2))*kappaMult
  statTable$SigmaTrue       <- ifelse(is.null(opMod$SigmaDiag),sqrt(opMod$pars$Sigma2[1:nS]),sqrt(opMod$SigmaDiag[1:nS]))
  statTable$kappaMult       <- kappaMult
  statTable$corr            <- ifelse(is.null(opMod$corrOffDiag),opMod$corrMult,opMod$corrOffDiag)
  statTable$nReps           <- apply(X = success, FUN = sum, MARGIN = 2)
  statTable$Umax            <- ifelse(is.null(opMod$Umax),opMod$Umult[2],opMod$Umax)
  statTable$tUpeak          <- opMod$tUpeak
  statTable$nS              <- opMod$nS
  statTable$UmsyOM          <- opMod$Umsy[1:nS]
  statTable$BmsyOM          <- opMod$Bmsy[1:nS]
  statTable$fYear           <- opMod$fYear[1:nS]
  statTable$initDep         <- opMod$initDep[1:nS]
  statTable$nSurv           <- nSurv
  statTable$fixProc         <- blob$ctrl$fixProc
  statTable$corrTargVar     <- blob$opMod$corrTargVar
  statTable$sigmaPriorCode  <- blob$assess$sigmaPriorCode
  statTable$initBioCode     <- blob$assess$initBioCode[1:nS]

  # Now values that change with the replicate
  for(gIdx in 1:nGood)
  {
    r <- goodReps[gIdx]
    for( survIdx in 1:nSurv )
    {
      statTable[ (gIdx-1)*nS+(1:nS), q_survOM[survIdx] ]      <- blob$om$q_os[r,survIdx,1:nS] 
    }
    
    statTable[(gIdx-1)*nS+(1:nS),"rep"]           <- r
    statTable[(gIdx-1)*nS+(1:nS),"ssBnT"]         <- (ss$predInt$BnT[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msBnT"]         <- (ms$predInt$BnT[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssUmsy"]        <- (ss$predInt$Umsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msUmsy"]        <- (ms$predInt$Umsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssBmsy"]        <- (ss$predInt$Bmsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msBmsy"]        <- (ms$predInt$Bmsy[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssDep"]         <- (ss$predInt$dep[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msDep"]         <- (ms$predInt$dep[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"ssHessPD"]      <- (ss$hesspd[r,] )
    statTable[(gIdx-1)*nS+(1:nS),"msHessPD"]      <-  ms$hesspd[r]
    statTable[(gIdx-1)*nS+(1:nS),"msTries"]       <-  ms$nTries[r]
    statTable[(gIdx-1)*nS+(1:nS),"ssTries"]       <-  ss$nTries[r,]
    statTable[(gIdx-1)*nS+(1:nS), q_survSS]       <- ss$predInt$q_os[r,,]
    statTable[(gIdx-1)*nS+(1:nS), q_survMS]       <- ms$predInt$q_os[r,,]

  }
  
  # return
  statTable
}

# .calcIntervalCoverage()
# This function takes a blob object produced by a sim-est procedure
# and produces distributions of absolute and relative errors
# inputs:   blob=list object output of simEstProc
# ouputs:   blob=list object with error distributions appended
# usage:    to create output that is saved for later analysis
.calcIntervalCoverage <- function ( blob )
{
  # recover control list, OM and AMs
  ctrl  <- blob$ctrl
  om    <- blob$om
  opMod <- blob$opMod
  ss    <- blob$am$ss
  ms    <- blob$am$ms

  # Get control constants
  nReps <- ctrl$nReps
  nS    <- opMod$nS
  nSurv <- opMod$nSurv
  fYear <- opMod$fYear
  lYear <- opMod$lYear
  nT    <- lYear - min(fYear) + 1

  # First get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  ssHess <- ss$hesspd
  msHess <- ms$hesspd
  hessPD <- ssHess
  for( s in 1:ncol(hessPD) )
    hessPD[,s] <- as.logical(ssHess[,s] * msHess)

  # Create a list to hold error values
  ssIC <- list ( 
                Bmsy  = matrix( NA, nrow = nReps, ncol = nS ),
                Umsy  = matrix( NA, nrow = nReps, ncol = nS ),
                q_os  = array(  NA, dim = c(nReps,nSurv,nS) ),
                dep   = matrix( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix( NA, nrow = nReps, ncol = nS )
              )
  # Slightly different structure for MS model
  msIC <- list ( 
                Bmsy    = matrix( NA, nrow = nReps, ncol = nS ),
                Umsy    = matrix( NA, nrow = nReps, ncol = nS ),
                q_os    = array(  NA, dim = c( nReps, nSurv, nS ) ),
                dep     = matrix( NA, nrow = nReps, ncol = nS ),
                BnT     = matrix( NA, nrow = nReps, ncol = nS )
              )

  # append IC lists to blob
  ss$intCov <- ssIC
  ms$intCov <- msIC

  # Fill in ss MLE relative errors
  for( rIdx in 1:nReps )
  {
    # Skip if no PD hessian
    if( !any( hessPD[ rIdx, ] ) ) next
    
    # Save MS sdreport
    sdrep.ms <- summary(ms$sdrep[[rIdx]])

    sdrep.ms <- data.frame( var = rownames(sdrep.ms),
                              est = sdrep.ms[,1],
                              sd  = sdrep.ms[,2] ) %>%
                  dplyr::mutate(  lbCI50 = est - .67*sd,
                                  ubCI50 = est + .67*sd)

    
    for ( s in 1:nS )
    {
      if( ! hessPD[ rIdx, s ] ) next

      sdrep.ss <- summary(ss$sdrep[[rIdx]][[s]]) 

      sdrep.ss <- data.frame( var = rownames(sdrep.ss),
                              est = sdrep.ss[,1],
                              sd  = sdrep.ss[,2] ) %>%
                  dplyr::mutate(  lbCI50 = est - .67*sd,
                                  ubCI50 = est + .67*sd)


      # SS interval coverage on management parameters
      ss$intCov$Bmsy[rIdx,s]    <- (log(opMod$Bmsy[s]) < sdrep.ss[sdrep.ss$var == "lnBmsy","ubCI50"]) & (log(opMod$Bmsy[s]) > sdrep.ss[sdrep.ss$var == "lnBmsy","lbCI50"])
      ss$intCov$Umsy[rIdx,s]    <- (log(opMod$Umsy[s]) < sdrep.ss[sdrep.ss$var == "lnUmsy","ubCI50"]) & (log(opMod$Umsy[s]) > sdrep.ss[sdrep.ss$var == "lnUmsy","lbCI50"])
      ss$intCov$q_os[rIdx,,s]   <- (log(om$q_os[rIdx,,s]) < sdrep.ss[sdrep.ss$var == "lnqhat_os","ubCI50"]) & (log(om$q_os[rIdx,,s]) > sdrep.ss[sdrep.ss$var == "lnqhat_os","lbCI50"])
      ss$intCov$dep[rIdx,s]     <- (log(om$dep[rIdx,s,nT]) < sdrep.ss[sdrep.ss$var == "lnDnT","ubCI50"]) & (log(om$dep[rIdx,s,nT]) > sdrep.ss[sdrep.ss$var == "lnDnT","lbCI50"])
      ss$intCov$BnT[rIdx,s]     <- (log(om$Bt[rIdx,s,nT]) < sdrep.ss[sdrep.ss$var == "lnBnT","ubCI50"]) & (log(om$Bt[rIdx,s,nT]) > sdrep.ss[sdrep.ss$var == "lnBnT","lbCI50"])

      # Now fill in ms interval coverage
      ms$intCov$Bmsy[rIdx,s]    <- (log(opMod$Bmsy[s]) < sdrep.ms[sdrep.ms$var == "lnBmsy","ubCI50"][s]) & (log(opMod$Bmsy[s]) > sdrep.ms[sdrep.ms$var == "lnBmsy","lbCI50"][s])
      ms$intCov$Umsy[rIdx,s]    <- (log(opMod$Umsy[s]) < sdrep.ms[sdrep.ms$var == "lnUmsy","ubCI50"][s]) & (log(opMod$Umsy[s]) > sdrep.ms[sdrep.ms$var == "lnUmsy","lbCI50"][s])
      ms$intCov$dep[rIdx,s]     <- (log(om$dep[rIdx,s,nT]) < sdrep.ms[sdrep.ms$var == "lnDnT","ubCI50"][s]) & (log(om$dep[rIdx,s,nT]) > sdrep.ms[sdrep.ms$var == "lnDnT","lbCI50"][s])
      ms$intCov$BnT[rIdx,s]     <- (log(om$Bt[rIdx,s,nT]) < sdrep.ms[sdrep.ms$var == "lnBnT","ubCI50"][s]) & (log(om$Bt[rIdx,s,nT]) > sdrep.ms[sdrep.ms$var == "lnBnT","lbCI50"][s])

      # Catchability takes some work cos we need to figure out what the order is
      repRows_qs <- 1:nSurv + nSurv * (s - 1) 
      ms$intCov$q_os[rIdx,,s]   <- (log(om$q_os[rIdx,,s]) < sdrep.ms[sdrep.ms$var == "lnqhat_os","ubCI50"][repRows_qs]) & (log(om$q_os[rIdx,,s]) > sdrep.ms[sdrep.ms$var == "lnqhat_os","lbCI50"][repRows_qs])
     
    }
  }
  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- hessPD

  blob
}

# .calcIntervalCoverage()
# This function takes a blob object produced by a sim-est procedure
# and produces distributions of absolute and relative errors
# inputs:   blob=list object output of simEstProc
# ouputs:   blob=list object with error distributions appended
# usage:    to create output that is saved for later analysis
.calcPredictiveInt <- function ( blob )
{
  # recover control list, OM and AMs
  ctrl  <- blob$ctrl
  om    <- blob$om
  opMod <- blob$opMod
  ss    <- blob$am$ss
  ms    <- blob$am$ms

  # Get control constants
  nReps <- ctrl$nReps
  nS    <- opMod$nS
  nSurv <- opMod$nSurv
  fYear <- opMod$fYear
  lYear <- opMod$lYear
  nT    <- lYear - min(fYear) + 1

  # First get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  ssHess <- ss$hesspd
  msHess <- ms$hesspd
  hessPD <- ssHess
  for( s in 1:ncol(hessPD) )
    hessPD[,s] <- as.logical(ssHess[,s] * msHess)

  # Create a list to hold error values
  ssPI <- list ( 
                Bmsy  = matrix( NA, nrow = nReps, ncol = nS ),
                Umsy  = matrix( NA, nrow = nReps, ncol = nS ),
                q_os  = array(  NA, dim = c(nReps,nSurv,nS) ),
                dep   = matrix( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix( NA, nrow = nReps, ncol = nS )
              )
  # Slightly different structure for MS model
  msPI <- list ( 
                Bmsy    = matrix( NA, nrow = nReps, ncol = nS ),
                Umsy    = matrix( NA, nrow = nReps, ncol = nS ),
                q_os    = array(  NA, dim = c( nReps, nSurv, nS ) ),
                dep     = matrix( NA, nrow = nReps, ncol = nS ),
                BnT     = matrix( NA, nrow = nReps, ncol = nS )
              )

  # append IC lists to blob
  ss$predInt <- ssPI
  ms$predInt <- msPI

  # Fill in ss MLE relative errors
  for( rIdx in 1:nReps )
  {
    # Skip if no PD hessian
    if( !any( hessPD[ rIdx, ] ) ) next
    
    # Save MS sdreport
    sdrep.ms <- summary(ms$sdrep[[rIdx]])

    sdrep.ms <- data.frame( var = rownames(sdrep.ms),
                              est = sdrep.ms[,1],
                              sd  = sdrep.ms[,2] )

    
    for ( s in 1:nS )
    {
      if( ! hessPD[ rIdx, s ] ) next

      sdrep.ss <- summary(ss$sdrep[[rIdx]][[s]]) 

      sdrep.ss <- data.frame( var = rownames(sdrep.ss),
                              est = sdrep.ss[,1],
                              sd  = sdrep.ss[,2] )


      # SS interval coverage on management parameters
      ss$predInt$Bmsy[rIdx,s]    <- pnorm(log(opMod$Bmsy[s]), mean = sdrep.ss[sdrep.ss$var == "lnBmsy","est"], sd = sdrep.ss[sdrep.ss$var == "lnBmsy","sd"])
      ss$predInt$Umsy[rIdx,s]    <- pnorm(log(opMod$Umsy[s]), mean = sdrep.ss[sdrep.ss$var == "lnUmsy","est"], sd = sdrep.ss[sdrep.ss$var == "lnUmsy","sd"])      
      ss$predInt$dep[rIdx,s]     <- pnorm(log(om$dep[rIdx,s,nT]), mean = sdrep.ss[sdrep.ss$var == "lnDnT","est"], sd = sdrep.ss[sdrep.ss$var == "lnDnT","sd"])
      ss$predInt$BnT[rIdx,s]     <- pnorm(log(om$Bt[rIdx,s,nT]), mean = sdrep.ss[sdrep.ss$var == "lnBnT","est"], sd = sdrep.ss[sdrep.ss$var == "lnBnT","sd"])
      for( oIdx in 1:nSurv)
        ss$predInt$q_os[rIdx,oIdx,s] <- pnorm(log(om$q_os[rIdx,oIdx,s]), mean = sdrep.ss[sdrep.ss$var == "lnqhat_os","est"][oIdx], sd =  sdrep.ss[sdrep.ss$var == "lnqhat_os","sd"][oIdx])

      # Now fill in ms interval coverage
      ms$predInt$Bmsy[rIdx,s]    <- pnorm(log(opMod$Bmsy[s]), mean = sdrep.ms[sdrep.ms$var == "lnBmsy","est"][s], sd = sdrep.ms[sdrep.ms$var == "lnBmsy","sd"][s])
      ms$predInt$Umsy[rIdx,s]    <- pnorm(log(opMod$Umsy[s]), mean = sdrep.ms[sdrep.ms$var == "lnUmsy","est"][s], sd = sdrep.ms[sdrep.ms$var == "lnUmsy","sd"][s])
      ms$predInt$dep[rIdx,s]     <- pnorm(log(om$dep[rIdx,s,nT]), mean = sdrep.ms[sdrep.ms$var == "lnDnT","est"][s], sd = sdrep.ms[sdrep.ms$var == "lnDnT","sd"][s])
      ms$predInt$BnT[rIdx,s]     <- pnorm(log(om$Bt[rIdx,s,nT]), mean = sdrep.ms[sdrep.ms$var == "lnBnT","est"][s], sd = sdrep.ms[sdrep.ms$var == "lnBnT","sd"][s])

      # Catchability takes some work cos we need to figure out what the order is
      repRows_qs <- 1:nSurv + nSurv * (s - 1) 
      for( oIdx in 1:nSurv )
        ms$predInt$q_os[rIdx,oIdx,s]   <- pnorm(log(om$q_os[rIdx,oIdx,s]), mean = sdrep.ms[sdrep.ms$var == "lnqhat_os","est"][repRows_qs[oIdx]], sd = sdrep.ms[sdrep.ms$var == "lnqhat_os","sd"][repRows_qs[oIdx]])
     
    }
  }
  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- hessPD

  blob
}





# makeErrorDists()
# This function takes a blob object produced by a sim-est procedure
# and produces distributions of absolute and relative errors
# inputs:   blob=list object output of simEstProc
# ouputs:   blob=list object with error distributions appended
# usage:    to create output that is saved for later analysis
.makeRelErrorDists <- function ( blob )
{
  # recover control list, OM and AMs
  ctrl  <- blob$ctrl
  om    <- blob$om
  opMod <- blob$opMod
  ss    <- blob$am$ss
  ms    <- blob$am$ms

  # Get control constants
  nReps <- ctrl$nReps
  nS    <- opMod$nS
  nSurv <- opMod$nSurv
  fYear <- opMod$fYear
  lYear <- opMod$lYear
  nT    <- lYear - min(fYear) + 1

  # First get the replicate numbers for succesful fits (MCMC runs) in BOTH models
  ssHess <- ss$hesspd
  msHess <- ms$hesspd
  hessPD <- ssHess
  for( s in 1:ncol(hessPD) )
    hessPD[,s] <- as.logical(ssHess[,s] * msHess)

  # Create a list to hold error values
  ssErr <- list ( 
                Bmsy  = matrix( NA, nrow = nReps, ncol = nS ),
                Umsy  = matrix( NA, nrow = nReps, ncol = nS ),
                kappa2= matrix( NA, nrow = nReps, ncol = nS ),
                tau2_o= array(  NA, dim = c(nReps,nSurv,nS) ),
                q_os  = array(  NA, dim = c(nReps,nSurv,nS) ),
                dep   = matrix( NA, nrow = nReps, ncol = nS ),
                BnT   = matrix( NA, nrow = nReps, ncol = nS ),
                totRE = matrix( NA, nrow = nReps, ncol = nS )
              )
  # Slightly different structure for MS model
  msErr <- list ( 
                Bmsy    = matrix( NA, nrow = nReps, ncol = nS ),
                Umsy    = matrix( NA, nrow = nReps, ncol = nS ),
                kappa2  = matrix( NA, nrow = nReps, ncol = 1 ),
                Sigma2  = matrix( NA, nrow = nReps, ncol = nS ),
                tau2_o  = array(  NA, dim = c( nReps, nSurv, nS ) ),
                tauq2_o = array(  NA, dim = c( nReps, nSurv ) ),
                q_os    = array(  NA, dim = c( nReps, nSurv, nS ) ),
                qbar_o  = array(  NA, dim = c( nReps, nSurv ) ),
                dep     = matrix( NA, nrow = nReps, ncol = nS ),
                BnT     = matrix( NA, nrow = nReps, ncol = nS ),
                totRE   = matrix( NA, nrow = nReps, ncol = nS )
              )

  # append error lists to blob
  # Single species only has one error term
  ss$err.mle <- ssErr

  # now append Sigma2 to the err list
  ms$err.mle <- msErr

  # Create matrices of the shared parameters
  # to avoid matrix-vector indexing errors
  tauSurv   <- matrix(opMod$tauSurv[1:nSurv], nrow = nReps, ncol = nSurv, byrow = T )
  qSurvM    <- matrix(opMod$qSurvM[1:nSurv], nrow = nReps, ncol = nSurv, byrow = T )
  qSurvSD   <- matrix(opMod$qSurvSD[1:nSurv], nrow = nReps, ncol = nSurv, byrow = T )

  # Fill in ss MLE relative errors
  for ( s in 1:nS )
  {
    # browser()
    ss$err.mle$Bmsy[hessPD[,s],s]    <- (ss$Bmsy[hessPD[,s],s] - opMod$Bmsy[s])/opMod$Bmsy[s]
    ss$err.mle$Umsy[hessPD[,s],s]    <- (ss$Umsy[hessPD[,s],s] - opMod$Umsy[s])/opMod$Umsy[s]
    ss$err.mle$kappa2[hessPD[,s],s]  <- (ss$kappa2[hessPD[,s],s] - (opMod$kappa2*(opMod$kappaMult^2)+opMod$SigmaDiag[s]))/(opMod$kappa2*(opMod$kappaMult^2)+opMod$SigmaDiag[s])
    ss$err.mle$tau2_o[hessPD[,s],,s] <- (ss$tau2_o[hessPD[,s],,s] - tauSurv[hessPD[,s],]^2)/(tauSurv[hessPD[,s],]^2)
    ss$err.mle$q_os[hessPD[,s],,s]   <- (ss$q_os[hessPD[,s],,s] - om$q_os[hessPD[,s],,s])/om$q_os[hessPD[,s],,s]
    ss$err.mle$dep[hessPD[,s],s]     <- (ss$dep[hessPD[,s],s,nT] - om$dep[hessPD[,s],s,nT])/om$dep[hessPD[,s],s,nT]
    ss$err.mle$BnT[hessPD[,s],s]     <- (ss$Bt[hessPD[,s],s,nT] - om$Bt[hessPD[,s],s,nT])/om$Bt[hessPD[,s],s,nT]
    ss$err.mle$totRE[hessPD[,s],s]   <- ss$err.mle$kappa2[hessPD[,s],s]

    # Now fill in ms MLE relative errors
    # some are only estimated once (instead of nS times)
    if (s == 1)
    {
      ms$err.mle$kappa2[hessPD[,s],]     <- (ms$kappa2[hessPD[,s],] - opMod$kappa2*(opMod$kappaMult^2))/opMod$kappa2/(opMod$kappaMult^2)
      ms$err.mle$tauq2_o[hessPD[,s],]    <- (ms$tauq2_o[hessPD[,s],] - qSurvSD[hessPD[,s],]^2)/(qSurvSD[hessPD[,s]]^2)
      ms$err.mle$qbar_o[hessPD[,s],]     <- (ms$qbar_o[hessPD[,s],] - qSurvM[hessPD[,s],])/qSurvM[hessPD[,s],]
    }
    # Now the rest of the pars
    ms$err.mle$tau2_o[hessPD[,s],,s] <- (ms$tau2_o[hessPD[,s],] - tauSurv[hessPD[,s],]^2)/(tauSurv[hessPD[,s],]^2)
    ms$err.mle$Sigma2[hessPD[,s],s]  <- (ms$Sigma2[hessPD[,s],s] - opMod$SigmaDiag[s])/opMod$SigmaDiag[s]    
    ms$err.mle$Bmsy[hessPD[,s],s]    <- (ms$Bmsy[hessPD[,s],s] - opMod$Bmsy[s])/opMod$Bmsy[s]
    ms$err.mle$Umsy[hessPD[,s],s]    <- (ms$Umsy[hessPD[,s],s] - opMod$Umsy[s])/opMod$Umsy[s]    
    ms$err.mle$q_os[hessPD[,s],,s]   <- (ms$q_os[hessPD[,s],,s] - om$q_os[hessPD[,s],,s])/om$q_os[hessPD[,s],,s]
    ms$err.mle$dep[hessPD[,s],s]     <- (ms$dep[hessPD[,s],s,nT] - om$dep[hessPD[,s],s,nT])/om$dep[hessPD[,s],s,nT]
    ms$err.mle$BnT[hessPD[,s],s]     <- (ms$Bt[hessPD[,s],s,nT] - om$Bt[hessPD[,s],s,nT])/om$Bt[hessPD[,s],s,nT]

    # calculate total RE variance
    totVarFit <- ms$kappa2 + ms$Sigma2[,s]
    totVarOM  <- opMod$kappa2*(opMod$kappaMult^2) + opMod$SigmaDiag[s]
    ms$err.mle$totRE[hessPD[,s],s]  <- (totVarFit[hessPD[,s]] - totVarOM) / totVarOM
  }

  # Append these to blob
  blob$am$ss <- ss
  blob$am$ms <- ms

  # Now save the good replicates
  blob$goodReps <- hessPD

  blob
}


