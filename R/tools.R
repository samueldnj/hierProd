# --------------------------------------------------------------------------
# tools.R
# 
# Script that contains tool functions for the coastwide multispecies 
# simulation-estimation procedure.
# 
# Author: Samuel D N Johnson
# Date: 24 August, 2016
#
#
# ToDo: 1. Complete list of tools
#       2. Write headers for saveSim()
#
#	List of tools:
#		lisread()			Function that reads ADMB style files as a list object (JS)
#   writeADMB()   Function that writes ADMB dat/pin files from list objects
#   saveSim()     Function that saves the output from a simulation run
#   read.admb()   Function that reads ADMB output (SJDM)
# 
# --------------------------------------------------------------------------


makeParEstTablePub <- function( tables = c("DoverAssessHigh","DoverAssessLow"),
                                models = c("qPriorOnly","qUpriors"),
                                pars = c("Umsy","BnT","DnT","U_Umsy","Bmsy","AICc"),
                                nStocks = 3,
                                saveName = "DoverStrongWeak.csv" )
{
  tabSkelParMod <- matrix( ncol = nStocks + 2, nrow = length(models) + 1 )

  tabList <- vector(mode = "list", length = length(tables) )
  names(tabList) = tables

  for( tIdx in 1:length(tables) )
  {
    tabName <- tables[tIdx]
    tabFile <- paste(tabName,".csv",sep = "")
    tabPath <- file.path(".","project","Statistics",tabFile)
    tab <- read.csv( tabPath, header = T, stringsAsFactors = F)
    parTabList <- vector( mode = "list", length = length(pars) )
    stocks <- unique( tab$stock )
    for( pIdx in 1:length( pars ) )
    {
      parTab <- tabSkelParMod
      parTab[,1] <- pars[pIdx]
      colnames(parTab) <- c("Parameter","Model",stocks)
      for( mIdx in 1:length(models) )
      {
        # filter down to the model
        tmpTab <- tab %>%
                  filter( mp == models[mIdx] )
        # If it's the first model then save the SS results
        if( mIdx == 1 )
        {
          ssPar <- paste( "ss", pars[pIdx], sep = "" )
          ssParSE <- paste(ssPar,"se",sep = "")
          parTab[1,2] <- "Single-Stock"
          if(pars[pIdx] != "AICc")
          {
            lCI     <- format(round(exp(tmpTab[,ssPar] - tmpTab[,ssParSE]),3),nsmall = 3)
            uCI     <- format(round(exp(tmpTab[,ssPar] + tmpTab[,ssParSE]),3),nsmall = 3)
            CV      <- format(round(sqrt(exp(tmpTab[,ssParSE]) - 1),2),nsmall = 2)
            parEst  <- format(as.numeric(round(exp(tmpTab[,ssPar]),3)),nsmall = 3)

            parTab[1,3:5] <- paste( parEst, " (", CV,")", sep = "")
          } else 
            parTab[1,3:5] <- format(round(as.numeric(tmpTab[,ssPar]),3),nsmall = 3)

        }
        # Now save the MS results
        msPar <- paste( "ms", pars[pIdx], sep = "" )
        msParSE <- paste( msPar, "se", sep = "" )
        parTab[1+mIdx,2] <- models[mIdx]
        if(pars[pIdx] != "AICc")
        {
          lCI     <- format(round(exp(tmpTab[,msPar] - tmpTab[,msParSE]),3),nsmall = 3)
          uCI     <- format(round(exp(tmpTab[,msPar] + tmpTab[,msParSE]),3),nsmall = 3)
          CV      <- format(round(sqrt(exp(tmpTab[,msParSE]) - 1),2),nsmall = 2)
          parEst  <- format(as.numeric(round(exp(tmpTab[,msPar]),3)),nsmall = 3)

          parTab[1+mIdx,3:5] <- paste( parEst, " (", CV ,")", sep = "")
        } else 
          parTab[1+mIdx,3:5] <- format(round(as.numeric(tmpTab[,msPar]),3),nsmall = 3)
      }
      parTab <- as.data.frame(parTab)
      parTabList[[pIdx]] <- parTab
    } 
    # Now bind the specific parameter tables together
    tabList[[tIdx]] <- do.call("rbind",parTabList)
    # tabList[[tIdx]]$tabName <- tabName
  }

  # And now the cbind each table
  aggTab <- do.call("cbind",tabList)

  # Now save
  savePath <- file.path(".","project","Statistics",saveName)
  write.csv(aggTab, savePath)

  aggTab
}

# Makes the metamodel ranking tables for
# factor effects
makeMetaModelTables <- function(  tabNameRoot, IC =TRUE,
                                  tabSavePath = "./Project/Statistics" )
{
  makeNewTabCols( tableRoot = tabNameRoot, IC = IC )
  tabName <- paste(tabNameRoot,"_MARE",sep = "")
  .makeDeltaCols( tabName = tabName)

  tabSavePath <- file.path( tabSavePath, tabNameRoot )
  if(!dir.exists(tabSavePath))
    dir.create(tabSavePath)

  # Complex level pars
  metaModels( tabName = paste(tabName,"_cplx",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = c("DeltaUmsy","DeltaBnT","DeltaBmsy","DeltaDep","Deltaq_1","Deltaq_2"),
              spec = c("Stock1"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = tabSavePath )

  # Then stock level pars for stock1, with Delta values
  metaModels( tabName = paste(tabName,"_Delta",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = c("DeltaUmsy","DeltaBnT","DeltaBmsy","DeltaDep","Deltaq_1","Deltaq_2"),
              spec = c("Stock1"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = tabSavePath )
  
  metaModels( tabName = paste(tabNameRoot,"_IC",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = NULL,
              spec = c("Stock1"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = tabSavePath )

  metaModels( tabName = paste(tabNameRoot,"_MRE",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = NULL,
              spec = c("Stock1"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = tabSavePath )

  # Then stock level pars for stock1, with Delta values
  metaModels( tabName = paste(tabName,"_Delta",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = c("DeltaUmsy","DeltaBnT","DeltaBmsy","DeltaDep","Deltaq_1","Deltaq_2"),
              spec = c("Stock4"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = file.path(tabSavePath,"Stock4") )
  
  metaModels( tabName = paste(tabNameRoot,"_IC",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = NULL,
              spec = c("Stock4"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = file.path(tabSavePath,"Stock4") )

  metaModels( tabName = paste(tabNameRoot,"_MRE",sep = ""),
              multiResp = c("Umsy","BnT","Bmsy","Dep", "q_1","q_2"),
              singleResp = NULL,
              spec = c("Stock4"),
              expVars = c("initDep","fYear","nDiff","Umax","nS","mp"),
              sig = .1, intercept = TRUE,
              scaled = TRUE, saveOut = TRUE, interactions = FALSE,
              tabSavePath = file.path(tabSavePath,"Stock4") )
}

# Converts the meta-model tables to the publication format
# Needs updating now that we're moving the SEs to a second row
makePubTableDF <- function( tableRoot = "test_epstSS",
                            tabPath = "./Project/Statistics" )
{
  tables <- list.files(file.path(tabPath,tableRoot), full.names = TRUE)
  tables <- tables[grepl(x = tables, pattern = ".csv") ]

  # OK, we have the two tables, one is cplx level and the other is
  # for stock 1
  stock1        <- tables[grepl(x = tables, pattern = paste(tableRoot,"_MARE_Delta",sep = ""))]
  stock1.int    <- stock1[grepl(x = stock1, pattern = "interleave")]
  stock1.inline <- stock1[!grepl(x = stock1, pattern = "interleave")]
  stock1.int    <- read.csv(stock1.int, header = T, stringsAsFactors =FALSE)
  stock1.inline <- read.csv(stock1.inline, header = T, stringsAsFactors =FALSE)

  cplx          <- tables[grepl(x = tables, paste(tableRoot,"_MARE_cplx",sep = ""))]
  cplx.int      <- cplx[grepl(x = cplx, pattern = "interleave")]
  cplx.inline   <- cplx[!grepl(x = cplx, pattern = "interleave")]
  cplx.int      <- read.csv(cplx.int, header = T, stringsAsFactors =FALSE)
  cplx.inline   <- read.csv(cplx.inline, header = T, stringsAsFactors =FALSE)

  # pull table row numbers for each section of table we're interested
  # in
  # Complex
  cplx.DeltaRows <- grepl(x = cplx.inline$par, pattern = "Delta")
  # Stock1
  Pars <- c("Umsy","BnT","Bmsy","Dep", "q_1","q_2")
  ssPars <- paste("ss",Pars,sep = "")
  msPars <- paste("ms",Pars,sep = "")
  stock1.DeltaRows <- grepl(x = stock1.inline$par, pattern = "Delta")
  stock1.ssRows <- which(stock1.inline$par %in% ssPars)
  stock1.msRows <- which(stock1.inline$par %in% msPars)

  # now rbind and tidy up
  outTable <- rbind(  stock1.inline[stock1.DeltaRows,], cplx.inline[cplx.DeltaRows,],
                      stock1.inline[c(stock1.ssRows,stock1.msRows),] ) %>%
              dplyr::select(  par,               # Response
                              beta0 = X.Intercept.,
                              q = mpqPriorOnly,      # 
                              r = mpUmsyPriorOnly,
                              qr = mpqUpriors,
                              D0 = initDep,
                              Tinit = fYear,
                              L = nDiff,
                              S = nS,
                              Fhist = Umax )

  outTableFileName  <- paste(tableRoot,"pubTable.csv",sep = "")
  outTablePath      <- file.path(tabPath,tableRoot)
  write.csv(outTable, file = file.path(outTablePath, outTableFileName) )

  return(outTable)
}


makeNewTabCols <- function( tableRoot = "allSame_infoScenarios", IC = FALSE )
{
  ## Path specific function - needs updating in mutate arguments ##
  # Cretes new column(s) in stats tables from the scenario label, 
  # used when factor levels for metamodels are not included in the 
  # control file (e.g. nDiff in infoScenarios)
  # inputs:   tableRoot = root file name of stats table
  #           ALSO REQUIRES UPDATING IN MUTATE ARGUMENTS
  # outputs:  NA
  # Side-eff: rewrites stat tables (MARE, MRE and RE versions) with
  #           new cols

  MAREtabRoot   <- paste(tableRoot,"_MARE", sep = "")
  MAREtabFile   <- paste(tableRoot,"_MARE.csv", sep = "")
  MAREtabPath   <- file.path(getwd(),"project","Statistics",MAREtabFile)
  MAREtab       <- read.csv( MAREtabPath, header=TRUE, stringsAsFactors=FALSE )

  MREtabRoot  <- paste(tableRoot,"_MRE", sep = "")
  MREtabFile  <- paste(tableRoot,"_MRE.csv", sep = "")
  MREtabPath  <- file.path(getwd(),"project","Statistics",MREtabFile)
  MREtab      <- read.csv( MREtabPath, header=TRUE, stringsAsFactors=FALSE )

  if(IC)
  {
    ICtabRoot   <- paste(tableRoot,"_IC", sep = "")
    ICtabFile   <- paste(tableRoot,"_IC.csv", sep = "")
    ICtabPath   <- file.path(getwd(),"project","Statistics",ICtabFile)
    ICtab       <- read.csv( ICtabPath, header=TRUE, stringsAsFactors=FALSE )
  }

  REtabRoot   <- paste(tableRoot,"_RE", sep = "")
  REtabFile   <- paste(tableRoot,"_RE.csv", sep = "")
  REtabPath   <- file.path(getwd(),"project","Statistics",REtabFile)
  REtab       <- read.csv( REtabPath, header=TRUE, stringsAsFactors=FALSE )


  splitScenLabel <- function ( scenLabel, match = "CVhi" )
  {
    splitLabel <- unlist(str_split(scenLabel,"_"))
    element <- splitLabel[grep( match, splitLabel, value = F)]
    splitElement <- unlist(str_split(element,match))
    num <- as.numeric(splitElement[2])
    num
  }

  ## UPDATE MUTATE ARGUMENTS ##
  MAREtab <-  MAREtab %>%
              group_by( scenario, mp, species ) %>%
              mutate( nDiff = splitScenLabel(scenario,match = "nDiff") )

  REtab <-  REtab %>%
            group_by( scenario, mp, species, rep ) %>%
            mutate( nDiff  = splitScenLabel(scenario,match = "nDiff") )

  if(IC)
    ICtab <-  ICtab %>%
              group_by( scenario, mp, species ) %>%
              mutate( nDiff  = splitScenLabel(scenario,match = "nDiff") )


  MREtab <- MREtab %>%
            group_by( scenario, mp, species ) %>%
            mutate( nDiff  = splitScenLabel(scenario,match = "nDiff") )


  write.csv( REtab, REtabPath )
  write.csv( MREtab, MREtabPath )
  write.csv( MAREtab, MAREtabPath )
  if(IC)
    write.csv( ICtab, ICtabPath )
}


# makeInfoScenarioDesign()
# Creates the .bch file for a historical fishing intensity experiment,
# without the MP section underneath.
makeInfoScenarioDesign <- function ( levels = list( Uhist = c("c(0.2,4,1)","c(1,1,1)"),
                                                    initYear = c(1984,2003),
                                                    nS = c(4,7,10),
                                                    initDep = c(.4,.7,1.0),
                                                    nDiff = c(0,1,2,3)
                                                  ),
                                    base = list(  initYear = 1984, initDep = 1,
                                                  nS = 10 ),
                                    nPoints = 1,
                                    bchName = "infoScenarios" )
{
  # First, make the LHR design
  lhrDesign <- makeLHRfromList( levels = levels, nPoints = nPoints )

  treatments <- lhrDesign$treatments

  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeInfoScenarioDesign() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(treatments) )
  {
    scenLabel <- ""
    for( cIdx in 1:ncol(treatments) )
    {
      scenLabel <- paste( scenLabel, colnames(treatments)[cIdx], 
                          levels[[ cIdx ]][ treatments[ rIdx, cIdx ] ], sep = "" )
      if( cIdx < ncol(treatments) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    nDiff <- levels$nDiff[ treatments[rIdx,"nDiff" ] ]
    nSame <- base$nS - nDiff
    initDep <- levels$initDep[ treatments[ rIdx, "initDep" ] ]
    initYear <- levels$initYear[ treatments[ rIdx, "initYear" ] ]
    if( (initDep != 1 | initYear != base$initYear) & nDiff > 0  ) initBioCode <- 1 else initBioCode <- 0
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Umult ", levels$Uhist[treatments[rIdx, "Uhist"] ], "\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$fYear c(rep(", levels$initYear[treatments[rIdx,"initYear"] ],
          ",", nDiff, "),rep(", base$initYear, ",", nSame, "))\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", levels$nS[treatments[rIdx,"nS"] ] ,"\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$initDep c(rep(", levels$initDep[treatments[rIdx,"initDep"] ],
          ",", nDiff, "),rep(", base$initDep, ",", nSame, "))\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$assess$initBioCode c(rep(", initBioCode,
          ",", nDiff, "),rep(0,", nSame, "))\n",
          sep = "", append = T, file = outFile  )

    cat( "#\n", file = outFile, append = T )
    
    cat( "#\n", file = outFile, append = T )
  }

  row.names(treatments) <- NULL

  treatments <- as.data.frame(treatments)
  for( k in 1:nrow(treatments) )
  {
    for( l in 1:length(names(levels)))
    {
      nom       <- names(treatments)[l]
      levelIdx  <- as.numeric(treatments[k, nom])
      treatments[ k, nom ] <- levels[nom][[1]][levelIdx]
    }
  }

  write.csv(x = treatments, file = paste(bchName,"LHrD.csv",sep="_"))

  treatments
}

# makeInfoScenarioDesign()
# Creates the .bch file for a historical fishing intensity experiment,
# without the MP section underneath.
makeFixedProcREDesign <- function ( levels = list(  m = c(0.5,1,1.5,2),
                                                    c = c(0.1,.3,0.5,.7,0.9),
                                                    Uhist = c("c(0.2,2,1)","c(1,1,1)"),
                                                    targVar = c("Ct","Bt"),
                                                    nS = c(4,7,10)
                                                  ),
                                    bchName = "FixedProcRE" ,
                                    nPoints = 1 )
{
  # First, make the LHR design
  lhrDesign <- makeLHRfromList( levels = levels, nPoints = nPoints )

  treatments <- lhrDesign$treatments

  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeFixedprocREDesign() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(treatments) )
  {
    scenLabel <- ""
    for( cIdx in 1:ncol(treatments) )
    {
      scenLabel <- paste( scenLabel, colnames(treatments)[cIdx], 
                          levels[[ cIdx ]][ treatments[ rIdx, cIdx ] ], sep = "" )
      if( cIdx < ncol(treatments) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Umult ", levels$Uhist[treatments[rIdx, "Uhist"] ], "\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$kappaMult ", levels$m[treatments[rIdx,"m"] ],"\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$ctrl$fixProc TRUE\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$corrOffDiag ", levels$c[treatments[rIdx,"c"] ],"\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", levels$nS[treatments[rIdx,"nS"] ],"\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$corrTargVar '", levels$targVar[treatments[rIdx,"targVar"] ],"'\n",
          sep = "", append = T, file = outFile  )
    cat( "#\n", file = outFile, append = T )
  }

  row.names(treatments) <- NULL

  treatments <- as.data.frame(treatments)
  for( k in 1:nrow(treatments) )
  {
    for( l in 1:length(names(levels)))
    {
      nom       <- names(treatments)[l]
      levelIdx  <- as.numeric(treatments[k, nom])
      treatments[ k, nom ] <- levels[nom][[1]][levelIdx]
    }
  }

  write.csv(x = treatments, file = paste(bchName,"LHrD.csv",sep="_"))

  treatments
}

# makeFhistDesign()
# Creates the .bch file for a historical fishing intensity experiment,
# without the MP section underneath.
makeFhistDesign <- function (  levels = list( nS  = seq(2,10,by=4),
                                              Umax  = c(0.8,1.4,2.0),
                                              tUpeak = c(2,6,10)
                                            ),
                                bchName = "Fhist" )
{
  # First, recover the number of factors
  k <- length(levels)
  
  # expand.grid the factor levels
  combos <- expand.grid(levels)
  
  totS <- max(levels$nS)
  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(combos) )
  {
    scenLabel <- ""
    repDiff <- combos[rIdx,"nDiff"]
    repBase <- totS - repDiff
    for( fIdx in 1:ncol(combos) )
    {
      scenLabel <- paste( scenLabel, colnames(combos)[fIdx], combos[rIdx,fIdx],sep = "" )
      if( fIdx < ncol(combos) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", combos[rIdx,"nS"], "\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$tUpeak", combos[rIdx,"tUpeak"] ,"\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$Umult c(0.2,", combos[rIdx,"Umax"], ",0.8)\n", 
          sep = "", append = T, file = outFile )    
    
    cat( "#\n", file = outFile, append = T )
  }
  cat( "# Management Procedures \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "# File Ends <not run>\n", append = T, file = outFile)

  combos
}

# makeREexpDesign()
# Creates the .bch file for the random effects and correlation experiment
makeREexpDesign <- function (  levels = list(   nS    = seq(2,10,by=4),
                                                m     = seq(0,2,by = 0.5),
                                                c     = seq(.2,.8,by=0.2)   
                                            ),
                                bchName = "REexp" )
{
  # First, recover the number of factors
  k <- length(levels)
  
  # expand.grid the factor levels
  combos <- expand.grid(levels)
  
  totS <- max(levels$nS)
  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(combos) )
  {
    scenLabel <- ""
    repDiff <- combos[rIdx,"nDiff"]
    repBase <- totS - repDiff
    for( fIdx in 1:ncol(combos) )
    {
      scenLabel <- paste( scenLabel, colnames(combos)[fIdx], combos[rIdx,fIdx],sep = "" )
      if( fIdx < ncol(combos) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", combos[rIdx,"nS"], "\n", 
          sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$Udevs 'corr'\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$kappaMult ", combos[rIdx,"m"] ,"\n",
          sep = "", append = T, file = outFile  )
    cat(  "scenario$scenario", rIdx, "$opMod$corrOffDiag ", combos[rIdx,"c"] ,"\n",
          sep = "", append = T, file = outFile  )
    cat( "#\n", file = outFile, append = T )
  }
  cat( "# Management Procedures \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "# File Ends <not run>\n", append = T, file = outFile)

  combos
}

# makeInitCondDesign()
# Creates the .bch file for the initial conditions experiment,
# without the MP section underneath.
makeInitCondDesign <- function (  levels = list(  nS        = seq(4,10,by=3),
                                                  nDiff     = 0:3,
                                                  sYear     = c(1995,2002),
                                                  initDep   = c(0.4,0.7),
                                                  survFreq  = c(1,2)
                                                ),
                                  bchName = "initConds" )
{
  # First, recover the number of factors
  k <- length(levels)
  
  # expand.grid the factor levels
  combos <- expand.grid(levels)
  
  totS <- max(levels$nS)
  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(combos) )
  {
    scenLabel <- ""
    repDiff <- combos[rIdx,"nDiff"]
    repBase <- totS - repDiff
    for( fIdx in 1:ncol(combos) )
    {
      scenLabel <- paste( scenLabel, colnames(combos)[fIdx], combos[rIdx,fIdx],sep = "" )
      if( fIdx < ncol(combos) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", combos[rIdx,"nS"], "\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$surveyFreq ", combos[rIdx,"survFreq"], "\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$sYear c(rep(", combos[rIdx,"sYear"], 
          ",", repDiff, "),rep(1988,",repBase, ")),\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$initDep c(rep(", combos[rIdx,"initDep"], 
          ",", repDiff, "),rep(1,",repBase, ")),\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$assess$initBioCode c(rep(1,",repDiff,"),rep(0,",repBase, ")),\n", sep = "", append = T, file = outFile )

    cat( "#\n", file = outFile, append = T )
  }

  # cat( "#\n", file = outFile, append = T )
  # cat( "# Management Procedures \n", file = outFile, append = T )
  # cat( "#\n", file = outFile, append = T )
  # cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  # cat( "# \n", append = T, file = outFile )
  # cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  # cat( "# \n", append = T, file = outFile )
  # cat( "# File Ends <not run>\n", append = T, file = outFile)

  combos
}

# makeExploreDesign()
# Creates the .bch batch design file for exploring the region where the 
# ms model did better, (see q prior only runs in 
# cwMSexperiments/.../1way_msInc_allSame experiment). 
makeExploreDesign <- function (  levels = list( Umax      = c(0.8,2),
                                                nS        = seq(4,10,by=3),
                                                nDiff     = 0:3,
                                                Umsy      = c(.1,.4),
                                                q         = c(.4,.8),
                                                obsCV     = c(.1,.3),
                                                sYear     = c(1994,2000)
                                              ),
                                bchName = "nonEqExplore" )
{
  # First, recover the number of factors
  k <- length(levels)
  # transform CVs to variances
  levels$obsCV <- round(log(levels$obsCV^2 + 1),2)
  
  # expand.grid the factor levels
  combos <- expand.grid(levels)
  
  totS <- max(levels$nS)
  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(combos) )
  {
    scenLabel <- ""
    repDiff <- combos[rIdx,"nDiff"]
    repBase <- totS - repDiff
    for( fIdx in 1:ncol(combos) )
    {
      scenLabel <- paste( scenLabel, colnames(combos)[fIdx], combos[rIdx,fIdx],sep = "" )
      if( fIdx < ncol(combos) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", combos[rIdx,"nS"], "\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Umult c(0.2,", combos[rIdx,"Umax"], ",0.8)\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Umsy c(rep(", combos[rIdx,"Umsy"], 
          ",", repDiff, "),rep(.25,",repBase, ")),\n", sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$q c(rep(", combos[rIdx,"q"], 
          ",", repDiff, "),rep(.6,",repBase, ")),\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$tau2 c(rep(", combos[rIdx,"obsCV"], 
          ",", repDiff, "),rep(.04,",repBase, ")),\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$sYear c(rep(", combos[rIdx,"sYear"], 
          ",", repDiff, "),rep(1988,",repBase, ")),\n", sep = "", append = T, file = outFile )
    cat( "#\n", file = outFile, append = T )
  }

  # cat( "#\n", file = outFile, append = T )
  # cat( "# Management Procedures \n", file = outFile, append = T )
  # cat( "#\n", file = outFile, append = T )
  # cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  # cat( "# \n", append = T, file = outFile )
  # cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  # cat( "# \n", append = T, file = outFile )
  # cat( "# File Ends <not run>\n", append = T, file = outFile)

  combos
}

# makeObsErrDesign()
# Creates a partial factorial design for the observation error
# simulation experiments.
# inputs:   levels = list of factor names and their levels
#           bchName = character root of batch file name
# ouputs:   table = design table data.frame
# side-eff: creates <bchName>.bch in working directory
makeObsErrDesign <- function (  levels = list(  nS        = seq(4,10,by=3),
                                                CVhi      = c(0.4,0.5),
                                                CVlo      = c(0.1,0.2),
                                                nHi       = c(0:4) 
                                              ),
                                bchName = "obsErr" )
{
  # First, recover the number of factors
  k <- length(levels)

  # transform CVs to variances
  levels$CVhi <- round(log(levels$CVhi^2 + 1),2)
  levels$CVlo <- round(log(levels$CVlo^2 + 1),2)
  # Count total nS
  totS <- max(levels$nS)
  
  # expand.grid the factor levels
  combos <- expand.grid(levels)
  

  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(combos) )
  {
    scenLabel <- ""
    for( fIdx in 1:ncol(combos) )
    {
      scenLabel <- paste( scenLabel, colnames(combos)[fIdx], combos[rIdx,fIdx],sep = "" )
      if( fIdx < ncol(combos) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    repHi <- combos[rIdx,"nHi"]
    repLo <- totS - repHi
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$nS ", combos[rIdx,"nS"], "\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$tau2 c(rep(", combos[rIdx,"CVhi"], 
          ",", repHi, "),rep(",combos[rIdx,"CVlo"],",",repLo, "))\n", sep = "", append = T, file = outFile )
    cat( "#\n", file = outFile, append = T )
  }

  cat( "#\n", file = outFile, append = T )
  cat( "# Management Procedures \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "# File Ends <not run>\n", append = T, file = outFile)

  combos
}

# makeFacDesignDover()
# Creates a fully factorial design for a simulation experiment,
# following DASE methodology (Kleijnen, 2008). Currently
# only modifies Dover Sole entries in DERPA complex parameters
# inputs:   levels = list of factor names and their levels
#           bchName = character root of batch file name
# ouputs:   table = design table data.frame
# side-eff: creates <bchName>.bch in working directory
makeFacDesignDover <- function (  levels = list(  Umsy  = seq(0.1,0.4, by = .1),
                                                  Bmsy  = c(seq(5,35,by = 10)),
                                                  q     = seq(0.2,0.8,by=.2)
                                                ),
                                  bchName = "rKqExp" )
{
  # First, recover the number of factors
  k <- length(levels)

  # expand.grid the factor levels
  combos <- expand.grid(levels)
  
  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(combos) )
  {
    scenLabel <- ""
    for( fIdx in 1:ncol(combos) )
    {
      scenLabel <- paste( scenLabel, colnames(combos)[fIdx], combos[rIdx,fIdx],sep = "" )
      if( fIdx < ncol(combos) ) scenLabel <- paste( scenLabel, "_", sep = "" )
    }
    cat( "# Scenario ", rIdx, " : ", scenLabel, "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$ctrl$scenarioName '", scenLabel, "'\n", 
          sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Umsy c(", combos[rIdx,"Umsy"], 
          ",rep(.25,9))\n", sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Bmsy c(", combos[rIdx,"Bmsy"], 
          ",rep(20,9))\n", sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$q c(", combos[rIdx,"q"], 
          ",rep(.6,9))\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$assess$mBmsy c(", combos[rIdx,"Bmsy"], 
          ",rep(20,9))\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$assess$s2Bmsy c(", combos[rIdx,"Bmsy"]^2, 
          ",rep(400,9))\n", sep = "", append = T, file = outFile )
    cat( "#\n", file = outFile, append = T )
  }

  cat( "#\n", file = outFile, append = T )
  cat( "# Management Procedures \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "# File Ends <not run>\n", append = T, file = outFile)

  combos
}

# makeDesignDover()
# Creates a l^(k-p) resolution III design for a simulation experiment,
# following DASE methodology (Kleijnen, 2008, Ch 2.4, ex 2.1). Currently
# only modifies Dover Sole entries in DERPA complex parameters
# inputs:   l = # of levels for each factor (make constant)
#           p = log2(fraction) of full factorial design
#           levels = list of factor names and their levels
#           bchName = character root of batch file name
# ouputs:   table = design table data.frame
# side-eff: creates <bchName>.bch in working directory
makeDesignDover <- function ( l = 2,
                              p = 4,
                              levels = list(  Umsy  = c( 0.292, 0.4 ),
                                              Bmsy  = c( 10, 15 ),
                                              q     = c( .6, .7 ),
                                              tau2  = c( .04, .06 ),
                                              kappa2= c( .025, .05),
                                              Sigma2= c( .025, .05),
                                              corrOD= c( 0, .5 )
                                            ),
                              generators = list( c(1,2), c(1,3), c(2,3),
                                                  c(1,2,3) ),
                              bchName = "DASEex" )
{
  # First, recover the number of factors
  k <- length(levels)
  # Squawk if there aren't enough combinations to saturate
  if( 2^(k-p) < k + 1 ) return("Not enough combinations to saturate, decrease p or k, or increase l")

  # Now, start building the design matrix
  desMtx <- matrix( 1, nrow = 2^(k-p), ncol = k, 
                    dimnames = list(  paste("exp",1:(2^(k-p)),sep=""),
                                      names(levels) ) )
  
  # Fill the first k-p columns with a fully factorial design
  for( c in 1:(k-p) )
  {
    pmPattern <- rep(x = c(+1,-1), times = c(2^(c-1),2^(c-1)))
    nReps     <- 2^(k-p) / length( pmPattern )
    desMtx[,c]<- rep( pmPattern, nReps )
  }
  
  # Now squawk if the list of generators is the wrong length 
  if( length(generators) != p ) return( "Wrong number of generators, revise and try again." )

  # Now apply generators to fill remaining columns
  for( gIdx in 1:p )
  {
    cIdx <- k - p + gIdx
    g <- generators[[gIdx]]
    for( genIdx in 1:length(g) ) desMtx[, cIdx] <- desMtx[,cIdx] * desMtx[,g[genIdx]]
  }

  # Create a copy of desMtx that has factor level list entry numbers
  entryMtx <- desMtx
  entryMtx[ entryMtx == 1 ] <- 2
  entryMtx[ entryMtx == -1 ] <- 1

  # Now start making the batch file for the simulation experiment
  outFile <- paste( bchName, ".bch", sep = "")
  cat(  "# Batch Control File, created ", date(), " by makeDesignDover() \n", 
        file = outFile, append = F, sep = "" )
  cat( "parameter value\n", sep = "", append = T, file = outFile)
  cat( "#\n", file = outFile, append = T )
  cat( "# Scenarios \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  # This will loop over the design matrix and create the scenario entry in the 
  # batch control file
  for( rIdx in 1:nrow(desMtx) )
  {
    cat( "# Scenario ", rIdx, " : ", rownames(desMtx)[rIdx], "\n", file = outFile, append = T, sep = "" )
    cat( "#\n", file = outFile, append = T )
    cat(  "scenario$scenario", rIdx, "$opMod$Umsy c(", levels$Umsy[entryMtx[rIdx,"Umsy"] ], 
          ",.29,.29,.29,.29)\n", sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$Bmsy c(", levels$Bmsy[entryMtx[rIdx,"Bmsy"] ], 
          ",10,10,10,10)\n", sep = "", append = T, file = outFile )  
    cat(  "scenario$scenario", rIdx, "$opMod$q c(", levels$q[entryMtx[rIdx,"q"] ], 
          ",.6,.6,.6,.6)\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$tau2 c(", levels$tau2[entryMtx[rIdx,"tau2"] ], 
          ",.04,.04,.04,.04)\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$kappa2 c(", levels$kappa2[entryMtx[rIdx,"kappa2"] ], 
          ")\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$SigmaDiag c(", levels$Sigma2[entryMtx[rIdx,"Sigma2"] ], 
          ",.025,.025,.025,.025)\n", sep = "", append = T, file = outFile )
    cat(  "scenario$scenario", rIdx, "$opMod$corrOffDiag c(", levels$corrOD[entryMtx[rIdx,"corrOD"] ], 
          ")\n", sep = "", append = T, file = outFile )
    cat( "#\n", file = outFile, append = T )
  }

  cat( "#\n", file = outFile, append = T )
  cat( "# Management Procedures \n", file = outFile, append = T )
  cat( "#\n", file = outFile, append = T )
  cat( "# MP 1 : baseAM \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "mp$mp1$ctrl$mpLabel 'baseAM' \n", append = T, file = outFile )
  cat( "# \n", append = T, file = outFile )
  cat( "# File Ends <not run>\n", append = T, file = outFile)

  desMtx
}

# Should only need to use a couple of times, to fix the naming in blob
fixBlob <- function ( sim = 1)
{
  # Check if blob is loaded, if not, load the blooob
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./project",full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  simList <- dirList[grep(pattern="sim",x=dirList)]
  folder <- simList[sim]

  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  assign( "blob",blob,pos=1 )
  cat("MSG (fixBlob) Simulation ", folder, " loaded from ./project/\n", sep="" )
  
  cat("MSG (fixBlob) Fixing blob names\n", sep="" )
  
  reps  <- blob$ctrl$nReps
  nT    <- length(blob$opMod$rep$yr)

  ss <- blob$am$ss
  ms <- blob$am$ms
  
  blob <- .makeRelErrorDists(blob)

  blob
  .saveSim(blob,folder,file.path(getwd(),"project",folder))
}

# makeBatch()
# Takes a batch control file and produces all the necessary structure
# to run the batch of sim-est experiments.
# inputs:     batchCtl = character naming the batch control file
# ouputs:     batchDesign = data.frame containing the batch design
# usage:      to run batch jobs for multiple scenario/mp combinations
# author:     S.D.N. Johnson
makeBatch <- function ( batchCtlFile = "batchControlFile.bch", prjFld = "project",
                        batchFld = "Batch", baseCtlFile = "simCtlFileBase.txt")
{
  .subChar <<- "__"
  # First, load the batch control file
  batchCtl <- .readParFile ( file.path(prjFld,batchFld,batchCtlFile) )
  # Now load the base control file
  baseCtl  <- .readParFile ( file.path(prjFld,batchFld,baseCtlFile) )

  # Set globals 
  #(THIS SHOULD BE MOVED TO ANOTHER FILE, OR DIFFERENT APPROACH FOUND)
  .PRJFLD <<- prjFld
  .DEFBATFLD <<- batchFld

  # Create Batch Design 
  batchDesign <- .createBatchDesign(  ctlPar = batchCtl, basePars = baseCtl )
  .FBATDES  <- file.path(getwd(),.PRJFLD,.DEFBATFLD,"batchDesign.txt")
  .CTLDES   <- file.path(getwd(),.PRJFLD,.DEFBATFLD,baseCtlFile)
  .BCTLDES  <- file.path(getwd(),.PRJFLD,.DEFBATFLD,batchCtlFile)
  .writeDesignFile (obj=batchDesign,out=.FBATDES)
  file.copy(baseCtlFile, .CTLDES)
  file.copy(batchCtlFile, .BCTLDES)

  return(batchDesign)
}

# Calls runSimEst in parallel inside a separate copy of the
# working directory, allowing for parallel system calls
doBatchRun <- function( arg )
{
  require(tools)
  cat("Running batchjob:", arg[1],"\n")
  # source control script to load DLL
  # source("control.r")
  
  # runMSE with the batch file
  # add random delay to offset simFolder names
  runSimEst(ctlFile=arg[1], folder=arg[2])
  return(NULL)
}


# .runBatchJob  (Runs the simulations specified in a design dataframe)
# Purpose:      Loops through the rows of the design dataframe, each of which
#               specifies an input parameter file and labels for a simulation.
#               The mseR function runMSE is called for each row which generates
#               the simulation results (e.g., blob) and simulation folders.
# Parameters:   batchDesign is a batch design dataframe created by the function
#                 .createBatchDesign.
#               prefix is the Simulation File Prefix/
# Returns:      NULL (invisibly)
# Side Effects: A simulation folder containing *.info and *.Rdata file (blob) and
#               for each row of the design dataframe, i.e., for each simulation.
# Source:       A.R. Kronlund
.runBatchJob <- function( batchDesign=NULL, par=FALSE,prefix=NULL, initPar = 1,
                          subset = NULL )
{
  # Runs simulations from the design data.frame specified in batchDesign object.
  # 1. Does the mseR input parameter file exist? If YES then read the file.
  # 2. Run the simulation using runMSE().
  # 3. Update the design data.frame somehow...
  if (is.null(batchDesign))
  {
    desPath <- file.path(getwd(),.PRJFLD,.DEFBATFLD,"batchDesign.txt")
    batchDesign <- read.csv(desPath, header=TRUE, skip=1, stringsAsFactors=FALSE)
  } 

  # Vector of simCtlFile names, like simCtlFile1.txt, simCtlFile2.txt, etc.
  batchParFile <- file.path( getwd(),.PRJFLD, .DEFBATFLD, basename( batchDesign$parFile ) )
  blobName     <- batchDesign$blobName  # Vector of blobNames.
  nSims        <- length(blobName)      # Number of blobs (i.e., simulations).

  # closeWin()                     # Close GUI to avoid TclTk ugliness.

  result <- data.frame( sim=character(nSims), scenarioLabel=character(nSims),
                        mpLabel=character(nSims), elapsed=numeric(nSims),
                        stringsAsFactors=FALSE )

  # This is part where shell script PPSS could be used to parallelize the batch job
  # Need to create a folder of inputParameters.par files to be processed via runMSE().
  # Might need to have the .par file be an argument to runMSE()...check what effects
  # that would have elsewhere.

  # if a parallel flag is set, run in parallel using multiple wds
  if (par)
  { 
    # First, load parallel package
    library(parallel)
    # turn off squawking
    options(warn=-1)
    # Get number of batch runs
    nBatchFiles   <- length(batchParFile)

    # Create folder names for batch running
    batchFolderNames <- paste("parBat",prefix,1:nBatchFiles,sep="")
    
    # combine folder and control file names
    if( !is.null( subset ) )
    {
      parBatchArgList <- vector(mode = "list", length = length(subset) )
      for( i in 1:length(subset) )
      {
        parBatchArgList[[i]] <- c( batchParFile[ subset[i] ],batchFolderNames[ subset[i] ])   
      }
    } else {
      parBatchArgList <- vector(mode = "list", length = length(batchParFile) - initPar + 1)
      for(i in initPar:nBatchFiles)
      {
        parBatchArgList[[i - initPar + 1]] <- c( batchParFile[i],batchFolderNames[i])
      }
    }

    nSims <- length(parBatchArgList)

    # Now set # of cores and make a cluster
    nCores  <- min(length(parBatchArgList),detectCores()-1)
    cl      <- makeCluster(nCores, outFile = "./parBatchMsg.txt")
    # Run parallel batch
    cat ("Running ", nSims, " simulations in parallel on ",
          nCores, " cores.\n", sep = "" )
    tBegin    <- proc.time()
    startDate <- date()
    clusterEvalQ( cl, expr = {source('control.R')} )
    tmp       <- clusterApplyLB(cl, x=parBatchArgList, fun=doBatchRun)
    # tmp <-lapply(X=parBatchArgList, FUN=doBatchRun)
    stopCluster(cl)

    # Now copy the contents of each batchFolderName to the project folder
    elapsed <- (proc.time() - tBegin)[ "elapsed" ]
    cat( "\nMSG (.runBatchJob): Elapsed time for parallel batch = ",
      round(elapsed/60.0,digits=2)," minutes.\n" )

  } else for ( i in initPar:nSims ) {
    
    if ( file.exists( batchParFile[i] ) )
    {
      fileName <- strsplit( batchParFile[i],"\\." )[[1]][1]
      cat( "\nRunning batch job: ",batchParFile[i],"...\n" )

      tBegin    <- proc.time()
      startDate <- date()

      cat( "\nMSG (.runBatchJob) Processing batchParFile = ", batchParFile[i], "\n" )
      
      # Overwrite the simCtlFile.txt file for mseR, the overwrite=TRUE is CRITICAL.
      
      cat( "\nWARNING (.runBatchJob) simCtlFile.txt overwritten by ",
           batchParFile[i],"\n" )
      file.copy(batchParFile[i],
                 "simCtlFile.txt", overwrite=TRUE )
     
      # runSimEst() assumes that input is simCtlFile.txt
      # browser()
      runSimEst()

      elapsed <- (proc.time() - tBegin)[ "elapsed" ]
      cat( "\nMSG (.runBatchJob): Elapsed time for simulation = ",
        round(elapsed/60.0,digits=2)," minutes.\n" )
    }
    
    result[ i,"sim" ]           <- batchDesign[ i,"parFile" ]
    result[ i,"scenarioLabel" ] <- batchDesign[ i,"scenarioLabel" ]
    result[ i,"mpLabel" ]       <- batchDesign[ i,"mpLabel" ]
    result[ i,"elapsed" ]       <- round( elapsed/60.0, digits=2)
    
    cat( "\nMSG (.runBatchJob) Progress as of ", date()," Simulation ",i," of ",nSims,"\n\n" )
    print( result[1:i,] )
  }

  cat( "\n" )
  invisible()
}     # END function .runBatchJob


# .writeDesignFile (Writes the batch job design dataframe to a file)
# Purpose:         Given a dataframe, assumed to contain a batch job design,
#                  created by .createBatchDesign, write to text file.
# Parameters:      obj is a dataframe.
#                  outFile is the desired filename.
# Returns:         NULL (invisibly)
# Side Effects:    Call to this function may generate a warning, not sure why.
# Source:          A.R. Kronlund
.writeDesignFile <- function( obj, outFile=.FBATDES )
{
  nRow <- nrow( obj )
  if ( is.null(nRow) )
    nRow <- 1

  cat( file=outFile, "# mseR batch design file, ",date(),"\n" )
  
  # This generates a warning, not sure why.
  write.table( obj, file=outFile,
               col.names=TRUE, row.names=FALSE, quote=TRUE, sep=",", append=TRUE )
  cat( "\nMSG (.writeDesignFile): Batch design file written to ",outFile,"\n" )
  invisible()
}     # END function .writeDesignFile


# .createBatchDesign  (Creates design dataframe and writes input par files)
# Purpose:      Given a batch control list, base parameter file, and simulation
#               file prefix, this function creates the batch job design by
#               crossing each management procedure node with each scenario node.
#               If the number of scenarios is nScenario and the number of
#               management procedures is nMP, then the result is a design
#               dataframe with nScenario*nMP rows and a *.par file corresponding
#               to each row (combination of scenario and management procedure).
#               The output design dataframe is passed to the .runBatchJob to
#               control completing the simulations using runMSE.
# Parameters:   ctlList is a list created by .createKeyList from the Batch File
#               basePars is a dataframe containing the base parameters that are
#                 modified by the ctlList values;
#               prefix is the Simulation File Prefix string concatenated to the
#               design filename, output parameter files, and simulation results
#               in .Rdata files (e.g., blobs).
# Returns:      result is the batch job design dataframe.
# Side Effects: The *.par files that specify each simulation are written.
# Source:       A.R. Kronlund, and probably under questionable conditions.
.createBatchDesign <- function( ctlPar, basePars, prefix="job" )
{
  # Input a control list object, usually the output from .createList.
  ctlList <- .createKeyList(ctlPar)
  
  scenario  <- ctlList$scenario              # List of scenarios.
  nScenario <- length( scenario )            # Number of scenarios.

  mp  <- ctlList$mp                          # List of management procedures.
  nMP <- length( mp )                        # Number of management procedures.

  nBatch <- nScenario * nMP                  # Number of mseR simulations.

  # Design dataframe - each row identifies the simulation and adds colors,
  # line types, line widths and symbols to be used for plotting outside of mseR.
  # This file is needed as input to runBatchJob to control the job, and for use
  # in plotting and performance calculations are to be done outside of mseR.
  
  parFile      <- character( nBatch )        # Name of each input parameter file.
  blobName     <- character( nBatch )        # Name of each blob.

  scenarioLabel <- character( nBatch )       # Name of each scenario.
  mpLabel       <- character( nBatch )       # Name of each management procedure.
  
  dataName     <- rep( "data",     nBatch )  # Name of data method.
  methodName   <- rep( "method",   nBatch )  # Name of assessment method.
  ruleName     <- rep( "rule",     nBatch )  # Name of HCR.
  scCol        <- rep( "black",    nBatch )  # Color for scenario.
  scLty        <- rep( 1,          nBatch )  # Line type for scenario.
  scLwd        <- rep( 1,          nBatch )  # Line width for scenario.
  scSym        <- rep( 1,          nBatch )  # Symbol for scenario.
  mpCol        <- rep( "black",    nBatch )  # Color for procedure.
  mpLty        <- rep( 1,          nBatch )  # Line type procedure.
  mpLwd        <- rep( 1,          nBatch )  # Line width for procedure.
  mpSym        <- rep( 1,          nBatch )  # Symbol for procedure.
  join         <- rep( 0,          nBatch )  # Join... ???

  iBatch <- 1

  # Loop over the scenarios.
  for ( i in 1:nScenario )
  {
    # Loop over the management procedures.
    for ( j in 1:nMP )
    {
      parFile[iBatch]   <- paste( prefix,iBatch,".txt",sep="" )
      
      # Create a unique blobname.
      # NOTE: This step is deferred until the simulation is launched to obtain
      #       a unique date-time stamp at time of execution.  For now simply
      #       provide a name indexed to the simulation number.
      blobName[iBatch] <- paste( "sim",prefix,iBatch,sep="" )
      
      # Set default scenarioLabel and mpLabel values in case not supplied.
      scenarioLabel[iBatch] <- paste( "scenario",i,sep="" )
      mpLabel[iBatch]       <- paste( "mp",j,sep="" )

      # Use the scenarioLabel if provided.
      #if ( !is.null(scenario[[i]]$scenarioLabel) )
      #  scenarioLabel[iBatch] <- scenario[[i]]$scenarioLabel

      # Scenario colour.
      if ( !is.null(scenario[[i]]$scCol ) )
        scCol[iBatch] <- scenario[[i]]$scCol

      # Scenario line type.
      if ( !is.null(scenario[[i]]$scLty ) )
        scLty[iBatch] <- scenario[[i]]$scLty

      # Scenario line width.
      if ( !is.null(scenario[[i]]$scLwd ) )
        scLwd[iBatch] <- scenario[[i]]$scLwd

      # Scenario symbol.
      if ( !is.null(scenario[[i]]$scSym ) )
        scSym[iBatch] <- scenario[[i]]$scSym

      # Use the mpLabel if provided.
      #if ( !is.null(mp[[j]]$mpLabel) )
      #  mpLabel[iBatch] <- mp[[j]]$mpLabel

      # Use the dataName if provided.
      if ( !is.null(mp[[j]]$dataName) )
        dataName[iBatch] <- mp[[j]]$dataName

      # Use the methodName if provided.
      if ( !is.null(mp[[j]]$methodName ) )
        methodName[iBatch] <- mp[[j]]$methodName

      # Use the ruleName if provided.
      if ( !is.null(mp[[j]]$ruleName) )
        ruleName[iBatch] <- mp[[j]]$ruleName

      # Management procedure colour.
      if ( !is.null(mp[[j]]$mpCol ) )
        mpCol[iBatch] <- mp[[j]]$mpCol

      # Management procedure line type.
      if ( !is.null(mp[[j]]$mpLty ) )
        mpLty[iBatch] <- mp[[j]]$mpLty

      # Management procedure line width.
      if ( !is.null(mp[[j]]$mpLwd ) )
        mpLwd[iBatch] <- mp[[j]]$mpLwd

      # Management procedure symbol.
      if ( !is.null(mp[[j]]$mpSym ) )
        mpSym[iBatch] <- mp[[j]]$mpSym

      # Join the procedures (groups share common integer value).
      if ( !is.null(mp[[j]]$join ) )
        join[iBatch] <- mp[[j]]$join

      # Replace values for any keywords that match those in the mseR input
      # parameters file:
      # 1. Find keywords shared between the batch job control list and the mseR
      #    parameter file.
      # 2. Replace the values in the mseR parameter file with those from the
      #    batch control list.
      
      newPars <- basePars
      
      # This step compares the names in the i-th scenario to the names in the
      # first column of newPars to find the common names using intersect.
      # browser()

      # Change all .subChar symbols to "$".
      names(scenario[[i]]) <- gsub( .subChar,"$",names(scenario[[i]]), fixed=TRUE )
       
      sharedNames <- intersect( names(scenario[[i]]),newPars[,1] )

      # Create a batch root, and copy new entries from the batchPar list
      scenarioRoot <- paste("scenario$scenario",i,"$",sep="")
      sharedNamesPar <- paste(scenarioRoot,sharedNames,sep="")
 
      if ( length(sharedNames) > 0 )
        for ( k in 1:length(sharedNames) )
        {
          val <- ctlPar[,2][ ctlPar[,1] == sharedNamesPar[k] ]
          # if ( is.character(val) )
          #   val <- paste("\"",val,"\"",sep="")
          newPars[,2][ newPars[,1]==sharedNames[k] ] <- val
        }

      # Change all .subChar symbols to "$".
      names( mp[[j]] ) <- gsub( .subChar,"$",names(mp[[j]]), fixed=TRUE )
       
      sharedNames <- intersect( names(mp[[j]]),newPars[,1] )
      # Create an MP root for copying from batchPars table
      mpRoot <- paste ("mp$mp",j,"$",sep="")
      sharedNamesPar <- paste(mpRoot,sharedNames,sep="")
      
      if ( length(sharedNames) > 0 )
        for ( k in 1:length(sharedNames) )
        {
          val <- ctlPar[,2][ ctlPar[,1] == sharedNamesPar[k] ]
          # if ( is.character(val) )
          #   val <- dQuote( val )
          newPars[,2][ newPars[,1]==sharedNames[k] ] <- val
        }

      # Check to see if scenarioLabel and mpLabel are updated.
      scenarioLabel[iBatch] <- newPars[,2][ newPars[,1]=="ctrl$scenarioName" ]
      mpLabel[iBatch] <- newPars[,2][ newPars[,1]=="ctrl$mpLabel" ]
      
      # Remove leading white space: gsub('^[[:space:]]+', '', " a test ")
      newPars[,2] <- gsub( '^[[:space:]]+','',newPars[,2] )
      
      fName <- parFile[iBatch]          # mseR input parameters file.
      
      fName <- file.path( .PRJFLD, .DEFBATFLD, fName )
      batchDate <- date()               # Current date and time.      
      
      # Open a new parameter file and write the title and date/time stamp.
      cat( file=fName,
        "# ",parFile[iBatch],": mseR parameter file written ",batchDate,".\n", sep="" )

      # NOTE: write.table wants a matrix or data.frame, not a vector.
 
      # Write the header field names for the parameter file.
      
      colNames <- names( basePars )
      #write.table( file=fName, matrix( colNames, nrow=1 ), quote=FALSE,
      #             col.names=FALSE, row.names=FALSE,
      #             sep=" ", append=TRUE )
                   
      write.table( file=fName, newPars, quote=FALSE,
                   col.names=TRUE, row.names=FALSE,
                   sep=" ", append=TRUE )          #RF changed append=FALSE to append=TRUE
   
      #options( warn=-1 )                        # Turn off whining.
      #for ( k in 1:nrow(newPars) )
      #{
      #  isNumericVal <- !is.na( as.numeric( newPars[k,2] ) )  # Coerce non-numeric to NA.
      #  if ( isNumericVal )
      #    cat( file=fName, newPars[k,1]," ",newPars[k,2],"\n", append=TRUE, sep="" )
      #  else
      #    cat( file=fName, newPars[k,1]," ",dQuote(newPars[k,2]),"\n", append=TRUE, sep="" )        
      #}
      #options( warn=0 )                 # Turn on whining.
      
      cat( "\nMSG(.createBatchDesign) mseR parameter file ",fName," written...\n" )

      iBatch <- iBatch + 1           # Increment the batch job counter.
    }
  }
  
  # Bind the vectors that make up the design dataframe.
  result <- data.frame( parFile, blobName, scenarioLabel,
              mpLabel, dataName, methodName, ruleName,
              scCol, scLty, scLwd, scSym,
              mpCol, mpLty, mpLwd, mpSym, join,
              stringsAsFactors=FALSE )
  result
}    # END function .createBatchDesign


# .readParFile   (reads an ASCII file with 1 comment line, header, data frame)
# Purpose:      Reads an ASCII file: 1 comment, 1 header, space-delimited
#               data frame usually containing columns "parameter" and "value".
# Parameters:   parFile is a character string indicating the input file.
# Returns:      result, a data frame.
# Source:       A.R. Kronlund
.readParFile <- function( parFile="inputFile.par" )
{
  # Read the file and store as a dataframe.
  result <- read.table( file=parFile, as.is=TRUE, header=TRUE, skip=1,
                        quote="",sep=" " )
  result
}

.createList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as parameter and values in "value".

  result <- list()

  # Shut off whining, then coerce to numeric to let NA indicate non-numerics.
  options( warn=-1 )
  numericVal <- as.numeric( obj[,"value"] )
  options( warn=0 )

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
      listText <- paste( "result$",obj[i,"parameter"],"=",
                    obj[i,"value"],sep="" )
    # Value is character, build the parse string.
    else
      listText <- paste( "result$",obj[i,"parameter"],"=",
                  obj[i,"value"], sep="" )

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )
    
    # Evaluate the parsed string.
    eval( parse( text=listText ) )
  }
  result
}

# .createKeyList (Creates a list from the Batch File for .createBatchDesign):
# Purpose:      Function to convert the Batch File dataframe (loaded by the
#               function .readParFile) into a list that structures parameters
#               into "scenario" and "mp" nodes, each of each contains the
#               parameters for unique scenarios and management procedures,
#               respectively.
# Parameters:   obj is a dataframe read by .readParFile with columns "parameter"
#                 and "value".
# Returns:      result, a list with the batch control structure required to form
#               the desired cross-combinations of scenarios and procedures.
# Source:       A.R. Kronlund
.createKeyList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as first level of parameter and the
  # balance as the key, with values in "value".

  result <- list()

  options( warn=-1 )                        # Turn off whining.
  numericVal <- as.numeric( obj[,"value"] ) # Coerce non-numeric to NA.
  options( warn=0 )                         # Turn on whining.

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
    {
      parName <- obj[i,"parameter"]
      # Replace all the "$" with "&" in parameter name.
      parName <- gsub( "$",.subChar,parName, fixed=TRUE )
      # Replace ONLY the first "&" with "$" in parameter name, TWICE.
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      
      listText <- paste( "result$",parName,"=",
                    obj[i,"value"],sep="" )
    }
    # Value is character, build the parse string.
    else
    {
      parName <- obj[i,"parameter"]
      # Replace all the "$" with "&" in parameter name.
      parName <- gsub( "$",.subChar,parName, fixed=TRUE )
      # Replace ONLY the first "&" with "$" in parameter name.
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      parName <- sub( .subChar,"$",parName, fixed=TRUE )
      
      listText <- paste( "result$",parName,"=",
                    obj[i,"value"],sep="" )
    }

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )
    # Evaluate the parse string.
    eval( parse( text=listText ) )
  }
  result
}    # END function .createKeyList

# panLegend   (Place legend in plot region)
# Purpose:    Place a legend in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "legend" to be
#             passed.
# Parameters: x, y are the coordinates of the legend
#             legTxt is the text associated with the legend
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLegend <- function( x, y, legTxt, ... )
{
  # Allows legend to be placed at 0<x<1, 0<y<1.
  usr  <- par( "usr" )
  yLog <- par("ylog")
  xLog <- par("xlog")
  # Check for log-transformed axes and adjust usr commands as needed
    # note: when a log scale is in use, 
    #           usr gives limits in the form 10 ^ par("usr")
  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE) {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
     usr[3:4]<-10 ^ par("usr")[3:4]
     par( usr=c(0,1,0,1), ylog=FALSE )
   } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
     usr[1:2]<-10 ^ par("usr")[1:2]
     par( usr=c(0,1,0,1), xlog=FALSE )
   } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
     usr[1:4]<-10 ^ par("usr")[1:4]
     par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
   } 
  legend( x, y, legend=legTxt, ... )
  par( usr=usr )
  return( NULL )
}

# panLab      (Place text labels in plot region)
# Purpose:    Place a text label in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "text" to be
#             passed.
# Parameters: x, y are the coordinates of the label
#             txt is the text
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLab <- function( x, y, txt, ... )
{
  # Allows text to be placed in plot panel at 0<x<1, 0<y<1.
  usr <- par( "usr" )
  
  yLog <- par("ylog")
  xLog <- par("xlog")
  
  # Check for log-transformed axes and adjust usr commands as needed
  # note: when a log scale is in use, 
  #           usr gives limits in the form 10 ^ par("usr")

  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE)
  {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
    usr[3:4]<-10 ^ par("usr")[3:4]
    par( usr=c(0,1,0,1), ylog=FALSE )
  } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
    usr[1:2]<-10 ^ par("usr")[1:2]
    par( usr=c(0,1,0,1), xlog=FALSE )
  } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
    usr[1:4]<-10 ^ par("usr")[1:4]
    par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
  } 
  text( x, y, txt, ... )
  par( usr=usr )
  return( NULL )
}

# loadSim()
# Loads the nominated simulation blob into memory, so that plot functions
# run faster.
# inputs:   sim=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadSim <- function (sim=1)
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./project",full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  simList <- dirList[grep(pattern="sim",x=dirList)]
  folder <- simList[sim]

  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  assign( "blob",blob,pos=1 )
  cat("MSG (loadSim) Simulation ", folder, " loaded from ./project/\n", sep="" )
}


# calcDetFit()
# For troubleshooting bad model fits, by stepping through the
# empirical covariance matrix and calculating determinants for sub-matrices.
# Hopefully, the sub-matrices without full rank can be identified and 
# the culprits can be found.
# inputs:   fit=list returned by read.fit
# outputs:  det=named vector of deteterminants for submatrices missing
#               the named variable
calcDetFit <- function (fit)
{
  npar  <- fit$npar
  cor   <- fit$cor

  subDet <- numeric ( length = npar )

  for(j in 1:npar)
  {
    subMat <- cor[-j,-j]
    subDet[j] <- det(subMat)
  }
  return(subDet)
}

# saveSim()
# Saves a simulation in the project folder
# inputs:   blob=list object containing simulation information
#           name=blob name (provided by runSimEst())
#           path=full path to save folder
.saveSim <- function(blob, name, path)
{
  rDataFile <- paste ( name, ".RData", sep = "" )
  dataFilePath <- file.path ( path, rDataFile )
  save ( blob, file = dataFilePath )
  cat( "\nMSG (saveSim) Saved ",rDataFile,"in ./project/",name,"/.\n",sep="" )
}


# lisread()
# lisread: Function to read a list of data objects from a file.
# The initial characters "##" denote a comment line (ignored).
# The initial characters "# " denote a variable name.
# All other lines must contain scalars or vectors of numbers.
# Furthermore, all rows in a matrix must contain the same number of
# columns. Row and column vectors are not converted to matrices.
#
# fname  : File name.
# quiet  : If true, shut up about reporting progress.
# result : List object with components named in the file.

# Original functions courtesy of Jon Schnute.
# Modifications by A.R. Kronlund.
lisread <- function( fname,quiet=TRUE )
{
  lis2var <- function( x )
  {
    # lis2var: Makes global variables from the components of a list
    # x      : list object with named components.
    # result : global variables with names and contents extracted from x.

    namx <- names( x )
    nx <- length( namx )
    if (nx > 0) for (i in 1:nx)
    {
      if (namx[i] != "")
      {
        cmd <- paste( namx[i],"<<- x[[i]]" )
        eval( parse(text=cmd) )
      }
    }
    namx[namx != ""]
  }

  # numvecX functions:
  #
  # Function to convert a single string with white spaces into a numeric
  # vector. The string is parsed into separate components and converted
  # to char and numeric. A direct conversion to numeric fails.

  numvec <- function( x )
  {
    # Deprecated.
    xp <- parse( text=x )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec2 <- function( x )
  {
    # Patch required for S6.0 where white=T option is defunct in parse.
    # Deprecated:  xp <- parse( text=x,white=T )
    # ARK 30-Oct-03: R text connections get used up, must open/close.
    tc <- textConnection( x )
    xp <- scan( tc )
    close( tc )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec3 <- function( x,quiet )
  {
    # ARK 12-Jan-06: Patch to read characters because Rashmi asked nicely.
    # This is a largely untested hack, no expressed or implied warrantee.
    tc <- textConnection( x )
    xp <- scan( tc, what="character",quiet=quiet )
    close( tc )
    xc <- as.character( xp )
    if ( !all(is.na(as.numeric(xc))) )
      xc <- as.numeric( xc )
    xc
  }

  #------------------------------------------------------------------#

  file <- scan( fname, what=character(), sep="\n", quiet=quiet )

  f2 <- file[ regexpr("##",file)==-1 ]           # remove comments
  nf2 <- length( f2 )                            # number of lines
  llab <- regexpr( "#",f2 )==1                   # identifies label lines
  vlab <- substring( f2[llab],3 )                # variable labels

  # ARK 30-Oct-03 R does not coerce logical to character for grep.
  ilab <- grep( "TRUE",as.character(llab) )      # label indices

  nvar <- length( vlab )                         # number of variables

  if( nvar==1 )
    nrow <- c( nf2 + 1 ) - ilab - 1
  else
    nrow <- c( ilab[2:nvar],nf2+1 ) - ilab - 1     # number of rows in var i
  zout <- list( NULL )

  for ( i in 1:nvar )
  {
    i1 <- ilab[i] + 1
    i2 <- i1 + nrow[i] - 1                       # range of lines for var i
    zstr <- paste(f2[i1:i2],collapse=" ")
#    zvec <- numvec2(zstr)                        # numeric vector
    zvec <- numvec3(zstr,quiet)                  # numeric or character vector

    nz <- length(zvec)
    zrow <- nrow[i]
    zcol <- nz / zrow                            # dimensions
    if ( (zrow>1) & (zcol>1) )                   # a true matrix
      zvec <- matrix( zvec,nrow=zrow,ncol=zcol,byrow=T )
    zout[[i]] <- zvec
    if ( !quiet )
      cat( "vlab = ", vlab[i], "\n" )
  }
  names(zout) <- vlab
  zout
}


# Write element i of list object x to the active file. Note that any
# elements must be numeric (as this is the only thing ADMB can read in),
# also sub-list structure is removed and elements of lists are written as 
# vectors, probably to be read in to ADMB as ragged arrays
# usage:
# cat ( "## writeAMDB output file, created ", format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", sep = "", file = activeFile )
# lapply ( X = seq_along(x), FUN = writeADMB, x, activeFile)
writeADMB <- function( i, x, activeFile ) {
  # Write element name
  cat( paste( "# ", names(x)[[i]], "\n", sep = "" ), file=activeFile, append=TRUE )
  # Write element value
  if( is.matrix(x[[i]]) ) {
    # Matrices
    write.table( x[[i]], file=activeFile, row.names=FALSE, col.names=FALSE, append=TRUE )
  } else if( is.list(x[[i]]) ) {
    # Lists
    lapply( x[[i]], "\n", FUN=cat, file=activeFile, append=TRUE )
  } else {
    # Vectors/Scalars
    cat( x[[i]], "\n", file=activeFile, append=TRUE )  
  }
}  # end writeADMB

# Read ragged arrays from ADMB rep object
readRagged <- function( x, var, ncols ) {
  outList <- list()
  for( i in 1:length(ncols) ) {
    outList[[i]] <- numeric(ncols[i])
    for( j in 1:ncols[i] ) {
      varName <- paste( var, "_", i, j, sep=""  )
      k <- which( names(x)==varName )
      outList[[i]][j] <- x[[k]]
    }
  }
  outList
}   # end readRagged()
read.admb <-
function(ifile)
{ 
  ret=read.fit(ifile)
  
  fn=paste(ifile,'.rep', sep='')
  A=read.rep(fn)
  A$fit=ret
  
  pfn=paste(ifile,'.psv',sep='')
  if(file.exists(pfn))
    A$post.samp=read.psv(pfn)
  
  return(A)
}

read.fit <-
function(ifile)
{
  # __Example:             
  # file <-("~/admb/simple")
  # A <- reptoRlist(file)
  # Note there is no extension on the file name.
  
  ## The following is a contribution from:
  ## Anders Nielsen that reads the par & cor files.
  ret<-list() 
  parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
   what='', n=16, quiet=TRUE)[c(6,11,16)]) 
  ret$nopar<-as.integer(parfile[1]) 
  ret$nlogl<-parfile[2] 
  ret$maxgrad<-parfile[3] 
  file<-paste(ifile,'.cor', sep='') 
  lin<-readLines(file) 
  ret$npar<-length(lin)-2 
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
  ret$names<-unlist(lapply(sublin,function(x)x[2])) 
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
  ret$cor<-matrix(NA, ret$npar, ret$npar) 
  corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
  ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
  ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  return(ret)
}

read.rep <- 
function(fn)
{
  # The following reads a report file
  # Then the 'A' object contains a list structure
  # with all the elemements in the report file.
  # In the REPORT_SECTION of the AMDB template use 
  # the following format to output objects:
  #   report<<"object \n"<<object<<endl;
  #
  # The part in quotations becomes the list name.
  # Created By Steven Martell
  options(warn=-1)  #Suppress the NA message in the coercion to double
  
  
  ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
    dum=NA
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
    if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))

    if(is.numeric(dum))#Logical test to ensure dealing with numbers
    {
      A[[vnam[i]]]=dum
    }
  }
  options(warn=0)
  
  return(A)
}

read.psv <-
function(fn, nsamples=10000)
{
  #This function reads the binary output from ADMB
  #-mcsave command line option.
  #fn = paste(ifile,'.psv',sep='')
  filen <- file(fn, "rb")
  nopar <- readBin(filen, what = integer())
  mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
  mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
  close(filen)
  return(mcmc)
}

gletter <-
function(i=1)
{
  usr=par("usr"); inset.x=0.05*(usr[2]-usr[1]); inset.y=0.05*(usr[4]-usr[3])
  text(usr[1]+inset.x,usr[4]-inset.y,paste("(",letters[i],")",sep=""),cex=1.,font=1)
}

makeLHRfromList <- function(  levels = list(  Uhist = c("c(0.2,2,1)","c(1,1,1)"),
                                              initYear = c(1984,2003),
                                              nS = c(4,7,10),
                                              initDep = c(0.4,0.7,1.0),
                                              nDiff = c(0,2,4) ),
                              nPoints = 1,
                              seed = 123 )
{
  # Creates a Latin Hyper Rectangle experimental design from a list of factor levels
  # inputs:     levels = list of factor levels
  #             nPoints = number of points to sample from the LHR design

  set.seed(seed)

  # First, create an array to hold the design
  # First, we need the number of factors and their levels to make dimensions
  nLevels     <- lapply ( X = levels, FUN = length )
  nLevels     <- unlist( nLevels )
  nFactors    <- length( nLevels )

  # Create the array, using the factor levels as dimnames
  LHRdesign   <- array( NA, dim = nLevels, dimnames = levels )

  # Choose the entries in the array as the size of the largest dimension
  maxEntry  <- max(nLevels)

  # I think adding the entry dimension indices mod maxEntry will
  # populate the matrix and preserve the latin property. To do this
  # we gotta expand.grid for all possible entries
  dimIndices <- vector(mode = "list", length = nFactors )
  for( lIdx in 1:nFactors )
  {
    dimIndices[[lIdx]] <- 1:nLevels[lIdx]
  }
  entryIndices <- expand.grid( dimIndices )

  # Loop over the entryIndices data.frame and pull the rows
  for( rIdx in 1:nrow(entryIndices) )
  {
    entryIdx <- as.numeric(entryIndices[rIdx,])
    LHRdesign[matrix(entryIdx,nrow=1)] <- (sum(entryIdx) %% maxEntry)
  }

  # Now randomly permute the dimensions
  idx <- lapply(dim(LHRdesign), sample)
  LHRdesign[] <- do.call(`[`, c( list(LHRdesign), idx) )

  # Now sample design space for treatments
  points      <- sample( x = 0:(maxEntry-1), size = nPoints )
  treatments  <- which( LHRdesign %in% points, arr.ind = F )
  treatments  <- arrayInd(  ind = treatments, .dim = dim(LHRdesign), 
                            .dimnames = dimnames(LHRdesign), useNames = T)

  out <- list(  levels = levels,
                designArray = LHRdesign,
                treatments =  treatments
              )

  out

}
