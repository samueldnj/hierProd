// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// msProd.cpp
// 
// A multi-stock surplus production (Schaefer) state-space model with 
// joint prior distributions applied to catchability (q) 
// and productivity (U_msy), and a process error year-effect component (eps_t) 
// shared between stocks
// 
// Author: Samuel Johnson
// Date: 1 November, 2016
//
// Purpose: This is the assessment model in a simulation study of multi-stock
// robin hood methods.
// 
// Updates:
//    May 13, 2017: Added multiple surveys
// 
// ><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>

#include <TMB.hpp>       // Links in the TMB libraries
#include <iostream>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// invLogit
template<class Type>
Type invLogit(Type x, Type scale, Type trans){
  return scale/(Type(1.0) + exp(-Type(1.0)*x)) - trans;
}

// invLogit
template<class Type>
Type square(Type x)
{
  return pow(x,2);
}



// objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*data section*/
  // Data Structures
  DATA_ARRAY(It);               // CPUE data
  DATA_ARRAY(Ct);               // Catch data
  // Model dimensions
  int nO = It.dim(0);           // No. of surveys
  int nS = It.dim(1);           // No. of species
  int nT = It.dim(2);           // No of time steps

  // Model switches
  DATA_INTEGER(kappaPriorCode); // 0 => OFF, 1 => IG on kappa2
  DATA_INTEGER(SigmaPriorCode); // 0 => IG on diagonal element, 1 => IW on cov matrix
  DATA_INTEGER(tauqPriorCode);  // 0 => IG on tauq2, 1 => normal
  DATA_INTEGER(sigUPriorCode);  // 0 => IG on sigU2, 1 => normal
  DATA_INTEGER(condMLEq);       // 0 => q leading par, 1 => q concentrated
  DATA_INTEGER(lnqPriorCode);   // 0 => hyperprior, 1 => multilevel
  DATA_INTEGER(lnUPriorCode);   // 0 => hyperprior, 1 => multilevel 
  DATA_INTEGER(BPriorCode);     // 0 => normal, 1 => Jeffreys 
  DATA_IVECTOR(initT);          // first year of reconstructed catch hist
  DATA_IVECTOR(initBioCode);    // initial biomass at 0 => unfished, 1=> fished
  DATA_SCALAR(posPenFactor);    // Positive-penalty multiplication factor


  /*parameter section*/
  // Leading Parameters
  PARAMETER_VECTOR(lnBmsy);             // Biomass at MSY
  PARAMETER_VECTOR(lnUmsy);             // Optimal exploitation rate            
  PARAMETER_ARRAY(lnq_os);              // Survey-species catchability
  PARAMETER_VECTOR(lntau2_o);           // survey obs error var
  PARAMETER_VECTOR(lnBinit);            // Non-equilibrium initial biomass
  // Priors
  PARAMETER_VECTOR(lnqbar_o);           // survey mean catchability
  PARAMETER_VECTOR(lntauq_o);           // survey catchability sd 
  PARAMETER(mq);                        // hyperprior mean catchability
  PARAMETER(sq);                        // hyperprior catchability sd
  PARAMETER(lnUmsybar);                 // shared prior mean Umsy
  PARAMETER(lnsigUmsy);                 // shared prior Umsy sd
  PARAMETER(mUmsy);                     // hyperprior mean Umsy
  PARAMETER(sUmsy);                     // hyperprior Umsy sd
  PARAMETER_VECTOR(mBmsy);              // prior mean eqbm biomass
  PARAMETER_VECTOR(sBmsy);              // prior eqbm biomass var
  PARAMETER_VECTOR(tau2IGa);            // Inverse Gamma Prior parameters for tau2 prior
  PARAMETER_VECTOR(tau2IGb);            // Inverse Gamma Prior parameters for tau2 prior
  PARAMETER_VECTOR(tauq2Prior);         // Hyperparameters for tauq2 prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(sigU2Prior);         // Hyperparameters for sigU2 prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(kappa2IG);           // Inverse Gamma Prior parameters for kappa2 prior
  PARAMETER_VECTOR(Sigma2IG);           // Inverse Gamma Prior parameters for Sigma2 prior
  PARAMETER_MATRIX(wishScale);          // IW scale matrix for Sigma prior
  PARAMETER(nu);                        // IW degrees of freedom for Sigma prior    
  PARAMETER(deltat);                    // Fractional time step used in pop dynamics to reduce chaotic behaviour
  // Random Effects
  PARAMETER(lnkappa2);                  // year effect uncorrelated variance
  PARAMETER(lnSigmaDiag);               // Species effect cov matrix diag
  PARAMETER_VECTOR(eps_t);              // year effect  
  PARAMETER_ARRAY(zeta_st);             // species effect
  PARAMETER_VECTOR(SigmaDiagMult);      // Sigma diagonal mults
  PARAMETER_VECTOR(logitSigmaOffDiag);  // Species effect corr chol factor off diag  
  PARAMETER(logit_gammaYr);             // AR1 auto-corr on year effect (eps)
  

  // State variables
  array<Type>       Bt(nS,nT);
  array<Type>       lnBt(nS,nT);
  // Leading parameters
  vector<Type>      Bmsy    = exp(lnBmsy);
  vector<Type>      Umsy    = exp(lnUmsy);
  vector<Type>      tau2_o  = exp(lntau2_o);
  vector<Type>      Binit   = exp(lnBinit);
  array<Type>       lnqhat_os(nO,nS);
  vector<Type>      tau2hat_o(nO);
  
  // Prior hyperpars
  vector<Type>      qbar_o  = exp(lnqbar_o);
  vector<Type>      tauq_o  = exp(lntauq_o);
  vector<Type>      tauq2_o = exp(Type(2) * lntauq_o);
  Type              Umsybar = exp(lnUmsybar);
  Type              sigUmsy2= exp(Type(2) * lnsigUmsy);
  // Random Effects
  Type              kappa2  = exp(lnkappa2);
  Type              gammaYr = Type(1.96) / (Type(1.) + exp(Type(-2.)*logit_gammaYr) ) - Type(0.98);
  vector<Type>      SigmaOffDiag(nS*(nS-1)/2);
                    SigmaOffDiag  = Type(1.96) / (Type(1.) + exp(Type(-2.)*logitSigmaOffDiag) ) - Type(0.98);     
  vector<Type>      SigmaDiag = exp(lnSigmaDiag) * SigmaDiagMult;
  matrix<Type>      SigmaCorr(nS,nS);
  matrix<Type>      SigmaD(nS,nS);
  matrix<Type>      Sigma(nS,nS);
  vector<Type>      omegat(nT);
  // Scalars
  Type              nll     = 0.0;   // objective function (neg log likelihood)
  Type              nllRE   = 0.0;   // process error likelihood
  Type              pospen  = 0.0;   // posfun penalty
  // Derived variables
  vector<Type>      MSY = Bmsy * Umsy;
  vector<Type>      lnMSY = log(MSY);
  array<Type>       Ut(nS,nT);
  vector<Type>      DnT(nS);
  vector<Type>      lnDnT(nS);
  vector<Type>      BnT(nS);
  vector<Type>      lnBnT(nS);
  array<Type>       U_Umsy(nS,nT);
  array<Type>       lnU_Umsy(nS,nT);
  vector<Type>      lnU_UmsyT(nS);

  // Procedure Section //
  // First create the correlation matrix for the species effects
  UNSTRUCTURED_CORR_t<Type> specEffCorr(SigmaOffDiag);
  SigmaCorr = specEffCorr.cov();
  // Standard deviation component
  SigmaD.fill(0.0);
  SigmaD.diagonal() = sqrt(SigmaDiag);
  // Now combine
  Sigma = (SigmaD * SigmaCorr);
  Sigma *= SigmaD;

  // Generate Population Dynamics
  // initialise first year effect at 0, then loop to fill in
  omegat(0) = 0;
  for( int t=1; t<nT; t++ )
  {
    omegat(t) = gammaYr * omegat(t-1) + (1 - gammaYr) * eps_t(t-1);
  }
  
  // Initialise biomass and log-biomass
  Bt.fill(-1);
  lnBt.fill(0.0);
  Type nSteps = 1 / deltat;
  // Now loop over species, reconstruct history from initT
  for( int s = 0; s < nS; s++ )
  {
    // initialise population, if initBioCode(s)==0 this will be eqbm
    if( initBioCode(s) == 0 ) Bt(s,initT(s)) = Type(2) * Bmsy(s);
    if( initBioCode(s) == 1 ) Bt(s,initT(s)) = Binit(s);
    lnBt(s,initT(s)) = log(Bt(s,initT(s)));
    for( int t = initT(s)+1; t < nT; t++ )
    {
      Type tmpBt = Bt(s,t-1);
      for( int dt = 0; dt < nSteps; dt ++ )
      {
        // Compute the deltat step of biomass (assuming catch is evenly distributed across the year)
        Type tmpBdt = 0;
        tmpBdt =  tmpBt + deltat * Umsy(s) * tmpBt * (Type(2.0) - tmpBt / Bmsy(s) ) - deltat*Ct(s,t-1);
        tmpBdt *= exp(deltat * (omegat(t) + zeta_st(s,t-1)));
        // Now update tmpBt
        tmpBt = posfun(tmpBdt, Type(1e-3), pospen);
      }
      Bt(s,t) = tmpBt;
      lnBt(s,t) = log(Bt(s,t));
    }

  }
  // Loop again to add species and year effects
  for( int t=1; t<nT; t++ )
  {
    // Add year effect contribution to objective function
    if( kappaPriorCode == 1 )
      nllRE -= dnorm( eps_t(t-1),Type(0.),sqrt(kappa2),1);
    // Add correlated species effects contribution to likelihood
    for( int s = 0; s < nS; s++ )
    {
      nllRE += 0.5*(lnSigmaDiag + pow(zeta_st(s,t-1),2)/SigmaDiag(s) );
    }
        
  }
  // add REs to joint nll
  nllRE += posPenFactor*pospen;
  nll += nllRE;

  // Concentrate species specific obs error likelihood?
  // Initialise arrays for concentrating conditional MLE qhat
  array<Type>   validObs(nO,nS);
  vector<Type>  totObs_o(nO);
  array<Type>   qhat_os(nO,nS);
  array<Type>   z_ost(nO,nS,nT);
  array<Type>   zSum_os(nO,nS);
  array<Type>   SS_os(nO,nS);
  vector<Type>  totSS_o(nO);
  // Fill with 0s
  Type nllObs = 0.0;
  validObs.fill(1e-6);
  totObs_o.fill(0);
  zSum_os.fill(0.0);
  z_ost.fill(0.0);
  qhat_os.fill(-1.0);
  SS_os.fill(0.);
  totSS_o.fill(0.);


  // Compute observation likelihood
  // Loop over surveys
  for( int o=0; o < nO; o++ )
  {
    // species
    for( int s = 0; s < nS; s++ )
    {
      // time
      for( int t = initT(s); t < nT; t++ )
      {
        // only add a contribution if the data exists (Iost < 0 is missing)
        if( ( It(o,s,t) > 0. ) ) 
        {
          validObs(o,s) += int(1);
          z_ost(o,s,t) = log( It( o, s, t ) ) - log( Bt( s, t ) );
          zSum_os(o,s) += z_ost(o,s,t);
        }       
      }
      if(condMLEq == 1)
      {
        // compute conditional MLE q from observation
        // SS model
        if( nS == 1 | lnqPriorCode == 0) 
          lnqhat_os(o,s) = zSum_os(o,s) / validObs(o,s);

        // MS model
        if( nS > 1 & lnqPriorCode == 1 ) 
          lnqhat_os(o,s) = ( zSum_os(o,s) / tau2_o(o) + lnqbar_o(o)/tauq2_o(o) ) / ( validObs(o,s) / tau2_o(o) + 1 / tauq2_o(o) );  
      } else
        lnqhat_os(o,s) = lnq_os(o,s);

      // Exponentiate
      qhat_os(o,s) = exp(lnqhat_os(o,s));

      // Subtract lnq_os from the resids
      for( int t = initT(s); t < nT; t++ )
        if( (It(o,s,t) > 0.0) )
          z_ost(o,s,t) -= lnqhat_os(o,s);

      // Calculate sum of squared resids
      for( int t = initT(s); t < nT; t++ ) 
        SS_os(o,s) += square(z_ost(o,s,t));
      
      // Add to likelihood
      nllObs += 0.5 * ( validObs(o,s)*lntau2_o(o) + SS_os(o,s)/tau2_o(o));

      // Add valid Obs and SS to a total for each survey
      totObs_o(o) += validObs(o,s);
      totSS_o(o) += SS_os(o,s);
    }
    tau2hat_o(o) = totSS_o(o) / totObs_o(o);
  }
  nll += nllObs;

  // Add priors
  Type nllBprior = 0.0;
  Type nllqPrior = 0.0;
  Type nllUprior = 0.0;

  // eqbm biomass
  for (int s=0; s<nS; s++ )
  {
    // Normal prior
    if(BPriorCode == 0)
    {
      // Bmsy
      nllBprior -= dnorm(Bmsy(s),mBmsy(s), sBmsy(s), 1); 

      // Initial biomass
      if(initBioCode(s) == 1) 
        nllBprior -= dnorm( Binit(s), mBmsy(s)/2, sBmsy(s)/2, 1);  

    }
    // Jeffreys prior
    if( BPriorCode == 1 )
      nllBprior += lnBmsy(s) + lnBinit(s);
    
  } // End biomass prior

  // multispecies shared priors
  // If not using shared priors, use the hyperpriors for all
  // species specific parameters
  if (nS > 1)
  {
    // 1st level priors
    for( int s = 0; s < nS; s++ )
    {
      // In here, loop over surveys
      for( int o = 0; o < nO; o++ )
      {
        // catchability
        // Shared Prior
        if( lnqPriorCode == 1 & condMLEq == 0 )
          nllqPrior -= dnorm( lnqhat_os(o,s), lnqbar_o(o), tauq_o(o), 1);

        // No shared prior (uses the same prior as the SS model)
        if( lnqPriorCode == 0 )
          nllqPrior += 0.5 * pow((qhat_os(0,s) - mq)/sq,2);

      }
      // productivity
      // Shared Prior
      if( lnUPriorCode == 1 )
        nllUprior -= dnorm( lnUmsy(s), lnUmsybar, sqrt(sigUmsy2), 1);

      // No shared prior (uses SS model prior)
      if( lnUPriorCode == 0 )
        nllUprior += pow( (Umsy(s) - mUmsy)/ sUmsy, 2);

    }  
    // Hyperpriors
    // catchability
    if( lnqPriorCode == 1 )
      for( int o = 0; o < nO; o++ ) 
        nllqPrior -= dnorm( qbar_o(o), mq, sq, 1); 

    // productivity
    if( lnUPriorCode == 1 )
      nllUprior -= dnorm(Umsybar, mUmsy, sUmsy, 1);
    
  } // End multispecies shared priors    
  
  // Now for single species model
  if( nS == 1 ) 
  { 
    // catchability
    for (int o=0; o<nO; o++)
      nllqPrior -= dnorm( qhat_os(o,0), mq, sq, 1); 
    
    // productivity
    nllUprior -= dnorm( Umsy(0), mUmsy, sUmsy, 1);
  }
  // Add all priors
  nll += nllBprior +  nllqPrior + nllUprior;
  
  // Variance IG priors
  // Obs error var
  Type nllVarPrior = 0.;
  Type nllSigPrior = 0.;
  for( int o = 0; o < nO; o++ )
    nllVarPrior += (tau2IGa(o)+Type(1))*lntau2_o(o)+tau2IGb(o)/tau2_o(o);  

  // year effect deviations var
  if( kappaPriorCode == 1 )  
    nllVarPrior += (kappa2IG(0)+Type(1))*lnkappa2 + kappa2IG(1)/kappa2;

  // Now multispecies priors
  if (nS > 1)
  {
    // shared q prior variance
    if( tauqPriorCode == 0 )
    {
      for( int o = 0; o < nO; o++)
      {
        nllVarPrior += (tauq2Prior(0)+Type(1))*Type(2)*lntauq_o(o)+tauq2Prior(1)/tauq2_o(o);
      }
    }
    if( tauqPriorCode == 1 )
    {
      for( int o = 0; o < nO; o++)
      {
        nllVarPrior += Type(0.5) * pow( tauq2_o(o) - tauq2Prior(0), 2) / tauq2Prior(1);
      }
    }
    // shared U prior variance
    // IG
    if( sigUPriorCode == 0 )
    {
      nllVarPrior += (sigU2Prior(0)+Type(1))*Type(2.)*lnsigUmsy+sigU2Prior(1)/sigUmsy2;
    }
    // Normal
    if( sigUPriorCode == 1 )
    {
      nllVarPrior += Type(0.5) * pow( sigUmsy2 - sigU2Prior(0), 2) / sigU2Prior(1);
    }
    
    if( SigmaPriorCode == 1 ) // Apply IW prior to Sigma matrix
    {
      matrix<Type> traceMat = wishScale * Sigma.inverse();
      Type trace = 0.0;
      for (int s=0;s<nS;s++) trace += traceMat(s,s);
      nllSigPrior += Type(0.5) *( (nu + nS + 1) * atomic::logdet(Sigma) + trace);
    }
  }

  // Apply Sigma Prior
  if( SigmaPriorCode == 0 ) // Apply IG to estimated SigmaDiag element
    nllSigPrior += (Sigma2IG(0)+Type(1))*lnSigmaDiag+Sigma2IG(1)/exp(lnSigmaDiag);

  nll += nllVarPrior + nllSigPrior;

  // Derive some output variables
  Ut      = Ct / Bt;
  DnT     = Bt.col(nT-1)/Bmsy/2;
  lnDnT   = log(DnT);
  lnBnT   = log(Bt.col(nT-1));
  for( int t = 0; t < nT; t ++ )
  {
    U_Umsy.col(t) = Ut.col(t) / Umsy;
    lnU_Umsy.col(t) = log(U_Umsy.col(t));
  }
  lnU_UmsyT = lnU_Umsy.col(nT-1);

  // Reporting Section //
  // Variables we want SEs for
  ADREPORT(lnBt);
  ADREPORT(lnqhat_os);
  ADREPORT(lnMSY);
  ADREPORT(lnBinit);
  ADREPORT(lntau2_o);
  ADREPORT(lnkappa2);
  ADREPORT(lnDnT);
  ADREPORT(lnBnT);
  ADREPORT(lnU_UmsyT);

  
  // Everything else //
  REPORT(Bt);
  REPORT(It);
  REPORT(Ut);
  REPORT(Ct);
  REPORT(eps_t);
  REPORT(omegat);
  REPORT(Binit);
  REPORT(DnT);
  REPORT(U_Umsy);
  REPORT(qhat_os);
  REPORT(MSY);
  REPORT(Bmsy);
  REPORT(Umsy);
  REPORT(tau2_o);
  REPORT(kappa2);
  REPORT(nT);
  REPORT(nO);
  REPORT(nS);
  REPORT(zSum_os);
  REPORT(z_ost);
  REPORT(validObs);
  REPORT(zeta_st);
  REPORT(Sigma);
  REPORT(SigmaDiag);
  if (nS > 1)
  {
    REPORT(SigmaCorr);
    REPORT(Umsybar);
    REPORT(sigUmsy2);
    REPORT(qbar_o);
    REPORT(tauq2_o);
  }
  REPORT(gammaYr);
  REPORT(nll);
  REPORT(nllRE);
  REPORT(nllObs);
  REPORT(nllqPrior);
  REPORT(nllUprior);
  REPORT(nllBprior);
  REPORT(nllSigPrior);
  REPORT(nllVarPrior);
  REPORT(pospen);
  
  return nll;
}
