/// \file PReMiuM.cpp
/// \author David Hastie
/// \brief Main file for running PReMiuM

/// \note (C) Copyright David Hastie and Silvia Liverani, 2012.

/// PReMiuM++ is free software; you can redistribute it and/or modify it under the
/// terms of the GNU Lesser General Public License as published by the Free Software
/// Foundation; either version 3 of the License, or (at your option) any later
/// version.

/// PReMiuM++ is distributed in the hope that it will be useful, but WITHOUT ANY
/// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
/// PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

/// You should have received a copy of the GNU Lesser General Public License
/// along with PReMiuM++ in the documentation directory. If not, see
/// <http://www.gnu.org/licenses/>.

/// The external linear algebra library Eigen, parts of which are included  in the
/// lib directory is released under the LGPL3+ licence. See comments in file headers
/// for details.

/// The Boost C++ header library, parts of which are included in the  lib directory
/// is released under the Boost Software Licence, Version 1.0, a copy of which is
/// included in the documentation directory.

/// Version 3.1.3 edited by Rob Johnson, March 2017


// Standard includes
#include<cmath>
#include<cstdio>
#include<iostream>
#include<numeric>
#include<vector>
#include<string>
#include<cstdlib>
#include<ctime>
#include<sstream>

// Custom includes
#include<PReMiuMProposals.h>

#include<MCMC/sampler.h>
#include<MCMC/model.h>
#include<MCMC/proposal.h>
#include<PReMiuMOptions.h>
#include<PReMiuMModel.h>
#include<PReMiuMData.h>

#include<PReMiuMIO.h>

#include "include/postProcess.h"
//#include <stdlib.h> srand(seed)
#include<Rcpp.h>



#include<unistd.h>
#include<stdio.h>
using std::vector;
using std::ostringstream;
using std::time;
using std::string;
using std::cout;

RcppExport SEXP profRegr(SEXP inputString) {

  string inputStr = Rcpp::as<string>(inputString);

  /* ---------- Start the timer ------------------*/
  time_t beginTime,currTime;
  beginTime = time(NULL);
  /* -----------Process the command line ---------*/
  pReMiuMOptions options = processCommandLine(inputStr);
  /* ---------- Set up the sampler object--------*/
  // Initialise the sampler object
  mcmcSampler<pReMiuMParams,pReMiuMOptions,
              pReMiuMPropParams,pReMiuMData> pReMiuMSampler;

  // Set the options
  pReMiuMSampler.options(options);
  // Set the model
  pReMiuMSampler.model(&importPReMiuMData,&initialisePReMiuM,
                       &pReMiuMLogPost,true);
  // Set the missing data function
  pReMiuMSampler.updateMissingDataFn(&updateMissingPReMiuMData);

  // Add the function for writing output
  pReMiuMSampler.userOutputFn(&writePReMiuMOutput);

  // Seed the random number generator
  pReMiuMSampler.seedGenerator(options.seed());
  //pReMiuMSampler.seedGeneratorMult(options.seedsMult());

  // Set the sampler specific variables
  pReMiuMSampler.nSweeps(options.nSweeps());
  pReMiuMSampler.nBurn(options.nBurn());
  pReMiuMSampler.nFilter(options.nFilter());
  pReMiuMSampler.nProgress(options.nProgress());
  pReMiuMSampler.reportBurnIn(options.reportBurnIn());

  /* ---------- Read in the data -------- */
  pReMiuMSampler.model().dataset().outcomeType(options.outcomeType());
  pReMiuMSampler.model().dataset().kernelType(options.kernelType());//AR
  pReMiuMSampler.model().dataset().covariateType(options.covariateType());
  pReMiuMSampler.model().dataset().includeCAR(options.includeCAR());
  pReMiuMSampler.importData(options.inFileName(),options.predictFileName(),options.neighbourFileName());
  pReMiuMData dataset = pReMiuMSampler.model().dataset();

  /* ---------- Add the proposals -------- */
  // Set the proposal parameters
  pReMiuMPropParams proposalParams(options.nSweeps(),dataset.nCovariates(),
                                   dataset.nFixedEffects(),dataset.nFixedEffects_mix(),dataset.nCategoriesY(),dataset.kernelType());
  pReMiuMSampler.proposalParams(proposalParams);

  // The gibbs update for the active V
  pReMiuMSampler.addProposal("gibbsForVActive",1.0,1,1,&gibbsForVActive);

  if(options.covariateType().compare("Discrete")==0){
    // For discrete X data we do a mixture of Categorical and ordinal updates
    //  Update for the active phi parameters
    pReMiuMSampler.addProposal("updateForPhiActive",1.0,1,1,&updateForPhiActive);

  }else if(options.covariateType().compare("Normal")==0){
    // Need to add the proposals for the normal case
    // Update for the active mu parameters
    if (options.useNormInvWishPrior()){
      // If Inverse Normal Inverse Wishart Prior is used
      pReMiuMSampler.addProposal("gibbsForMuActiveNIWP",1.0,1,1,&gibbsForMuActiveNIWP);
    }else{
      // If independant Normal and Inverse Wishart priors are used
      pReMiuMSampler.addProposal("gibbsForMuActive",1.0,1,1,&gibbsForMuActive);
    }

    // Update for the active Sigma parameters
    pReMiuMSampler.addProposal("gibbsForTauActive",1.0,1,1,&gibbsForTauActive);

  }else if(options.covariateType().compare("Mixed")==0){
    // For discrete X data we do a mixture of Categorical and ordinal updates
    //  Update for the active phi parameters
    pReMiuMSampler.addProposal("updateForPhiActive",1.0,1,1,&updateForPhiActive);

    // Need to add the proposals for the normal case
    // Update for the active mu parameters
    if (options.useNormInvWishPrior()){
      // If Inverse Normal Inverse Wishart Prior is used
      pReMiuMSampler.addProposal("gibbsForMuActiveNIWP",1.0,1,1,&gibbsForMuActiveNIWP);
    }else{
      // If independant Normal and Inverse Wishart priors are used
      pReMiuMSampler.addProposal("gibbsForMuActive",1.0,1,1,&gibbsForMuActive);
    }

    // Update for the active Sigma parameters
    pReMiuMSampler.addProposal("gibbsForTauActive",1.0,1,1,&gibbsForTauActive);

  }

  if(options.varSelectType().compare("None")!=0){
    // Add the variable selection moves
    unsigned int firstSweep;
    firstSweep=1+(unsigned int)(options.nBurn()/10);
    if(options.varSelectType().compare("Continuous")!=0){
      // Gibbs update for gamma
      pReMiuMSampler.addProposal("gibbsForGammaActive",1.0,1,firstSweep,&gibbsForGammaActive);
    }
  }

  if(options.includeResponse()){
    // The Metropolis Hastings update for the active theta
    if(options.outcomeType().compare("Longitudinal")!=0 && options.outcomeType().compare("LME")!=0) //AR
      pReMiuMSampler.addProposal("metropolisHastingsForThetaActive",1.0,1,1,&metropolisHastingsForThetaActive);

    // Adaptive MH for beta
    if(dataset.nFixedEffects()>0 || dataset.nFixedEffects_mix()>0){
      if(options.outcomeType().compare("LME")!=0){
        pReMiuMSampler.addProposal("metropolisHastingsForBeta",1.0,1,1,&metropolisHastingsForBeta);
      }else{
        pReMiuMSampler.addProposal("gibbsForBeta",1.0,1,1,&GibbsForBeta);
      }
    }

    if(options.responseExtraVar()){
      // Adaptive MH for lambda
      pReMiuMSampler.addProposal("metropolisHastingsForLambda",1.0,1,1,&metropolisHastingsForLambda);

      // Gibbs for tauEpsilon
      pReMiuMSampler.addProposal("gibbsForTauEpsilon",1.0,1,1,&gibbsForTauEpsilon);
    }

    //if spatial random term
    if (options.includeCAR()){
      //Adaptive rejection sampling for uCAR
      pReMiuMSampler.addProposal("gibbsforUCAR", 1.0,1,1,&gibbsForUCAR);

      //Gibbs for TauCAR
      pReMiuMSampler.addProposal("gibbsForTauCAR", 1.0,1,1,&gibbsForTauCAR);
    }

    if(options.outcomeType().compare("Survival")==0){
      if (!options.weibullFixedShape()) {
        // Gibbs for shape parameter (cluster specific) of Weibull for survival response model
        pReMiuMSampler.addProposal("gibbsForNu",1.0,1,1,&gibbsForNu);
      }
    }
    //RJ add sampler for L & MVN
    if(options.outcomeType().compare("Longitudinal")==0){
      pReMiuMSampler.addProposal("metropolisHastingsForL",1.0,1,1,&metropolisHastingsForL);
    }
    if(options.outcomeType().compare("MVN")==0){
      //pReMiuMSampler.addProposal("gibbsForMVNMuActive",1.0,1,1,&gibbsForMVNMuActive);
      // Update for the active Sigma parameters
      pReMiuMSampler.addProposal("gibbsForMVNTauActive",1.0,1,1,&gibbsForMVNTauActive);
    }
    //AR
    if(options.outcomeType().compare("LME")==0){
      // Update covariance matrices and random effects
      pReMiuMSampler.addProposal("gibbsForForCovRELMEActive",1.0,1,1,&gibbsForForCovRELMEActive);
      // Update variance error
      pReMiuMSampler.addProposal("gibbsForSigmaEpsilonLME",1.0,1,1,&gibbsForSigmaEpsilonLME);
    }
  }

  // The Metropolis Hastings update for labels
  if(options.whichLabelSwitch().compare("123")==0){
    pReMiuMSampler.addProposal("metropolisHastingsForLabels123",1.0,1,1,&metropolisHastingsForLabels123);
  } else if (options.whichLabelSwitch().compare("12")==0){
    pReMiuMSampler.addProposal("metropolisHastingsForLabels12",1.0,1,1,&metropolisHastingsForLabels12);
  } else if(options.whichLabelSwitch().compare("3")==0){
    pReMiuMSampler.addProposal("metropolisHastingsForLabels3",1.0,1,1,&metropolisHastingsForLabels3);
  }
  // Gibbs for U
  if(options.samplerType().compare("Truncated")!=0){
    pReMiuMSampler.addProposal("gibbsForU",1.0,1,1,&gibbsForU);
  }

  // The Metropolis Hastings update for alpha
  if(options.fixedAlpha()<=-1){
    pReMiuMSampler.addProposal("metropolisHastingsForAlpha",1.0,1,1,&metropolisHastingsForAlpha);
  }

  // The gibbs update for the inactive V
  pReMiuMSampler.addProposal("gibbsForVInActive",1.0,1,1,&gibbsForVInActive);

  if(options.covariateType().compare("Discrete")==0){
    // For discrete X data we do a mixture of Categorical and ordinal updates
    //  Update for the inactive phi parameters
    pReMiuMSampler.addProposal("gibbsForPhiInActive",1.0,1,1,&gibbsForPhiInActive);

  }else if(options.covariateType().compare("Normal")==0){
    // Need to add the proposals for the normal case
    // Update for the active mu parameters
    if (options.useNormInvWishPrior()){
      // If Inverse Normal Inverse Wishart Prior is used
      pReMiuMSampler.addProposal("gibbsForMuInActiveNIWP",1.0,1,1,&gibbsForMuInActiveNIWP);
    }else{
      // If independant Normal and Inverse Wishart priors are used
      pReMiuMSampler.addProposal("gibbsForMuInActive",1.0,1,1,&gibbsForMuInActive);
    }
    // Update for the active Sigma parameters
    pReMiuMSampler.addProposal("gibbsForTauInActive",1.0,1,1,&gibbsForTauInActive);

  }else if(options.covariateType().compare("Mixed")==0){

    // For discrete X data we do a mixture of Categorical and ordinal updates
    //  Update for the inactive phi parameters
    pReMiuMSampler.addProposal("gibbsForPhiInActive",1.0,1,1,&gibbsForPhiInActive);

    // Need to add the proposals for the normal case
    // Update for the active mu parameters
    if (options.useNormInvWishPrior()){
      // If Inverse Normal Inverse Wishart Prior is used
      pReMiuMSampler.addProposal("gibbsForMuInActiveNIWP",1.0,1,1,&gibbsForMuInActiveNIWP);
    }else{
      // If independant Normal and Inverse Wishart priors are used
      pReMiuMSampler.addProposal("gibbsForMuInActive",1.0,1,1,&gibbsForMuInActive);
    }
    // Update for the active Sigma parameters
    pReMiuMSampler.addProposal("gibbsForTauInActive",1.0,1,1,&gibbsForTauInActive);
  }

  if(dataset.nFixedEffects_mix()>0 && options.outcomeType().compare("LME")!=0){
    pReMiuMSampler.addProposal("gibbsForBetaInActive",1.0,1,1,&gibbsForBetaInActive);
  }

  if(options.varSelectType().compare("None")!=0){
    // Add the variable selection moves
    unsigned int firstSweep;
    firstSweep=1+(unsigned int)(options.nBurn()/10);
    if(options.varSelectType().compare("Continuous")!=0){
      // Gibbs update for gamma
      pReMiuMSampler.addProposal("gibbsForGammaInActive",1.0,1,firstSweep,&gibbsForGammaInActive);
    }
  }

  if(options.includeResponse()){
    // The Metropolis Hastings update for the inactive theta
    if(options.outcomeType().compare("Longitudinal")!=0 && options.outcomeType().compare("LME")!=0) //AR
      pReMiuMSampler.addProposal("gibbsForThetaInActive",1.0,1,1,&gibbsForThetaInActive);

    if(options.outcomeType().compare("Survival")==0&&!options.weibullFixedShape()) {
      pReMiuMSampler.addProposal("gibbsForNuInActive",1.0,1,1,&gibbsForNuInActive);
    }
    //RJ add sampler for L & MVN inactive
    if(options.outcomeType().compare("Longitudinal")==0){
      // Gibbs for L for longitudinal response
      pReMiuMSampler.addProposal("gibbsForLInActive",1.0,1,1,&gibbsForLInActive);
    }else if(options.outcomeType().compare("MVN")==0){
      pReMiuMSampler.addProposal("gibbsForMVNMuInActive",1.0,1,1,&gibbsForMVNMuInActive);
      pReMiuMSampler.addProposal("gibbsForMVNTauInActive",1.0,1,1,&gibbsForMVNTauInActive);
    }
  }

  if(options.varSelectType().compare("None")!=0){
    // Add the variable selection moves
    // Metropolis Hastings for joint update of rho and omega
    unsigned int firstSweep;
    firstSweep=1+(unsigned int)(options.nBurn()/10);

    pReMiuMSampler.addProposal("metropolisHastingsForRhoOmega",1.0,1,firstSweep,&metropolisHastingsForRhoOmega);
  }

  if(options.includeResponse()){
    if(options.outcomeType().compare("Normal")==0){
      // Gibbs for sigmaSqY for Normal response model
      pReMiuMSampler.addProposal("gibbsForSigmaSqY",1.0,1,1,&gibbsForSigmaSqY);
    }

    if(options.outcomeType().compare("Survival")==0){
      if (options.weibullFixedShape()) {
        // Gibbs for shape parameter of Weibull for survival response model
        pReMiuMSampler.addProposal("gibbsForNu",1.0,1,1,&gibbsForNu);
      }
    }
  }
  // Gibbs update for the allocation parameters
  pReMiuMSampler.addProposal("gibbsForZ",1.0,1,1,&gibbsForZ);
  /* ---------- Initialise the output files -----*/
  pReMiuMSampler.initialiseOutputFiles(options.outFileStem());
  /* ---------- Write the log file ------------- */
  // The standard log file
  pReMiuMSampler.writeLogFile();
  /* ---------- Initialise the chain ---- */
  pReMiuMSampler.initialiseChain();

  pReMiuMHyperParams hyperParams = pReMiuMSampler.chain().currentState().parameters().hyperParams();
  unsigned int nClusInit = pReMiuMSampler.chain().currentState().parameters().workNClusInit();

  // The following is only used if the sampler type is truncated
  unsigned int maxNClusters = pReMiuMSampler.chain().currentState().parameters().maxNClusters();

  /* ---------- Run the sampler --------- */
  // Note: in this function the output gets written

  pReMiuMSampler.run();


  /* -- End the clock time and write the full run details to log file --*/
  currTime = time(NULL);
  double timeInSecs=(double)currTime-(double)beginTime;
  string tmpStr = storeLogFileData(options,dataset,hyperParams,nClusInit,maxNClusters,timeInSecs);
  pReMiuMSampler.appendToLogFile(tmpStr);

  /* ---------- Clean Up ---------------- */
  pReMiuMSampler.closeOutputFiles();

  //int err = 0;
  return Rcpp::wrap(0);
  // alternative output
  // return Rcpp::List::create(Rcpp::Named("yModel")=options.outcomeType());

}

