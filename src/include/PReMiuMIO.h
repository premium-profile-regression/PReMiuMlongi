/// \file PReMiuMIO.h
/// \author David Hastie
/// \brief Header file for handling input and output for PReMiuM

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
/// is released under the Boost Software Licence, Version 1.0, a copy  of which is
/// included in the documentation directory.

/// Version 3.1.3 edited by Rob Johnson, March 2017


#ifndef DIPBACIO_H_
#define DIPBACIO_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<numeric>
#include<limits>
#include<cstdlib>
#include<sstream>
#include<algorithm>
#include<iterator>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>

#include<Rcpp.h>

// Custom includes
#include<MCMC/chain.h>
#include<MCMC/model.h>
#include<MCMC/sampler.h>
#include<PReMiuMOptions.h>
#include<PReMiuMData.h>
#include<PReMiuMModel.h>

using namespace Eigen;

using std::vector;
using std::ostringstream;
using std::istringstream;
using std::string;
using std::endl;
using std::stringstream;

// Process the command line run time options
pReMiuMOptions processCommandLine(string inputStr){

  /* ---------- Handle run time options ----------*/
  pReMiuMOptions options;

  vector<string> inputStrings;
  istringstream iss(inputStr);
  copy(std::istream_iterator<string>(iss),
       std::istream_iterator<string>(),
       std::back_inserter<vector<string> >(inputStrings));

  int argc = inputStrings.size();
  int currArg=1;
  bool wasError=false;
  if(argc>1){
    string inString = inputStrings[currArg];

    if(inString.compare("--help")==0){
      // Output help if requested
      Rprintf("\n\n ### PReMiuMpp Help Page. ### \n\n");
      Rprintf("Possible arguments (defaults in parentheses):\n\n");
      Rprintf("--help\n\tShow this help page\n");
      Rprintf("--input=<string>\n\tThe full file path to the input data (./input.txt)\n");
      Rprintf("--output=<string>\n\tThe file stem (including full file path) where the data\n\tshould be written (./output)\n");
      Rprintf("--hyper=<string>\n\tThe full file path to the file containing hyper\n\t parameters. (Hyper parameter file not used)\n");
      Rprintf("--predict=<string>\n\tThe full file path to the file containing prediction\n\t covariates. (Prediction file not used)\n");
      Rprintf("--nSweeps=<unsigned int>\n\tThe number of sweeps (after burn in) to run the\n\t sampler for (10000)\n");
      Rprintf("--nBurn=<unsigned int>\n\tThe number of sweeps in the burn in period (1000)\n");
      Rprintf("--reportBurnIn=<bool>\n\tIt enables reporting in the output files of the burn-in period (true).\n");
      Rprintf("--nProgress=<unsigned int>\n\tThe number of sweeps at which to print a\nprogress update (500)\n");
      Rprintf("--nFilter=<unsigned int>\n\tThe frequency (in sweeps) with which to write\n\tthe output to file (1)\n");
      Rprintf("--nClusInit=<unsigned int>\n\tThe number of clusters individuals should be\n\tinitially randomly assigned to (Unif[50,60])\n");
      Rprintf("--seed=<unsigned int>\n\tThe value for the seed for the random number\n\tgenerator (current time)\n");
      //RJ add Longitudinal and MVN to yModel printed options
      Rprintf("--yModel=<string>\n\tThe model type for the outcome variable. Options are\n\tcurrently 'Bernoulli','Poisson','Binomial', 'Categorical', 'Survival', 'Normal', 'MVN', 'Longitudinal' and 'LME'  (Bernoulli)\n");
      Rprintf("--xModel=<string>\n\tThe model type for the covariates. Options are\n\tcurrently 'Discrete', 'Normal' and 'Mixed' (Discrete)\n");
      Rprintf("--sampler=<string>\n\tThe sampler type to be used. Options are\n\tcurrently 'SliceDependent', 'SliceIndependent' and 'Truncated' (SliceDependent)\n");
      Rprintf("--alpha=<double>\n\tThe value to be used if alpha is to remain fixed.\n\tIf a negative value is used then alpha is updated (-2)\n");
      Rprintf("--dPitmanYor=<double>\n\tThe value to be used for the discount parameter of the Pitman-Yor process prior.\n\tThe default corresponds to the Dirichlet process prior (0)\n");
      Rprintf("--excludeY\n\tIf included only the covariate data X is modelled (not included)\n");
      Rprintf("--extraYVar\n\tIf included extra Gaussian variance is included in the\n\tresponse model (not included).\n");
      Rprintf("--varSelect=<string>\n\tThe type of variable selection to be used 'None',\n\t'BinaryCluster' or 'Continuous' (None)\n");
      Rprintf("--entropy\n\tIf included then we compute allocation entropy (not included)\n");
      Rprintf("--predictType=<string>\n\tThe type of predictions to be used 'RaoBlackwell' or 'random' (RaoBlackwell)\n");
      Rprintf("--weibullFixedShape=<bool>\n\tWhether the shape parameter of the Weibull distribution is fixed.\n");
      Rprintf("--kernel=<string>\n\tThe kernel type for the covariance function of the GP if yModel == Longitudinal. Options are\n\tcurrently 'SQexponential' and 'Quadratic' (SQexponential)\n");
      Rprintf("--sampleGPmean=<string>\n\tIndicator if the underlying mean function is sampled, if yModel == Longitudinal.\n");
      Rprintf("--estim_ratio=<string>\n\t Indicator of estimation of Ratio between L1k and L3k defining the variance of the GP, if yModel == Longitudinal.\n");

    }else{
      while(currArg < argc){
        inString.assign(inputStrings[currArg]);
        if(inString.find("--input")!=string::npos){
          size_t pos = inString.find("=")+1;
          string inFileName = inString.substr(pos,inString.size()-pos);
          options.inFileName(inFileName);
        }else if(inString.find("--output")!=string::npos){
          size_t pos = inString.find("=")+1;
          string outFileStem = inString.substr(pos,inString.size()-pos);
          options.outFileStem(outFileStem);
        }else if(inString.find("--hyper")!=string::npos){
          size_t pos = inString.find("=")+1;
          string hyperParamFileName = inString.substr(pos,inString.size()-pos);
          options.hyperParamFileName(hyperParamFileName);
        }else if(inString.find("--predict")!=string::npos){
          size_t pos = inString.find("=")+1;
          string predictFileName = inString.substr(pos,inString.size()-pos);
          options.predictFileName(predictFileName);
        }else if(inString.find("--nSweeps")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          unsigned int nSweeps = (unsigned int)atoi(tmpStr.c_str());
          options.nSweeps(nSweeps);
        }else if(inString.find("--nBurn")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          unsigned int nBurn=(unsigned int)atoi(tmpStr.c_str());
          options.nBurn(nBurn);
        }else if(inString.find("--reportBurnIn")!=string::npos){
          options.reportBurnIn(true);
        }else if(inString.find("--nProgress")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          unsigned int nProgress=(unsigned int)atoi(tmpStr.c_str());
          options.nProgress(nProgress);
        }else if(inString.find("--nFilter")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          unsigned int nFilter=(unsigned int)atoi(tmpStr.c_str());
          options.nFilter(nFilter);
        }else if(inString.find("--nClusInit")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          unsigned int nClusInit = (unsigned int)atoi(tmpStr.c_str());
          options.nClusInit(nClusInit);
        }else if(inString.find("--seed")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          long rndSeed=(long)atoi(tmpStr.c_str());
          options.seed(rndSeed);
        }else if(inString.find("--yModel")!=string::npos){
          size_t pos = inString.find("=")+1;
          string outcomeType = inString.substr(pos,inString.size()-pos);
          if(outcomeType.compare("Poisson")!=0&&outcomeType.compare("Bernoulli")!=0&&
             outcomeType.compare("Categorical")!=0&&outcomeType.compare("Survival")!=0&&
             outcomeType.compare("Binomial")!=0&&outcomeType.compare("Normal")!=0&&
             outcomeType.compare("MVN")!=0&&outcomeType.compare("LME")!=0&&
             outcomeType.compare("Longitudinal")!=0){//RJ add Longitudinal and MVN to outcome checklist
            // Illegal outcome model entered
            Rprintf("This yModel type is not accounted for\n");
            wasError=true;
            break;
          }
          options.outcomeType(outcomeType);
          if(outcomeType.compare("Normal")==0&&options.responseExtraVar()){
            Rprintf("Response extra variation not permitted with Normal response\n");
            options.responseExtraVar(false);
          }
          if(outcomeType.compare("Survival")==0&&options.responseExtraVar()){
            Rprintf("Response extra variation not permitted with Survival response\n");
            options.responseExtraVar(false);
          }
          //RJ cannot have Longitudinal or MVN and Response Extra Variation
          if((outcomeType.compare("MVN")==0||outcomeType.compare("Longitudinal")==0||outcomeType.compare("LME")==0)&&options.responseExtraVar()){
            Rprintf("Response extra variation not permitted with Longitudinal response\n");
            options.responseExtraVar(false);
          }
        }else if(inString.find("--kernel")!=string::npos){ //AR
          size_t pos = inString.find("=")+1;
          string kernelType = inString.substr(pos,inString.size()-pos);
          if(kernelType.compare("SQexponential")!=0&&kernelType.compare("Quadratic")!=0){
            // Illegal covariate type entered
            Rprintf("Problem with kernelType\n");
            wasError=true;
            break;
          }
          options.kernelType(kernelType);
        }else if(inString.find("--sampleGPmean")!=string::npos){ //AR
          options.sampleGPmean(true);
        }else if(inString.find("--estim_ratio")!=string::npos){ //AR
          options.estim_ratio(true);
        }else if(inString.find("--xModel")!=string::npos){
          size_t pos = inString.find("=")+1;
          string covariateType = inString.substr(pos,inString.size()-pos);
          if(covariateType.compare("Discrete")!=0&&covariateType.compare("Normal")!=0&&covariateType.compare("Mixed")!=0){
            // Illegal covariate type entered
            Rprintf("Problem with covariateType\n");
            wasError=true;
            break;
          }
          options.covariateType(covariateType);
        }else if(inString.find("--whichLabelSwitch")!=string::npos){
          size_t pos = inString.find("=")+1;
          string whichLabelSwitch = inString.substr(pos,inString.size()-pos);
          if(whichLabelSwitch.compare("123")!=0&&whichLabelSwitch.compare("12")!=0&&whichLabelSwitch.compare("3")!=0){
            // Illegal covariate type entered
            Rprintf("Problem with whichLabelSwitch\n");
            wasError=true;
            break;
          }
          options.whichLabelSwitch(whichLabelSwitch);
        }else if(inString.find("--sampler")!=string::npos){
          size_t pos = inString.find("=")+1;
          string samplerType = inString.substr(pos,inString.size()-pos);
          if(samplerType.compare("SliceDependent")!=0&&samplerType.compare("SliceIndependent")!=0
               &&samplerType.compare("Truncated")!=0){
            // Illegal sampler type entered
            Rprintf("Problem with samplerType\n");
            wasError=true;
            break;
          }
          options.samplerType(samplerType);
        }else if(inString.find("--alpha")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          double alpha=(double)atof(tmpStr.c_str());
          options.fixedAlpha(alpha);
        }else if(inString.find("--dPitmanYor")!=string::npos){
          size_t pos = inString.find("=")+1;
          string tmpStr = inString.substr(pos,inString.size()-pos);
          double dPitmanYor=(double)atof(tmpStr.c_str());
          options.dPitmanYor(dPitmanYor);
        }else if(inString.find("--excludeY")!=string::npos){
          options.includeResponse(false);
        }else if(inString.find("--includeCAR")!=string::npos){
          options.includeCAR(true);
        }else if(inString.find("--neighbours")!=string::npos){
          size_t pos = inString.find("=")+1;
          string neighboursFile = inString.substr(pos,inString.size()-pos);
          options.neighbourFileName(neighboursFile);
        }else if(inString.find("--extraYVar")!=string::npos){
          if(options.outcomeType().compare("Normal")!=0){
            options.responseExtraVar(true);
          }else{
            Rprintf("Response extra variation not permitted with Normal response\n");
          }
          if(options.outcomeType().compare("Survival")==0) Rprintf("Response extra variation not permitted with Survival response\n");
        }else if(inString.find("--varSelect")!=string::npos){
          size_t pos = inString.find("=")+1;
          string varSelectType = inString.substr(pos,inString.size()-pos);
          if(varSelectType.compare("None")!=0&&
             varSelectType.compare("BinaryCluster")!=0&&varSelectType.compare("Continuous")!=0){
            // Illegal type for variable selection entered
            Rprintf("Problem with varSelectType\n");
            wasError=true;
            break;
          }
          options.varSelectType(varSelectType);
        }else if(inString.find("--entropy")!=string::npos){
          options.computeEntropy(true);
        }else if(inString.find("--predType")!=string::npos){
          size_t pos = inString.find("=")+1;
          string predictType = inString.substr(pos,inString.size()-pos);
          if(predictType.compare("RaoBlackwell")!=0&&predictType.compare("random")!=0){
            // Illegal predictType type entered
            Rprintf("Problem with predictType\n");
            wasError=true;
            break;
          }
          options.predictType(predictType);
        }else if(inString.find("--weibullFixedShape")!=string::npos){
          options.weibullFixedShape(true);
        }else if(inString.find("--useNormInvWishPrior")!=string::npos){
          options.useNormInvWishPrior(true);

        }else{
          Rprintf("Unknown command line option.\n");
          wasError=true;
          break;
        }
        currArg++;
      }
    }
  }

  // Return if there was an error
  if(wasError){
    Rprintf("There is a mistake in the arguments provided in profRegr.\n");
    Rprintf("The code will be run with default values.\n");
    //	Rprintf("Please use:\n");
    //	Rprintf("\t profileRegression --help\n");
    //	Rprintf("to get help on correct usage.\n");
    //	exit(-1);
  }

  return options;

}

// Read the PReMiuM data set
void importPReMiuMData(const string& fitFilename,const string& predictFilename, const string& neighboursFilename, pReMiuMData& dataset){

  ifstream inputFile,predictFile;
  inputFile.open(fitFilename.c_str());
  if(!inputFile.is_open()){
    Rprintf("Input file not found\n");
    //	exit(-1);
  }
  if(predictFilename.compare("")!=0){
    predictFile.open(predictFilename.c_str());
    if(!predictFile.is_open()){
      Rprintf("Prediction covariate file not found\n");
      //		exit(-1);
    }
  }


  unsigned int& nSubjects=dataset.nSubjects();
  unsigned int& nOutcomes=dataset.nOutcomes();
  unsigned int& nTimes=dataset.nTimes();
  unsigned int& nTimes_unique=dataset.nTimes_unique(); //AR
  unsigned int& nCovariates=dataset.nCovariates();
  unsigned int& nDiscreteCovs=dataset.nDiscreteCovs();
  unsigned int& nContinuousCovs=dataset.nContinuousCovs();
  unsigned int& nCategoriesY=dataset.nCategoriesY();
  unsigned int& nPredictSubjects=dataset.nPredictSubjects();
  //RJ initialise pointer to equalTimes
  unsigned int& equalTimes=dataset.equalTimes();
  vector<unsigned int>& nTimes_m=dataset.nTimes_m();
  vector<unsigned int>& nCategories=dataset.nCategories();
  vector<unsigned int>& discreteY=dataset.discreteY();
  vector<double>& continuousY=dataset.continuousY();
  vector<vector<int> >& discreteX=dataset.discreteX();
  vector<vector<double> >& continuousX=dataset.continuousX();
  vector<string>& covNames=dataset.covariateNames();
  vector<vector<bool> >& missingX=dataset.missingX();
  vector<unsigned int>& nContinuousCovariatesNotMissing=dataset.nContinuousCovariatesNotMissing();
  vector<vector<double> >& W=dataset.W();
  vector<vector<double> >& W_mix=dataset.W_mix(); // Cluster-specific fixed efects if yModel != LME
  vector<MatrixXd>& W_RE=dataset.W_RE();
  vector<MatrixXd>& W_LME=dataset.W_LME();
  vector<MatrixXd>& W_LME_mix=dataset.W_LME_mix(); // Cluster-specific fixed efects if yModel == LME


  vector<unsigned int>& nFixedEffects=dataset.nFixedEffects();
  vector<unsigned int>& nFixedEffects_mix=dataset.nFixedEffects_mix();
  vector<vector<string>>& confNames=dataset.fixedEffectNames();
  vector<vector<string>>& confNames_mix=dataset.fixedEffectNames_mix();
  vector<unsigned int>& nRandomEffects=dataset.nRandomEffects();
  vector<vector<string>>& confNames_RE=dataset.randomEffectNames();
  vector<string>& outcNames=dataset.outcNames();


  string outcomeType = dataset.outcomeType();
  string covariateType = dataset.covariateType();
  //RJ add vectors for times, tStart, tStop
  vector<double>& times=dataset.times();
  vector<double>& times_corr=dataset.times_corr(); //AR
  vector<double>& times_unique=dataset.times_unique(); //AR
  vector<int>& tStart=dataset.tStart();
  vector<int>& tStop=dataset.tStop();
  vector<double>& logOffset=dataset.logOffset();
  vector<unsigned int>& nTrials=dataset.nTrials();
  vector<unsigned int>& censoring=dataset.censoring();
  vector<vector<unsigned int> >& neighbours=dataset.neighbours();
  vector<unsigned int>& nNeighbours=dataset.nNeighbours();
  bool& includeCAR=dataset.includeCAR();
  bool wasError=false;

  // Get the number of subjects
  inputFile >> nSubjects;

  // Get the number of outcomes
  inputFile >> nOutcomes;

  outcNames.resize(nOutcomes);
  for(unsigned int m=0;m<nOutcomes;m++){
    inputFile >> outcNames[m];
  }

  // Get the number of covariates
  inputFile >> nCovariates;
  if(outcomeType.compare("LME")!=0){
    //RJ get nTimes from inputFile
    inputFile >> nTimes;
  }else{
    nTimes_m.resize(nOutcomes);
    for(unsigned int i=0;i<nOutcomes;i++){
      inputFile >> nTimes_m[i];
      nTimes += nTimes_m[i];
    }
  }

  covNames.resize(nCovariates);
  if(covariateType.compare("Mixed")==0){
    inputFile >> nDiscreteCovs;
    inputFile >> nContinuousCovs;
    if(nDiscreteCovs+nContinuousCovs!=nCovariates){
      Rprintf("Illegal number of covariates, discrete covariates or continuous covariates\n");
      // Illegal number of covariates, discrete covariates or continuous covariates
      wasError=true;
    }
    if(nDiscreteCovs==0 || nContinuousCovs==0){
      Rprintf("If xModel=Mixed a positive number of discrete and continuous covariates must be provided\n");
      // Illegal number of discrete covariates or continuous covariates
      wasError=true;
    }
  } else {
    nDiscreteCovs = 0;
    nContinuousCovs = 0;
  }
  for(unsigned int i=0;i<nCovariates;i++){
    inputFile >> covNames[i];
  }

  // Get the number of fixed effects
  nFixedEffects.resize(nOutcomes);
  for(unsigned int m=0;m<nOutcomes;m++)
    inputFile >> nFixedEffects[m];

  confNames.resize(nOutcomes);
  for(unsigned int m=0;m<nOutcomes;m++){
    confNames[m].resize(nFixedEffects[m]);
    for(unsigned int i=0;i<nFixedEffects[m];i++)
      inputFile >> confNames[m][i];
  }

  nFixedEffects_mix.resize(nOutcomes);
  confNames_mix.resize(nOutcomes);

  if(outcomeType.compare("LME")!=0){

    for(unsigned int m=0;m<nOutcomes;m++)
      nFixedEffects_mix[m] = 0;

    // for(unsigned int m=0;m<nOutcomes;m++){
    //   confNames_mix[m].resize(nFixedEffects_mix[m]);
    //   for(unsigned int i=0;i<nFixedEffects_mix[m];i++)
    //     confNames_mix[m][i]=0;
    // }
  }else{

    for(unsigned int m=0;m<nOutcomes;m++)
      inputFile >> nFixedEffects_mix[m];


    for(unsigned int m=0;m<nOutcomes;m++){
      confNames_mix[m].resize(nFixedEffects_mix[m]);
      for(unsigned int i=0;i<nFixedEffects_mix[m];i++)
        inputFile >> confNames_mix[m][i];
    }

    nRandomEffects.resize(nOutcomes);
    for(unsigned int m=0;m<nOutcomes;m++)
      inputFile >> nRandomEffects[m];

    confNames_RE.resize(nOutcomes);
    for(unsigned int m=0;m<nOutcomes;m++){
      confNames_RE[m].resize(nRandomEffects[m]);
      for(unsigned int i=0;i<nRandomEffects[m];i++)
        inputFile >> confNames_RE[m][i];
    }
  }

  if(nOutcomes==0){
    nFixedEffects.resize(1);
    nFixedEffects[0]=0;
    nFixedEffects_mix.resize(1);
    nFixedEffects_mix[0]=0;
    nRandomEffects.resize(1);
    nRandomEffects[0]=0;
  }


  // Get the number of categories of outcome Y
  if(outcomeType.compare("Categorical")==0){
    inputFile >> nCategoriesY;
    nCategoriesY--;
  } else {
    nCategoriesY=1;
  }

  nCategories.resize(nCovariates);
  if(covariateType.compare("Discrete")==0){
    // Get the number of categories for each covariate
    for(unsigned int j=0;j<nCovariates;j++){
      inputFile >> nCategories[j];
    }
  }else if(covariateType.compare("Normal")==0){
    for(unsigned int j=0;j<nCovariates;j++){
      nCategories[j]=0;
    }
  }else if(covariateType.compare("Mixed")==0){
    for(unsigned int j=0;j<nDiscreteCovs;j++){
      inputFile >> nCategories[j];
    }
    for(unsigned int j=nDiscreteCovs;j<nCovariates;j++){
      nCategories[j]=0;
    }
  }

  if(predictFile.is_open()){
    predictFile >> nPredictSubjects;
  }

  // Get the data
  discreteY.resize(nSubjects);
  //RJ resize Y (nTimes), and times, tStart, and tStop
  if(outcomeType.compare("MVN")==0)
    nOutcomes = nTimes/nSubjects;
  continuousY.resize(nTimes);
  times.resize(nTimes);

  if(outcomeType.compare("LME")!=0){
    tStart.resize(nSubjects);
    tStop.resize(nSubjects);
  }else{
    tStart.resize(nSubjects*nOutcomes);
    tStop.resize(nSubjects*nOutcomes);
  }
  times_corr.resize(nTimes);

  discreteX.resize(nSubjects+nPredictSubjects);
  continuousX.resize(nSubjects+nPredictSubjects);
  W.resize(nSubjects);
  W_mix.resize(nSubjects);

  if(outcomeType.compare("Poisson")==0){
    logOffset.resize(nSubjects);
  }
  if(outcomeType.compare("Binomial")==0){
    nTrials.resize(nSubjects);
  }
  if(outcomeType.compare("Survival")==0){
    censoring.resize(nSubjects);
  }
  missingX.resize(nSubjects+nPredictSubjects);
  nContinuousCovariatesNotMissing.resize(nSubjects+nPredictSubjects);
  vector<double> meanX(nCovariates,0);
  vector<unsigned int> nXNotMissing(nCovariates,0);


  for(unsigned int i=0;i<nSubjects;i++){
    if(outcomeType.compare("LME")!=0){
      if(outcomeType.compare("Normal")==0||outcomeType.compare("Survival")==0||outcomeType.compare("MVN")==0){
        for(unsigned int j=0; j<nOutcomes; j++){//RJ read in MVN data
          inputFile >> continuousY[i*nOutcomes+j];
        }
      }else{
        inputFile >> discreteY[i];
      }
    }
    if(covariateType.compare("Discrete")==0 || covariateType.compare("Normal")==0){
      discreteX[i].resize(nCovariates);
      continuousX[i].resize(nCovariates);
    } else if(covariateType.compare("Mixed")==0){
      discreteX[i].resize(nDiscreteCovs);
      continuousX[i].resize(nContinuousCovs);
    }
    missingX[i].resize(nCovariates);

    for(unsigned int j=0;j<nCovariates;j++){
      missingX[i][j]=true;
      if(covariateType.compare("Discrete")==0){
        inputFile >> discreteX[i][j];
        // -999 is missing data indicator
        if(discreteX[i][j]!=-999){
          meanX[j]+=(double)discreteX[i][j];
          nXNotMissing[j]+=1;
          missingX[i][j]=false;
        }
      }else if(covariateType.compare("Normal")==0){
        inputFile >> continuousX[i][j];
        // -999 is missing data indicator
        if(fabs(continuousX[i][j]+999)>0.00000000001){
          nContinuousCovariatesNotMissing[i]++;
          meanX[j]+=continuousX[i][j];
          nXNotMissing[j]+=1;
          missingX[i][j]=false;
        }
      }else if(covariateType.compare("Mixed")==0){
        if (j < nDiscreteCovs) {
          inputFile >> discreteX[i][j];
          if(discreteX[i][j]!=-999){
            meanX[j]+=(double)discreteX[i][j];
            nXNotMissing[j]+=1;
            missingX[i][j]=false;
          }
        } else {
          inputFile >> continuousX[i][j-nDiscreteCovs];
          if(fabs(continuousX[i][j-nDiscreteCovs]+999)>0.00000000001){
            nContinuousCovariatesNotMissing[i]++;
            meanX[j]+=continuousX[i][j-nDiscreteCovs];
            nXNotMissing[j]+=1;
            missingX[i][j]=false;
          }
        }
      }
    }

    if(outcomeType.compare("LME")!=0 ){
      //Read outcomes
      if(nOutcomes>0){


        W[i].resize(nFixedEffects[0]);
        //W_mix[i].resize(nFixedEffects_mix);

        for(unsigned int j=0;j<nFixedEffects[0];j++){
          inputFile >> W[i][j]; //used only if Ymodel != LME
        }
        // for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //   inputFile >> W_mix[i][j];//used only if Ymodel != LME
        // }

        if(outcomeType.compare("Poisson")==0){
          double tmp;
          inputFile >> tmp;
          logOffset[i]=log(tmp);
        }
        if(outcomeType.compare("Binomial")==0){
          inputFile >> nTrials[i];
        }
        if(outcomeType.compare("Survival")==0){
          inputFile >> censoring[i];
        }
      }
    }
  }

  //RJ Read in timepoints and longitudinal data
  if(outcomeType.compare("Longitudinal")==0 || outcomeType.compare("LME")==0){

    int ind=0;
    for(unsigned int m=0; m<nOutcomes; m++){
      for(unsigned int i=0; i<nSubjects; i++){
        inputFile >> tStart[ind];
        inputFile >> tStop[ind];
        ind++;
      }
    }

    if(outcomeType.compare("Longitudinal")==0){
      for(unsigned int i=0; i<nTimes; i++){
        inputFile >> times[i];
        inputFile >> continuousY[i];
      }
    }else{//LME
      int dd=0;
      for(unsigned int m=0; m<nOutcomes; m++){
        for(unsigned int i=0; i<nTimes_m[m]; i++){
          inputFile >> times[dd];
          inputFile >> continuousY[dd];
          dd++;
        }
      }
    }

    int length = tStop[0] - tStart[0] + 1;
    equalTimes = length;
    for(unsigned int i=1; i<nSubjects; i++){
      if(!std::equal(times.begin(),times.begin()+length,times.begin()+length*i)){
        equalTimes = 0;
      }
    }

    if(outcomeType.compare("LME")==0){
      W_LME.resize(nOutcomes);
      for(unsigned int m=0;m<nOutcomes;m++){
        W_LME[m].setZero(nTimes_m[m],nFixedEffects[m]);
        for(unsigned int i=0;i<nTimes_m[m];i++){
          for(unsigned int k=0;k<nFixedEffects[m];k++){
            inputFile >> W_LME[m](i,k);
          }
        }
      }

      W_LME_mix.resize(nOutcomes);
      for(unsigned int m=0;m<nOutcomes;m++){
        W_LME_mix[m].setZero(nTimes_m[m],nFixedEffects_mix[m]);
        for(unsigned int i=0;i<nTimes_m[m];i++){
          for(unsigned int k=0;k<nFixedEffects_mix[m];k++){
            inputFile >> W_LME_mix[m](i,k);
          }
        }
      }

      W_RE.resize(nOutcomes);
      for(unsigned int m=0;m<nOutcomes;m++){
        W_RE[m].setZero(nTimes_m[m],nRandomEffects[m]);
        for(unsigned int i=0;i<nTimes_m[m];i++){
          for(unsigned int k=0;k<nRandomEffects[m];k++){
            inputFile >> W_RE[m](i,k);
          }
        }
      }
    }
  }

  if(outcomeType.compare("Longitudinal")==0){
    //AR Read in timepoints and longitudinal data
    inputFile >> nTimes_unique;
    times_unique.resize(nTimes_unique);
    for(unsigned int i=0; i<nTimes_unique; i++){
      inputFile >> times_unique[i];
    }
    for(unsigned int i=0; i<nTimes; i++){
      inputFile >> times_corr[i];
    }
  }

  for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
    if(covariateType.compare("Discrete")==0 || covariateType.compare("Normal")==0){
      discreteX[i].resize(nCovariates);
      continuousX[i].resize(nCovariates);
    } else if(covariateType.compare("Mixed")==0){
      discreteX[i].resize(nDiscreteCovs);
      continuousX[i].resize(nContinuousCovs);
    }
    missingX[i].resize(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      missingX[i][j]=true;
      if(covariateType.compare("Discrete")==0){
        predictFile >> discreteX[i][j];
        // -999 is missing data indicator
        if(discreteX[i][j]!=-999){
          missingX[i][j]=false;
        }
      }else if(covariateType.compare("Normal")==0){
        predictFile >> continuousX[i][j];
        // -999 is missing data indicator
        if(fabs(continuousX[i][j]+999)>0.00000000001){
          nContinuousCovariatesNotMissing[i]++;
          missingX[i][j]=false;
        }
      }else if(covariateType.compare("Mixed")==0){
        if (j < nDiscreteCovs) {
          predictFile >> discreteX[i][j];
          if(discreteX[i][j]!=-999){
            missingX[i][j]=false;
          }
        } else {
          predictFile >> continuousX[i][j-nDiscreteCovs];
          if(fabs(continuousX[i][j-nDiscreteCovs]+999)>0.00000000001){
            nContinuousCovariatesNotMissing[i]++;
            missingX[i][j]=false;
          }
        }
      }
    }
  }

  /// Initially we just replace missing values by their means
  for(unsigned int j=0;j<nCovariates;j++){
    meanX[j]=meanX[j]/(double)nXNotMissing[j];
  }


  for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
    for(unsigned int j=0;j<nCovariates;j++){
      if(missingX[i][j]){
        if(covariateType.compare("Discrete")==0){
          discreteX[i][j]=(int)meanX[j];
        }else if(covariateType.compare("Normal")==0){
          continuousX[i][j]=meanX[j];
        }else if(covariateType.compare("Mixed")==0){
          if (j < nDiscreteCovs) {
            discreteX[i][j]=(int)meanX[j];
          } else {
            continuousX[i][j-nDiscreteCovs]=(double)meanX[j-nDiscreteCovs];
          }
        }
      }
    }
  }

    inputFile.close();
    if(predictFile.is_open()){
      predictFile.close();
    }

    //Fill nNeighbours and Neighbours
    if (includeCAR){
      ifstream neighFile;
      neighFile.open(neighboursFilename.c_str());

      if (!neighFile.is_open()){
        Rprintf("Neighbourhood structure file not found\n");
        wasError = true;
      }
      if (neighFile.good()){
        string line;
        getline(neighFile, line);
        stringstream streamline(line);
        unsigned int nsub;
        streamline >> nsub;
        nNeighbours.resize(nSubjects);
        neighbours.resize(nSubjects);
      }
      int i=0;
      while (neighFile.good()){
        string line;
        getline(neighFile, line);
        stringstream streamline(line);
        int j;
        streamline>>j;
        streamline>>nNeighbours[j-1];
        neighbours[j-1].resize(nNeighbours[j-1]);
        int k=0;
        while (streamline.good()){
          streamline>>neighbours[j-1][k];
          k++;
        }
        i++;
      }
      neighFile.close();
    }

  // Return if there was an error
  if(wasError){
    Rprintf("Please use:\n");
    Rprintf("\t profileRegression --help\n");
    Rprintf("to get help on correct usage.\n");
    //	exit(-1);
  }
}

// Function to read the hyper parameters from file
void readHyperParamsFromFile(const string& filename,pReMiuMHyperParams& hyperParams, unsigned int nOutcomes){

  ifstream inputFile;
  inputFile.open(filename.c_str());
  if(!inputFile.is_open()){
    Rprintf("Parameter file not found\n");
    //	exit(-1);
  }

  string inString;

  bool wasError=false;

  while(!inputFile.eof()){
    getline(inputFile,inString);
    if(inString.find("shapeAlpha")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double shapeAlpha = (double)atof(tmpStr.c_str());
      hyperParams.shapeAlpha(shapeAlpha);
    }else if(inString.find("rateAlpha")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double rateAlpha = (double)atof(tmpStr.c_str());
      hyperParams.rateAlpha(rateAlpha);
      //		}else if(inString.find("useReciprocalNCatsPhi")==0){
      //			size_t pos = inString.find("=")+1;
      //			string tmpStr = inString.substr(pos,inString.size()-pos);
      //			bool useRecip = false;
      //			if(tmpStr.compare("true")==0){
      //				useRecip = true;
      //			}
      //			hyperParams.useReciprocalNCatsPhi(useRecip);
    }else if(inString.find("aPhi")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> aVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          aVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        aVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size());
      }
      hyperParams.aPhi(aVec);
    }else if(inString.find("mu0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> muVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          muVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        muVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size());
      }
      VectorXd mu0=VectorXd::Zero(muVec.size());
      for(unsigned int j=0;j<muVec.size();j++){
        mu0(j)=muVec[j];
      }
      hyperParams.mu0(mu0);
    }else if(inString.find("Tau0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> TauVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          TauVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        TauVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
      }
      unsigned int dim = (unsigned int)sqrt((double)TauVec.size());
      MatrixXd Tau0=MatrixXd::Zero(dim,dim);
      for(unsigned int j1=0;j1<dim;j1++){
        for(unsigned int j2=0;j2<dim;j2++){
          Tau0(j1,j2)=TauVec[j1*dim+j2];
        }
      }
      hyperParams.Tau0(Tau0);
    }else if(inString.find("R0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> RVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          RVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        RVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
      }
      unsigned int dim = (unsigned int)sqrt((double)RVec.size());
      MatrixXd R0=MatrixXd::Zero(dim,dim);
      for(unsigned int j1=0;j1<dim;j1++){
        for(unsigned int j2=0;j2<dim;j2++){
          R0(j1,j2)=RVec[j1*dim+j2];
        }
      }
      hyperParams.R0(R0);
    }else if(inString.find("kappa0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      unsigned int kappa0 = (unsigned int)atoi(tmpStr.c_str());
      hyperParams.kappa0(kappa0);
    }else if(inString.find("nu0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double nu0 = (double)atof(tmpStr.c_str());
      hyperParams.nu0(nu0);
    }else if(inString.find("muTheta")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double muTheta = (double)atof(tmpStr.c_str());
      hyperParams.muTheta(muTheta);
    }else if(inString.find("sigmaTheta")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double sigmaTheta = (double)atof(tmpStr.c_str());
      hyperParams.sigmaTheta(sigmaTheta);
    }else if(inString.find("dofTheta")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      unsigned int dofTheta = (unsigned int)atoi(tmpStr.c_str());
      hyperParams.dofTheta(dofTheta);
    }else if(inString.find("muBeta")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double muBeta = (double)atof(tmpStr.c_str());
      hyperParams.muBeta(muBeta);
    }else if(inString.find("sigmaBeta")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double sigmaBeta = (double)atof(tmpStr.c_str());
      hyperParams.sigmaBeta(sigmaBeta);
    }else if(inString.find("dofBeta")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      unsigned int dofBeta = (unsigned int)atoi(tmpStr.c_str());
      hyperParams.dofBeta(dofBeta);
      //RJ Read muL and SigmaL from hyper file
    }else if(inString.find("muL_signal")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double muL = (double)atof(tmpStr.c_str());
      hyperParams.muL(0,muL);
    }else if(inString.find("muL_lengthscale")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double muL = (double)atof(tmpStr.c_str());
      hyperParams.muL(1,muL);
    }else if(inString.find("muL_noise")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double muL = (double)atof(tmpStr.c_str());
      hyperParams.muL(2,muL);
    }else if(inString.find("muL_time")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double muL = (double)atof(tmpStr.c_str());
      hyperParams.muL(3,muL);
    }else if(inString.find("sigmaL_signal")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double sigmaL = (double)atof(tmpStr.c_str());
      hyperParams.sigmaL(0,sigmaL);
    }else if(inString.find("sigmaL_lengthscale")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double sigmaL = (double)atof(tmpStr.c_str());
      hyperParams.sigmaL(1,sigmaL);
    }else if(inString.find("sigmaL_noise")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double sigmaL = (double)atof(tmpStr.c_str());
      hyperParams.sigmaL(2,sigmaL);
    }else if(inString.find("sigmaL_time")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double sigmaL = (double)atof(tmpStr.c_str());
      hyperParams.sigmaL(3,sigmaL);
    }else if(inString.find("aRatio")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double aRatio = (double)atof(tmpStr.c_str());
      hyperParams.aRatio(aRatio);
    }else if(inString.find("bRatio")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double bRatio = (double)atof(tmpStr.c_str());
      hyperParams.bRatio(bRatio);
    }else if(inString.find("eps_vu")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double eps_vu = (double)atof(tmpStr.c_str());
      hyperParams.eps_vu(eps_vu);
    }else if(inString.find("eps_sigma2_0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double eps_sigma2_0 = (double)atof(tmpStr.c_str());
      hyperParams.eps_sigma2_0(eps_sigma2_0);
    }else if(inString.find("eps_shape")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double eps_shape = (double)atof(tmpStr.c_str());
      hyperParams.eps_shape(eps_shape);
    }else if(inString.find("eps_scale")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double eps_scale = (double)atof(tmpStr.c_str());
      hyperParams.eps_scale(eps_scale);
    }else if(inString.find("SigmaLME_R0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> SigmaVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          SigmaVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        SigmaVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
      }
      unsigned int dim = (unsigned int)sqrt((double)SigmaVec.size());
      MatrixXd SigmaLME_R0=MatrixXd::Zero(dim,dim);
      for(unsigned int j1=0;j1<dim;j1++){
        for(unsigned int j2=0;j2<dim;j2++){
          SigmaLME_R0(j1,j2)=SigmaVec[j1*dim+j2];
        }
      }
      for(unsigned int m=0;m<nOutcomes;m++){
        hyperParams.SigmaLME_R0(m,SigmaLME_R0);
      }

    }else if(inString.find("MVNmu0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> muVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          muVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        muVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size());
      }
      VectorXd mu0=VectorXd::Zero(muVec.size());
      for(unsigned int j=0;j<muVec.size();j++){
        mu0(j)=muVec[j];
      }
      hyperParams.MVNmu0(mu0);
    }else if(inString.find("MVNTau0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> TauVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          TauVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        TauVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
      }
      unsigned int dim = (unsigned int)sqrt((double)TauVec.size());
      MatrixXd Tau0=MatrixXd::Zero(dim,dim);
      for(unsigned int j1=0;j1<dim;j1++){
        for(unsigned int j2=0;j2<dim;j2++){
          Tau0(j1,j2)=TauVec[j1*dim+j2];
        }
      }
      hyperParams.MVNTau0(Tau0);
    }else if(inString.find("MVNR0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> RVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          RVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        RVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
      }
      unsigned int dim = (unsigned int)sqrt((double)RVec.size());
      MatrixXd R0=MatrixXd::Zero(dim,dim);
      for(unsigned int j1=0;j1<dim;j1++){
        for(unsigned int j2=0;j2<dim;j2++){
          R0(j1,j2)=RVec[j1*dim+j2];
        }
      }
      hyperParams.MVNR0(R0);
    }else if(inString.find("MVNkappa0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double kappa0 = (double)atof(tmpStr.c_str());
      hyperParams.MVNkappa0(kappa0);
    }else if(inString.find("MVNnu0")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      unsigned int nu0 = (unsigned int)atoi(tmpStr.c_str());
      hyperParams.MVNnu0(nu0);
    }else if(inString.find("shapeTauEpsilon")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double shapeTauEpsilon = (double)atof(tmpStr.c_str());
      hyperParams.shapeTauEpsilon(shapeTauEpsilon);
    }else if(inString.find("rateTauEpsilon")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double rateTauEpsilon = (double)atof(tmpStr.c_str());
      hyperParams.rateTauEpsilon(rateTauEpsilon);
    }else if(inString.find("aRho")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double aRho = (double)atof(tmpStr.c_str());
      hyperParams.aRho(aRho);
    }else if(inString.find("bRho")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double bRho = (double)atof(tmpStr.c_str());
      hyperParams.bRho(bRho);
    }else if(inString.find("atomRho")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double atomRho = (double)atof(tmpStr.c_str());
      hyperParams.atomRho(atomRho);
      if(hyperParams.atomRho()<=0 || hyperParams.atomRho()>1){
        // Illegal atomRho value entered - it must be in (0,1] where 1 corresponds to the non-sparsity inducing var selection
        wasError=true;
        break;
      }
    }else if(inString.find("shapeSigmaSqY")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double shapeSigmaSqY = (double)atof(tmpStr.c_str());
      hyperParams.shapeSigmaSqY(shapeSigmaSqY);
    }else if(inString.find("scaleSigmaSqY")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double scaleSigmaSqY = (double)atof(tmpStr.c_str());
      hyperParams.scaleSigmaSqY(scaleSigmaSqY);
    }else if(inString.find("shapeNu")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double shapeNu = (double)atof(tmpStr.c_str());
      hyperParams.shapeNu(shapeNu);
    }else if(inString.find("scaleNu")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double scaleNu = (double)atof(tmpStr.c_str());
      hyperParams.scaleNu(scaleNu);
    }else if(inString.find("rSlice")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double rSlice = (double)atof(tmpStr.c_str());
      hyperParams.rSlice(rSlice);
    }else if(inString.find("truncationEps")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double truncationEps = (double)atof(tmpStr.c_str());
      hyperParams.truncationEps(truncationEps);
    }else if(inString.find("shapeTauCAR")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double shapeTauCAR = (double)atof(tmpStr.c_str());
      hyperParams.shapeTauCAR(shapeTauCAR);
    }else if(inString.find("rateTauCAR")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      double rateTauCAR = (double)atof(tmpStr.c_str());
      hyperParams.rateTauCAR(rateTauCAR);
    }else if(inString.find("initAlloc")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> initAl;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          initAl.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        initAl.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size());
      }
      hyperParams.initAlloc(initAl);
    }
    else if(inString.find("initL")==0){
      size_t pos = inString.find("=")+1;
      string tmpStr = inString.substr(pos,inString.size()-pos);
      vector<double> LVec;
      while(tmpStr.find(" ")!=string::npos){
        pos = tmpStr.find(" ");
        if(pos==(tmpStr.size()-1)){
          string elem = tmpStr.substr(0,pos);
          LVec.push_back((double)atof(elem.c_str()));
          tmpStr = tmpStr.substr(pos+1,tmpStr.size());
          break;
        }
        string elem = tmpStr.substr(0,pos);
        LVec.push_back((double)atof(elem.c_str()));
        tmpStr = tmpStr.substr(pos+1,tmpStr.size()-pos-1);
      }

      unsigned int nL=4;
      if(!hyperParams.muL().empty()){ //AR
        nL=3;
      }
      unsigned int dim = (unsigned int)(double)LVec.size()/nL;
      MatrixXd L0=MatrixXd::Zero(dim,nL);
      for(unsigned int j1=0;j1<dim;j1++){
        for(unsigned int j2=0;j2<nL;j2++){
          L0(j1,j2)=LVec[j1*nL+j2];
        }
      }
      hyperParams.initL(L0);
    }

  }

  // Return if there was an error
  if(wasError){
    Rprintf("There is a mistake in the arguments provided in profRegr.\n");
    Rprintf("The code will be run with default values.\n");
    //	Rprintf("Please use:\n");
    //	Rprintf("\t profileRegression --help\n");
    //	Rprintf("to get help on correct usage.\n");
    //	exit(-1);
  }
}

// Initialise the PReMiuM object (needed in this file as it calls
// function to read hyper parameters)
void initialisePReMiuM(baseGeneratorType& rndGenerator,
                       const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                       pReMiuMParams& params){

  const pReMiuMData& dataset = model.dataset();
  const string kernelType=model.options().kernelType(); //AR
  const pReMiuMOptions& options = model.options();
  pReMiuMHyperParams& hyperParams = params.hyperParams();
  unsigned int nTimes_unique = dataset.nTimes_unique();

  unsigned int nSubjects=dataset.nSubjects();
  //RJ add int nTimes
  unsigned int nTimes=dataset.nTimes();
  unsigned int nOutcomes=dataset.nOutcomes();
  unsigned int nCovariates=dataset.nCovariates();
  unsigned int nDiscreteCovs=dataset.nDiscreteCovs();
  unsigned int nContinuousCovs=dataset.nContinuousCovs();
  vector<unsigned int> nRandomEffects=dataset.nRandomEffects();
  string outcomeType = options.outcomeType();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY=dataset.nCategoriesY();
  unsigned int nPredictSubjects=dataset.nPredictSubjects();
  unsigned int nClusInit = options.nClusInit();
  string covariateType = options.covariateType();
  string hyperParamFileName = options.hyperParamFileName();
  string varSelectType = options.varSelectType();
  string samplerType = options.samplerType();
  bool includeResponse = options.includeResponse();
  bool responseExtraVar = options.responseExtraVar();
  bool includeCAR=options.includeCAR();
  string predictType = options.predictType();
  bool weibullFixedShape = options.weibullFixedShape();
  vector<int> tStart=dataset.tStart();
  vector<int> tStop=dataset.tStop();
  vector<double> y=dataset.continuousY();
  vector<unsigned int> nCategories;
  nCategories = dataset.nCategories();

  // Set the hyper parameters to their default values
  hyperParams.setSizes(nCovariates,nDiscreteCovs,
                       nContinuousCovs,nOutcomes,covariateType,outcomeType, nRandomEffects);

  hyperParams.setDefaults(dataset,options);

  // Read the parameters from file if file provided
  if(hyperParamFileName.compare("")!=0){
    readHyperParamsFromFile(hyperParamFileName,hyperParams, nOutcomes);
  }

  // Allocate the right sizes for each of the parameter variables
  // This also switches "on" all variable indicators (gamma)
  // This gets changed below if variable selection is being done
  params.setSizes(nSubjects,nCovariates,nDiscreteCovs,nContinuousCovs,nFixedEffects,nFixedEffects_mix,nCategoriesY,	nPredictSubjects,nTimes,nTimes_unique,nOutcomes, nCategories,nClusInit,covariateType,outcomeType,weibullFixedShape,kernelType, nRandomEffects); //AR

  unsigned int maxNClusters=params.maxNClusters();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  if(nClusInit==0){
    nClusInit=50+(unsigned int)11*unifRand(rndGenerator);
  }

  // Fix the number of clusters if we are using the truncated sampler
  if(samplerType.compare("Truncated")==0){
    maxNClusters=20;
    if((nClusInit+10)>maxNClusters){
      maxNClusters=nClusInit+10;
    }
    // Now compute the bound recommended in Ishwaran and James 2001
    double multiplier=0.0;
    if(options.fixedAlpha()>-1){
      multiplier=options.fixedAlpha();
    }else{
      // Use the expected value of alpha as the multiplier
      multiplier=hyperParams.shapeAlpha()/hyperParams.rateAlpha();
    }
    double computedBound=1+multiplier*(log(4.0*nSubjects)-log(hyperParams.truncationEps()));
    if(computedBound>maxNClusters){
      maxNClusters=computedBound;
    }
    params.maxNClusters(maxNClusters,covariateType,outcomeType,kernelType, nTimes_unique, nRandomEffects);//AR
  }

  // Copy the dataset X matrix to a working object in params
  params.workDiscreteX(dataset.discreteX());
  params.workContinuousX(dataset.continuousX());

  // Now initialise the actual parameters
  randomGamma gammaRand(hyperParams.shapeAlpha(),1.0/hyperParams.rateAlpha());

  double alpha=gammaRand(rndGenerator);
  if(options.fixedAlpha()>-1){
    alpha=options.fixedAlpha();
  }
  params.alpha(alpha);

  double dPitmanYor = options.dPitmanYor();
  params.dPitmanYor(dPitmanYor);

  vector<unsigned int> nXInCluster(maxNClusters,0);
  unsigned int maxZ=0;
  params.workNClusInit(nClusInit);
  if (hyperParams.initAlloc().empty()){
    for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
      int c=(int) nClusInit*unifRand(rndGenerator);
      params.z(i,c,covariateType);
      if(c>(int)maxZ){
        maxZ=c;
      }
      if(i<nSubjects){
        nXInCluster[c]++;
      }
    }
  } else {
    for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
      int c = hyperParams.initAlloc(i);
      params.z(i,c,covariateType);
      if(c>(int)maxZ){
        maxZ=c;
      }
      if(i<nSubjects){
        nXInCluster[c]++;
      }
    }
  }
  params.workNXInCluster(nXInCluster);
  params.workMaxZi(maxZ);

  // Sample v (for logPsi)
  // This is sampled from the posterior given the z vector above
  // Prior comes from the conjugacy of the dirichlet and multinomial
  // See Ishwaran and James 2001

  // Sample active V
  vector<unsigned int> sumCPlus1ToMaxMembers(maxZ+1,0);
  for(int c=maxZ-1;c>=0;c--){
    sumCPlus1ToMaxMembers[c]=sumCPlus1ToMaxMembers[c+1]+params.workNXInCluster(c+1);
  }

  double tmp=0.0;
  for(unsigned int c=0;c<=maxZ;c++){
    double vVal = betaRand(rndGenerator,1.0+params.workNXInCluster(c)-dPitmanYor,alpha+sumCPlus1ToMaxMembers[c]+dPitmanYor*(c+1));
    params.v(c,vVal);
    // Set logPsi
    params.logPsi(c,tmp+log(vVal));
    tmp += log(1-vVal);
  }

  if(samplerType.compare("Truncated")==0){
    // Just sample the remaining V from the prior
    vector<double> vNew=params.v();
    vector<double> logPsiNew=params.logPsi();

    for(unsigned int c=maxZ+1;c<maxNClusters;c++){
      double v=betaRand(rndGenerator,1.0-dPitmanYor,alpha+dPitmanYor*c);
      double logPsi=log(v)+log(1-vNew[c-1])-log(vNew[c-1])+logPsiNew[c-1];
      if(c>=vNew.size()){
        vNew.push_back(v);
        logPsiNew.push_back(logPsi);
      }else{
        vNew[c]=v;
        logPsiNew[c]=logPsi;
      }
    }
    params.v(vNew);
    params.logPsi(logPsiNew);

  }else{

    // Sample u (auxilliary variables). This will determine the maximum number of clusters
    double minU=1.0;
    for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
      int zi=params.z(i);
      double ui=0.0;
      if(samplerType.compare("SliceDependent")==0){
        ui = exp(params.logPsi(zi))*unifRand(rndGenerator);
      }else if(samplerType.compare("SliceIndependent")==0){
        ui = hyperParams.workXiSlice(zi)*unifRand(rndGenerator);
      }
      if(ui<minU){
        minU=ui;
      }
      params.u(i,ui);
    }
    params.workMinUi(minU);

    // Sample V
    vector<double> cumPsi(maxZ+1,0.0);
    cumPsi[0] = exp(params.logPsi(0));
    for(unsigned int c=1;c<=maxZ;c++){
      cumPsi[c]=cumPsi[c-1]+exp(params.logPsi(c));
    }

    vector<double> vNew=params.v();
    vector<double> logPsiNew=params.logPsi();

    maxNClusters = maxZ+1;
    if(samplerType.compare("SliceIndependent")==0){

      maxNClusters=2+(int)((log(params.workMinUi())-log(1.0-hyperParams.rSlice()))/log(hyperParams.rSlice()));
    }

    bool continueLoop=true;
    unsigned int c=maxZ;
    while(continueLoop){
      if(samplerType.compare("SliceDependent")==0&&cumPsi[c]>1-minU){
        // We can stop
        maxNClusters=c+1;
        continueLoop=false;
      }else if(samplerType.compare("SliceIndependent")==0&&c>=maxNClusters){
        continueLoop=false;
      }else{
        c++;
        // We need a new sampled value of v
        double v=betaRand(rndGenerator,1.0-dPitmanYor,alpha+dPitmanYor*c);
        double logPsi=log(v)+log(1-vNew[c-1])-log(vNew[c-1])+logPsiNew[c-1];
        if(c>=vNew.size()){
          vNew.push_back(v);
          logPsiNew.push_back(logPsi);
        }else{
          vNew[c]=v;
          logPsiNew[c]=logPsi;
        }
        cumPsi.push_back(cumPsi[c-1]+exp(logPsi));
      }
    }

    params.maxNClusters(maxNClusters,covariateType,outcomeType,kernelType, nTimes_unique, nRandomEffects);
    params.v(vNew);
    params.logPsi(logPsiNew);
  }

  if(covariateType.compare("Discrete")==0){
    // Sample logPhi
    // Need to count the number of X[i][j]==p for each covariate and category of p
    vector<vector<unsigned int> > nXpMembers(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      nXpMembers[j].resize(nCategories[j]);
      for(unsigned int p=0;p<nCategories[j];p++){
        nXpMembers[j][p]=0;
        for(unsigned int i=0;i<nSubjects;i++){
          if(dataset.discreteX(i,j)==(int)p&&!dataset.missingX(i,j)){
            nXpMembers[j][p]++;
          }
        }
      }
    }

    // Now we can sample. We don't use the priors, but instead look at the number
    // of people in each category and do a dirichlet sample that takes account of that
    boost::math::normal_distribution<double> norm01(0.0,1.0);
    for(unsigned int c=0;c<maxNClusters;c++){
      for(unsigned int j=0;j<nCovariates;j++){
        vector<double> dirichParams(nCategories[j]);
        for(unsigned int p=0;p<nCategories[j];p++){
          dirichParams[p]=(double)nXpMembers[j][p]+hyperParams.aPhi(j);
        }
        vector<double> logDirichSample(nCategories[j]);
        vector<double> dirichSample(nCategories[j]);
        dirichSample=dirichletRand(rndGenerator,dirichParams);
        for(unsigned int p=0;p<nCategories[j];p++){
          logDirichSample[p]=log(dirichSample[p]);
        }
        params.logPhi(c,j,logDirichSample);
      }
    }
    // Initialise the null parameters for the variable selection case
    // In all cases, initialise it at the value it will be fixed at for
    // the continuous indicator case
    if(varSelectType.compare("None")!=0){
      for(unsigned int j=0;j<nCovariates;j++){
        double sumVec=0.0;
        vector<double> probVec(nCategories[j],0.0000001);
        for(unsigned int p=0;p<nCategories[j];p++){
          probVec[p]+=(double)(nXpMembers[j][p]);
          sumVec+=(double)nXpMembers[j][p];
        }
        vector<double> logProbVec(nCategories[j]);
        for(unsigned int p=0;p<nCategories[j];p++){
          logProbVec[p]=log(probVec[p]/sumVec);
        }
        params.logNullPhi(j,logProbVec);
      }
    }

  }else if(covariateType.compare("Normal")==0){
    // In the following it is useful to have the rows of X as
    // Eigen dynamic vectors
    vector<VectorXd> xi(nSubjects);
    for(unsigned int i=0;i<nSubjects;i++){
      xi[i].setZero(nCovariates);
      for(unsigned int j=0;j<nCovariates;j++){
        xi[i](j)=dataset.continuousX(i,j);
      }
    }

    // Now we can sample from the conditionals (using Sigma_c=Sigma_0 and
    // mu_c=mu_0 for all c) to get mu_c and Sigma_c for each cluster

    // First we sample mu_c for each cluster

    // We begin by computing the mean X for individuals in each cluster
    vector<VectorXd> meanX(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      meanX[c].setZero(nCovariates);
    }
    for(unsigned int i=0;i<nSubjects;i++){
      meanX[params.z(i)]=meanX[params.z(i)]+xi[i];
    }

    for(unsigned int c=0;c<maxNClusters;c++){
      // Having computed this we can calcuate the posterior mean
      // and posterior covariance for each mu_c
      if(params.workNXInCluster(c)>0){
        meanX[c]=meanX[c]/(double)params.workNXInCluster(c);
      }else{
        meanX[c].setZero(nCovariates);
      }
      MatrixXd covMat(nCovariates,nCovariates);
      covMat = (hyperParams.Tau0()+params.workNXInCluster(c)*hyperParams.Tau0()).inverse();
      VectorXd meanVec(nCovariates);
      meanVec = hyperParams.Tau0()*hyperParams.mu0()+params.workNXInCluster(c)*hyperParams.Tau0()*meanX[c];
      meanVec = covMat*meanVec;

      VectorXd mu(nCovariates);
      // We sample from this posterior
      mu = multivarNormalRand(rndGenerator,meanVec,covMat);

      // We store our sample
      params.mu(c,mu);

    }

    // Now we can sample Tau_c for each cluster
    vector<MatrixXd> Rc(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      Rc[c].setZero(nCovariates,nCovariates);
    }

    for(unsigned int i=0;i<nSubjects;i++){
      unsigned int zi = params.z(i);
      Rc[zi]=Rc[zi]+(xi[i]-params.mu(zi))*((xi[i]-params.mu(zi)).transpose());
    }

    for(unsigned int c=0;c<maxNClusters;c++){
      Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
      MatrixXd Tau = wishartRand(rndGenerator,Rc[c],params.workNXInCluster(c)+hyperParams.kappa0());
      params.Tau(c,Tau);
    }

    // Now do the null mu for variable selection
    // In all cases, initialise it at the value it will be fixed at for
    // the continuous indicator case
    if(varSelectType.compare("None")!=0){
      vector<double> meanXVec(nCovariates,0.0);
      vector<unsigned int> countXVec(nCovariates,0);
      for(unsigned int i=0;i<nSubjects;i++){
        for(unsigned int j=0;j<nCovariates;j++){
          if(!dataset.missingX(i,j)){
            meanXVec[j]+=dataset.continuousX(i,j);
            countXVec[j]+=1;
          }
        }
      }
      VectorXd nullMu=VectorXd::Zero(nCovariates);
      for(unsigned int j=0;j<nCovariates;j++){
        nullMu(j)=meanXVec[j]/(double)countXVec[j];
      }
      params.nullMu(nullMu);
    }

  }else if(covariateType.compare("Mixed")==0){

    // Sample logPhi
    // Need to count the number of X[i][j]==p for each covariate and category of p
    vector<vector<unsigned int> > nXpMembers(nDiscreteCovs);
    for(unsigned int j=0;j<nDiscreteCovs;j++){
      nXpMembers[j].resize(nCategories[j]);
      for(unsigned int p=0;p<nCategories[j];p++){
        nXpMembers[j][p]=0;
        for(unsigned int i=0;i<nSubjects;i++){
          if(dataset.discreteX(i,j)==(int)p&&!dataset.missingX(i,j)){
            nXpMembers[j][p]++;
          }
        }
      }
    }

    // Now we can sample. We don't use the priors, but instead look at the number
    // of people in each category and do a dirichlet sample that takes account of that
    boost::math::normal_distribution<double> norm01(0.0,1.0);
    for(unsigned int c=0;c<maxNClusters;c++){
      for(unsigned int j=0;j<nDiscreteCovs;j++){
        vector<double> dirichParams(nCategories[j]);
        for(unsigned int p=0;p<nCategories[j];p++){
          dirichParams[p]=(double)nXpMembers[j][p]+hyperParams.aPhi(j);
        }
        vector<double> logDirichSample(nCategories[j]);
        vector<double> dirichSample(nCategories[j]);
        dirichSample=dirichletRand(rndGenerator,dirichParams);
        for(unsigned int p=0;p<nCategories[j];p++){
          logDirichSample[p]=log(dirichSample[p]);
        }
        params.logPhi(c,j,logDirichSample);
      }
    }
    // Initialise the null parameters for the variable selection case
    // In all cases, initialise it at the value it will be fixed at for
    // the continuous indicator case
    if(varSelectType.compare("None")!=0){
      for(unsigned int j=0;j<nDiscreteCovs;j++){
        double sumVec=0.0;
        vector<double> probVec(nCategories[j],0.0000001);
        for(unsigned int p=0;p<nCategories[j];p++){
          probVec[p]+=(double)(nXpMembers[j][p]);
          sumVec+=(double)nXpMembers[j][p];
        }
        vector<double> logProbVec(nCategories[j]);
        for(unsigned int p=0;p<nCategories[j];p++){
          logProbVec[p]=log(probVec[p]/sumVec);
        }
        params.logNullPhi(j,logProbVec);
      }
    }

    // In the following it is useful to have the rows of X as
    // Eigen dynamic vectors
    vector<VectorXd> xi(nSubjects);
    for(unsigned int i=0;i<nSubjects;i++){
      xi[i].setZero(nContinuousCovs);
      for(unsigned int j=0;j<nContinuousCovs;j++){
        xi[i](j)=dataset.continuousX(i,j);
      }
    }

    // Now we can sample from the conditionals (using Sigma_c=Sigma_0 and
    // mu_c=mu_0 for all c) to get mu_c and Sigma_c for each cluster

    // First we sample mu_c for each cluster

    // We begin by computing the mean X for individuals in each cluster
    vector<VectorXd> meanX(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      meanX[c].setZero(nContinuousCovs);
    }
    for(unsigned int i=0;i<nSubjects;i++){
      meanX[params.z(i)]=meanX[params.z(i)]+xi[i];
    }

    for(unsigned int c=0;c<maxNClusters;c++){
      // Having computed this we can calcuate the posterior mean
      // and posterior covariance for each mu_c
      if(params.workNXInCluster(c)>0){
        meanX[c]=meanX[c]/(double)params.workNXInCluster(c);
      }else{
        meanX[c].setZero(nContinuousCovs);
      }
      MatrixXd covMat(nContinuousCovs,nContinuousCovs);
      covMat = (hyperParams.Tau0()+params.workNXInCluster(c)*hyperParams.Tau0()).inverse();
      VectorXd meanVec(nContinuousCovs);
      meanVec = hyperParams.Tau0()*hyperParams.mu0()+params.workNXInCluster(c)*hyperParams.Tau0()*meanX[c];
      meanVec = covMat*meanVec;

      VectorXd mu(nContinuousCovs);
      // We sample from this posterior
      mu = multivarNormalRand(rndGenerator,meanVec,covMat);

      // We store our sample
      params.mu(c,mu);

    }

    // Now we can sample Tau_c for each cluster
    vector<MatrixXd> Rc(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      Rc[c].setZero(nContinuousCovs,nContinuousCovs);
    }

    for(unsigned int i=0;i<nSubjects;i++){
      unsigned int zi = params.z(i);
      Rc[zi]=Rc[zi]+(xi[i]-params.mu(zi))*((xi[i]-params.mu(zi)).transpose());
    }

    for(unsigned int c=0;c<maxNClusters;c++){
      Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
      MatrixXd Tau = wishartRand(rndGenerator,Rc[c],params.workNXInCluster(c)+hyperParams.kappa0());
      params.Tau(c,Tau);
    }

    // Now do the null mu for variable selection
    // In all cases, initialise it at the value it will be fixed at for
    // the continuous indicator case
    if(varSelectType.compare("None")!=0){
      vector<double> meanXVec(nContinuousCovs,0.0);
      vector<unsigned int> countXVec(nContinuousCovs,0);
      for(unsigned int i=0;i<nSubjects;i++){
        for(unsigned int j=0;j<nContinuousCovs;j++){
          if(!dataset.missingX(i,j)){
            meanXVec[j]+=dataset.continuousX(i,j);
            countXVec[j]+=1;
          }
        }
      }
      VectorXd nullMu=VectorXd::Zero(nContinuousCovs);
      for(unsigned int j=0;j<nContinuousCovs;j++){
        nullMu(j)=meanXVec[j]/(double)countXVec[j];
      }
      params.nullMu(nullMu);
    }
  }

  // Initialise the variable selection variables if appropriate
  // Bias towards having variables in
  if(varSelectType.compare("None")!=0){
    vector<vector<double> > gamma(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      gamma[c].resize(nCovariates);
    }
    vector<unsigned int> omega(nCovariates);
    vector<double> rho(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      if((unifRand(rndGenerator)<0.01) && (hyperParams.atomRho()!=1)){
        // We are in the point mass at 0 case - variable is switched off
        omega[j]=0;
        rho[j]=0;
        if(varSelectType.compare("BinaryCluster")==0){
          for(unsigned int c=0;c<maxNClusters;c++){
            gamma[c][j]=0;
          }
        }else{
          gamma[0][j]=0;
        }
      }else{
        omega[j]=1;
        rho[j]=0.75+0.25*unifRand(rndGenerator);
        if(varSelectType.compare("BinaryCluster")==0){
          for(unsigned int c=0;c<maxNClusters;c++){
            if(unifRand(rndGenerator)<rho[j]){
              gamma[c][j]=1;
            }else{
              gamma[c][j]=0;
            }
          }
        }
      }

      params.omega(j,omega[j]);
      params.rho(j,rho[j],covariateType,varSelectType);

      if(varSelectType.compare("BinaryCluster")==0){
        for(unsigned int c=0;c<maxNClusters;c++){
          params.gamma(c,j,gamma[c][j],covariateType);
        }
      }
      // Note in the case of the continuous variable selection indicators
      // gamma is deterministically equal to rho, and so is set in the method
      // for rho so we do nothing here.
    }
  }


  if(includeResponse){
    // Finally we sample the theta and beta values from uniform distributions
    for(unsigned int c=0;c<maxNClusters;c++){
      for (unsigned int k=0;k<nCategoriesY;k++){
        // Thetas are randomly between -2 and 2
        params.theta(c,k,-2.0+4.0*unifRand(rndGenerator));
      }
    }

    for(unsigned int m=0;m<nOutcomes;m++){
      for(unsigned int j=0;j<nFixedEffects[m];j++){
        for (unsigned int k=0;k<nCategoriesY;k++){
          // Betas are randomly between -2 and 2
          if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){ //AR
            params.beta(m,j,k,nCategoriesY,-2.0+4.0*unifRand(rndGenerator));
          }else{
            params.beta(m,j,k,nCategoriesY,0);
          }
        }
      }

      if(outcomeType.compare("LME")==0){ //AR
        for(unsigned int j=0;j<nFixedEffects_mix[m];j++){
          for(unsigned int c=0;c<maxNClusters;c++){
            for (unsigned int k=0;k<nCategoriesY;k++)
              params.beta_mix(m,c,j,k,nCategoriesY,0);
          }
        }
      }
    }


    if(outcomeType.compare("Normal")==0){
      randomGamma gammaRand(hyperParams.shapeSigmaSqY(),1.0/hyperParams.scaleSigmaSqY());
      double sigmaSqY=1.0/(gammaRand(rndGenerator));
      params.sigmaSqY(sigmaSqY);
    }

    if(outcomeType.compare("Survival")==0){
      randomGamma gammaRand(hyperParams.shapeNu(),hyperParams.scaleNu());
      if(weibullFixedShape){
        double nu=gammaRand(rndGenerator);
        params.nu(0,nu);

      } else {
        for (unsigned int c=0;c<maxNClusters;c++){
          double nu=gammaRand(rndGenerator);
          params.nu(c,nu);
        }
      }
    }else if(outcomeType.compare("MVN")==0){
      vector<VectorXd> xi(nSubjects);
      for(unsigned int i=0;i<nSubjects;i++){
        xi[i].setZero(nOutcomes);
        for(unsigned int j=0;j<nOutcomes;j++){
          xi[i](j)=dataset.continuousY(i*nOutcomes+j);
        }
      }
      vector<VectorXd> meanX(maxNClusters);
      for(unsigned int c=0;c<maxNClusters;c++){
        meanX[c].setZero(nOutcomes);
      }
      for(unsigned int i=0;i<nSubjects;i++){
        meanX[params.z(i)]=meanX[params.z(i)]+xi[i];
      }
      for(unsigned int c=0;c<maxNClusters;c++){
        if(params.workNXInCluster(c)>0){
          meanX[c]=meanX[c]/(double)params.workNXInCluster(c);
        }else{
          meanX[c].setZero(nOutcomes);
        }
        MatrixXd covMat(nOutcomes,nOutcomes);
        covMat = (hyperParams.MVNTau0()+params.workNXInCluster(c)*hyperParams.MVNTau0()).inverse();
        VectorXd meanVec(nOutcomes);
        meanVec = hyperParams.MVNTau0()*hyperParams.MVNmu0()+params.workNXInCluster(c)*hyperParams.MVNTau0()*meanX[c];
        meanVec = covMat*meanVec;
        VectorXd mu(nOutcomes);
        mu = multivarNormalRand(rndGenerator,meanVec,covMat);
        params.MVNmu(c,mu);
      }
      vector<MatrixXd> Rc(maxNClusters);
      for(unsigned int c=0;c<maxNClusters;c++){
        Rc[c].setZero(nOutcomes,nOutcomes);
      }
      for(unsigned int i=0;i<nSubjects;i++){
        unsigned int zi = params.z(i);
        Rc[zi]=Rc[zi]+(xi[i]-params.MVNmu(zi))*((xi[i]-params.MVNmu(zi)).transpose());
      }
      for(unsigned int c=0;c<maxNClusters;c++){
        Rc[c]=(hyperParams.MVNR0().inverse()+Rc[c]).inverse();
        //MatrixXd Tau = wishartRand(rndGenerator,Rc[c],params.workNXInCluster(c)+hyperParams.MVNkappa0());
        MatrixXd Tau ;
        Tau = wishartRand(rndGenerator,Rc[c],params.workNXInCluster(c)+hyperParams.MVNnu0());
        params.MVNTau(c,Tau);
      }
    }

    if(outcomeType.compare("LME")==0){ //AR

      int ind=0;
      int ind_y=0;
      for(unsigned int m=0;m<nOutcomes;m++){
        MatrixXd Tau = wishartRand(rndGenerator,hyperParams.workTauLME_R0(m),hyperParams.SigmaLME_kappa0(m));

        for(unsigned int c=0;c<maxNClusters;c++){
          MatrixXd Tauinv = Tau.inverse();
          params.covRE(m,c, Tauinv);
          //MatrixXd covRE=params.covRE(m,c);
          //_workLogDetTauLME(m,c)=-log(cov.determinant());
          //_workSqrtTauLME[m][c]=(llt.compute(cov.inverse())).matrixU();

          params.workLogDetTauLME(m,c,log(Tau.determinant()));
          LLT<MatrixXd> llt;
          params.workSqrtTauLME(m,c,(llt.compute(Tau)).matrixU());
        }



        // Initialise random effects
        for(unsigned int i=0;i<nSubjects;i++){
          VectorXd yi;
          VectorXd ui(nRandomEffects[m]);
          unsigned int zi= params.z(i);

          unsigned int ni =  (tStop[ind] - tStart[ind] + 1);
          yi.resize(ni);

          for(unsigned int j=0;j<(tStop[ind]-tStart[ind]+1);j++){

            yi(j) = y[ind_y];//yi(j) = y[tStart[ind]-1+j];

            for(unsigned int b=0;b<nFixedEffects[m];b++){
              yi(j)-=params.beta(m,b,0,nCategoriesY)*dataset.W_LME(m,tStart[ind]-1+j,b);
            }
            for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
              yi(j)-=params.beta_mix(m,zi,b,0,nCategoriesY)*dataset.W_LME_mix(m,tStart[ind]-1+j,b);
            }

            ind_y++;
          }

          MatrixXd block=dataset.W_RE(m,tStart[ind]-1, 0, ni, nRandomEffects[m]);
          MatrixXd sigmae=MatrixXd::Identity(ni, ni) * params.SigmaE(m);

          MatrixXd V = block *params.covRE(m,0)* block.transpose() + sigmae;
          LLT<MatrixXd> lltOfA(V); // compute the Cholesky decomposition of A
          MatrixXd L = lltOfA.matrixL();
          //double logDetPrecMat=  2*log(L.determinant());
          MatrixXd Vi_inv = L.inverse().transpose()*L.inverse();
          VectorXd mu = params.covRE(m,0)*block.transpose()*Vi_inv*yi;

          //B - B*Zi^T*Vi^{-1}* (Zi*B^T)
          MatrixXd cov = params.covRE(m,0) - params.covRE(m,0)*block.transpose()*Vi_inv*block*params.covRE(m,0);

          ui = multivarNormalRand(rndGenerator,mu,cov);
          params.RandomEffects(m,i,ui);

          ind ++;
        }

        // for(unsigned int i=0;i<nSubjects;i++){
        //   //int zi = params.z(i);
        //   VectorXd mu(nRandomEffects[m]);
        //   VectorXd meanRE(nRandomEffects[m]);
        //   meanRE.setZero();
        //   MatrixXd covRE=params.covRE(m,0);
        //   mu = multivarNormalRand(rndGenerator,meanRE,covRE);
        //   //mu.setZero();
        //   params.RandomEffects(m,i,mu);
        // }

        //randomChiSquare chisQRand(hyperParams.eps_vu());
        //double temp = chisQRand(rndGenerator);
        //double epsilon = hyperParams.eps_vu()*hyperParams.eps_sigma2_0()/temp;
        //params.SigmaE(epsilon);
        // Define a inverse gamma random number generator

        randomGamma gammaRand(hyperParams.eps_shape(),hyperParams.eps_scale()); //k, 1/theta


        double temp = gammaRand(rndGenerator);
        double epsilon = 1.0/temp;

        params.SigmaE(m,epsilon);
      }
    }
    //RJ populate params.L
    if(outcomeType.compare("Longitudinal")==0){
      unsigned int nL;
      if(kernelType.compare("SQexponential")==0){ //AR
        nL=3;
      }else{
        nL=4;
      }

      if (hyperParams.initL().size()==0){
        for (unsigned int c=0;c<maxNClusters;c++){
          for(unsigned int l=0;l<nL;l++){
            randomNormal normalRand(hyperParams.muL(l),hyperParams.sigmaL(l));
            params.L(c,l,normalRand(rndGenerator));
          }

          if(options.sampleGPmean()){ //AR
            VectorXd Fval(dataset.nTimes_unique());

            Fval =  Sample_GPmean(params, dataset, c, rndGenerator,1);
            for(unsigned int j=0;j<dataset.nTimes_unique();j++){
              params.meanGP(c,j, Fval(j));
            }

            if(options.estim_ratio()){
              double vVal = betaRand(rndGenerator,hyperParams.aRatio(),hyperParams.bRatio());
              params.ratio(c, vVal);
            }
          }
        }
      }else{
        for (unsigned int c=0;c<maxNClusters;c++){
          VectorXd lval =hyperParams.initL(c);
          for(unsigned int l=0;l<nL;l++){
            params.L(c,l,lval(l));
          }
        }
      }
    }
  }

  // And also the extra variation values if necessary
  if(responseExtraVar){
    // Shape and rate parameters
    double a=5,b=2;
    // Boost parameterised in terms of shape and scale
    randomGamma gammaRand(a,1.0/b);
    // Tau is now a Gamma(a,b)
    double tau = gammaRand(rndGenerator);
    params.tauEpsilon(tau);

    randomNormal normalRand(0,1.0/sqrt(tau));

    if(outcomeType.compare("LME")!=0){
      for(unsigned int i=0;i<nSubjects;i++){
        double eps = normalRand(rndGenerator);
        int zi = params.z(i);
        double meanVal = 0;

        meanVal+= params.theta(zi,0);

        if(outcomeType.compare("Categorical")==0){
          for(unsigned int j=0;j<nFixedEffects[0];j++){
            meanVal+=params.beta(0,j,dataset.discreteY(i),nCategoriesY)*dataset.W(i,j);
          }
        } else{ // Bernoulli, Poisson, Binomial, Survivak, LVN, Longitudinal
          for(unsigned int j=0;j<nFixedEffects[0];j++){
            meanVal+=params.beta(0,j,0,nCategoriesY)*dataset.W(i,j);
            // }
            // for(unsigned int j=0;j<nFixedEffects_mix[0];j++){
            //   meanVal+=params.beta_mix(0,zi,j, 0, nCategoriesY)*dataset.W_mix(i,j);
            // }
          }
          if(outcomeType.compare("Poisson")==0){
            meanVal+=dataset.logOffset(i);
          }
          params.lambda(i,meanVal+eps);
        }
      }
    }else{//LME
      for(unsigned int i=0;i<nSubjects;i++)
        params.lambda(i,0);
    }

    // And also _uCAR and _TauCAR if includeCAR==TRUE
    if (includeCAR){
      // Boost parameterised in terms of shape and scale
      randomGamma gammaRand(5.0,0.5);
      // Tau is now a Gamma(shape,rate)
      double tau = gammaRand(rndGenerator);
      params.TauCAR(tau);

      double mean_w=0;
      for (unsigned int i=0; i<nSubjects; i++ ) mean_w+=dataset.nNeighbours(i);
      mean_w /= nSubjects;
      randomNormal normalRand(0,sqrt(mean_w/tau));
      for(unsigned int i=0;i<nSubjects;i++){
        double eps = normalRand(rndGenerator);
        params.uCAR(i,eps);
      }
    }
  }


  // std::cout << " parameters "<<endl;
  // std::cout << " beta "<<endl;
  // for(unsigned int m=0;m<nOutcomes;m++){
  //   for(unsigned int b=0;b<nFixedEffects[m];b++)
  //     std::cout << m << " b "<< b <<  " beta "<< params.beta(m,b,0,nCategoriesY)<<endl;
  //   for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
  //     for(unsigned int c=0;c<maxNClusters;c++)
  //       std::cout << m << " b "<< b <<  " c "<<c<< " betamix "<< params.beta_mix(m,c,b,0, nCategoriesY)<<endl;
  //   }
  //   std::cout<<endl << m << " sigmaE "<< params.SigmaE(m)<<endl<<endl;
  //    for(unsigned int c=0;c<maxNClusters;c++){
  //      std::cout<<m << " covRE "<< params.covRE(m,c)<<endl;
  //      std::cout<<m << " workLogDetTauLME "<<params.workLogDetTauLME(m,c)<<endl;
  //      std::cout<<m << " workSqrtTauLME "<<params.workSqrtTauLME(m,c)<<endl;
  //      }
  //   //params.RandomEffects(m,i)
  // }

  if(2<1 & nOutcomes ==1 ){
    int zi;
    ifstream inputFile;

    string fitFilename = "/Users/naisr/Documents/2022_MCF/code/Applications/Plongi_3C_AXE/Simu/ui2_500.txt";


    inputFile.open(fitFilename.c_str());
    for (unsigned int i=0; i<nSubjects; i++ ){
      if(i<251){
        zi=0;
      }else{
        zi=1;
      }

      params.z(i,zi,covariateType);

      VectorXd ui(2);
      inputFile >>ui(0);
      inputFile >>ui(1);
      params.RandomEffects(0,i,ui);
    }
    inputFile.close();
    params.beta(0,0,0,nCategoriesY,2.45);
    params.beta_mix(0,0,0,0, nCategoriesY, 20.71);
    params.beta_mix(0,0,1,0, nCategoriesY, -1.00 );
    params.beta_mix(0,1,0,0, nCategoriesY, 29.39);
    params.beta_mix(0,1,1,0, nCategoriesY, 0.13);


    MatrixXd Rfixed(nRandomEffects[0],nRandomEffects[0]);
    Rfixed.setZero();
    Rfixed(0,0)=4.50;
    Rfixed(1,0)=0.36;
    Rfixed(0,1)=0.36;
    Rfixed(1,1)=0.79;

    for(unsigned int c=0;c<=maxZ;c++)
      params.covRE(0,c, Rfixed);

    double epsilon=1.0;
    // Define a inverse gamma random number generator
    //randomGamma gammaRand(shape_post, 1/scale_rate); //gammaRand(shape, scale=1/rate)

    params.SigmaE(0,epsilon);// variance

  }else if(2<1 & nOutcomes == 2){

    std::cout<< " init "<<endl;
    int zi;
    ifstream inputFile;

    string fitFilename = "/Users/naisr/Documents/2022_MCF/code/Applications/Plongi_3C_AXE/Simu/ui_M2_G2_R3_500.txt";


    inputFile.open(fitFilename.c_str());
    for (unsigned int i=0; i<nSubjects; i++ ){
      if(i<253){
        zi=0;
      }else{
        zi=1;
      }

      params.z(i,zi,covariateType);

      for (unsigned int m=0; m<nOutcomes; m++ ){
        VectorXd ui(3);
        inputFile >>ui(0);
        inputFile >>ui(1);
        inputFile >>ui(2);
        params.RandomEffects(m,i,ui);
      }
    }
    inputFile.close();
    params.beta(0,0,0,nCategoriesY,1);
    params.beta_mix(0,0,0,0, nCategoriesY, 3);
    params.beta_mix(0,0,1,0, nCategoriesY, -0.2);
    params.beta_mix(0,0,2,0, nCategoriesY, 0.2);

    params.beta_mix(0,1,0,0, nCategoriesY, -2);
    params.beta_mix(0,1,1,0, nCategoriesY, 0.3);
    params.beta_mix(0,1,2,0, nCategoriesY, -0.2);

    params.beta(1,0,0,nCategoriesY,-1);
    params.beta_mix(1,0,0,0, nCategoriesY, -3);
    params.beta_mix(1,0,1,0, nCategoriesY, 0.3 );
    params.beta_mix(1,0,2,0, nCategoriesY, -0.2);

    params.beta_mix(1,1,0,0, nCategoriesY, 2);
    params.beta_mix(1,1,1,0, nCategoriesY, -0.3);
    params.beta_mix(1,1,2,0, nCategoriesY, 0.3);


    MatrixXd Rfixed(nRandomEffects[0],nRandomEffects[0]);
    Rfixed.setZero();
    Rfixed(0,0)=2.5;
    Rfixed(0,1)=0.8;
    Rfixed(0,2)=-0.29;
    Rfixed(1,0)=0.8;
    Rfixed(1,1)=1.1;
    Rfixed(1,2)=0.4;
    Rfixed(2,0)=-0.29;
    Rfixed(2,1)=0.4;
    Rfixed(2,2)=0.5;

    for(unsigned int c=0;c<=maxZ;c++)
      params.covRE(0,c, Rfixed);

    Rfixed.setZero();
    Rfixed(0,0)=1.5;
    Rfixed(0,1)=-0.8;
    Rfixed(0,2)=-0.24;
    Rfixed(1,0)=-0.8;
    Rfixed(1,1)=1.12;
    Rfixed(1,2)=0.01;
    Rfixed(2,0)=-0.24;
    Rfixed(2,1)=0.01;
    Rfixed(2,2)=0.1;

    for(unsigned int c=0;c<=maxZ;c++)
      params.covRE(1,c, Rfixed);

    double epsilon=0.3;
    params.SigmaE(0,epsilon);// variance
    epsilon=0.5;
    params.SigmaE(1,epsilon);// variance

    VectorXd mu(1);
    mu(0)=0;
    params.mu(0,mu);
    mu(0)=3;
    params.mu(1,mu);
    MatrixXd Tau(1,1);
    Tau(0,0)=1.0;
    params.Tau(0,Tau);
    params.Tau(1,Tau);

    //cout << " params.logPhi "<<params.logPhi().size()<<endl;
    vector<double> phi(2);
    phi[0]=log(0.3);
    phi[1]=log(0.7);
    params.logPhi(0,0,phi);
    phi[0]=log(0.7);
    phi[1]=log(0.3);
    params.logPhi(1,0,phi);
  }
}


// Write the sampler output
void writePReMiuMOutput(mcmcSampler<pReMiuMParams,pReMiuMOptions, pReMiuMPropParams,pReMiuMData>& sampler,
                        const unsigned int& sweep){

  bool reportBurnIn = sampler.reportBurnIn();
  unsigned int nBurn = sampler.nBurn();
  unsigned int nFilter = sampler.nFilter();
  vector<ofstream*>& outFiles = sampler.outFiles();

  // Check if we need to do anything
  if((reportBurnIn||((!reportBurnIn)&&sweep>nBurn))&&(sweep%nFilter==0)){
    const pReMiuMParams& params = sampler.chain().currentState().parameters();

    unsigned int nSubjects = sampler.model().dataset().nSubjects();
    unsigned int nOutcomes = sampler.model().dataset().nOutcomes();
    unsigned int nPredictSubjects = params.nPredictSubjects();
    unsigned int maxNClusters = params.maxNClusters();
    unsigned int nCovariates = params.nCovariates();
    unsigned int nDiscreteCovs=params.nDiscreteCovs();
    unsigned int nContinuousCovs=params.nContinuousCovs();
    unsigned int nCategoriesY = params.nCategoriesY();
    string covariateType = sampler.model().dataset().covariateType();
    bool includeResponse = sampler.model().options().includeResponse();
    bool includeCAR = sampler.model().options().includeCAR();
    bool responseExtraVar = sampler.model().options().responseExtraVar();
    double fixedAlpha = sampler.model().options().fixedAlpha();
    string outcomeType = sampler.model().options().outcomeType();
    string kernelType = sampler.model().dataset().kernelType();
    bool computeEntropy = sampler.model().options().computeEntropy();
    vector<unsigned int> nFixedEffects = sampler.model().dataset().nFixedEffects();
    vector<unsigned int> nFixedEffects_mix = sampler.model().dataset().nFixedEffects_mix();
    string varSelectType = sampler.model().options().varSelectType();
    string predictType = sampler.model().options().predictType();
    bool weibullFixedShape = sampler.model().options().weibullFixedShape();

    const pReMiuMData& dataset = sampler.model().dataset();
    pReMiuMPropParams& proposalParams = sampler.proposalParams();

    vector<unsigned int> nCategories;
    if(covariateType.compare("Discrete")==0||covariateType.compare("Mixed")==0){
      nCategories = params.nCategories();
    }

    // Check if the files are already open
    if(outFiles.size()==0){
      string fileStem =sampler.outFileStem();
      string fileName = fileStem + "_nClusters.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));
      fileName = fileStem + "_psi.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));
      if(covariateType.compare("Discrete")==0){
        fileName = fileStem + "_phi.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
      }else if(covariateType.compare("Normal")==0){
        fileName = fileStem + "_mu.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
        fileName = fileStem + "_Sigma.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
      }else if(covariateType.compare("Mixed")==0){
        fileName = fileStem + "_phi.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
        fileName = fileStem + "_mu.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
        fileName = fileStem + "_Sigma.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
      }
      fileName = fileStem + "_z.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));
      fileName = fileStem + "_entropy.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));
      fileName = fileStem + "_alpha.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));
      fileName = fileStem + "_logPost.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));
      fileName = fileStem + "_nMembers.txt";
      outFiles.push_back(new ofstream(fileName.c_str()));

      if(fixedAlpha<=-1){
        fileName = fileStem + "_alphaProp.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
      }
      if(includeResponse){
        if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
          fileName = fileStem + "_theta.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(nFixedEffects[0]>0){
          fileName = fileStem + "_beta.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(nFixedEffects_mix[0]>0){
          fileName = fileStem + "_beta_cluster-specific.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));

        }
        if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
          fileName = fileStem + "_thetaProp.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(nFixedEffects[0]>0){
          fileName = fileStem + "_betaProp.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }

        //fileName = fileStem + "_betamixProp.txt";
        //outFiles.push_back(new ofstream(fileName.c_str()));
        if(outcomeType.compare("Normal")==0){
          fileName = fileStem + "_sigmaSqY.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(outcomeType.compare("Survival")==0){
          fileName = fileStem + "_nu.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        //RJ open _L.txt
        if(outcomeType.compare("Longitudinal")==0){
          fileName = fileStem + "_L.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          if(sampler.model().options().sampleGPmean()){//AR
            fileName = fileStem + "_meanGP.txt";
            outFiles.push_back(new ofstream(fileName.c_str()));
            if(sampler.model().options().estim_ratio()){//AR
              fileName = fileStem + "_ratio_variances.txt";
              outFiles.push_back(new ofstream(fileName.c_str()));
            }
          }
        }
        if(outcomeType.compare("MVN")==0){
          fileName = fileStem + "_MVNmu.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_MVNSigma.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(outcomeType.compare("LME")==0){
          fileName = fileStem + "_cov_RandomEffects_LME.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_epsilonLME.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_RandomEffectsLME.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(responseExtraVar){
          fileName = fileStem + "_epsilon.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_sigmaEpsilon.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_epsilonProp.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(nPredictSubjects>0){
          fileName = fileStem + "_predictThetaRaoBlackwell.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if (includeCAR){
          fileName = fileStem + "_TauCAR.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_uCAR.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
      }
      if(varSelectType.compare("None")!=0){
        fileName = fileStem + "_omega.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
        fileName = fileStem + "_rho.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
        fileName = fileStem + "_rhoOmegaProp.txt";
        outFiles.push_back(new ofstream(fileName.c_str()));
        if(varSelectType.compare("Continuous")!=0){
          fileName = fileStem + "_gamma.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
        if(covariateType.compare("Discrete")==0){
          fileName = fileStem + "_nullPhi.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }else if(covariateType.compare("Normal")==0){
          fileName = fileStem + "_nullMu.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }else if(covariateType.compare("Mixed")==0){
          fileName = fileStem + "_nullPhi.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
          fileName = fileStem + "_nullMu.txt";
          outFiles.push_back(new ofstream(fileName.c_str()));
        }
      }
    }

    // File indices
    int nClustersInd=-1,psiInd=-1,phiInd=-1,muInd=-1,SigmaInd=-1,zInd=-1,entropyInd=-1,alphaInd=-1;
    int logPostInd=-1,nMembersInd=-1,alphaPropInd=-1;
    //LME and Longitudinal
    int thetaInd=-1,betaInd=-1,betamixInd=-1,thetaPropInd=-1,betaPropInd=-1,sigmaSqYInd=-1,nuInd=-1,LInd=-1,
      meanGPInd=-1,estim_ratioInd=-1,MVNmuInd=-1,MVNSigmaInd=-1,epsilonInd=-1,
      CovRELMEInd=-1, EpsilonLMEInd=-1, RandomEffectsLMEInd=-1 ;//LME

    //betamixPropInd=-1,
    int sigmaEpsilonInd=-1,epsilonPropInd=-1,omegaInd=-1,rhoInd=-1;
    int rhoOmegaPropInd=-1,gammaInd=-1,nullPhiInd=-1,nullMuInd=-1;
    int predictThetaRaoBlackwellInd=-1;
    int TauCARInd=-1,uCARInd=-1;

    int r=0;
    nClustersInd=r++;
    psiInd=r++;
    if(covariateType.compare("Discrete")==0){
      phiInd=r++;
    }else if(covariateType.compare("Normal")==0){
      muInd=r++;
      SigmaInd=r++;
    }else if(covariateType.compare("Mixed")==0){
      phiInd=r++;
      muInd=r++;
      SigmaInd=r++;
    }
    zInd=r++;
    entropyInd=r++;
    alphaInd=r++;
    logPostInd=r++;
    nMembersInd=r++;

    if(fixedAlpha<=-1){
      alphaPropInd=r++;
    }
    if(includeResponse){
      if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0)
        thetaInd=r++;
      if(nFixedEffects[0]>0)
        betaInd=r++;
      if(nFixedEffects_mix[0]>0)
        betamixInd=r++;
      if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0)
        thetaPropInd=r++;
      if(nFixedEffects[0]>0)
        betaPropInd=r++;
      if(outcomeType.compare("Normal")==0){
        sigmaSqYInd=r++;
      }
      if(outcomeType.compare("Survival")==0){
        nuInd=r++;
      }
      if(outcomeType.compare("Longitudinal")==0){//RJ
        LInd=r++;
        if(sampler.model().options().sampleGPmean()){
          meanGPInd=r++;//AR
          if(sampler.model().options().estim_ratio())
            estim_ratioInd=r++;//AR
        }
      }
      if(outcomeType.compare("MVN")==0){//RJ
        MVNmuInd=r++;
        MVNSigmaInd=r++;
      }
      if(outcomeType.compare("LME")==0){//AR
        CovRELMEInd=r++;
        EpsilonLMEInd=r++;
        RandomEffectsLMEInd=r++;
      }

      if(responseExtraVar){
        epsilonInd=r++;
        sigmaEpsilonInd=r++;
        epsilonPropInd=r++;
      }
      if(nPredictSubjects>0){
        predictThetaRaoBlackwellInd=r++;
      }
      if (includeCAR){
        TauCARInd=r++;
        uCARInd=r++;
      }
    }

    if(varSelectType.compare("None")!=0){
      omegaInd=r++;
      rhoInd=r++;
      rhoOmegaPropInd=r++;
      if(varSelectType.compare("Continuous")!=0){
        gammaInd=r++;
      }
      if(covariateType.compare("Discrete")==0){
        nullPhiInd=r++;
      }else if (covariateType.compare("Normal")==0){
        nullMuInd=r++;
      }else if (covariateType.compare("Mixed")==0){
        nullPhiInd=r++;
        nullMuInd=r++;
      }
    }


    *(outFiles[nClustersInd]) << maxNClusters << endl;

    unsigned int sumMembers=0;
    for(unsigned int c=0;c<maxNClusters;c++){
      // Print logPsi
      *(outFiles[psiInd]) << exp(params.logPsi(c));
      if(includeResponse){
        // Print theta
        if(outcomeType.compare("Categorical")==0){
          for (unsigned int k=0;k<nCategoriesY;k++){
            *(outFiles[thetaInd]) << params.theta(c,k);
            if (k<(nCategoriesY-1)) {
              *(outFiles[thetaInd]) <<" ";
            }
          }
        } else if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
          //if(outcomeType.compare("Longitudinal")!=0)
          *(outFiles[thetaInd]) << params.theta(c,0);
        }
      }
      // Print number of members of each cluster
      unsigned int nXinC = params.workNXInCluster(c);
      *(outFiles[nMembersInd]) << nXinC;
      sumMembers +=nXinC;
      if(c<maxNClusters-1){
        *(outFiles[psiInd]) << " ";
        if(includeResponse){
          if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0 )
            *(outFiles[thetaInd]) << " ";
        }
        *(outFiles[nMembersInd]) << " ";
      }else{
        *(outFiles[psiInd]) << endl;
        if(includeResponse){
          if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0 )
            *(outFiles[thetaInd]) << endl;
        }
        *(outFiles[nMembersInd]) << " " << sumMembers << endl;
      }
    }


    unsigned int maxNCategories=0;
    if(covariateType.compare("Discrete")==0){
      for(unsigned int j=0;j<nCovariates;j++){
        if(nCategories[j]>maxNCategories){
          maxNCategories=nCategories[j];
        }
      }

      for(unsigned int j=0;j<nCovariates;j++){
        for(unsigned int p=0;p<maxNCategories;p++){
          for(unsigned int c=0;c<maxNClusters;c++){
            // Print Phi
            if(p<nCategories[j]){
              *(outFiles[phiInd]) << exp(params.logPhi(c,j,p));
            }else{
              // pad the output with dummy variables
              // to make reading in R easier
              *(outFiles[phiInd]) << -999;
            }
            if(c<(maxNClusters-1)||p<(maxNCategories-1)||j<(nCovariates-1)){
              *(outFiles[phiInd]) << " ";
            }

          }
        }
      }
      *(outFiles[phiInd]) << endl;
    }else if(covariateType.compare("Normal")==0){
      // To make the output comparable with discrete, we will write the
      // output grouped by covariate (for each cluster)
      for(unsigned int j=0;j<nCovariates;j++){
        for(unsigned int c=0;c<maxNClusters;c++){
          *(outFiles[muInd]) << params.mu(c,j);

          if(c<(maxNClusters-1)||j<(nCovariates-1)){
            *(outFiles[muInd]) << " ";
          }
        }
      }
      *(outFiles[muInd]) << endl;


      // For the covariance matrices we write by covariate x covariate (for each cluster)

      for(unsigned int j1=0;j1<nCovariates;j1++){
        for(unsigned int j2=0;j2<nCovariates;j2++){
          for(unsigned int c=0;c<maxNClusters;c++){
            *(outFiles[SigmaInd]) << params.Sigma(c,j1,j2);
            if(c<(maxNClusters-1)||j1<(nCovariates-1)||j2<(nCovariates-1)){
              *(outFiles[SigmaInd]) << " ";
            }
          }
        }
      }
      *(outFiles[SigmaInd]) << endl;

    }else if(covariateType.compare("Mixed")==0){
      for(unsigned int j=0;j<nDiscreteCovs;j++){
        if(nCategories[j]>maxNCategories){
          maxNCategories=nCategories[j];
        }
      }
      for(unsigned int j=0;j<nDiscreteCovs;j++){
        for(unsigned int p=0;p<maxNCategories;p++){
          for(unsigned int c=0;c<maxNClusters;c++){
            // Print Phi
            if(p<nCategories[j]){
              *(outFiles[phiInd]) << exp(params.logPhi(c,j,p));
            }else{
              // pad the output with dummy variables
              // to make reading in R easier
              *(outFiles[phiInd]) << -999;
            }
            if(c<(maxNClusters-1)||p<(maxNCategories-1)||j<(nDiscreteCovs-1)){
              *(outFiles[phiInd]) << " ";
            }

          }
        }
      }
      *(outFiles[phiInd]) << endl;

      // To make the output comparable with discrete, we will write the
      // output grouped by covariate (for each cluster)
      for(unsigned int j=0;j<nContinuousCovs;j++){
        for(unsigned int c=0;c<maxNClusters;c++){
          *(outFiles[muInd]) << params.mu(c,j);
          if(c<(maxNClusters-1)||j<(nContinuousCovs-1)){
            *(outFiles[muInd]) << " ";
          }
        }
      }
      *(outFiles[muInd]) << endl;


      // For the covariance matrices we write by covariate x covariate (for each cluster)

      for(unsigned int j1=0;j1<nContinuousCovs;j1++){
        for(unsigned int j2=0;j2<nContinuousCovs;j2++){
          for(unsigned int c=0;c<maxNClusters;c++){
            *(outFiles[SigmaInd]) << params.Sigma(c,j1,j2);
            if(c<(maxNClusters-1)||j1<(nContinuousCovs-1)||j2<(nContinuousCovs-1)){
              *(outFiles[SigmaInd]) << " ";
            }
          }
        }
      }
      *(outFiles[SigmaInd]) << endl;
    }


    if(includeResponse){
      if(outcomeType.compare("Categorical")==0){
        // Print beta
        for(unsigned int j=0;j<nFixedEffects[0];j++){
          for (unsigned int k=0;k<nCategoriesY;k++){
            *(outFiles[betaInd]) << params.beta(0,j,k,nCategoriesY) <<" ";
          }
          if(j==(nFixedEffects[0]-1)){
            *(outFiles[betaInd]) << endl;
          }
        }
        // for(unsigned int c=0;c<maxNClusters;c++){
        //   for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //     for (unsigned int k=0;k<nCategoriesY;k++){
        //       *(outFiles[betamixInd]) << params.beta_mix(c, j, k, nCategoriesY) <<" ";
        //     }
        //     if(j==(nFixedEffects_mix-1) && c == (maxNClusters-1)){
        //       *(outFiles[betamixInd]) << endl;
        //     }
        //   }
        // }
      } else {
        if(outcomeType.compare("Normal")==0){
          *(outFiles[sigmaSqYInd]) << params.sigmaSqY() << endl;
        }
        if(outcomeType.compare("Survival")==0){
          // Print parameter nu for each cluster
          if (weibullFixedShape){
            *(outFiles[nuInd]) << params.nu(0) << endl;
          } else {
            for(unsigned int c=0;c<maxNClusters;c++){
              *(outFiles[nuInd]) << params.nu(c);
              if(c<maxNClusters-1){
                *(outFiles[nuInd]) << " ";
              }else{
                *(outFiles[nuInd]) << endl;
              }
            }
          }
        }
        if(outcomeType.compare("Longitudinal")==0){
          //RJ Print parameter L for each cluster
          for(unsigned int c=0;c<maxNClusters;c++){
            unsigned int nL;
            if(kernelType.compare("SQexponential")==0){
              nL=3;
            }else{
              nL=4;
            }
            for(unsigned int l=0;l<nL;l++){
              *(outFiles[LInd]) << params.L(c,l);
              if(c<(maxNClusters-1) || l<(nL-1)){
                *(outFiles[LInd]) << " ";
              }else{
                *(outFiles[LInd]) << endl;
              }
            }

            if(sampler.model().options().sampleGPmean()){ //AR
              for(unsigned int l=0;l<dataset.nTimes_unique();l++){
                *(outFiles[meanGPInd]) << params.meanGP(c,l)<< " ";
                if(c<(maxNClusters-1) || l<(dataset.nTimes_unique()-1)){
                  *(outFiles[meanGPInd]) << " ";
                }else{
                  *(outFiles[meanGPInd]) << endl;
                }
              }

              if(sampler.model().options().estim_ratio()){ //AR
                *(outFiles[estim_ratioInd]) << params.ratio(c);
                if(c<(maxNClusters-1)){
                  *(outFiles[estim_ratioInd]) << " ";
                }else{
                  *(outFiles[estim_ratioInd]) << endl;
                }
              }
            }
          }
        }
        if(outcomeType.compare("MVN")==0){
          //RJ Print MVN parameters for each cluster
          for(unsigned int c=0;c< maxNClusters;c++){
            for(unsigned int l=0;l<nOutcomes;l++){
              *(outFiles[MVNmuInd]) << " "<< params.MVNmu(c,l);
              //*(outFiles[MVNmuInd]) << " "<<sweep << " c "<< c << " l "<< l << " maxNClusters " << maxNClusters << " nOutcomes " << nOutcomes << " p "<< params.MVNmu(c,l)<<endl;//RJ!!params.MVNmu(c,l);
              if(c<maxNClusters-1 || l<(nOutcomes-1)){
                *(outFiles[MVNmuInd]) << " ";
              }else{
                *(outFiles[MVNmuInd]) << endl;
              }
              for(unsigned int l2=0;l2<=l;l2++){
                *(outFiles[MVNSigmaInd]) << ""<<params.MVNSigma(c,l,l2);//RJ!!params.MVNSigma(c,l,l2);
                if(c<maxNClusters-1 || l2<(nOutcomes-1)){
                  *(outFiles[MVNSigmaInd]) << " ";
                }else{
                  *(outFiles[MVNSigmaInd]) << endl;
                }
              }
            }
          }
        }

        if(outcomeType.compare("LME")==0){
          vector<unsigned int> nRandomEffects = dataset.nRandomEffects();

          for(unsigned int m=0;m<nOutcomes;m++){
            //for(unsigned int c=0;c< maxNClusters;c++){
              for(unsigned int l=0;l<nRandomEffects[m];l++){
                for(unsigned int l2=0;l2<=l;l2++){
                  *(outFiles[CovRELMEInd]) << ""<< params.covRE(m,0,l,l2);
                  if( l2<(nRandomEffects[m]-1)){
                    *(outFiles[CovRELMEInd]) << " ";
                  }else{
                    *(outFiles[CovRELMEInd]) << endl;
                  }
                }
              }
            //}

            *(outFiles[EpsilonLMEInd]) << params.SigmaE(m) << endl; //params.sigmakInd(c);
          }


          for(unsigned int m=0;m<nOutcomes;m++){
            for(unsigned int i=0;i<nSubjects;i++){
              *(outFiles[RandomEffectsLMEInd]) << " "<< params.RandomEffects(m,i).transpose()<< " "; //params.sigmakInd(c);
            }
            *(outFiles[RandomEffectsLMEInd]) << endl;
          }

          for(unsigned int m=0;m<nOutcomes;m++){
            for(unsigned int j=0;j<nFixedEffects[m];j++)
              *(outFiles[betaInd]) << params.beta(m,j,0,nCategoriesY)<< " ";
            *(outFiles[betaInd]) << endl;
          }

          for(unsigned int m=0;m<nOutcomes;m++){
            for(unsigned int c=0;c< maxNClusters;c++){
              for(unsigned int j=0;j<nFixedEffects_mix[m];j++)
                *(outFiles[betamixInd]) << params.beta_mix(m,c,j, 0, nCategoriesY)<< " "; //j*1 (Ycategory)
             }
            *(outFiles[betamixInd]) << endl;
          }
        }

        if (includeCAR){
          for(unsigned int i=0;i<nSubjects;i++){
            double uCARi = params.uCAR(i);
            *(outFiles[uCARInd]) << uCARi;
            if(i<nSubjects-1){
              *(outFiles[uCARInd]) << " ";
            }else{
              *(outFiles[uCARInd]) << endl;
            }
          }
          double tau=params.TauCAR();
          *(outFiles[TauCARInd]) << tau;
          *(outFiles[TauCARInd]) << endl;
        }
      }
    }

    // Print the allocations
    for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
      *(outFiles[zInd]) << params.z(i);
      if(computeEntropy){
        *(outFiles[entropyInd]) << params.workEntropy(i);
      }
      if(i<nSubjects+nPredictSubjects-1){
        if(computeEntropy){
          *(outFiles[entropyInd]) << " ";
        }
        *(outFiles[zInd]) << " ";
      }else{
        if(computeEntropy){
          *(outFiles[entropyInd]) << endl;
        }
        *(outFiles[zInd]) << endl;
      }
      // And print the expected theta for the prediction subjects
      if(i>=nSubjects){
        if(includeResponse){
          for (unsigned int k=0;k<nCategoriesY;k++){
            *(outFiles[predictThetaRaoBlackwellInd]) << params.workPredictExpectedTheta(i-nSubjects,k)<<" ";
          }
          if(i<nSubjects+nPredictSubjects-1){
            *(outFiles[predictThetaRaoBlackwellInd]) << " ";
          }else{
            *(outFiles[predictThetaRaoBlackwellInd]) << endl;
          }
        }
      }
    }

    // Print alpha
    if(fixedAlpha<=-1||sweep==0){
      *(outFiles[alphaInd]) << params.alpha() << endl;
    }

    // Print the log posterior
    *(outFiles[logPostInd]) << sampler.chain().currentState().logPosterior() << " ";
    *(outFiles[logPostInd]) << sampler.chain().currentState().logLikelihood() << " ";
    *(outFiles[logPostInd]) << sampler.chain().currentState().logPrior() << endl;

    bool anyUpdates;
    if(includeResponse){
      // Print the acceptance rates for theta
      anyUpdates = proposalParams.thetaAnyUpdates();
      if(anyUpdates && outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
        *(outFiles[thetaPropInd]) << sampler.proposalParams().thetaAcceptRate() <<
          " " << sampler.proposalParams().thetaStdDev() << endl;
        proposalParams.thetaAnyUpdates(false);
      }


      // Print the acceptance rates for beta
      anyUpdates = proposalParams.betaAnyUpdates();
      if(anyUpdates){//yModel-=LME
        for(unsigned int j=0;j<nFixedEffects[0];j++){
          *(outFiles[betaPropInd]) << sampler.proposalParams().betaAcceptRate(j) <<
            " " << sampler.proposalParams().betaStdDev(j);
          if(j<(nFixedEffects[0]-1)){
            *(outFiles[betaPropInd]) << endl;
          }
          *(outFiles[betaPropInd]) << endl;
          proposalParams.betaAnyUpdates(false);
        }
        // for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //   //			    for(unsigned int c=0;c< maxNClusters;c++){
        //
        //   *(outFiles[betamixPropInd]) << sampler.proposalParams().betamixAcceptRate(j) <<
        //     " " << sampler.proposalParams().betamixStdDev(j);
        //   if(j<(nFixedEffects_mix-1)){
        //     *(outFiles[betamixPropInd]) << endl;
        //   }
        //   *(outFiles[betaPropInd]) << endl;
        //   proposalParams.betaAnyUpdates(false);
        // }

        if(responseExtraVar){

          if(outcomeType.compare("LME")!=0){
            vector<double> meanVec(nSubjects,0.0);
            if(outcomeType.compare("Poisson")==0){
              meanVec=dataset.logOffset();
            }
            for(unsigned int i=0;i<nSubjects;i++){
              int zi = params.z(i);
              double meanVal = meanVec[i]+params.theta(zi,0);

              for(unsigned int j=0;j<nFixedEffects[0];j++)
                meanVal+=params.beta(0,j,0,nCategoriesY)*dataset.W(i,j);


              double eps=params.lambda(i)-meanVal;
              *(outFiles[epsilonInd]) << eps;
              if(i<nSubjects-1){
                *(outFiles[epsilonInd]) << " ";
              }else{
                *(outFiles[epsilonInd]) << endl;
              }
            }
            anyUpdates = proposalParams.lambdaAnyUpdates();
            if(anyUpdates){
              *(outFiles[epsilonPropInd]) << sampler.proposalParams().lambdaAcceptRate() <<
                " " << sampler.proposalParams().lambdaStdDev() << endl;
              proposalParams.lambdaAnyUpdates(false);
            }
            *(outFiles[sigmaEpsilonInd]) << 1.0/sqrt(params.tauEpsilon()) << endl;
          }else{
            printf("responseExtraVar not possible with yModel = LME.");
          }
        }
      }
    }

    // Print the acceptance rates for alpha
    if(fixedAlpha<=-1){
      anyUpdates = proposalParams.alphaAnyUpdates();
      if(anyUpdates){
        *(outFiles[alphaPropInd]) << sampler.proposalParams().alphaAcceptRate() <<
          " " << sampler.proposalParams().alphaStdDev() << endl;
        proposalParams.alphaAnyUpdates(false);
      }
    }
    if(varSelectType.compare("None")!=0){
      // Print variable selection related quantities
      for(unsigned int j=0;j<nCovariates;j++){
        *(outFiles[omegaInd]) << params.omega(j);
        *(outFiles[rhoInd]) << params.rho(j);
        if(j<nCovariates-1){
          *(outFiles[omegaInd]) << " ";
          *(outFiles[rhoInd]) << " ";

        }else{
          *(outFiles[omegaInd]) << endl;
          *(outFiles[rhoInd]) << endl;
        }
        if(sweep!=0){ // this was "==0", it might be worth double checking
          if(covariateType.compare("Discrete")==0){
            for(unsigned int p=0;p<maxNCategories;p++){
              if(p<nCategories[j]){
                *(outFiles[nullPhiInd]) << exp(params.logNullPhi(j,p));
              }else{
                *(outFiles[nullPhiInd]) << -999;
              }

              if(p<(maxNCategories-1)||j<(nCovariates-1)){
                *(outFiles[nullPhiInd]) << " ";
              }else{
                *(outFiles[nullPhiInd]) << endl;
              }

            }
          }else if(covariateType.compare("Normal")==0){
            *(outFiles[nullMuInd]) << params.nullMu(j);
            if(j<nCovariates-1){
              *(outFiles[nullMuInd]) << " ";
            }else{
              *(outFiles[nullMuInd]) << endl;
            }

          }else if(covariateType.compare("Mixed")==0){
            if (j < nDiscreteCovs){
              for(unsigned int p=0;p<maxNCategories;p++){
                if(p<nCategories[j]){
                  *(outFiles[nullPhiInd]) << exp(params.logNullPhi(j,p));
                }else{
                  *(outFiles[nullPhiInd]) << -999;
                }
                if(p<(maxNCategories-1)||j<(nDiscreteCovs-1)){
                  *(outFiles[nullPhiInd]) << " ";
                }else{
                  *(outFiles[nullPhiInd]) << endl;
                }
              }
            } else {
              *(outFiles[nullMuInd]) << params.nullMu(j-nDiscreteCovs);
              if(j<nCovariates-1){
                *(outFiles[nullMuInd]) << " ";
              }else{
                *(outFiles[nullMuInd]) << endl;
              }
            }
          }
        }
        if(varSelectType.compare("BinaryCluster")==0){
          for(unsigned int c=0;c<maxNClusters;c++){
            *(outFiles[gammaInd]) << params.gamma(c,j);
            if(c<maxNClusters-1||j<nCovariates-1){
              *(outFiles[gammaInd]) << " ";
            }else{
              *(outFiles[gammaInd]) << endl;
            }
          }
        }
      }

      anyUpdates = proposalParams.rhoAnyUpdates();
      if(anyUpdates){
        for(unsigned int j=0;j<nCovariates;j++){
          *(outFiles[rhoOmegaPropInd]) << sampler.proposalParams().rhoAcceptRate(j) <<
            " " << sampler.proposalParams().rhoStdDev(j);
          if(j<(nCovariates-1)){
            *(outFiles[rhoOmegaPropInd]) << " ";
          }else{
            *(outFiles[rhoOmegaPropInd]) << endl;
          }
        }
        proposalParams.rhoAnyUpdates(false);

      }
    }
  }
}

string storeLogFileData(const pReMiuMOptions& options,
                        const pReMiuMData& dataset,
                        const pReMiuMHyperParams& hyperParams,
                        const unsigned int& nClusInit,
                        const unsigned int& maxNClusters,
                        const double& timeInSecs){

  ostringstream tmpStr;
  tmpStr << "Number of subjects: " << dataset.nSubjects() << endl;
  tmpStr << "Number of prediction subjects: " << dataset.nPredictSubjects() << endl;
  tmpStr << "Prediction type: " << options.predictType() << endl;
  tmpStr << "Sampler type: " << options.samplerType();
  if(options.samplerType().compare("Truncated")==0){
    tmpStr << " " << maxNClusters << " clusters" << endl;
  }else{
    tmpStr << endl;
  }
  tmpStr << "Number of initial clusters: " << nClusInit;
  if(options.nClusInit()==0){
    tmpStr << " (Random, Unif[50,60])" << endl;
  }else{
    tmpStr << endl;
  }
  tmpStr << "Covariates: " ;//<< endl;
  if(options.covariateType().compare("Mixed")==0){
    tmpStr << "Number of discrete covariates: " << dataset.nDiscreteCovs() << endl;
    tmpStr << "Number of continuous covariates: " << dataset.nContinuousCovs() << endl;
  }
  for(unsigned int j=0;j<dataset.nCovariates();j++){
    tmpStr << "\t" << dataset.covariateNames(j);
    if(options.covariateType().compare("Discrete")==0){
      tmpStr << " (categorical)";
    } else if(options.covariateType().compare("Mixed")==0){
      if (j < dataset.nDiscreteCovs()){
        tmpStr << " (categorical)";
      }
    }
    tmpStr << endl;
  }
  if(dataset.nFixedEffects(0)>0){
    tmpStr << "FixedEffects: " ;//<< endl;
    for(unsigned int m=0;m<dataset.nOutcomes();m++){
      for(unsigned int j=0;j<dataset.nFixedEffects(m);j++){
        tmpStr << "\t" << dataset.fixedEffectNames(m,j) << " ";
      }
    }
    tmpStr << endl;
  }else{
    tmpStr<< "No fixed effects" << endl;
  }
  if(dataset.nFixedEffects_mix(0)>0){
    tmpStr << "Cluster-specific FixedEffects: " ;//<< endl;
    for(unsigned int m=0;m<dataset.nOutcomes();m++){
      for(unsigned int j=0;j<dataset.nFixedEffects_mix(m);j++){
        tmpStr << "\t" << dataset.fixedEffectNames_mix(m,j) << " ";
      }
    }
    tmpStr << endl;
  }else{
    tmpStr<< "No cluster-specific fixed effects" << endl;
  }
  if(dataset.nCategoriesY()>1){
    tmpStr << "NumberOfCategoriesY: " <<  dataset.nCategoriesY() << endl;
  }
  tmpStr << "Model for Y: " << options.outcomeType() << endl;
  if(options.responseExtraVar()){
    tmpStr << "Extra Y variance: True" << endl;
  }else{
    tmpStr << "Extra Y variance: False" << endl;
  }
  if(options.includeResponse()){
    tmpStr << "Include response: True" << endl;
  }else{
    tmpStr << "Include response: False" << endl;
  }
  if(dataset.outcomeType().compare("Longitudinal")==0){//AR
    tmpStr << "GP Kernel type:" << dataset.kernelType() << endl;
    if(options.sampleGPmean()){
      tmpStr << "GP mean sampling: True" << endl;
      if(options.estim_ratio())
        tmpStr << "Estimation ratio between L1k and L3k: True" << endl;
    }
  }

  if(options.outcomeType().compare("LME")==0){
    tmpStr << "RandomEffects: " ;
    for(unsigned int m=0;m<dataset.nOutcomes();m++){
      for(unsigned int j=0;j<dataset.nRandomEffects(m);j++){
        tmpStr << "\t" << dataset.randomEffectNames(m,j) << " ";
      }
      tmpStr << endl;
    }
  }
  if(options.fixedAlpha()<=-1){
    tmpStr << "Update alpha: True" << endl;
  }else{
    tmpStr << "Update alpha: False" << endl;
    tmpStr << "Fixed alpha: " << options.fixedAlpha() << endl;
    if(options.dPitmanYor()==0) {
      tmpStr << "Dirichlet process prior, so dPitmanYor: " << options.dPitmanYor() << endl;
    } else {
      tmpStr << "dPitmanYor: " << options.dPitmanYor() << endl;
    }
  }
  if(options.computeEntropy()){
    tmpStr << "Compute allocation entropy: True" << endl;
  }else{
    tmpStr << "Compute allocation entropy: False" << endl;
  }

  tmpStr << "Model for X: " << options.covariateType() << endl;
  tmpStr << "Variable selection: " << options.varSelectType() << endl;

  tmpStr << endl << "Hyperparameters:" << endl;
  if(options.fixedAlpha()<=-1){
    tmpStr << "shapeAlpha: " << hyperParams.shapeAlpha() << endl;
    tmpStr << "rateAlpha: " << hyperParams.rateAlpha() << endl;
  }
  if(options.covariateType().compare("Discrete")==0 ||options.covariateType().compare("Mixed")==0 ){
    tmpStr << "aPhi[j]: ";
    if(options.covariateType().compare("Discrete")==0){
      for(unsigned int j=0;j<dataset.nCovariates();j++){
        tmpStr << hyperParams.aPhi(j) << " ";
      }
    }
    if(options.covariateType().compare("Mixed")==0){
      for(unsigned int j=0;j<dataset.nDiscreteCovs();j++){
        tmpStr << hyperParams.aPhi(j) << " ";
      }
    }
    tmpStr << endl;
  }

  if(options.covariateType().compare("Normal")==0 ||options.covariateType().compare("Mixed")==0 ){
    tmpStr << "mu0: " << endl;
    tmpStr << hyperParams.mu0() << endl;
    tmpStr << "Tau0:" << endl;
    tmpStr << hyperParams.Tau0() << endl;
    tmpStr << "R0: "  << endl;
    tmpStr << hyperParams.R0() << endl;
    tmpStr << "kappa0: " << hyperParams.kappa0() << endl;
    tmpStr << "nu0: " << hyperParams.nu0() << endl;
  }

  tmpStr << "muTheta: " << hyperParams.muTheta() << endl;
  tmpStr << "sigmaTheta: " << hyperParams.sigmaTheta() << endl;
  tmpStr << "dofTheta: " << hyperParams.dofTheta() << endl;

  if(dataset.nFixedEffects(0)>0|| dataset.nFixedEffects_mix(0)>0){
    tmpStr << "muBeta: " << hyperParams.muBeta() << endl;
    tmpStr << "sigmaBeta: " << hyperParams.sigmaBeta() << endl;
    tmpStr << "dofBeta: " << hyperParams.dofBeta() << endl;
  }

  if(options.outcomeType().compare("LME")==0){
    tmpStr << " cov RE, Sigma R0: ";
    if(dataset.nRandomEffects(0)>1)
      tmpStr<<endl;
    tmpStr << hyperParams.SigmaLME_R0(0)<<endl;
    tmpStr << " cov RE, kappa0: "<< hyperParams.SigmaLME_kappa0(0)<< endl;

    //tmpStr << "eps_vu: "<< hyperParams.eps_vu() << endl;
    //tmpStr << "eps_sigma2_0: "<< hyperParams.eps_sigma2_0() <<endl;
    tmpStr << "eps_shape: "<< hyperParams.eps_shape() << endl;
    tmpStr << "eps_scale: "<< hyperParams.eps_scale() << endl;
  }

  if(options.responseExtraVar()){
    tmpStr << "shapetauEpsilon: " << hyperParams.shapeTauEpsilon() << endl;
    tmpStr << "ratetauEpsilon: " << hyperParams.rateTauEpsilon() << endl;
  }

  if(options.varSelectType().compare("None")!=0){
    tmpStr << "aRho: " << hyperParams.aRho() << endl;
    tmpStr << "bRho: " << hyperParams.bRho() << endl;
    tmpStr << "atomRho: " << hyperParams.atomRho() << endl;
  }

  if(dataset.outcomeType().compare("Normal")==0){
    tmpStr << "shapeSigmaSqY: " << hyperParams.shapeSigmaSqY() << endl;
    tmpStr << "scaleSigmaSqY: " << hyperParams.scaleSigmaSqY() << endl;
  }
  if(dataset.outcomeType().compare("Survival")==0){
    tmpStr << "Weibull with fixed shape parameter: " << options.weibullFixedShape() << endl;
    tmpStr << "shapeNu: " << hyperParams.shapeNu() << endl;
    tmpStr << "scaleNu: " << hyperParams.scaleNu() << endl;
  }
  //RJ write hypers to file
  if(dataset.outcomeType().compare("Longitudinal")==0){
    tmpStr << "muL: " ;
    unsigned int nL;
    if(dataset.kernelType().compare("SQexponential")==0){
      nL=3;
    }else{
      nL=4;
    }
    for(int l=0; l<nL; l++)
      tmpStr << hyperParams.muL(l)  << " ";
    tmpStr << endl;
    tmpStr << "sigmaL: " ;
    for(int l=0; l<nL; l++)
      tmpStr <<  hyperParams.sigmaL(l) << " ";
    tmpStr << endl;
    if(options.estim_ratio()){
      tmpStr << "ratio a: " << hyperParams.aRatio()<< " ";
      tmpStr << "ratio b: "<< hyperParams.bRatio()<< endl;
    }
  }
  if(options.outcomeType().compare("MVN")==0 ){
    tmpStr << "MVNmu0: " << endl;
    tmpStr << hyperParams.MVNmu0() << endl;
    tmpStr << "MVNTau0:" << endl;
    tmpStr << hyperParams.MVNTau0() << endl;
    tmpStr << "MVNR0: "  << endl;
    tmpStr << hyperParams.MVNR0() << endl;
    tmpStr << "MVNkappa0: " << hyperParams.MVNkappa0() << endl;
    tmpStr << "MVNnu0: " << hyperParams.MVNnu0() << endl;
  }

  if(options.includeCAR()){
    tmpStr << "shapeTauCAR: " << endl;
    tmpStr << hyperParams.shapeTauCAR() << endl;
    tmpStr << "rateTauCAR:" << endl;
    tmpStr << hyperParams.rateTauCAR() << endl;
  }

  tmpStr << endl << options.nSweeps()+options.nBurn() << " sweeps done in " <<
    timeInSecs << " seconds" << endl;

  return tmpStr.str();
}

#endif // DIPBACIO_H_
