/// \file PReMiuMProposals.h
/// \author David Hastie
/// \brief Header file for model specification for PReMiuMpp

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

#ifndef DIPBACPROPOSALS_H_
#define DIPBACPROPOSALS_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<numeric>
#include<limits>

#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions/beta.hpp>
#include<boost/math/distributions/students_t.hpp>
#include<boost/math/special_functions/gamma.hpp>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>


// Custom includes
#include<MCMC/chain.h>
#include<MCMC/model.h>
#include<Math/random.h>
#include<Math/distribution.h>
#include<PReMiuMOptions.h>
#include<PReMiuMModel.h>
#include<PReMiuMData.h>
#include<PReMiuMArs.h>

using namespace Eigen;

using std::vector;
using std::accumulate;
using std::numeric_limits;
using boost::math::normal_distribution;
using boost::math::students_t_distribution;
using boost::math::lgamma;


class pReMiuMPropParams{

public:
  // Default constructor
  pReMiuMPropParams() {};

  pReMiuMPropParams(const unsigned int& nSweeps, const unsigned int& nCovariates,//const unsigned int& nOutcomes,
                    const vector<unsigned int>& nFixedEffects, //const vector<unsigned int>& nFixedEffects_mix,
                    const unsigned int& nCategoriesY, const string& kernelType){
    _thetaStdDev=1.0;
    _thetaStdDevLower=0.1;
    _thetaStdDevUpper=99.9;
    _nTryTheta=0;
    _nAcceptTheta=0;
    _nLocalAcceptTheta=0;
    _nResetTheta=0;
    _thetaAcceptTarget=0.44;
    _thetaUpdateFreq=25;
    _thetaAnyUpdates=true;

    _nTryBeta.resize(nFixedEffects[0]);
    _nAcceptBeta.resize(nFixedEffects[0]);
    _nLocalAcceptBeta.resize(nFixedEffects[0]);
    _nResetBeta.resize(nFixedEffects[0]);
    _betaStdDev.resize(nFixedEffects[0]);
    _betaStdDevLower.resize(nFixedEffects[0]);
    _betaStdDevUpper.resize(nFixedEffects[0]);
    for(unsigned int j=0;j<nFixedEffects[0];j++){
      _betaStdDev[j]=1.0;
      _betaStdDevLower[j]=0.1;
      _betaStdDevUpper[j]=99.9;
      _nTryBeta[j]=0;
      _nAcceptBeta[j]=0;
      _nLocalAcceptBeta[j]=0;
      _nResetBeta[j]=0;
    }
    _betaAcceptTarget = 0.44;
    _betaUpdateFreq = 25;
    _betaAnyUpdates=true;

    // _nTryBetamix.resize(nFixedEffects_mix);
    // _nAcceptBetamix.resize(nFixedEffects_mix);
    // _nLocalAcceptBetamix.resize(nFixedEffects_mix);
    // _nResetBetamix.resize(nFixedEffects_mix);
    // _betamixStdDev.resize(nFixedEffects_mix);
    // _betamixStdDevLower.resize(nFixedEffects_mix);
    // _betamixStdDevUpper.resize(nFixedEffects_mix);
    // for(unsigned int j=0;j<nFixedEffects_mix;j++){
    //   _betamixStdDev[j]=1.0;
    //   _betamixStdDevLower[j]=0.1;
    //   _betamixStdDevUpper[j]=99.9;
    //   _nTryBetamix[j]=0;
    //   _nAcceptBetamix[j]=0;
    //   _nLocalAcceptBetamix[j]=0;
    //   _nResetBetamix[j]=0;
    // }
    _betamixAcceptTarget = 0.44;
    _betamixUpdateFreq = 25;
    _betamixAnyUpdates=true;

    //AR for LME outcome
    _SigmaEStdDev=1.0;
    _SigmaEStdDevLower=0.1;
    _SigmaEStdDevUpper=99.9;
    _nTrySigmaE=0;
    _nAcceptSigmaE=0;
    _nLocalAcceptSigmaE=0;
    _nResetSigmaE=0;
    _SigmaEAcceptTarget=0.44;
    _SigmaEUpdateFreq=25;
    _SigmaEAnyUpdates=true;



    //RJ resize L step values
    unsigned int nL;

    if(kernelType.compare("SQexponential")==0){
      nL=3;
    }else{
      nL=4;
    }
    _nTryL.resize(nL);
    _nAcceptL.resize(nL);
    _nLocalAcceptL.resize(nL);
    _nResetL.resize(nL);
    _LStdDev.resize(nL);
    _LStdDevLower.resize(nL);
    _LStdDevUpper.resize(nL);

    for(unsigned int j=0;j<nL;j++){
      _LStdDev[j]=1.0;
      _LStdDevLower[j]=0.1;
      _LStdDevUpper[j]=99.9;
      _nTryL[j]=0;
      _nAcceptL[j]=0;
      _nLocalAcceptL[j]=0;
      _nResetL[j]=0;
    }
    _LAcceptTarget = 0.44;
    _LUpdateFreq = 25;
    _LAnyUpdates=true;

    _alphaStdDev=1.0;
    _alphaStdDevLower=0.1;
    _alphaStdDevUpper=99.9;
    _nTryAlpha=0;
    _nAcceptAlpha=0;
    _nLocalAcceptAlpha=0;
    _nResetAlpha=0;
    _alphaAcceptTarget = 0.44;
    _alphaUpdateFreq = 25;
    _alphaAnyUpdates=true;


    _nTryRho.resize(nCovariates);
    _nAcceptRho.resize(nCovariates);
    _nLocalAcceptRho.resize(nCovariates);
    _nResetRho.resize(nCovariates);
    _rhoStdDev.resize(nCovariates);
    _rhoStdDevLower.resize(nCovariates);
    _rhoStdDevUpper.resize(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      _rhoStdDev[j]=0.5;
      _rhoStdDevLower[j]=0.0001;
      _rhoStdDevUpper[j]=9.9999;
      _nTryRho[j]=0;
      _nAcceptRho[j]=0;
      _nLocalAcceptRho[j]=0;
      _nResetRho[j]=0;
    }
    _rhoAcceptTarget = 0.44;
    _rhoUpdateFreq = 10;
    _rhoAnyUpdates=true;


    _lambdaStdDev=1.0;
    _lambdaStdDevLower=0.1;
    _lambdaStdDevUpper=99.9;
    _nTryLambda=0;
    _nAcceptLambda=0;
    _nLocalAcceptLambda=0;
    _nResetLambda=0;
    _lambdaAcceptTarget = 0.44;
    _lambdaUpdateFreq = 500;
    _lambdaAnyUpdates=true;

  };

  ~pReMiuMPropParams(){};


  unsigned int nTryTheta() const{
    return _nTryTheta;
  }

  unsigned int nAcceptTheta() const{
    return _nAcceptTheta;
  }

  double thetaAcceptRate() const{
    if(_nTryTheta>0){
      return (double)_nAcceptTheta/(double)_nTryTheta;
    }else{
      return 0.0;
    }
  }

  unsigned int thetaUpdateFreq() const{
    return _thetaUpdateFreq;
  }

  unsigned int nLocalAcceptTheta() const{
    return _nLocalAcceptTheta;
  }

  double thetaLocalAcceptRate() const{
    return (double)_nLocalAcceptTheta/(double)_thetaUpdateFreq;
  }


  double thetaAcceptTarget() const{
    return _thetaAcceptTarget;
  }

  void thetaAddTry(){
    _nTryTheta++;
  }

  void thetaAddAccept(){
    _nAcceptTheta++;
    _nLocalAcceptTheta++;
  }

  void thetaLocalReset(){
    _nLocalAcceptTheta=0;
  }

  unsigned int nResetTheta() const{
    return _nResetTheta;
  }

  void thetaStdDevReset(){
    _thetaStdDev = 1.0;
    _nResetTheta++;
    _thetaStdDevLower = pow(10.0,-((double)_nResetTheta+1.0));
    _thetaStdDevUpper = 100.0-pow(10.0,-((double)_nResetTheta+1.0));
  }

  double& thetaStdDev(){
    return _thetaStdDev;
  }

  double thetaStdDev() const{
    return _thetaStdDev;
  }

  // Member function for setting the standard deviation for
  // proposal for theta
  void thetaStdDev(const double& sd){
    _thetaStdDev=sd;
  }

  double thetaStdDevLower() const{
    return _thetaStdDevLower;
  }

  double thetaStdDevUpper() const{
    return _thetaStdDevUpper;
  }

  bool thetaAnyUpdates() const{
    return _thetaAnyUpdates;
  }

  void thetaAnyUpdates(const bool& newStatus){
    _thetaAnyUpdates = newStatus;
  }


  //AR LME outcome option
  unsigned int nTrySigmaE() const{
    return _nTrySigmaE;
  }

  unsigned int nAcceptSigmaE() const{
    return _nAcceptSigmaE;
  }

  double SigmaEAcceptRate() const{
    if(_nTrySigmaE>0){
      return (double)_nAcceptSigmaE/(double)_nTrySigmaE;
    }else{
      return 0.0;
    }
  }

  unsigned int SigmaEUpdateFreq() const{
    return _SigmaEUpdateFreq;
  }

  unsigned int nLocalAcceptSigmaE() const{
    return _nLocalAcceptSigmaE;
  }

  double SigmaELocalAcceptRate() const{
    return (double)_nLocalAcceptSigmaE/(double)_SigmaEUpdateFreq;
  }


  double SigmaEAcceptTarget() const{
    return _SigmaEAcceptTarget;
  }

  void SigmaEAddTry(){
    _nTrySigmaE++;
  }

  void SigmaEAddAccept(){
    _nAcceptSigmaE++;
    _nLocalAcceptSigmaE++;
  }

  void SigmaELocalReset(){
    _nLocalAcceptSigmaE=0;
  }

  unsigned int nResetSigmaE() const{
    return _nResetSigmaE;
  }

  void SigmaEStdDevReset(){
    _SigmaEStdDev = 1.0;
    _nResetSigmaE++;
    _SigmaEStdDevLower = pow(10.0,-((double)_nResetSigmaE+1.0));
    _SigmaEStdDevUpper = 100.0-pow(10.0,-((double)_nResetSigmaE+1.0));
  }

  double& SigmaEStdDev(){
    return _SigmaEStdDev;
  }

  double SigmaEStdDev() const{
    return _SigmaEStdDev;
  }

  // Member function for setting the standard deviation for
  // proposal for SigmaE
  void SigmaEStdDev(const double& sd){
    _SigmaEStdDev=sd;
  }

  double SigmaEStdDevLower() const{
    return _SigmaEStdDevLower;
  }

  double SigmaEStdDevUpper() const{
    return _SigmaEStdDevUpper;
  }

  bool SigmaEAnyUpdates() const{
    return _SigmaEAnyUpdates;
  }

  void SigmaEAnyUpdates(const bool& newStatus){
    _SigmaEAnyUpdates = newStatus;
  }


  vector<unsigned int> nTryBeta() const{
    return _nTryBeta;
  }

  unsigned int nTryBeta(const unsigned int& j) const{
    return _nTryBeta[j];
  }

  vector<unsigned int> nAcceptBeta() const{
    return _nAcceptBeta;
  }

  double betaAcceptRate(const unsigned int& j) const{
    if(_nTryBeta[j]>0){
      return (double)_nAcceptBeta[j]/(double)_nTryBeta[j];
    }else{
      return 0.0;
    }
  }

  unsigned int betaUpdateFreq() const{
    return _betaUpdateFreq;
  }

  vector<unsigned int> nLocalAcceptBeta() const{
    return _nLocalAcceptBeta;
  }

  double betaLocalAcceptRate(const unsigned int& j) const{
    return (double)_nLocalAcceptBeta[j]/(double)_betaUpdateFreq;
  }

  double betaAcceptTarget() const{
    return _betaAcceptTarget;
  }

  void betaAddTry(const unsigned int& j){
    _nTryBeta[j]++;
  }

  void betaAddAccept(const unsigned int& j){
    _nAcceptBeta[j]++;
    _nLocalAcceptBeta[j]++;
  }

  void betaLocalReset(const unsigned int& j){
    _nLocalAcceptBeta[j]=0;
  }

  vector<unsigned int> nResetBeta() const{
    return _nResetBeta;
  }

  void betaStdDevReset(const unsigned int& j){
    _betaStdDev[j] = 1.0;
    _nResetBeta[j]++;
    _betaStdDevLower[j] = pow(10.0,-((double)_nResetBeta[j]+1.0));
    _betaStdDevUpper[j] = 100.0-pow(10.0,-((double)_nResetBeta[j]+1.0));
  }

  vector<double> betaStdDev() const{
    return _betaStdDev;
  }

  vector<double>& betaStdDev(){
    return _betaStdDev;
  }

  double& betaStdDev(const unsigned int& j){
    return _betaStdDev[j];
  }

  const double& betaStdDev(const unsigned int& j) const{
    return _betaStdDev[j];
  }

  vector<double> betaStdDevLower() const{
    return _betaStdDevLower;
  }

  double betaStdDevLower(const unsigned int& j) const{
    return _betaStdDevLower[j];
  }

  vector<double> betaStdDevUpper() const{
    return _betaStdDevUpper;
  }

  double betaStdDevUpper(const unsigned int& j) const{
    return _betaStdDevUpper[j];
  }

  // Member function for setting the standard deviation for
  // proposal for beta for fixed effect j
  void betaStdDev(const unsigned int& j,const double& sd){
    _betaStdDev[j]=sd;
  }

  bool betaAnyUpdates() const{
    return _betaAnyUpdates;
  }

  void betaAnyUpdates(const bool& newStatus){
    _betaAnyUpdates = newStatus;
  }

  // AR Cluster-specific Beta
  vector<unsigned int> nTryBetamix() const{
    return _nTryBetamix;
  }

  unsigned int nTryBetamix(const unsigned int& j) const{
    return _nTryBetamix[j];
  }

  vector<unsigned int> nAcceptBetamix() const{
    return _nAcceptBetamix;
  }

  double betamixAcceptRate(const unsigned int& j) const{
    if(_nTryBetamix[j]>0){
      return (double)_nAcceptBetamix[j]/(double)_nTryBetamix[j];
    }else{
      return 0.0;
    }
  }

  unsigned int betamixUpdateFreq() const{
    return _betamixUpdateFreq;
  }

  vector<unsigned int> nLocalAcceptBetamix() const{
    return _nLocalAcceptBetamix;
  }

  double betamixLocalAcceptRate(const unsigned int& j) const{
    return (double)_nLocalAcceptBetamix[j]/(double)_betamixUpdateFreq;
  }

  double betamixAcceptTarget() const{
    return _betamixAcceptTarget;
  }

  void betamixAddTry(const unsigned int& j){
    _nTryBetamix[j]++;
  }

  void betamixAddAccept(const unsigned int& j){
    _nAcceptBetamix[j]++;
    _nLocalAcceptBetamix[j]++;
  }

  void betamixLocalReset(const unsigned int& j){
    _nLocalAcceptBetamix[j]=0;
  }

  vector<unsigned int> nResetBetamix() const{
    return _nResetBetamix;
  }

  void betamixStdDevReset(const unsigned int& j){
    _betamixStdDev[j] = 1.0;
    _nResetBetamix[j]++;
    _betamixStdDevLower[j] = pow(10.0,-((double)_nResetBetamix[j]+1.0));
    _betamixStdDevUpper[j] = 100.0-pow(10.0,-((double)_nResetBetamix[j]+1.0));
  }

  vector<double> betamixStdDev() const{
    return _betamixStdDev;
  }

  vector<double>& betamixStdDev(){
    return _betamixStdDev;
  }

  double& betamixStdDev(const unsigned int& j){
    return _betamixStdDev[j];
  }

  const double& betamixStdDev(const unsigned int& j) const{
    return _betamixStdDev[j];
  }

  vector<double> betamixStdDevLower() const{
    return _betamixStdDevLower;
  }

  double betamixStdDevLower(const unsigned int& j) const{
    return _betamixStdDevLower[j];
  }

  vector<double> betamixStdDevUpper() const{
    return _betamixStdDevUpper;
  }

  double betamixStdDevUpper(const unsigned int& j) const{
    return _betamixStdDevUpper[j];
  }

  // Member function for setting the standard deviation for
  // proposal for beta for fixed effect j
  void betamixStdDev(const unsigned int& j,const double& sd){
    _betamixStdDev[j]=sd;
  }

  bool betamixAnyUpdates() const{
    return _betamixAnyUpdates;
  }

  void betamixAnyUpdates(const bool& newStatus){
    _betamixAnyUpdates = newStatus;
  }

  //RJ L step handling functions
  vector<unsigned int> nTryL() const{
    return _nTryL;
  }
  unsigned int nTryL(const unsigned int& j) const{
    return _nTryL[j];
  }
  vector<unsigned int> nAcceptL() const{
    return _nAcceptL;
  }
  unsigned int nAcceptL(const unsigned int& j) const{
    return _nAcceptL[j];
  }
  double LAcceptRate(const unsigned int& j) const{
    if(_nTryL[j]>0){
      return (double)_nAcceptL[j]/(double)_nTryL[j];
    }else{
      return 0.0;
    }
  }
  unsigned int LUpdateFreq() const{
    return _LUpdateFreq;
  }
  vector<unsigned int> nLocalAcceptL() const{
    return _nLocalAcceptL;
  }
  unsigned int nLocalAcceptL(const unsigned int& j) const{
    return _nLocalAcceptL[j];
  }
  double LLocalAcceptRate(const unsigned int& j) const{
    return (double)_nLocalAcceptL[j]/(double)_LUpdateFreq;
  }

  double LAcceptTarget() const{
    return _LAcceptTarget;
  }

  void LAddTry(const unsigned int& j){
    _nTryL[j]++;
  }

  void LAddAccept(const unsigned int& j){
    _nAcceptL[j]++;
    _nLocalAcceptL[j]++;
  }
  void LLocalReset(const unsigned int& j){
    _nLocalAcceptL[j]=0;
  }
  vector<unsigned int> nResetL() const{
    return _nResetL;
  }
  void LStdDevReset(const unsigned int& j){
    _LStdDev[j] = 1.0;
    _nResetL[j]++;
    _LStdDevLower[j] = pow(10.0,-((double)_nResetL[j]+1.0));
    _LStdDevUpper[j] = 100.0-pow(10.0,-((double)_nResetL[j]+1.0));
  }
  vector<double> LStdDev() const{
    return _LStdDev;
  }
  vector<double>& LStdDev(){
    return _LStdDev;
  }
  double& LStdDev(const unsigned int& j){
    return _LStdDev[j];
  }
  const double& LStdDev(const unsigned int& j) const{
    return _LStdDev[j];
  }
  vector<double> LStdDevLower() const{
    return _LStdDevLower;
  }
  double LStdDevLower(const unsigned int& j) const{
    return _LStdDevLower[j];
  }
  vector<double> LStdDevUpper() const{
    return _LStdDevUpper;
  }
  double LStdDevUpper(const unsigned int& j) const{
    return _LStdDevUpper[j];
  }
  // Member function for setting the standard deviation for
  // proposal for L
  void LStdDev(const unsigned int& j,const double& sd){
    _LStdDev[j]=sd;
  }

  bool LAnyUpdates() const{
    return _LAnyUpdates;
  }

  void LAnyUpdates(const bool& newStatus){
    _LAnyUpdates = newStatus;
  }

  unsigned int nTryAlpha() const{
    return _nTryAlpha;
  }

  unsigned int nAcceptAlpha() const{
    return _nAcceptAlpha;
  }

  double alphaAcceptRate() const{
    if(_nTryAlpha>0){
      return (double)_nAcceptAlpha/(double)_nTryAlpha;
    }else{
      return 0.0;
    }
  }

  unsigned int alphaUpdateFreq() const{
    return _alphaUpdateFreq;
  }

  unsigned int nLocalAcceptAlpha() const{
    return _nLocalAcceptAlpha;
  }

  double alphaLocalAcceptRate() const{
    return (double)_nLocalAcceptAlpha/(double)_alphaUpdateFreq;
  }

  double alphaAcceptTarget() const{
    return _alphaAcceptTarget;
  }

  void alphaAddTry(){
    _nTryAlpha++;
  }

  void alphaAddAccept(){
    _nAcceptAlpha++;
    _nLocalAcceptAlpha++;
  }

  void alphaLocalReset(){
    _nLocalAcceptAlpha=0;
  }

  unsigned int nResetAlpha() const{
    return _nResetAlpha;
  }

  void alphaStdDevReset(){
    _alphaStdDev = 1.0;
    _nResetAlpha++;
    _alphaStdDevLower = pow(10.0,-((double)_nResetAlpha+1.0));
    _alphaStdDevUpper = 100.0-pow(10.0,-((double)_nResetAlpha+1.0));
  }


  const double alphaStdDev() const{
    return _alphaStdDev;
  }

  double& alphaStdDev(){
    return _alphaStdDev;
  }

  double alphaStdDevLower() const{
    return _alphaStdDevLower;
  }

  double alphaStdDevUpper() const{
    return _alphaStdDevUpper;
  }

  // Member function for setting the standard deviation for
  // proposal for beta for fixed effect j
  void alphaStdDev(const double& sd){
    _alphaStdDev=sd;
  }

  bool alphaAnyUpdates() const{
    return _alphaAnyUpdates;
  }

  void alphaAnyUpdates(const bool& newStatus){
    _alphaAnyUpdates = newStatus;
  }

  vector<unsigned int> nTryRho() const{
    return _nTryRho;
  }

  unsigned int nTryRho(const unsigned int& j) const{
    return _nTryRho[j];
  }

  vector<unsigned int> nAcceptRho() const{
    return _nAcceptRho;
  }

  double rhoAcceptRate(const unsigned int& j) const{
    if(_nTryRho[j]>0){
      return (double)_nAcceptRho[j]/(double)_nTryRho[j];
    }else{
      return 0.0;
    }
  }

  unsigned int rhoUpdateFreq() const{
    return _rhoUpdateFreq;
  }

  vector<unsigned int> nLocalAcceptRho() const{
    return _nLocalAcceptRho;
  }


  double rhoLocalAcceptRate(const unsigned int& j) const{
    return (double)_nLocalAcceptRho[j]/(double)_rhoUpdateFreq;
  }


  double rhoAcceptTarget() const{
    return _rhoAcceptTarget;
  }

  void rhoAddTry(const unsigned int& j){
    _nTryRho[j]++;
  }

  void rhoAddAccept(const unsigned int& j){
    _nAcceptRho[j]++;
    _nLocalAcceptRho[j]++;
  }

  void rhoLocalReset(const unsigned int& j){
    _nLocalAcceptRho[j]=0;
  }

  vector<unsigned int> nResetRho() const{
    return _nResetRho;
  }

  void rhoStdDevReset(const unsigned int& j){
    _rhoStdDev[j] = 0.5;
    _nResetRho[j]++;
    _rhoStdDevLower[j] = pow(10.0,-((double)_nResetRho[j]+4.0));
    _rhoStdDevUpper[j] = 10.0-pow(10.0,-((double)_nResetRho[j]+4.0));
  }

  vector<double>& rhoStdDev(){
    return _rhoStdDev;
  }

  vector<double> rhoStdDev() const{
    return _rhoStdDev;
  }

  double& rhoStdDev(const unsigned int& j){
    return _rhoStdDev[j];
  }

  const double& rhoStdDev(const unsigned int& j) const{
    return _rhoStdDev[j];
  }

  vector<double> rhoStdDevLower() const{
    return _rhoStdDevLower;
  }

  double rhoStdDevLower(const unsigned int& j) const{
    return _rhoStdDevLower[j];
  }

  vector<double> rhoStdDevUpper() const{
    return _rhoStdDevUpper;
  }

  double rhoStdDevUpper(const unsigned int& j) const{
    return _rhoStdDevUpper[j];
  }

  void rhoStdDev(const unsigned int& j,const double& sd){
    _rhoStdDev[j]=sd;
  }

  bool rhoAnyUpdates() const{
    return _rhoAnyUpdates;
  }

  void rhoAnyUpdates(const bool& newStatus){
    _rhoAnyUpdates = newStatus;
  }


  unsigned int nTryLambda() const{
    return _nTryLambda;
  }

  unsigned int nAcceptLambda() const{
    return _nAcceptLambda;
  }

  double lambdaAcceptRate() const{
    if(_nTryLambda>0){
      return (double)_nAcceptLambda/(double)_nTryLambda;
    }else{
      return 0.0;
    }
  }

  unsigned int nLocalAcceptLambda() const{
    return _nLocalAcceptLambda;
  }

  unsigned int lambdaUpdateFreq() const{
    return _lambdaUpdateFreq;
  }

  double lambdaLocalAcceptRate() const{
    return (double)_nLocalAcceptLambda/(double)_lambdaUpdateFreq;
  }

  double lambdaAcceptTarget() const{
    return _lambdaAcceptTarget;
  }

  void lambdaAddTry(){
    _nTryLambda++;
  }

  void lambdaAddAccept(){
    _nAcceptLambda++;
    _nLocalAcceptLambda++;
  }

  void lambdaLocalReset(){
    _nLocalAcceptLambda=0;
  }

  unsigned int nResetLambda() const{
    return _nResetLambda;
  }

  void lambdaStdDevReset(){
    _lambdaStdDev = 1.0;
    _nResetLambda++;
    _lambdaStdDevLower = pow(10.0,-((double)_nResetLambda+1.0));
    _lambdaStdDevUpper = 100.0-pow(10.0,-((double)_nResetLambda+1.0));
  }

  double lambdaStdDev() const{
    return _lambdaStdDev;
  }

  double& lambdaStdDev(){
    return _lambdaStdDev;
  }

  double lambdaStdDevLower() const{
    return _lambdaStdDevLower;
  }

  double lambdaStdDevUpper() const{
    return _lambdaStdDevUpper;
  }

  // Member function for setting the standard deviation for
  // proposal for beta for fixed effect j
  void lambdaStdDev(const double& sd){
    _lambdaStdDev=sd;
  }

  bool lambdaAnyUpdates() const{
    return _lambdaAnyUpdates;
  }

  void lambdaAnyUpdates(const bool& newStatus){
    _lambdaAnyUpdates = newStatus;
  }


  // Need to define a copy iterator
  pReMiuMPropParams& operator=(const pReMiuMPropParams& propParams){
    _nTryTheta=propParams.nTryTheta();
    _nAcceptTheta=propParams.nAcceptTheta();
    _nLocalAcceptTheta=propParams.nLocalAcceptTheta();
    _nResetTheta=propParams.nResetTheta();
    _thetaStdDev=propParams.thetaStdDev();
    _thetaStdDevLower=propParams.thetaStdDevLower();
    _thetaStdDevUpper=propParams.thetaStdDevUpper();
    _thetaAcceptTarget=propParams.thetaAcceptTarget();
    _thetaUpdateFreq=propParams.thetaUpdateFreq();
    _thetaAnyUpdates=propParams.thetaAnyUpdates();
    _nTryBeta=propParams.nTryBeta();
    _nAcceptBeta=propParams.nAcceptBeta();
    _nLocalAcceptBeta=propParams.nLocalAcceptBeta();
    _nResetBeta=propParams.nResetBeta();
    _betaStdDev=propParams.betaStdDev();
    _betaStdDevLower=propParams.betaStdDevLower();
    _betaStdDevUpper=propParams.betaStdDevUpper();
    _betaAcceptTarget=propParams.betaAcceptTarget();
    _betaUpdateFreq=propParams.betaUpdateFreq();
    _betaAnyUpdates=propParams.betaAnyUpdates();

    // AR cluster-specific betas
    _nTryBetamix=propParams.nTryBetamix();
    _nAcceptBetamix=propParams.nAcceptBetamix();
    _nLocalAcceptBetamix=propParams.nLocalAcceptBetamix();
    _nResetBetamix=propParams.nResetBetamix();
    _betamixStdDev=propParams.betamixStdDev();
    _betamixStdDevLower=propParams.betamixStdDevLower();
    _betamixStdDevUpper=propParams.betamixStdDevUpper();
    _betamixAcceptTarget=propParams.betamixAcceptTarget();
    _betamixUpdateFreq=propParams.betamixUpdateFreq();
    _betamixAnyUpdates=propParams.betamixAnyUpdates();

    //AR LME outcome option
    _nTrySigmaE=propParams.nTrySigmaE();
    _nAcceptSigmaE=propParams.nAcceptSigmaE();
    _nLocalAcceptSigmaE=propParams.nLocalAcceptSigmaE();
    _nResetSigmaE=propParams.nResetSigmaE();
    _SigmaEStdDev=propParams.SigmaEStdDev();
    _SigmaEStdDevLower=propParams.SigmaEStdDevLower();
    _SigmaEStdDevUpper=propParams.SigmaEStdDevUpper();
    _SigmaEAcceptTarget=propParams.SigmaEAcceptTarget();
    _SigmaEUpdateFreq=propParams.SigmaEUpdateFreq();
    _SigmaEAnyUpdates=propParams.SigmaEAnyUpdates();

    //RJ set L step values
    _nTryL=propParams.nTryL();
    _nAcceptL=propParams.nAcceptL();
    _nLocalAcceptL=propParams.nLocalAcceptL();
    _nResetL=propParams.nResetL();
    _LStdDev=propParams.LStdDev();
    _LStdDevLower=propParams.LStdDevLower();
    _LStdDevUpper=propParams.LStdDevUpper();
    _LAcceptTarget=propParams.LAcceptTarget();
    _LUpdateFreq=propParams.LUpdateFreq();
    _LAnyUpdates=propParams.LAnyUpdates();

    _nTryAlpha=propParams.nTryAlpha();
    _nAcceptAlpha=propParams.nAcceptAlpha();
    _nLocalAcceptAlpha=propParams.nLocalAcceptAlpha();
    _nResetAlpha=propParams.nResetAlpha();
    _alphaStdDev=propParams.alphaStdDev();
    _alphaStdDevLower=propParams.alphaStdDevLower();
    _alphaStdDevUpper=propParams.alphaStdDevUpper();
    _alphaAcceptTarget=propParams.alphaAcceptTarget();
    _alphaUpdateFreq=propParams.alphaUpdateFreq();
    _alphaAnyUpdates=propParams.alphaAnyUpdates();
    _nTryRho=propParams.nTryRho();
    _nAcceptRho=propParams.nAcceptRho();
    _nLocalAcceptRho=propParams.nLocalAcceptRho();
    _nResetRho=propParams.nResetRho();
    _rhoStdDev=propParams.rhoStdDev();
    _rhoStdDevLower=propParams.rhoStdDevLower();
    _rhoStdDevUpper=propParams.rhoStdDevUpper();
    _rhoAcceptTarget=propParams.rhoAcceptTarget();
    _rhoUpdateFreq=propParams.rhoUpdateFreq();
    _rhoAnyUpdates=propParams.rhoAnyUpdates();
    _nTryLambda=propParams.nTryLambda();
    _nAcceptLambda=propParams.nAcceptLambda();
    _nLocalAcceptLambda=propParams.nLocalAcceptLambda();
    _nResetLambda=propParams.nResetLambda();
    _lambdaStdDev=propParams.lambdaStdDev();
    _lambdaStdDevLower=propParams.lambdaStdDevLower();
    _lambdaStdDevUpper=propParams.lambdaStdDevUpper();
    _lambdaAcceptTarget=propParams.lambdaAcceptTarget();
    _lambdaUpdateFreq=propParams.lambdaUpdateFreq();
    _lambdaAnyUpdates=propParams.lambdaAnyUpdates();
    return *this;

  }

private:
  unsigned int _nTryTheta;
  unsigned int _nAcceptTheta;
  unsigned int _nLocalAcceptTheta;
  unsigned int _nResetTheta;
  double _thetaStdDev;
  double _thetaStdDevLower;
  double _thetaStdDevUpper;
  double _thetaAcceptTarget;
  unsigned int _thetaUpdateFreq;
  bool _thetaAnyUpdates;

  unsigned int _nTrySigmaE;
  unsigned int _nAcceptSigmaE;
  unsigned int _nLocalAcceptSigmaE;
  unsigned int _nResetSigmaE;
  double _SigmaEStdDev;
  double _SigmaEStdDevLower;
  double _SigmaEStdDevUpper;
  double _SigmaEAcceptTarget;
  unsigned int _SigmaEUpdateFreq;
  bool _SigmaEAnyUpdates;


  vector<unsigned int> _nTryBeta;
  vector<unsigned int> _nAcceptBeta;
  vector<unsigned int> _nLocalAcceptBeta;
  vector<unsigned int> _nResetBeta;
  vector<double> _betaStdDev;
  vector<double> _betaStdDevLower;
  vector<double> _betaStdDevUpper;
  double _betaAcceptTarget;
  unsigned int _betaUpdateFreq;
  bool _betaAnyUpdates;

  vector<unsigned int> _nTryBetamix;
  vector<unsigned int> _nAcceptBetamix;
  vector<unsigned int> _nLocalAcceptBetamix;
  vector<unsigned int> _nResetBetamix;
  vector<double> _betamixStdDev;
  vector<double> _betamixStdDevLower;
  vector<double> _betamixStdDevUpper;
  double _betamixAcceptTarget;
  unsigned int _betamixUpdateFreq;
  bool _betamixAnyUpdates;
  //RJ declare L step values
  vector<unsigned int> _nTryL;
  vector<unsigned int> _nAcceptL;
  vector<unsigned int> _nLocalAcceptL;
  vector<unsigned int> _nResetL;
  vector<double> _LStdDev;
  vector<double> _LStdDevLower;
  vector<double> _LStdDevUpper;
  double _LAcceptTarget;
  unsigned int _LUpdateFreq;
  bool _LAnyUpdates;
  unsigned int _nTryAlpha;
  unsigned int _nAcceptAlpha;
  unsigned int _nLocalAcceptAlpha;
  unsigned int _nResetAlpha;
  double _alphaStdDev;
  double _alphaStdDevLower;
  double _alphaStdDevUpper;
  double _alphaAcceptTarget;
  unsigned int _alphaUpdateFreq;
  bool _alphaAnyUpdates;
  vector<unsigned int> _nTryRho;
  vector<unsigned int> _nAcceptRho;
  vector<unsigned int> _nLocalAcceptRho;
  vector<unsigned int> _nResetRho;
  vector<double> _rhoStdDev;
  vector<double> _rhoStdDevLower;
  vector<double> _rhoStdDevUpper;
  double _rhoAcceptTarget;
  unsigned int _rhoUpdateFreq;
  bool _rhoAnyUpdates;
  unsigned int _nTryLambda;
  unsigned int _nAcceptLambda;
  unsigned int _nLocalAcceptLambda;
  unsigned int _nResetLambda;
  double _lambdaStdDev;
  double _lambdaStdDevLower;
  double _lambdaStdDevUpper;
  double _lambdaAcceptTarget;
  unsigned int _lambdaUpdateFreq;
  bool _lambdaAnyUpdates;

};

/*********** BLOCK 1 p(v^A,Theta^A,u|.) **********************************/
// A=Active, and Theta contains: phi, mu, Tau, gamma, theta
// We proceed by sampling from p(v^A,Theta^A|.) i.e. marginalising out the u
// and then sampling from p(u|v^A,Theta^A,.). The first of these steps, is achieved
// in a number of stages, updating the various components in turn, and then
// performing label switching.

// Gibbs move for updating the v (for psi) which are active i.e. v_c where c<=Z_max
// This is done by using Gibbs to sample from p(v^A|.), which is the conditional
// for this block with u integrated out (the marginal over u). This can be thought of as
// a MH sample from p(v^A,Theta^A|.), where we sample v^A from conditional, and
// leave Theta^A unchanged (thus having an acceptance rate of 1).
void gibbsForVActive(mcmcChain<pReMiuMParams>& chain,
                     unsigned int& nTry,unsigned int& nAccept,
                     const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                     pReMiuMPropParams& propParams,
                     baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();

  nTry++;
  nAccept++;

  // Find the active clusters
  unsigned int maxZ = currentParams.workMaxZi();

  // This is sampled from the posterior given the z vector above
  // Prior comes from the conjugacy of the dirichlet and multinomial
  vector<unsigned int> sumCPlus1ToMaxMembers(maxZ+1);
  sumCPlus1ToMaxMembers[maxZ]=0;
  for(int c=maxZ-1;c>=0;c--){
    sumCPlus1ToMaxMembers[c]=sumCPlus1ToMaxMembers[c+1]+currentParams.workNXInCluster(c+1);
  }

  double tmp=0.0;
  double alpha = currentParams.alpha();
  double dPitmanYor = currentParams.dPitmanYor();

  for(unsigned int c=0;c<=maxZ;c++){//
    double vVal = betaRand(rndGenerator,1.0+currentParams.workNXInCluster(c)-dPitmanYor,alpha+sumCPlus1ToMaxMembers[c]+dPitmanYor*(c+1));
    currentParams.v(c,vVal);
    // Set psi
    currentParams.logPsi(c,tmp+log(vVal));
    tmp += log(1-vVal);

  }
}

// Moves for updating the Theta which are active i.e. Theta_c where c<=Z_max
// Several different moves, which each update specific elements of Theta

// Move for updating phi
// Gibbs used if no variable selection, or binary switch based variable selection
// Otherwise Metropolis Hastings is used
void updateForPhiActive(mcmcChain<pReMiuMParams>& chain,
                        unsigned int& nTry,unsigned int& nAccept,
                        const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                        pReMiuMPropParams& propParams,
                        baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  const pReMiuMData& dataset = model.dataset();
  string varSelectType = model.options().varSelectType();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nDiscreteCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }
  // Find the number of subjects
  unsigned int nSubjects = dataset.nSubjects();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  double currentLogPost=0.0;

  for(unsigned int c=0;c<=maxZ;c++){
    // Loop over the covariates
    for(unsigned int j=0;j<nCovariates;j++){
      nTry++;

      if(varSelectType.compare("Continuous")==0){
        currentLogPost=logCondPostPhicj(currentParams,model,c,j);
      }

      unsigned int nCategories = currentParams.nCategories(j);
      // We are updating phis
      // First we must count how many individuals have Xij in each of the
      // possible categories for covariate j.
      vector<double> dirichParams(nCategories,hyperParams.aPhi(j));
      double gammacj = currentParams.gamma(c,j);
      for(unsigned int i=0;i<nSubjects;i++){
        int zi = currentParams.z(i);
        if(zi==(int)c){
          int Xij = dataset.discreteX(i,j);
          // When no variable selection will always add 1.
          // In the binary variable selection case
          // this will add a 1 only when the
          // variable is switched on (as required)
          // In the continuous case this seems like a sensible proposal
          dirichParams[Xij]=dirichParams[Xij]+gammacj;
        }
      }
      vector<double> currentLogPhi(nCategories);
      currentLogPhi=currentParams.logPhi(c,j);

      vector<double> proposedLogPhi(nCategories);
      proposedLogPhi=dirichletRand(rndGenerator,dirichParams);
      for(unsigned int p=0;p<nCategories;p++){
        proposedLogPhi[p]=log(proposedLogPhi[p]);
      }
      currentParams.logPhi(c,j,proposedLogPhi);
      // If no variable selection or binary variable selection
      // this is a sample from full conditional so no accept reject
      // step need. If it is continuous variable selection we
      // need an accept reject decision
      if(varSelectType.compare("Continuous")==0){
        double proposedLogPost = logCondPostPhicj(currentParams,model,c,j);
        double logAcceptRatio=0.0;
        logAcceptRatio=proposedLogPost-currentLogPost;
        logAcceptRatio+=logPdfDirichlet(currentLogPhi,dirichParams,true);
        logAcceptRatio-=logPdfDirichlet(proposedLogPhi,dirichParams,true);
        if(unifRand(rndGenerator)<exp(logAcceptRatio)){
          // Move accepted
          nAccept++;
        }else{
          // Move rejected
          // Reset phi
          currentParams.logPhi(c,j,currentLogPhi);
        }
      }else{
        nAccept++;
      }
    }
  }
}

// Gibbs update for mu in Normal covariate case
void gibbsForMuActive(mcmcChain<pReMiuMParams>& chain,
                      unsigned int& nTry,unsigned int& nAccept,
                      const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                      pReMiuMPropParams& propParams,
                      baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  const pReMiuMData& dataset = model.dataset();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nContinuousCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }
  // Find the number of subjects
  unsigned int nSubjects = dataset.nSubjects();

  nTry++;
  nAccept++;

  // In the following it is useful to have the rows of X as
  // Eigen dynamic vectors
  vector<VectorXd> xi(nSubjects);
  for(unsigned int i=0;i<nSubjects;i++){
    xi[i].setZero(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      xi[i](j)=dataset.continuousX(i,j);
    }
  }

  // We begin by computing the mean X for individuals in each cluster
  vector<VectorXd> meanX(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    meanX[c].setZero(nCovariates);
  }

  for(unsigned int i=0;i<nSubjects;i++){
    meanX[currentParams.z(i)]=meanX[currentParams.z(i)]+xi[i];
  }

  vector<MatrixXd> gammaMat(maxZ+1);
  vector<MatrixXd> oneMinusGammaMat(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    gammaMat[c].setZero(nCovariates,nCovariates);
    oneMinusGammaMat[c].setZero(nCovariates,nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      gammaMat[c](j,j)=currentParams.gamma(c,currentParams.nDiscreteCovs()+j);//j);
      oneMinusGammaMat[c](j,j)=1-gammaMat[c](j,j);
    }
  }

  for(unsigned int c=0;c<=maxZ;c++){
    // Having computed this we can calcuate the posterior mean
    // and posterior covariance for each mu_c
    int nXInC = currentParams.workNXInCluster(c);
    if(nXInC>0){
      meanX[c]=meanX[c]/(double)nXInC;
    }else{
      meanX[c].setZero(nCovariates);
    }
    MatrixXd covMat(nCovariates,nCovariates);
    covMat = (hyperParams.Tau0()+nXInC*gammaMat[c]*currentParams.Tau(c)*gammaMat[c]).inverse();
    VectorXd meanVec(nCovariates);
    meanVec = hyperParams.Tau0()*hyperParams.mu0()+
      nXInC*gammaMat[c]*currentParams.Tau(c)*(meanX[c]-oneMinusGammaMat[c]*currentParams.nullMu());
    meanVec = covMat*meanVec;
    VectorXd mu(nCovariates);
    // We sample from this posterior
    mu = multivarNormalRand(rndGenerator,meanVec,covMat);

    // We store our sample
    currentParams.mu(c,mu);
  }
}

// Gibbs update for mu in Normal covariate case and use of normal inverse prior
void gibbsForMuActiveNIWP(mcmcChain<pReMiuMParams>& chain,
                          unsigned int& nTry,unsigned int& nAccept,
                          const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                          pReMiuMPropParams& propParams,
                          baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  const pReMiuMData& dataset = model.dataset();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nContinuousCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }
  // Find the number of subjects
  unsigned int nSubjects = dataset.nSubjects();

  nTry++;
  nAccept++;

  // In the following it is useful to have the rows of X as
  // Eigen dynamic vectors
  vector<VectorXd> xi(nSubjects);
  for(unsigned int i=0;i<nSubjects;i++){
    xi[i].setZero(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      xi[i](j)=dataset.continuousX(i,j);
    }
  }

  // We begin by computing the mean X for individuals in each cluster
  vector<VectorXd> meanX(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    meanX[c].setZero(nCovariates);
  }

  for(unsigned int i=0;i<nSubjects;i++){
    meanX[currentParams.z(i)]=meanX[currentParams.z(i)]+xi[i];
  }

  vector<MatrixXd> gammaMat(maxZ+1);
  vector<MatrixXd> oneMinusGammaMat(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    gammaMat[c].setZero(nCovariates,nCovariates);
    oneMinusGammaMat[c].setZero(nCovariates,nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      gammaMat[c](j,j)=currentParams.gamma(c,currentParams.nDiscreteCovs()+j);
      oneMinusGammaMat[c](j,j)=1-gammaMat[c](j,j);
    }
  }

  for(unsigned int c=0;c<=maxZ;c++){
    // Having computed this we can calcuate the posterior mean
    // and posterior covariance for each mu_c
    int nXInC = currentParams.workNXInCluster(c);
    if(nXInC>0){
      meanX[c]=meanX[c]/(double)nXInC;
    }else{
      meanX[c].setZero(nCovariates);
    }
    MatrixXd covMat(nCovariates,nCovariates);
    covMat=(gammaMat[c]*currentParams.Sigma(c)*gammaMat[c])/(hyperParams.nu0()+nXInC);
    //There are 0 on the diagonal elements of the covariance matrix when the covariate j is not selected
    //we replace them by 0.1, the value does not care, since the sampling values will be replaced by \bar{x}_j
    for (unsigned int j=0; j<=nCovariates; j++){
      if (covMat(j,j)==0) covMat(j,j)=0.1;
    }
    VectorXd meanVec(nCovariates);
    meanVec = hyperParams.nu0()*hyperParams.mu0()+
      nXInC*(gammaMat[c]*meanX[c]-oneMinusGammaMat[c]*currentParams.nullMu());
    meanVec /= (hyperParams.nu0()+nXInC);
    VectorXd mu(nCovariates);
    // We sample from this posterior
    mu = multivarNormalRand(rndGenerator,meanVec,covMat);

    // We store our sample
    currentParams.mu(c,mu);

  }

}


// Gibbs update for Tau in the Normal covariate case
void gibbsForTauActive(mcmcChain<pReMiuMParams>& chain,
                       unsigned int& nTry,unsigned int& nAccept,
                       const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                       pReMiuMPropParams& propParams,
                       baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  const pReMiuMData& dataset = model.dataset();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nContinuousCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }
  // Find the number of subjects
  unsigned int nSubjects = dataset.nSubjects();

  nTry++;
  nAccept++;

  // In the following it is useful to have the rows of X as
  // Eigen dynamic vectors
  vector<VectorXd> xi(nSubjects);
  for(unsigned int i=0;i<nSubjects;i++){
    xi[i].setZero(nCovariates);
    for(unsigned int j=0;j<nCovariates;j++){
      xi[i](j)=dataset.continuousX(i,j);
    }
  }

  vector<MatrixXd> Rc(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    Rc[c].setZero(nCovariates,nCovariates);
  }

  for(unsigned int i=0;i<nSubjects;i++){
    unsigned int zi = currentParams.z(i);
    Rc[zi]=Rc[zi]+(xi[i]-currentParams.workMuStar(zi))*((xi[i]-currentParams.workMuStar(zi)).transpose());
  }

  for(unsigned int c=0;c<=maxZ;c++){
    Rc[c]=(hyperParams.R0().inverse()+Rc[c]).inverse();
    MatrixXd Tau = wishartRand(rndGenerator,Rc[c],currentParams.workNXInCluster(c)+hyperParams.kappa0());
    currentParams.Tau(c,Tau);
  }
}

//AR update the cluster-specific covariance matrix of the random effects for the LME outcome option
void gibbsForCovRELMEActive(mcmcChain<pReMiuMParams>& chain,
                               unsigned int& nTry,unsigned int& nAccept,
                               const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                               pReMiuMPropParams& propParams,
                               baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  // Find the number of subjects
  unsigned int nSubjects = dataset.nSubjects();
  unsigned int nOutcomes = dataset.nOutcomes();
  vector<unsigned int> nRandomEffects = dataset.nRandomEffects();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> y = dataset.continuousY();
  vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  vector<MatrixXd> RandomEffects = currentParams.RandomEffects();
  unsigned int nCategoriesY = dataset.nCategoriesY();


  //Sampling of the hyperprior Lambda: covRE~ IW(n,Lambda) and Lambda~W(delta+nRandomEffects+1,workTauLME_R0^-1)
  // Lambda|theta ~ W(SigmaLME_kappa0 + delta, ((workTauLME_R0^-1)^-1 + covRE^-1)^-1)
  // delta = SigmaLME_kappa0
  MatrixXd Lambda;

  // In the following it is useful to have the rows of X as
  // Get individual random effects


  int ind=0;//tstart
  int ind_y=0;//y and t
  for(unsigned int m=0;m<nOutcomes;m++){
    vector<VectorXd> bi(nSubjects);

    for(unsigned int i=0;i<nSubjects;i++){
      bi[i].setZero(nRandomEffects[m]);
      for(unsigned int j=0;j<nRandomEffects[m];j++){
        bi[i](j)=currentParams.RandomEffects(m,i,j);//dataset.continuousX(i,j);
      }
    }

    MatrixXd R;

    // Generate variance covariance matrix
    //R=(Lambda+R);
    //R = R.inverse();
    MatrixXd Tau;
    VectorXd nS_c;
    nS_c.setZero(maxZ+1);
    vector<MatrixXd> Rc(maxZ+1);
    for(unsigned int c=0;c<=maxZ;c++){
      Rc[c].setZero(nRandomEffects[m],nRandomEffects[m]);
    }

    for(unsigned int i=0;i<nSubjects;i++){
      unsigned int zi = currentParams.z(i);
      Rc[zi]=Rc[zi]+bi[i]*(bi[i].transpose());
      nS_c[zi]++;
    }
    //Tau = wishartRand(rndGenerator,R,nSubjects+hyperParams.SigmaLME_kappa0());
    //Tau(0,0)=1;
    // MatrixXd covv = Tau.inverse();
    // covv(0,0)=0.3;covv(0,1)=-0.1;covv(0,2)=0.1;
    // covv(1,0)=-0.1;covv(1,1)=0.1;covv(1,2)=0.2;
    // covv(2,0)=0.1;covv(2,1)=0.2;covv(2,2)=1;
    // Tau = covv.inverse();
    // currentParams.covRE(covv);

    //currentParams.covRE(Tau.inverse());

    for(unsigned int c=0;c<=maxZ;c++){

      Rc[c]=(hyperParams.workTauLME_R0(m).inverse()+Rc[c]).inverse();
      MatrixXd Tau = wishartRand(rndGenerator,Rc[c],nS_c[c]+hyperParams.SigmaLME_kappa0(m));
      currentParams.covRE(m,c, Tau.inverse());
    }

    // Update Random effects
    // Mean B*Z^T V^{-1}(Yi-Xi beta)
    // Variance B - B*Zi^T*Vi^{-1}* (Zi*B^T)

    for(unsigned int i=0;i<nSubjects;i++){
      VectorXd yi;
      VectorXd ui(nRandomEffects[0]);
      unsigned int zi= currentParams.z(i);

      unsigned int ni =  (tStop[ind] - tStart[ind] + 1);
      yi.resize(ni);

      for(unsigned int j=0;j<tStop[ind]-tStart[ind]+1;j++){

        yi(j) = y[ind_y];//yi(j) = y[tStart[ind]-1+j];

        for(unsigned int b=0;b<nFixedEffects[m];b++){
          yi(j)-=currentParams.beta(m,b,0,nCategoriesY)*dataset.W_LME(m,tStart[ind]-1+j,b);
        }
        for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
          yi(j)-=currentParams.beta_mix(m,zi,b,0,nCategoriesY)*dataset.W_LME_mix(m,tStart[ind]-1+j,b);
        }

        ind_y++;
      }

      MatrixXd block=dataset.W_RE(m,tStart[ind]-1, 0, ni, nRandomEffects[m]);
      MatrixXd sigmae=MatrixXd::Identity(ni, ni) * currentParams.SigmaE(m);

      MatrixXd V = block *currentParams.covRE(m,zi)* block.transpose() + sigmae;
      LLT<MatrixXd> lltOfA(V); // compute the Cholesky decomposition of A
      MatrixXd L = lltOfA.matrixL();
      //double logDetPrecMat=  2*log(L.determinant());
      MatrixXd Vi_inv = L.inverse().transpose()*L.inverse();
      VectorXd mu = currentParams.covRE(m,zi)*block.transpose()*Vi_inv*yi;

      //B - B*Zi^T*Vi^{-1}* (Zi*B^T)
      MatrixXd cov = currentParams.covRE(m,zi) - currentParams.covRE(m,zi)*block.transpose()*Vi_inv*block*currentParams.covRE(m,zi);

      ui = multivarNormalRand(rndGenerator,mu,cov);
      currentParams.RandomEffects(m,i,ui);
      ind ++;
    }

  }
}

//AR update the cluster-specific covariance matrix of the random effects for the LME outcome option
void gibbsForCovRELMEInActive(mcmcChain<pReMiuMParams>& chain,
                            unsigned int& nTry,unsigned int& nAccept,
                            const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                            pReMiuMPropParams& propParams,
                            baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();


  // Find the number of subjects
  const pReMiuMData& dataset = model.dataset();
  unsigned int nOutcomes = dataset.nOutcomes();
  for(unsigned int m=0;m<nOutcomes;m++){
    for(unsigned int c=maxZ+1;c<currentParams.maxNClusters();c++){
      MatrixXd Tau = wishartRand(rndGenerator,hyperParams.workTauLME_R0(m),hyperParams.SigmaLME_kappa0(m)); //added
      currentParams.covRE(m,c, Tau.inverse());
    }
  }
}


void gibbsForMVNTauActive(mcmcChain<pReMiuMParams>& chain,
                          unsigned int& nTry,unsigned int& nAccept,
                          const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                          pReMiuMPropParams& propParams,
                          baseGeneratorType& rndGenerator){
  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  // Find the number of subjects
  unsigned int nSubjects = dataset.nSubjects();
  // Find the number of outcomes
  unsigned int nOutcomes = dataset.nTimes()/nSubjects;
  nTry++;
  nAccept++;
  // In the following it is useful to have the rows of X as
  // Eigen dynamic vectors
  vector<VectorXd> yi(nSubjects);
  for(unsigned int i=0;i<nSubjects;i++){
    yi[i].setZero(nOutcomes);
    for(unsigned int j=0;j<nOutcomes;j++){
      yi[i](j)=dataset.continuousY(i*nOutcomes+j);
    }
  }
  vector<VectorXd> meanY(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    meanY[c].setZero(nOutcomes);
  }

  for(unsigned int i=0;i<nSubjects;i++){
    meanY[currentParams.z(i)]=meanY[currentParams.z(i)]+yi[i];
  }
  for(unsigned int c=0;c<=maxZ;c++){
    int nXInC = currentParams.workNXInCluster(c);
    if(nXInC>0){
      meanY[c]=meanY[c]/(double)nXInC;
    }
  }
  vector<MatrixXd> Rc(maxZ+1);
  for(unsigned int c=0;c<=maxZ;c++){
    Rc[c].setZero(nOutcomes,nOutcomes);
    int nXInC = currentParams.workNXInCluster(c);
    Rc[c] = (meanY[c]-hyperParams.MVNmu0())*((meanY[c]-hyperParams.MVNmu0()).transpose())*hyperParams.MVNkappa0()*nXInC/(hyperParams.MVNkappa0()+nXInC);
  }
  for(unsigned int i=0;i<nSubjects;i++){
    unsigned int zi = currentParams.z(i);
    Rc[zi] = Rc[zi]+(yi[i]-meanY[zi])*((yi[i]-meanY[zi]).transpose());
    //RJ from normal covariate case:
    //Rc[zi] = Rc[zi]+(yi[i]-currentParams.MVNmu(zi))*((yi[i]-currentParams.MVNmu(zi)).transpose());
  }
  for(unsigned int c=0;c<=maxZ;c++){

    Rc[c] = hyperParams.MVNR0()+Rc[c];
    MatrixXd Sigma = invWishartRand(rndGenerator,(Rc[c]),currentParams.workNXInCluster(c)+hyperParams.MVNnu0());
    MatrixXd Tau = Sigma.inverse();
    currentParams.MVNTau(c,Tau);

    int nXInC = currentParams.workNXInCluster(c);
    VectorXd meanVec(nOutcomes);
    meanVec = hyperParams.MVNkappa0()*hyperParams.MVNmu0()+nXInC*meanY[c];
    meanVec /= (hyperParams.MVNkappa0()+nXInC);
    VectorXd mu(nOutcomes);
    MatrixXd covMat(nOutcomes,nOutcomes);
    covMat=(currentParams.MVNSigma(c))/(hyperParams.MVNkappa0()+nXInC);
    mu = multivarNormalRand(rndGenerator,meanVec,covMat);
    currentParams.MVNmu(c,mu);
  }
}


// Gibbs update for update of gamma (only used in the binary variable selection case)
void gibbsForGammaActive(mcmcChain<pReMiuMParams>& chain,
                         unsigned int& nTry,unsigned int& nAccept,
                         const mcmcModel<pReMiuMParams,
                                         pReMiuMOptions,
                                         pReMiuMData>& model,
                                         pReMiuMPropParams& propParams,
                                         baseGeneratorType& rndGenerator){


  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of subjects
  unsigned int nCovariates = currentParams.nCovariates();
  unsigned int nSubjects = currentParams.nSubjects();
  unsigned int maxZ = currentParams.workMaxZi();
  string covariateType = model.options().covariateType();
  string varSelectType = model.options().varSelectType();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  nTry++;
  nAccept++;

  for(unsigned int j=0;j<nCovariates;j++){

    for(unsigned int c=0;c<=maxZ;c++){

      vector<double> currentGamma=currentParams.gamma(c);
      // work in terms of prob of sticking with current value and switching value
      // Compute probability of sticking
      double logProbStick=0;
      double logProbSwitch=0;
      double probStick=0;
      if(currentParams.omega(j)==0){
        // Nothing to do - not allowed to change
        continue;
      }else{
        for(unsigned int i=0;i<nSubjects;i++){
          unsigned int zi = currentParams.z(i);
          if(zi==c){
            logProbStick+=currentParams.workLogPXiGivenZi(i);
          }
        }
        logProbStick+=(currentGamma[j]*log(currentParams.rho(j))+
          (1-currentGamma[j])*log(1-currentParams.rho(j)));

        // Now compute probability of switching
        currentGamma[j]=1-currentGamma[j];
        currentParams.gamma(c,j,currentGamma[j],covariateType);
        for(unsigned int i=0;i<nSubjects;i++){
          unsigned int zi = currentParams.z(i);
          if(zi==c){
            logProbSwitch+=currentParams.workLogPXiGivenZi(i);
          }
        }
        logProbSwitch+=(currentGamma[j]*log(currentParams.rho(j))+
          (1-currentGamma[j])*log(1-currentParams.rho(j)));

        double maxLogProb;
        if(logProbSwitch<logProbStick){
          maxLogProb=logProbStick;
        }else{
          maxLogProb=logProbSwitch;
        }


        probStick=exp(logProbStick-maxLogProb)/(exp(logProbStick-maxLogProb)+exp(logProbSwitch-maxLogProb));
      }
      if(unifRand(rndGenerator)<probStick){
        // Sticking (we actually need to revert back to what we were
        // before doing the calculations)
        currentGamma[j]=1-currentGamma[j];
        currentParams.gamma(c,j,currentGamma[j],covariateType);
      }
      // Otherwise switching but nothing to do here as we had already done the
      // switch in the calculations
    }
  }
}


// Adaptive Metropolis-Hastings for theta
void metropolisHastingsForThetaActive(mcmcChain<pReMiuMParams>& chain,
                                      unsigned int& nTry,unsigned int& nAccept,
                                      const mcmcModel<pReMiuMParams,
                                                      pReMiuMOptions,
                                                      pReMiuMData>& model,
                                                      pReMiuMPropParams& propParams,
                                                      baseGeneratorType& rndGenerator){

  const string outcomeType = model.dataset().outcomeType();

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  unsigned int nCategoriesY = currentParams.nCategoriesY();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);
  // Define a normal random number generator
  randomNormal normRand(0,1);

  double thetaTargetRate = propParams.thetaAcceptTarget();
  unsigned int thetaUpdateFreq = propParams.thetaUpdateFreq();

  double currentCondLogPost = logCondPostThetaBeta(currentParams,model);

  for(unsigned int c=0;c<=maxZ;c++){
    for (unsigned int k=0;k<nCategoriesY;k++){
      nTry++;
      propParams.thetaAddTry();
      double& stdDev = propParams.thetaStdDev();
      double thetaOrig = currentParams.theta(c,k);
      double thetaProp = thetaOrig +stdDev*normRand(rndGenerator);

      currentParams.theta(c,k,thetaProp);
      double propCondLogPost = logCondPostThetaBeta(currentParams,model);
      double logAcceptRatio = propCondLogPost - currentCondLogPost;
      if(unifRand(rndGenerator)<exp(logAcceptRatio)){
        nAccept++;
        propParams.thetaAddAccept();
        currentCondLogPost = propCondLogPost;
        // Update the std dev of the proposal
        if(propParams.nTryTheta()%thetaUpdateFreq==0){
          stdDev += 10*(propParams.thetaLocalAcceptRate()-thetaTargetRate)/
            pow((double)(propParams.nTryTheta()/thetaUpdateFreq)+2.0,0.75);
          propParams.thetaAnyUpdates(true);
          if(stdDev>propParams.thetaStdDevUpper()||stdDev<propParams.thetaStdDevLower()){
            propParams.thetaStdDevReset();
          }
          propParams.thetaLocalReset();
        }
      }else{
        currentParams.theta(c,k,thetaOrig);
        // Update the std dev of the proposal
        if(propParams.nTryTheta()%thetaUpdateFreq==0){
          stdDev += 10*(propParams.thetaLocalAcceptRate()-thetaTargetRate)/
            pow((double)(propParams.nTryTheta()/thetaUpdateFreq)+2.0,0.75);
          propParams.thetaAnyUpdates(true);
          if(stdDev<propParams.thetaStdDevLower()||stdDev>propParams.thetaStdDevUpper()){
            propParams.thetaStdDevReset();
          }
          propParams.thetaLocalReset();
        }
      }
    }
  }
}


// Label switching moves (as recommended in Papaspiliopoulos and Roberts, 2008)
void metropolisHastingsForLabels123(mcmcChain<pReMiuMParams>& chain,
                                    unsigned int& nTry,unsigned int& nAccept,
                                    const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                                    pReMiuMPropParams& propParams,
                                    baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();

  unsigned int maxZ = currentParams.workMaxZi();
  if(maxZ==0){
    // If there is only one cluster with individuals in, don't do anything
    return;
  }
  string varSelectType = model.options().varSelectType();
  string covariateType = model.options().covariateType();
  string outcomeType = model.options().outcomeType();
  vector<unsigned int> nFixedEffects_mix=model.dataset().nFixedEffects_mix();
  vector<unsigned int> nFixedEffects=model.dataset().nFixedEffects();
  unsigned int nCategoriesY=model.dataset().nCategoriesY();
  unsigned int nOutcomes=model.dataset().nOutcomes();

  randomUniform unifRand(0,1);

  // Move 1 - swap labels of 2 randomly selected non-empty clusters,
  //          leaving psi_c^prop = psi_c for all c

  // Compute how many non-empty clusters
  unsigned int nNotEmpty=0;
  vector<unsigned int> nonEmptyIndices;
  for(unsigned int c=0;c<=maxZ;c++){
    if(currentParams.workNXInCluster(c)>0){
      nNotEmpty++;
      nonEmptyIndices.push_back(c);
    }
  }

  // Select two non-empty clusters at random
  nTry++;
  unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
  unsigned int c1=nonEmptyIndices[i1];
  nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);
  unsigned int i2=(unsigned int)(nNotEmpty-1)*unifRand(rndGenerator);
  unsigned int c2=nonEmptyIndices[i2];

  // Check whether we accept the move
  double logAcceptRatio = ((double)currentParams.workNXInCluster(c2)-
                           (double)currentParams.workNXInCluster(c1))
    *(currentParams.logPsi(c1)-currentParams.logPsi(c2));

  double uii=unifRand(rndGenerator);

  if(uii<exp(logAcceptRatio)){
    //		nAccept++;
    // Switch the labels
    currentParams.switchLabels(c1,c2,covariateType,outcomeType,varSelectType, nFixedEffects, nFixedEffects_mix, nCategoriesY, nOutcomes);
  }
  // Move 2 - swap labels of 2 randomly selected neighbouring clusters,
  //			also swapping the v at the same time
  //	nTry++;
  c1=(unsigned int)maxZ*unifRand(rndGenerator);

  logAcceptRatio=(double)currentParams.workNXInCluster(c1)*
    log(1-currentParams.v(c1+1))
    - (double)currentParams.workNXInCluster(c1+1)*log(1-currentParams.v(c1));

  uii=unifRand(rndGenerator);

  if(uii<exp(logAcceptRatio)){
    nAccept++;

    // Switch the labels
    currentParams.switchLabels(c1,c1+1,covariateType,outcomeType,varSelectType, nFixedEffects, nFixedEffects_mix, nCategoriesY, nOutcomes);

    // Also switch the v's
    double v1=currentParams.v(c1);
    double v2=currentParams.v(c1+1);
    double logPsi1=currentParams.logPsi(c1);
    double logPsi2=currentParams.logPsi(c1+1);

    currentParams.logPsi(c1,log(v2)+logPsi1-log(v1));
    currentParams.logPsi(c1+1,log(v1)+log(1-v2)+logPsi2-log(v2)-log(1-v1));
    currentParams.v(c1,v2);
    currentParams.v(c1+1,v1);

    if(c1==maxZ-1){
      // Just check if the maximum Z has now changed to be one less
      if(currentParams.workNXInCluster(c1+1)==0){
        currentParams.workMaxZi(c1);
        maxZ=c1;
      }
    }
  }

  // Move 3

  //	nTry++;
  c1=(unsigned int)maxZ*unifRand(rndGenerator);

  // Compute the acceptance ratio
  unsigned int sumNAfterC1Plus1=0;
  for(unsigned int c=c1+2;c<=maxZ;c++){
    sumNAfterC1Plus1+=currentParams.workNXInCluster(c);
  }
  double const1=0.0,const2=0.0;
  double alpha=currentParams.alpha();
  const1=(1.0+alpha+(double)currentParams.workNXInCluster(c1+1)+(double)sumNAfterC1Plus1)/
    (alpha+(double)currentParams.workNXInCluster(c1+1)+(double)sumNAfterC1Plus1);
  const2=(alpha+(double)currentParams.workNXInCluster(c1)+(double)sumNAfterC1Plus1)/
    (1.0+alpha+(double)currentParams.workNXInCluster(c1)+(double)sumNAfterC1Plus1);
  logAcceptRatio=(double)(currentParams.workNXInCluster(c1)+currentParams.workNXInCluster(c1+1))*
    log(exp(currentParams.logPsi(c1))+exp(currentParams.logPsi(c1+1)));
  logAcceptRatio-=(double)(currentParams.workNXInCluster(c1)+currentParams.workNXInCluster(c1+1))*
    log(exp(currentParams.logPsi(c1))*const1+exp(currentParams.logPsi(c1+1))*const2);
  logAcceptRatio+=(double)(currentParams.workNXInCluster(c1+1))*log(const1);
  logAcceptRatio+=(double)(currentParams.workNXInCluster(c1))*log(const2);

  uii=unifRand(rndGenerator);

  if(uii<exp(logAcceptRatio)){
    //		nAccept++;
    currentParams.switchLabels(c1,c1+1,covariateType,outcomeType,varSelectType, nFixedEffects, nFixedEffects_mix, nCategoriesY, nOutcomes);
    double currPsiC1 = exp(currentParams.logPsi(c1));
    double currPsiC1Plus1 = exp(currentParams.logPsi(c1+1));
    double sumCurrPsi = currPsiC1+currPsiC1Plus1;
    double normConst = sumCurrPsi/(const1*currPsiC1Plus1+const2*currPsiC1);
    double propPsiC1 = normConst*const1*currPsiC1Plus1;
    double propPsiC1Plus1 = normConst*const2*currPsiC1;

    double productPrev1MinusV = 1.0;
    if(c1>0){
      productPrev1MinusV = exp(currentParams.logPsi(c1-1))*(1-currentParams.v(c1-1))/
        currentParams.v(c1-1);
    }

    double propVC1=propPsiC1/productPrev1MinusV;
    double propVC1Plus1=propPsiC1Plus1/(productPrev1MinusV*(1-propVC1));

    currentParams.logPsi(c1,log(propPsiC1));
    currentParams.logPsi(c1+1,log(propPsiC1Plus1));
    currentParams.v(c1,propVC1);
    currentParams.v(c1+1,propVC1Plus1);
    if(c1==maxZ-1){
      // Just check if the maximum Z has now changed to be one less
      if(currentParams.workNXInCluster(c1+1)==0){
        currentParams.workMaxZi(c1);
      }
    }
  }
}

void metropolisHastingsForLabels12(mcmcChain<pReMiuMParams>& chain,
                                   unsigned int& nTry,unsigned int& nAccept,
                                   const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                                   pReMiuMPropParams& propParams,
                                   baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  vector<unsigned int> nFixedEffects_mix=model.dataset().nFixedEffects_mix();
  vector<unsigned int> nFixedEffects=model.dataset().nFixedEffects();
  unsigned int nOutcomes=model.dataset().nOutcomes();
  unsigned int nCategoriesY=model.dataset().nCategoriesY();

  unsigned int maxZ = currentParams.workMaxZi();
  if(maxZ==0){
    // If there is only one cluster with individuals in, don't do anything
    return;
  }
  string varSelectType = model.options().varSelectType();
  string covariateType = model.options().covariateType();
  string outcomeType = model.options().outcomeType();

  randomUniform unifRand(0,1);

  // Move 1 - swap labels of 2 randomly selected non-empty clusters,
  //          leaving psi_c^prop = psi_c for all c

  // Compute how many non-empty clusters
  unsigned int nNotEmpty=0;
  vector<unsigned int> nonEmptyIndices;
  for(unsigned int c=0;c<=maxZ;c++){
    if(currentParams.workNXInCluster(c)>0){
      nNotEmpty++;
      nonEmptyIndices.push_back(c);
    }
  }

  // Select two non-empty clusters at random
  nTry++;
  unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
  unsigned int c1=nonEmptyIndices[i1];
  nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);
  unsigned int i2=(unsigned int)(nNotEmpty-1)*unifRand(rndGenerator);
  unsigned int c2=nonEmptyIndices[i2];

  // Check whether we accept the move
  double logAcceptRatio = ((double)currentParams.workNXInCluster(c2)-(double)currentParams.workNXInCluster(c1))
    *(currentParams.logPsi(c1)-currentParams.logPsi(c2));

  if(unifRand(rndGenerator)<exp(logAcceptRatio)){
    nAccept++;
    // Switch the labels
    currentParams.switchLabels(c1,c2,covariateType,outcomeType,varSelectType, nFixedEffects, nFixedEffects_mix, nCategoriesY, nOutcomes);
  }

  // Move 2 - swap labels of 2 randomly selected neighbouring clusters,
  //			also swapping the v at the same time
  //	nTry++;
  c1=(unsigned int)maxZ*unifRand(rndGenerator);

  logAcceptRatio=(double)currentParams.workNXInCluster(c1)*log(1-currentParams.v(c1+1))
    - (double)currentParams.workNXInCluster(c1+1)*log(1-currentParams.v(c1));

  if(unifRand(rndGenerator)<exp(logAcceptRatio)){
    //		nAccept++;
    // Switch the labels
    currentParams.switchLabels(c1,c1+1,covariateType,outcomeType,varSelectType, nFixedEffects, nFixedEffects_mix, nCategoriesY, nOutcomes);

    // Also switch the v's
    double v1=currentParams.v(c1);
    double v2=currentParams.v(c1+1);
    double logPsi1=currentParams.logPsi(c1);
    double logPsi2=currentParams.logPsi(c1+1);

    currentParams.logPsi(c1,log(v2)+logPsi1-log(v1));
    currentParams.logPsi(c1+1,log(v1)+log(1-v2)+logPsi2-log(v2)-log(1-v1));
    currentParams.v(c1,v2);
    currentParams.v(c1+1,v1);

    if(c1==maxZ-1){
      // Just check if the maximum Z has now changed to be one less
      if(currentParams.workNXInCluster(c1+1)==0){
        currentParams.workMaxZi(c1);
      }
    }


  }

}

void metropolisHastingsForLabels3(mcmcChain<pReMiuMParams>& chain,
                                  unsigned int& nTry,unsigned int& nAccept,
                                  const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                                  pReMiuMPropParams& propParams,
                                  baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  vector<unsigned int> nFixedEffects_mix=model.dataset().nFixedEffects_mix();
  vector<unsigned int> nFixedEffects=model.dataset().nFixedEffects();
  unsigned int nCategoriesY=model.dataset().nCategoriesY();
  unsigned int nOutcomes=model.dataset().nOutcomes();

  unsigned int maxZ = currentParams.workMaxZi();
  if(maxZ==0){
    // If there is only one cluster with individuals in, don't do anything
    return;
  }
  string varSelectType = model.options().varSelectType();
  string covariateType = model.options().covariateType();
  string outcomeType = model.options().outcomeType();

  randomUniform unifRand(0,1);

  // Move 3

  // Compute how many non-empty clusters
  unsigned int nNotEmpty=0;
  vector<unsigned int> nonEmptyIndices;
  for(unsigned int c=0;c<=maxZ;c++){
    if(currentParams.workNXInCluster(c)>0){
      nNotEmpty++;
      nonEmptyIndices.push_back(c);
    }
  }

  // Select two non-empty clusters at random
  nTry++;
  unsigned int i1=(unsigned int)nNotEmpty*unifRand(rndGenerator);
  unsigned int c1=nonEmptyIndices[i1];
  nonEmptyIndices.erase(nonEmptyIndices.begin()+i1);

  // Check whether we accept the move
  double logAcceptRatio=0;

  c1=(unsigned int)maxZ*unifRand(rndGenerator);
  // Compute the acceptance ratio
  unsigned int sumNAfterC1Plus1=0;
  for(unsigned int c=c1+2;c<=maxZ;c++){
    sumNAfterC1Plus1+=currentParams.workNXInCluster(c);
  }
  double const1=0.0,const2=0.0;
  double alpha=currentParams.alpha();
  const1=(1.0+alpha+(double)currentParams.workNXInCluster(c1+1)+(double)sumNAfterC1Plus1)/
    (alpha+(double)currentParams.workNXInCluster(c1+1)+(double)sumNAfterC1Plus1);
  const2=(alpha+(double)currentParams.workNXInCluster(c1)+(double)sumNAfterC1Plus1)/
    (1.0+alpha+(double)currentParams.workNXInCluster(c1)+(double)sumNAfterC1Plus1);
  logAcceptRatio=(double)(currentParams.workNXInCluster(c1)+currentParams.workNXInCluster(c1+1))*
    log(exp(currentParams.logPsi(c1))+exp(currentParams.logPsi(c1+1)));
  logAcceptRatio-=(double)(currentParams.workNXInCluster(c1)+currentParams.workNXInCluster(c1+1))*
    log(exp(currentParams.logPsi(c1+1))*const1+exp(currentParams.logPsi(c1))*const2);
  logAcceptRatio+=(double)(currentParams.workNXInCluster(c1+1))*log(const1);
  logAcceptRatio+=(double)(currentParams.workNXInCluster(c1))*log(const2);

  if(unifRand(rndGenerator)<exp(logAcceptRatio)){
    nAccept++;
    currentParams.switchLabels(c1,c1+1,covariateType,outcomeType,varSelectType, nFixedEffects, nFixedEffects_mix, nCategoriesY, nOutcomes);
    double currPsiC1 = exp(currentParams.logPsi(c1));
    double currPsiC1Plus1 = exp(currentParams.logPsi(c1+1));
    double sumCurrPsi = currPsiC1+currPsiC1Plus1;
    double normConst = sumCurrPsi/(const1*currPsiC1Plus1+const2*currPsiC1);
    double propPsiC1 = normConst*const1*currPsiC1Plus1;
    double propPsiC1Plus1 = normConst*const2*currPsiC1;

    double productPrev1MinusV = 1.0;
    if(c1>0){
      productPrev1MinusV = exp(currentParams.logPsi(c1-1))*(1-currentParams.v(c1-1))/
        currentParams.v(c1-1);
    }

    double propVC1=propPsiC1/productPrev1MinusV;
    double propVC1Plus1=propPsiC1Plus1/(productPrev1MinusV*(1-propVC1));

    currentParams.logPsi(c1,log(propPsiC1));
    currentParams.logPsi(c1+1,log(propPsiC1Plus1));
    currentParams.v(c1,propVC1);
    currentParams.v(c1+1,propVC1Plus1);

    if(c1==maxZ-1){
      // Just check if the maximum Z has now changed to be one less
      if(currentParams.workNXInCluster(c1+1)==0){
        currentParams.workMaxZi(c1);
      }
    }
  }

}


// Gibbs move for updating the auxiliary variables u
// This is the second part of block 1. The first part used the marginal
// distribution with u integrated out, we now use the conditional distribution
// for u, conditional on the v^A,Theta^A parameters generated above
// This is done by using Gibbs to sample from p(u|v^A.), which is the conditional

void gibbsForU(mcmcChain<pReMiuMParams>& chain,
               unsigned int& nTry,unsigned int& nAccept,
               const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
               pReMiuMPropParams& propParams,
               baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  string samplerType = model.options().samplerType();

  nTry++;
  nAccept++;

  unsigned int nSubjects = currentParams.nSubjects();
  unsigned int nPredictSubjects = currentParams.nPredictSubjects();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  double minUi = 1.0;
  for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
    int zi = currentParams.z(i);
    double ui=unifRand(rndGenerator);
    if(samplerType.compare("SliceDependent")==0){
      ui*=exp(currentParams.logPsi((unsigned int)zi));
    }else if(samplerType.compare("SliceIndependent")==0){
      ui*=hyperParams.workXiSlice((unsigned int)zi);
    }

    // This is to avoid numerical errors be
    if(ui<0.0000000001){
      ui=0.0000000001;
    }

    // We only take the minimum over the fitting subjects, because
    // we don't allow predictiction subjects to be allocated to clusters
    // where there are no fitting members, so no need to calc probabilities
    // for cluster which are potential only for predictions subjects
    if(ui<minUi&&i<nSubjects){
      minUi=ui;
    }
    currentParams.u(i,ui);
  }
  currentParams.workMinUi(minUi);
}

/*********** BLOCK 2 p(alpha,v^I|.) **********************************/
// I=Inactive
// We proceed by sampling alpha from p(alpha|.) i.e. marginalising out
// v^I. Then we sample from p(v^I|alpha,.)

// Adaptive Metropolis Hastings move for alpha
void metropolisHastingsForAlpha(mcmcChain<pReMiuMParams>& chain,
                                unsigned int& nTry,unsigned int& nAccept,
                                const mcmcModel<pReMiuMParams,
                                                pReMiuMOptions,
                                                pReMiuMData>& model,
                                                pReMiuMPropParams& propParams,
                                                baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  double& stdDev = propParams.alphaStdDev();
  double alphaCurrent = currentParams.alpha();

  double alphaProp;
  alphaProp = truncNormalRand(rndGenerator,alphaCurrent,stdDev,"L",0,0);

  double logAcceptRatio = 0.0;
  for(unsigned int c=0;c<=maxZ;c++){
    double v=currentParams.v(c);
    logAcceptRatio += logPdfBeta(v,1.0,alphaProp)-logPdfBeta(v,1.0,
                                 alphaCurrent);
  }

  // Add in the prior contribution (no contribution if uniform)
  logAcceptRatio += logPdfGamma(alphaProp,hyperParams.shapeAlpha(),hyperParams.rateAlpha());
  logAcceptRatio -= logPdfGamma(alphaCurrent,hyperParams.shapeAlpha(),hyperParams.rateAlpha());


  // Add the proposal contribution
  logAcceptRatio += logPdfTruncatedNormal(alphaCurrent,alphaProp,stdDev,"L",0,0);
  logAcceptRatio -= logPdfTruncatedNormal(alphaProp,alphaCurrent,stdDev,"L",0,0);

  propParams.alphaAddTry();
  nTry++;
  double uni = unifRand(rndGenerator);
  if(uni <exp(logAcceptRatio)){
    nAccept++;
    propParams.alphaAddAccept();
    // If the move was accepted update the state
    currentParams.alpha(alphaProp);

    // Also update the proposal standard deviation
    if(propParams.nTryAlpha()%propParams.alphaUpdateFreq()==0){
      stdDev += 10*(propParams.alphaLocalAcceptRate()-propParams.alphaAcceptTarget())/
        pow((double)(propParams.nTryAlpha()/propParams.alphaUpdateFreq())+2.0,0.75);
      propParams.alphaAnyUpdates(true);
      if(stdDev>propParams.alphaStdDevUpper()||stdDev<propParams.alphaStdDevLower()){
        propParams.alphaStdDevReset();
      }
      propParams.alphaLocalReset();
    }
  }else{
    // Otherwise update the proposal standard deviation
    if(propParams.nTryAlpha()%propParams.alphaUpdateFreq()==0){
      stdDev += 10*(propParams.alphaLocalAcceptRate()-propParams.alphaAcceptTarget())/
        pow((double)(propParams.nTryAlpha()/propParams.alphaUpdateFreq())+2.0,0.75);
      propParams.alphaAnyUpdates(true);
      if(stdDev>propParams.alphaStdDevUpper()||stdDev<propParams.alphaStdDevLower()){
        propParams.alphaStdDevReset();
      }
      propParams.alphaLocalReset();
    }
  }
}

// Gibbs move for v which are inactive. Only update to maxNClusters, which
// needs to be computed here
void gibbsForVInActive(mcmcChain<pReMiuMParams>& chain,
                       unsigned int& nTry,unsigned int& nAccept,
                       const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                       pReMiuMPropParams& propParams,
                       baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  string samplerType = model.options().samplerType();
  string covariateType = model.options().covariateType();
  string outcomeType = model.options().outcomeType();
  string kernelType = model.options().kernelType(); //AR
  unsigned int  nTimes_unique = model.dataset().nTimes_unique(); //AR
  vector<unsigned int>  nRandomEffects = model.dataset().nRandomEffects();
  nTry++;
  nAccept++;

  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();

  double minUi = currentParams.workMinUi();

  vector<double> vNew=currentParams.v();
  vector<double> logPsiNew=currentParams.logPsi();

  double alpha = currentParams.alpha();
  double dPitmanYor = currentParams.dPitmanYor();

  if(samplerType.compare("Truncated")==0){
    // Just sample from the prior

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
  }else{
    if(samplerType.compare("SliceIndependent")==0){
      maxNClusters=2+(int)((log(minUi)-log(1.0-hyperParams.rSlice()))/log(hyperParams.rSlice()));
    }

    // Sample V
    vector<double> cumPsi(maxZ+1,0.0);
    cumPsi[0] = exp(currentParams.logPsi(0));
    for(unsigned int c=1;c<=maxZ;c++){
      cumPsi[c]=cumPsi[c-1]+exp(currentParams.logPsi(c));
    }
    bool continueLoop=true;
    unsigned int c=maxZ;

    while(continueLoop){
      if(samplerType.compare("SliceDependent")==0&&cumPsi[c]>(1-minUi)){
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
    currentParams.maxNClusters(maxNClusters,covariateType,outcomeType,kernelType,nTimes_unique, nRandomEffects);
  }

  currentParams.v(vNew);
  currentParams.logPsi(logPsiNew);
}

/*********** BLOCK 3 p(Theta^I|.) **********************************/
// I=Inactive. Sample the inactive cluster variables from the prior
// Theta contains phi, mu, Tau, gamma, theta. Only need to sample
// up to maxNClusters = max_i{Ci}. Several different routines here for
// each of the variables

// Gibbs move for updating phi
void gibbsForPhiInActive(mcmcChain<pReMiuMParams>& chain,
                         unsigned int& nTry,unsigned int& nAccept,
                         const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                         pReMiuMPropParams& propParams,
                         baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  string varSelectType = model.options().varSelectType();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();

  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nDiscreteCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }

  nTry++;
  nAccept++;

  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    // Loop over the covariates
    for(unsigned int j=0;j<nCovariates;j++){
      unsigned int nCategories = currentParams.nCategories(j);
      vector<double> dirichParams(nCategories,hyperParams.aPhi(j));
      vector<double> proposedLogPhi(nCategories);
      proposedLogPhi=dirichletRand(rndGenerator,dirichParams);

      for(unsigned int p=0;p<nCategories;p++){
        proposedLogPhi[p]=log(proposedLogPhi[p]);
      }
      currentParams.logPhi(c,j,proposedLogPhi);
    }
  }
}


// Gibbs update for mu in Normal covariate case
void gibbsForMuInActive(mcmcChain<pReMiuMParams>& chain,
                        unsigned int& nTry,unsigned int& nAccept,
                        const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                        pReMiuMPropParams& propParams,
                        baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nContinuousCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }

  nTry++;
  nAccept++;

  MatrixXd covMat(nCovariates,nCovariates);
  covMat = hyperParams.Tau0().inverse();
  VectorXd meanVec(nCovariates);
  meanVec = hyperParams.mu0();

  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    VectorXd mu(nCovariates);
    // We sample from this posterior
    mu = multivarNormalRand(rndGenerator,meanVec,covMat);
    // We store our sample
    currentParams.mu(c,mu);
  }

}

// Gibbs update for mu in Normal covariate case when the normal inverse Wishart prior is used
void gibbsForMuInActiveNIWP(mcmcChain<pReMiuMParams>& chain,
                            unsigned int& nTry,unsigned int& nAccept,
                            const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                            pReMiuMPropParams& propParams,
                            baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  // Find the number of covariates
  unsigned int nCovariates = 0;
  if(model.options().covariateType().compare("Mixed")==0){
    nCovariates = currentParams.nContinuousCovs();
  } else {
    nCovariates = currentParams.nCovariates();
  }

  nTry++;
  nAccept++;


  VectorXd meanVec(nCovariates);
  meanVec = hyperParams.mu0();

  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    MatrixXd covMat(nCovariates,nCovariates);
    covMat = currentParams.Sigma(c)/hyperParams.nu0();

    VectorXd mu(nCovariates);
    // We sample from this posterior
    mu = multivarNormalRand(rndGenerator,meanVec,covMat);
    // We store our sample
    currentParams.mu(c,mu);
  }

}

// Gibbs update for Tau in the Normal covariate case
void gibbsForTauInActive(mcmcChain<pReMiuMParams>& chain,
                         unsigned int& nTry,unsigned int& nAccept,
                         const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                         pReMiuMPropParams& propParams,
                         baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();

  nTry++;
  nAccept++;

  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    MatrixXd Tau = wishartRand(rndGenerator,hyperParams.R0(),hyperParams.kappa0());
    currentParams.Tau(c,Tau);
  }

}

void gibbsForMVNMuInActive(mcmcChain<pReMiuMParams>& chain,
                           unsigned int& nTry,unsigned int& nAccept,
                           const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                           pReMiuMPropParams& propParams,
                           baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  unsigned int nOutcomes = currentParams.nOutcomes();
  nTry++;
  nAccept++;
  VectorXd meanVec(nOutcomes);
  meanVec = hyperParams.MVNmu0();
  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    MatrixXd covMat(nOutcomes,nOutcomes);
    covMat = currentParams.MVNSigma(c)/hyperParams.MVNkappa0();
    VectorXd mu(nOutcomes);
    mu = multivarNormalRand(rndGenerator,meanVec,covMat);
    currentParams.MVNmu(c,mu);
  }
}
void gibbsForMVNTauInActive(mcmcChain<pReMiuMParams>& chain,
                            unsigned int& nTry,unsigned int& nAccept,
                            const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                            pReMiuMPropParams& propParams,
                            baseGeneratorType& rndGenerator){
  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  nTry++;
  nAccept++;
  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    MatrixXd Sigma = invWishartRand(rndGenerator,hyperParams.MVNR0(),hyperParams.MVNnu0());
    MatrixXd Tau = Sigma.inverse();
    //MatrixXd Tau = wishartRand(rndGenerator,hyperParams.MVNR0(),hyperParams.MVNnu0());
    currentParams.MVNTau(c,Tau);
  }
}


// Gibbs update for update of gamma (only used in the binary variable selection case)
void gibbsForGammaInActive(mcmcChain<pReMiuMParams>& chain,
                           unsigned int& nTry,unsigned int& nAccept,
                           const mcmcModel<pReMiuMParams,
                                           pReMiuMOptions,
                                           pReMiuMData>& model,
                                           pReMiuMPropParams& propParams,
                                           baseGeneratorType& rndGenerator){


  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of subjects
  unsigned int nCovariates = currentParams.nCovariates();
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  string covariateType = model.options().covariateType();
  string varSelectType = model.options().varSelectType();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  nTry++;
  nAccept++;

  for(unsigned int j=0;j<nCovariates;j++){

    for(unsigned int c=maxZ+1;c<maxNClusters;c++){

      double currentGamma=currentParams.gamma(c,j);
      double proposedGamma=0.0;
      // work in terms of prob of sticking with current value and switching value
      // Compute probability of sticking
      double logProbStick=0;
      double logProbSwitch=0;
      double probSwitch=0;
      if(currentParams.omega(j)==0){
        // Nothing to do - not allowed to change
        continue;
      }else{
        logProbStick+=(currentGamma*log(currentParams.rho(j))+
          (1-currentGamma)*log(1-currentParams.rho(j)));

        // Now compute probability of switching
        proposedGamma=1-currentGamma;

        logProbSwitch+=(proposedGamma*log(currentParams.rho(j))+
          (1-proposedGamma)*log(1-currentParams.rho(j)));

        double maxLogProb;
        if(logProbSwitch<logProbStick){
          maxLogProb=logProbStick;
        }else{
          maxLogProb=logProbSwitch;
        }


        probSwitch=exp(logProbSwitch-maxLogProb)/(exp(logProbStick-maxLogProb)+exp(logProbSwitch-maxLogProb));
      }
      if(unifRand(rndGenerator)<probSwitch){
        // Switching
        currentParams.gamma(c,j,proposedGamma,covariateType);
      }
    }
  }
}

// Gibbs for theta
void gibbsForThetaInActive(mcmcChain<pReMiuMParams>& chain,
                           unsigned int& nTry,unsigned int& nAccept,
                           const mcmcModel<pReMiuMParams,
                                           pReMiuMOptions,
                                           pReMiuMData>& model,
                                           pReMiuMPropParams& propParams,
                                           baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();
  unsigned int nCategoriesY=dataset.nCategoriesY();
  const string outcomeType = model.dataset().outcomeType();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();

  nTry++;
  nAccept++;

  double location = hyperParams.muTheta();
  double scale = hyperParams.sigmaTheta();
  unsigned int dof = hyperParams.dofTheta();
  randomStudentsT studentsTRand(dof);
  for (unsigned int k=0;k<nCategoriesY;k++){
    for(unsigned int c=maxZ+1;c<maxNClusters;c++){
      double theta=location+scale*studentsTRand(rndGenerator);
      currentParams.theta(c,k,theta);
    }
  }
}


void gibbsForBetaInActive(mcmcChain<pReMiuMParams>& chain,
                                  unsigned int& nTry,unsigned int& nAccept,
                                  const mcmcModel<pReMiuMParams,
                                                  pReMiuMOptions,
                                                  pReMiuMData>& model,
                                                  pReMiuMPropParams& propParams,
                                                  baseGeneratorType& rndGenerator){
  // Define a normal random number generator
  randomNormal normRand(0,1);

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();
  unsigned int nCategoriesY=dataset.nCategoriesY();
  unsigned int nOutcomes=dataset.nOutcomes();
  vector<unsigned int> nFixedEffects_mix=model.dataset().nFixedEffects_mix();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  double mu_b  = currentParams.hyperParams().muBeta();
  double sigma2beta = currentParams.hyperParams().sigmaBeta(); //variance

  for(unsigned int m=0;m<nOutcomes;m++){
    for(unsigned int j=0;j<nFixedEffects_mix[m];j++){
      for (unsigned int k=0;k<nCategoriesY;k++){
        for(unsigned int c=maxZ+1;c<maxNClusters;c++){
          double beta=mu_b+pow(sigma2beta,0.5)*normRand(rndGenerator);
          currentParams.beta_mix(c,j,k,m,nCategoriesY,beta);
        }
      }
    }
  }
}

// void gibbsForBetaInActive_student(mcmcChain<pReMiuMParams>& chain,
//                           unsigned int& nTry,unsigned int& nAccept,
//                           const mcmcModel<pReMiuMParams,
//                                           pReMiuMOptions,
//                                           pReMiuMData>& model,
//                                           pReMiuMPropParams& propParams,
//                                           baseGeneratorType& rndGenerator){
//   randomUniform unifRand(0,1);
//
//   mcmcState<pReMiuMParams>& currentState = chain.currentState();
//   pReMiuMParams& currentParams = currentState.parameters();
//   pReMiuMHyperParams hyperParams = currentParams.hyperParams();
//   const pReMiuMData& dataset = model.dataset();
//   unsigned int nCategoriesY=dataset.nCategoriesY();
//   const string outcomeType = model.dataset().outcomeType();
//   vector<unsigned int> nFixedEffects_mix=model.dataset().nFixedEffects_mix();
//
//   // Find the number of clusters
//   unsigned int maxZ = currentParams.workMaxZi();
//   unsigned int maxNClusters = currentParams.maxNClusters();
//   unsigned int nOutcomes = model.dataset().nOutcomes();
//
//   nTry++;
//   nAccept++;
//
//   double location = hyperParams.muBeta();
//   double scale = hyperParams.sigmaBeta();
//   unsigned int dof = hyperParams.dofBeta();
//   randomStudentsT studentsTRand(dof);
//   for(unsigned int m=0;m< nOutcomes;m++){
//     for(unsigned int j=0;j<nFixedEffects_mix[m];j++){
//       for (unsigned int k=0;k<nCategoriesY;k++){
//         for(unsigned int c=maxZ+1;c<maxNClusters;c++){
//           double beta=location+scale*studentsTRand(rndGenerator);
//           beta=-2.0+4.0*unifRand(rndGenerator);
//           currentParams.beta_mix(m, c,j,k,nCategoriesY,beta);
//           //currentParams.beta_mix(c,j+k*nFixedEffects_mix,0);
//           //        params.beta(j,k,-2.0+4.0*unifRand(rndGenerator));
//
//         }
//       }
//     }
//   }
// }

//RJ Gibbs for L
void gibbsForLInActive(mcmcChain<pReMiuMParams>& chain,
                       unsigned int& nTry,unsigned int& nAccept,
                       const mcmcModel<pReMiuMParams,
                                       pReMiuMOptions,
                                       pReMiuMData>& model,
                                       pReMiuMPropParams& propParams,
                                       baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  string kernelType =model.dataset().kernelType();
  const string outcomeType = model.dataset().outcomeType();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();
  bool ratio_estim = model.options().estim_ratio();

  nTry++;
  nAccept++;

  unsigned int nL;
  if(kernelType.compare("SQexponential")==0){
    nL=3;
  }else{
    nL=4;
  }

  if(model.options().estim_ratio()){
    for(unsigned int c=maxZ+1;c<maxNClusters;c++){
      double vVal = betaRand(rndGenerator,hyperParams.aRatio(),hyperParams.bRatio());
      currentParams.ratio(c, vVal);
    }
  }


  for (unsigned int l=0;l<nL;l++){
    for(unsigned int c=maxZ+1;c<maxNClusters;c++){
      if(!ratio_estim || l != 2){
        randomNormal normalRand(hyperParams.muL(l),hyperParams.sigmaL(l));
        double L = normalRand(rndGenerator);
        currentParams.L(c,l,L);
        if(ratio_estim && l == 0)
          currentParams.L(c,2,currentParams.L(c,0)+log(currentParams.ratio(c)));
      }
    }
  }

  if(model.options().sampleGPmean()){ //AR model.options().sampleGPmean()

     unsigned int nTimes_unique = model.dataset().nTimes_unique();
    // int threads = omp_get_max_threads();
    // vector<baseGeneratorType> rngArray(threads);

    // Seed by taking random numbers from the existing generator
    // boost::random::uniform_int_distribution<> seeder;
    // for (auto &rng : rngArray)
    //   rng.seed(seeder(rndGenerator));


    //#pragma omp parallel for
    for(unsigned int c=maxZ+1;c<maxNClusters;c++){
      //int threadNum = omp_get_thread_num();
      VectorXd Fval(nTimes_unique);
      //Fval =  Sample_GPmean(currentParams, model.dataset(), c, rngArray[threadNum],0);
      Fval =  Sample_GPmean(currentParams, model.dataset(), c, rndGenerator,0);
      for(unsigned int j=0;j<nTimes_unique;j++){
        currentParams.meanGP(c,j, Fval(j));
      }
    }
  }
}

// Gibbs for nu inactive (for survival case)
void gibbsForNuInActive(mcmcChain<pReMiuMParams>& chain,
                        unsigned int& nTry,unsigned int& nAccept,
                        const mcmcModel<pReMiuMParams,
                                        pReMiuMOptions,
                                        pReMiuMData>& model,
                                        pReMiuMPropParams& propParams,
                                        baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const string outcomeType = model.dataset().outcomeType();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int maxNClusters = currentParams.maxNClusters();

  nTry++;
  nAccept++;

  randomGamma gammaRand(hyperParams.shapeNu(),hyperParams.scaleNu());
  for(unsigned int c=maxZ+1;c<maxNClusters;c++){
    double nu=gammaRand(rndGenerator);
    currentParams.nu(c,nu);
  }
}

/*********** BLOCK 4 p(Theta^N|.) **********************************/
// N=Non-cluster, and Theta contains: beta, rho, omega, lambda, tau_epsilon, uCAR and TauCAR

// void GibbsForBeta_multiv(mcmcChain<pReMiuMParams>& chain,
//                   unsigned int& nTry,unsigned int& nAccept,
//                   const mcmcModel<pReMiuMParams,
//                                   pReMiuMOptions,
//                                   pReMiuMData>& model,
//                                   pReMiuMPropParams& propParams,
//                                   baseGeneratorType& rndGenerator){
//
//
//   mcmcState<pReMiuMParams>& currentState = chain.currentState();
//   pReMiuMParams& currentParams = currentState.parameters();
//   const pReMiuMData& dataset = model.dataset();
//   const string outcomeType = dataset.outcomeType();
//
//   // Find the number of clusters
//   vector<unsigned int> nFixedEffects = dataset.nFixedEffects();
//   vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
//   unsigned int nCategoriesY = dataset.nCategoriesY();
//   unsigned int nOutcomes = dataset.nOutcomes();
//
//   // Define a normal random number generator
//   randomNormal normRand(0,1);
//   int k=0; // Works if only 1 category for Y
//
//   unsigned int nSubjects=dataset.nSubjects();
//   vector<double> y = dataset.continuousY();
//   const string kernelType = dataset.kernelType(); //AR
//   vector<double> times = dataset.times();
//   vector<int> tStart = dataset.tStart();
//   vector<int> tStop = dataset.tStop();
//
//   unsigned int ind =0;
//   unsigned int ind_y=0;
//   for(unsigned int m=0;m<nOutcomes;m++){
//     if(nFixedEffects[m]>0){
//
//       VectorXd mu_b(nFixedEffects[m]);
//       MatrixXd V_b;
//       V_b.setZero(nFixedEffects[m],nFixedEffects[m]);
//
//       MatrixXd S;
//       S.setZero(nFixedEffects[m],nFixedEffects[m]);
//       VectorXd S2(nFixedEffects[m]);
//
//       for(unsigned int i=0;i<nFixedEffects[m];i++){
//         mu_b(i)  = currentParams.hyperParams().muBeta();
//         V_b(i,i) = currentParams.hyperParams().sigmaBeta();
//         S2(i)    = 0;
//       }
//
//
//       for(unsigned int i=0;i<nSubjects;i++){
//
//         int nmes = (tStop[ind] - tStart[ind] + 1);
//         VectorXd   Yi(nmes);
//         int zi   = currentParams.z(i);
//         MatrixXd   Xi(nmes, nFixedEffects[m]);
//
//         for(unsigned int j=0;j<tStop[ind]-tStart[ind]+1;j++){
//           Yi(j) = y[ind_y + tStart[ind]-1+j];
//
//           for(unsigned int b=0;b<nFixedEffects[m];b++){
//             Xi(j,b) = dataset.W_LME(m,tStart[ind]-1+j,b);
//           }
//
//           for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
//             Yi(j) -= currentParams.beta_mix(m,zi,b,0,nCategoriesY)*dataset.W_LME_mix(m,tStart[ind]-1+j,b);
//           }
//         }
//
//         MatrixXd block=dataset.W_RE(m,tStart[ind]-1, 0, tStop[ind]-tStart[ind]+1 , dataset.nRandomEffects(m));
//         Yi -= block*currentParams.RandomEffects(m,i);
//
//         MatrixXd Sigmae_inv = MatrixXd::Identity(nmes, nmes) * 1/currentParams.SigmaE(m);
//         S += Xi.transpose()*Sigmae_inv*Xi;
//         S2 += Xi.transpose()*Sigmae_inv*Yi;
//         ind++;
//         if(i==(nSubjects-1))
//           ind_y += nSubjects;
//       }
//
//       MatrixXd temp = V_b.inverse() + S;
//       LLT<MatrixXd> lltOfA(temp); // compute the Cholesky decomposition of A
//       MatrixXd L = lltOfA.matrixL();
//       MatrixXd V_b_post = L.inverse().transpose()*L.inverse();
//
//       VectorXd mu_b_post = V_b_post * (V_b.inverse()*mu_b + S2);
//
//       VectorXd betaProp=multivarNormalRand(rndGenerator, mu_b_post, V_b_post);
//       for(unsigned int b=0;b<nFixedEffects[m];b++)
//         currentParams.beta(m,b,k,nCategoriesY,betaProp(b));
//     }
//   }
//
//   ind = 0;
//   ind_y = 0;
//   for(unsigned int m=0;m<nOutcomes;m++){
//     if(nFixedEffects_mix[m]>0){
//       unsigned int maxZ = currentParams.workMaxZi();
//       VectorXd mu_b(nFixedEffects_mix);
//       MatrixXd V_b;
//       V_b.setZero(nFixedEffects_mix[m],nFixedEffects_mix[m]);
//
//       for(unsigned int i=0;i<nFixedEffects_mix[m];i++){
//         mu_b(i)  = currentParams.hyperParams().muBeta();
//         V_b(i,i) = currentParams.hyperParams().sigmaBeta();
//       }
//
//       for(unsigned int c=0;c<=maxZ;c++){
//
//         MatrixXd S;
//         S.setZero(nFixedEffects_mix[m], nFixedEffects_mix[m]);
//         VectorXd S2(nFixedEffects_mix[m]);
//
//         for(unsigned int i=0;i<nFixedEffects_mix[m];i++)
//           S2(i)    = 0;
//
//         for(unsigned int i=0;i<nSubjects;i++){
//           int zi   = currentParams.z(i);
//
//           if(zi == c){
//             int nmes = (tStop[ind] - tStart[ind] + 1);
//             VectorXd   Yi(nmes);
//             MatrixXd   Xi(nmes, nFixedEffects_mix[m]);
//
//             for(unsigned int j=0;j<tStop[ind]-tStart[ind]+1;j++){
//               Yi(j) = y[tStart[ind]-1+j];
//
//               for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
//                 Xi(j,b) = dataset.W_LME_mix(m,tStart[ind]-1+j,b);
//               }
//
//               for(unsigned int b=0;b<nFixedEffects[m];b++){
//                 Yi(j) -= currentParams.beta(m,b,k, nCategoriesY)*dataset.W_LME(m,tStart[ind]-1+j,b);
//               }
//             }
//
//             MatrixXd block=dataset.W_RE(m,tStart[ind]-1, 0, tStop[ind]-tStart[ind]+1 , dataset.nRandomEffects(m));
//             Yi -= block*currentParams.RandomEffects(m,i);
//
//             MatrixXd Sigmae_inv = MatrixXd::Identity(nmes, nmes) * 1/currentParams.SigmaE(m);
//             S += Xi.transpose()*Sigmae_inv*Xi;
//             S2 += Xi.transpose()*Sigmae_inv*Yi;
//           }
//           ind++;
//           if(i==(nSubjects-1))
//             ind_y += nSubjects;
//         }
//         MatrixXd temp = V_b.inverse() + S;
//         LLT<MatrixXd> lltOfA(temp); // compute the Cholesky decomposition of A
//         MatrixXd L = lltOfA.matrixL();
//         MatrixXd V_b_post = L.inverse().transpose()*L.inverse();
//
//         VectorXd mu_b_post = V_b_post * (V_b.inverse()*mu_b + S2);
//
//
//         VectorXd beta_mixProp=multivarNormalRand(rndGenerator, mu_b_post, V_b_post);
//         for(unsigned int b=0;b<nFixedEffects_mix;b++)
//           currentParams.beta_mix(m,c,b,0,nCategoriesY, beta_mixProp(b));
//       }
//
//       unsigned int maxNClusters=currentParams.maxNClusters();
//       for(unsigned int c=maxZ+1;c<maxNClusters;c++){
//         VectorXd beta_mixProp=multivarNormalRand(rndGenerator, mu_b, V_b);
//         for(unsigned int b=0;b<nFixedEffects_mix;b++)
//           currentParams.beta_mix(m,c,b,0,nCategoriesY,beta_mixProp(b));
//       }
//     }
//   }
// }

void GibbsForBeta(mcmcChain<pReMiuMParams>& chain,
                  unsigned int& nTry,unsigned int& nAccept,
                  const mcmcModel<pReMiuMParams,
                                  pReMiuMOptions,
                                  pReMiuMData>& model,
                                  pReMiuMPropParams& propParams,
                                  baseGeneratorType& rndGenerator){


  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  const pReMiuMData& dataset = model.dataset();
  const string outcomeType = dataset.outcomeType();

  // Find the number of clusters
  vector<unsigned int> nFixedEffects = dataset.nFixedEffects();
  vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();

  // Define a normal random number generator
  randomNormal normRand(0,1);

  unsigned int nSubjects=dataset.nSubjects();
  vector<double> y = dataset.continuousY();
  const string kernelType = dataset.kernelType(); //AR
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  unsigned int nCategoriesY=dataset.nCategoriesY();

  unsigned int maxNClusters=currentParams.maxNClusters();
  unsigned int maxZ = currentParams.workMaxZi();
  unsigned int nOutcomes = dataset.nOutcomes();

  double mu_b  = currentParams.hyperParams().muBeta();
  double sigma2beta = currentParams.hyperParams().sigmaBeta(); //variance

  unsigned int ind_y=0;
  unsigned int ind=0;

  for (unsigned int m=0;m<nOutcomes;m++){
    for(unsigned int b=0;b<(nFixedEffects[m] + nFixedEffects_mix[m]);b++){


      MatrixXd S2;
      S2.setZero(1,1);
      MatrixXd S;
      S.setZero(1,1);

      for(unsigned int i=0;i<nSubjects;i++){

        int zi   = currentParams.z(i);
        int nmes = (tStop[ind] - tStart[ind] + 1);
        MatrixXd   Xi(nmes, nFixedEffects[m]+nFixedEffects_mix[m]-1); // Xi\{\beta}
        VectorXd Xib(nmes); // Xi{\beta}
        VectorXd Yi(nmes); // XY{\beta}

        for(unsigned int j=0;j<tStop[ind]-tStart[ind]+1;j++){
          Yi(j) = y[ind_y + tStart[ind]-1+j];

          if(i==0 & j==0 & b ==0)
            std::cout << " Yi(j) "<<Yi(j)<<endl;


          for(unsigned int bb=0;bb<nFixedEffects[m];bb++){
            if(bb!=b)
              Yi(j) -= currentParams.beta(m,bb,0,nCategoriesY)*dataset.W_LME(m,tStart[ind]-1+j,bb);
          }

          for(unsigned int bb=0;bb<nFixedEffects_mix[m];bb++){
            if((bb+nFixedEffects[m])!=b)
              Yi(j) -= currentParams.beta_mix(m,zi,bb,0,nCategoriesY)*dataset.W_LME_mix(m,tStart[ind]-1+j,bb);

            if(i==0 & j==0 & b ==0 & bb==0){
              std::cout << bb<<" beta_mix(c,bb) "<<currentParams.beta_mix(m,zi,bb,0,nCategoriesY)<< " ";
              std::cout << bb<<" W_LME_mix(ij,bb) "<<dataset.W_LME_mix(m,tStart[ind]-1+j,bb)<<endl;
            }
          }

          for(unsigned int bb=0;bb<nFixedEffects[m];bb++){
            if(bb==b)
              Xib(j) = dataset.W_LME(m,tStart[ind]-1+j,bb);
          }
          for(unsigned int bb=0;bb<nFixedEffects_mix[m];bb++){
            if((nFixedEffects[m]+bb)==b)
              Xib(j) = dataset.W_LME_mix(m,tStart[ind]-1+j,bb);
            if(i==0 & j==0 & b ==0 & (nFixedEffects[m]+bb)==b)
              std::cout << " W_LME_mix "<<dataset.W_LME_mix(m,tStart[ind]-1+j,bb)<<endl;

          }
        }

        MatrixXd block=dataset.W_RE(m,tStart[ind]-1, 0, tStop[ind]-tStart[ind]+1 , dataset.nRandomEffects(m));
        Yi -= block*currentParams.RandomEffects(m,i);

        S += 1/pow(currentParams.SigmaE(m),nmes)*Xib.transpose()*Yi;
        S2 += 1/pow(currentParams.SigmaE(m),nmes)*Xib.transpose()*Xib;

        if(i==0 & b ==0 )
          std::cout << " Yi(0) "<<Yi(0)<< endl<<" S(0,0) "<< S(0,0) << " S2(0,0) "<< S2(0,0) <<endl;

        if(i==(nSubjects-1))
          ind_y += nSubjects;
        ind++;
      }

      S2(0,0) += 1/sigma2beta;
      S2(0,0) = 1 / S2(0,0);

      if(b ==0 )
        std::cout <<" S(0,0) "<< S(0,0) << " S2(0,0) "<< S2(0,0)  << "S*S2 "<<S(0,0) *S2(0,0)<<endl;


      S(0,0) *=S2(0,0);

      double betaProp = S(0,0)+pow(S2(0,0),0.5)*normRand(rndGenerator);

      if(b<nFixedEffects[m]){
        currentParams.beta(m,b,0,nCategoriesY,betaProp);
      }else{
        for(unsigned int c=0;c<=maxZ;c++){
          betaProp = S(0,0)+pow(S2(0,0),0.5)*normRand(rndGenerator);
          currentParams.beta_mix(m,c,b-nFixedEffects[m],0, nCategoriesY, betaProp);
          std::cout << b<< " c "<<c<<" beta_mix " <<currentParams.beta_mix(m,c,b-nFixedEffects[m],0, nCategoriesY)<<endl;
        }
      }

      if(b>=nFixedEffects[m]){
        for(unsigned int c=maxZ+1;c<maxNClusters;c++){
          double betaProp = mu_b+pow(sigma2beta,0.5)*normRand(rndGenerator);
          currentParams.beta_mix(m,c,b-nFixedEffects[m],0, nCategoriesY, betaProp);
        }
      }
    }
  }
}

// Adaptive Metropolis-Hastings for beta - fixed effects
void metropolisHastingsForBeta(mcmcChain<pReMiuMParams>& chain,
                               unsigned int& nTry,unsigned int& nAccept,
                               const mcmcModel<pReMiuMParams,
                                               pReMiuMOptions,
                                               pReMiuMData>& model,
                                               pReMiuMPropParams& propParams,
                                               baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  const string outcomeType = model.dataset().outcomeType();

  // Find the number of clusters
  unsigned int nFixedEffects = model.dataset().nFixedEffects(0);
  //unsigned int nFixedEffects_mix = model.dataset().nFixedEffects_mix(0);

  // Find the number of categories of Y
  unsigned int nCategoriesY = currentParams.nCategoriesY();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);
  // Define a normal random number generator
  randomNormal normRand(0,1);

  double betaTargetRate = propParams.betaAcceptTarget();
  unsigned int betaUpdateFreq = propParams.betaUpdateFreq();

  double currentCondLogPost = logCondPostThetaBeta(currentParams,model);
  currentCondLogPost = logCondPostThetaBeta(currentParams,model);

  for(unsigned int j=0;j<nFixedEffects;j++){
    for (unsigned int k=0;k<nCategoriesY;k++){
      nTry++;
      propParams.betaAddTry(j);
      double& stdDev = propParams.betaStdDev(j);
      double betaOrig = currentParams.beta(0,j,k,nCategoriesY);
      double betaProp = betaOrig+stdDev*normRand(rndGenerator);
      currentParams.beta(0,j,k,nCategoriesY,betaProp);

      double propCondLogPost = logCondPostThetaBeta(currentParams,model);
      double logAcceptRatio = propCondLogPost - currentCondLogPost;
      if(unifRand(rndGenerator)<exp(logAcceptRatio)){
        nAccept++;
        propParams.betaAddAccept(j);
        currentCondLogPost = propCondLogPost;
        // Update the std dev of the proposal
        if(propParams.nTryBeta(j)%betaUpdateFreq==0){
          stdDev += 10*(propParams.betaLocalAcceptRate(j)-betaTargetRate)/
            pow((double)(propParams.nTryBeta(j)/betaUpdateFreq)+2.0,0.75);
          propParams.betaAnyUpdates(true);
          if(stdDev>propParams.betaStdDevUpper(j)||stdDev<propParams.betaStdDevLower(j)){
            propParams.betaStdDevReset(j);
          }
          propParams.betaLocalReset(j);
        }
      }else{
        currentParams.beta(0,j,k,nCategoriesY,betaOrig);
        // Update the std dev of the proposal
        if(propParams.nTryBeta(j)%betaUpdateFreq==0){
          stdDev += 10*(propParams.betaLocalAcceptRate(j)-betaTargetRate)/
            pow((double)(propParams.nTryBeta(j)/betaUpdateFreq)+2.0,0.75);
          propParams.betaAnyUpdates(true);
          if(stdDev<propParams.betaStdDevLower(j)||stdDev>propParams.betaStdDevUpper(j)){
            propParams.betaStdDevReset(j);
          }
          propParams.betaLocalReset(j);
        }
      }
    }
  }

  // std::cout << " nTry nFixedEffects_mix "<<nTry<<std::endl;
  // for(unsigned int j=0;j<nFixedEffects_mix;j++)
  //   std::cout << j<<" propParams.nTryBetamix "<<propParams.nTryBetamix(j)<<std::endl;
  //
  // if(nFixedEffects_mix>0){
  //   unsigned int maxZ = currentParams.workMaxZi();
  //   double betamixTargetRate = propParams.betamixAcceptTarget();
  //   unsigned int betamixUpdateFreq = propParams.betamixUpdateFreq();
  //
  //   currentCondLogPost = logCondPostThetaBeta(currentParams,model);
  //   double currentCondLogPost = logCondPostThetaBeta(currentParams,model);
  //   for(unsigned int c=0;c<=maxZ;c++){
  //     for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //
  //       //int jjj=0;
  //       //while(jjj<200){
  //         for (unsigned int k=0;k<nCategoriesY;k++){
  //           nTry++;
  //           propParams.betamixAddTry(j);
  //           double& stdDev = propParams.betamixStdDev(j);
  //           double betaOrig = currentParams.beta_mix(c,j,k,nCategoriesY);
  //           double betaProp = betaOrig+stdDev*normRand(rndGenerator);
  //
  //           currentParams.beta_mix(c,j+k,nCategoriesY,betaProp);
  //           double propCondLogPost = logCondPostThetaBeta(currentParams,model);
  //           double logAcceptRatio = propCondLogPost - currentCondLogPost;
  //
  //           double uii=unifRand(rndGenerator);
  //           if(uii<exp(logAcceptRatio)){
  //
  //             nAccept++;
  //             propParams.betamixAddAccept(j);
  //             currentCondLogPost = propCondLogPost;
  //             // Update the std dev of the proposal
  //             if(propParams.nTryBetamix(j)%betaUpdateFreq==0){
  //               stdDev += 10*(propParams.betamixLocalAcceptRate(j)-betamixTargetRate)/
  //                 pow((double)(propParams.nTryBetamix(j)/betamixUpdateFreq)+2.0,0.75);
  //               propParams.betamixAnyUpdates(true);
  //               if(stdDev>propParams.betamixStdDevUpper(j)||stdDev<propParams.betamixStdDevLower(j)){
  //                 propParams.betamixStdDevReset(j);
  //               }
  //               propParams.betamixLocalReset(j);
  //             }
  //             //jjj++;
  //           }else{
  //             currentParams.beta_mix(c,j,k,nCategoriesY,betaOrig);
  //             // Update the std dev of the proposal
  //             if(propParams.nTryBetamix(j)%betaUpdateFreq==0){
  //               stdDev += 10*(propParams.betamixLocalAcceptRate(j)-betamixTargetRate)/
  //                 pow((double)(propParams.nTryBetamix(j)/betamixUpdateFreq)+2.0,0.75);
  //               propParams.betamixAnyUpdates(true);
  //               if(stdDev<propParams.betamixStdDevLower(j)||stdDev>propParams.betamixStdDevUpper(j)){
  //                 propParams.betamixStdDevReset(j);
  //               }
  //               propParams.betamixLocalReset(j);
  //             }
  //           }
  //         }
  //       //}
  //     }
  //   }
  //   std::cout << " nTry nFixedEffects_mix "<<nTry<<std::endl;
  //   for(unsigned int j=0;j<nFixedEffects_mix;j++)
  //     std::cout << j<<" propParams.nTryBetamix "<<propParams.nTryBetamix(j)<<std::endl;
  // }
}


// Adaptive Metropolis-Hastings for lambda
void metropolisHastingsForLambda(mcmcChain<pReMiuMParams>& chain,
                                 unsigned int& nTry,unsigned int& nAccept,
                                 const mcmcModel<pReMiuMParams,
                                                 pReMiuMOptions,
                                                 pReMiuMData>& model,
                                                 pReMiuMPropParams& propParams,
                                                 baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  const string outcomeType = model.dataset().outcomeType();

  // Find the number of subjects
  unsigned int nSubjects = currentParams.nSubjects();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);
  // Define a normal random number generator
  randomNormal normRand(0,1);


  double lambdaTargetRate = propParams.lambdaAcceptTarget();
  unsigned int lambdaUpdateFreq = propParams.lambdaUpdateFreq();

  double (*logCondPostLambdai)(const pReMiuMParams&,
          const mcmcModel<pReMiuMParams,
                          pReMiuMOptions,
                          pReMiuMData>&,
                          const unsigned int&) = NULL;

  if(outcomeType.compare("Bernoulli")==0){
    logCondPostLambdai = &logCondPostLambdaiBernoulli;
  }else if(outcomeType.compare("Binomial")==0){
    logCondPostLambdai = &logCondPostLambdaiBinomial;
  }else if(outcomeType.compare("Poisson")==0){
    logCondPostLambdai = &logCondPostLambdaiPoisson;
  }

  for(unsigned int i=0;i<nSubjects;i++){
    // Only update each lambda i with probability 0.1
    if(unifRand(rndGenerator)>0.1){
      continue;
    }

    nTry++;
    propParams.lambdaAddTry();

    double currentCondLogPost = logCondPostLambdai(currentParams,model,i);
    double& stdDev = propParams.lambdaStdDev();
    double lambdaOrig = currentParams.lambda(i);
    double lambdaProp = lambdaOrig+stdDev*normRand(rndGenerator);
    currentParams.lambda(i,lambdaProp);
    double propCondLogPost = logCondPostLambdai(currentParams,model,i);
    double logAcceptRatio = propCondLogPost - currentCondLogPost;
    if(unifRand(rndGenerator)<exp(logAcceptRatio)){
      nAccept++;
      propParams.lambdaAddAccept();
      // Update the std dev of the proposal
      if(propParams.nTryLambda()%lambdaUpdateFreq==0){
        stdDev += 10*(propParams.lambdaLocalAcceptRate()-lambdaTargetRate)/
          pow((double)(propParams.nTryLambda()/lambdaUpdateFreq)+2.0,0.75);
        propParams.lambdaAnyUpdates(true);
        if(stdDev>propParams.lambdaStdDevUpper()||stdDev<propParams.lambdaStdDevLower()){
          propParams.lambdaStdDevReset();
        }
        propParams.lambdaLocalReset();
      }
    }else{
      currentParams.lambda(i,lambdaOrig);
      // Update the std dev of the proposal
      if(propParams.nTryLambda()%lambdaUpdateFreq==0){
        stdDev += 10*(propParams.lambdaLocalAcceptRate()-lambdaTargetRate)/
          pow((double)(propParams.nTryLambda()/lambdaUpdateFreq)+2.0,0.75);
        propParams.lambdaAnyUpdates(true);
        if(stdDev<propParams.lambdaStdDevLower()||stdDev>propParams.lambdaStdDevUpper()){
          propParams.lambdaStdDevReset();
        }
        propParams.lambdaLocalReset();
      }
    }
  }

}

// Gibbs update for the precision of extra variation epsilon
void gibbsForTauEpsilon(mcmcChain<pReMiuMParams>& chain,
                        unsigned int& nTry,unsigned int& nAccept,
                        const mcmcModel<pReMiuMParams,
                                        pReMiuMOptions,
                                        pReMiuMData>& model,
                                        pReMiuMPropParams& propParams,
                                        baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();
  const string& outcomeType = model.dataset().outcomeType();

  unsigned int nSubjects=dataset.nSubjects();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  // unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix(0);
  unsigned int nCategoriesY = dataset.nCategoriesY();
  Rprintf("Verify that nCategoriesY=1 in this case");

  nTry++;
  nAccept++;

  double a=hyperParams.shapeTauEpsilon(),b=hyperParams.rateTauEpsilon();

  double sumEpsilon = 0.0;
  vector<double> meanVec(nSubjects,0.0);
  if(outcomeType.compare("Poisson")==0){
    meanVec=dataset.logOffset();
  }
  for(unsigned int i=0;i<nSubjects;i++){
    int zi=currentParams.z(i);
    double meanVal=meanVec[i]+currentParams.theta(zi,0);
    for(unsigned int j=0;j<nFixedEffects;j++){
      meanVal+=currentParams.beta(0,j,0,nCategoriesY)*dataset.W(i,j);
    }
    // for(unsigned int j=0;j<nFixedEffects_mix;j++){
    //   meanVal+=currentParams.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
    // }
    double eps = currentParams.lambda(i)-meanVal;
    sumEpsilon+=pow(eps,2.0);
  }
  a+=(double)nSubjects/2.0;
  b+=sumEpsilon/2.0;

  // Boost uses shape and scale parameterisation
  randomGamma gammaRand(a,1.0/b);
  double tau = gammaRand(rndGenerator);
  currentParams.tauEpsilon(tau);

}

// Metropolis-Hastings for joint uptdate of rho and omega
void metropolisHastingsForRhoOmega(mcmcChain<pReMiuMParams>& chain,
                                   unsigned int& nTry,unsigned int& nAccept,
                                   const mcmcModel<pReMiuMParams,
                                                   pReMiuMOptions,
                                                   pReMiuMData>& model,
                                                   pReMiuMPropParams& propParams,
                                                   baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  // Find the number of subjects
  unsigned int nCovariates = currentParams.nCovariates();
  string covariateType = model.options().covariateType();
  string varSelectType = model.options().varSelectType();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);
  // Define a normal random number generator
  randomNormal normRand(0,1);

  double currentLogPost = 0;
  if(varSelectType.compare("Continuous")==0){
    currentLogPost = logCondPostRhoOmegaj(currentParams,model,0);
  }


  double proposedLogPost = currentLogPost;
  vector<unsigned int> currentOmega = currentParams.omega();
  unsigned int proposedOmega;
  vector<double> currentRho = currentParams.rho();
  double proposedRho;

  for(unsigned int j=0;j<nCovariates;j++){

    if(varSelectType.compare("Continuous")!=0){
      currentLogPost = logCondPostRhoOmegaj(currentParams,model,j);
    }

    nTry++;

    // Propose from the priors
    double& stdDev = propParams.rhoStdDev(j);

    if(unifRand(rndGenerator)>hyperParams.atomRho()){
      // Proposing an omega 0
      if(currentOmega[j]==0){
        // Nothing to do as move to the same place
        nAccept++;
        continue;
      }
      proposedOmega=0;
      proposedRho=0.0;

      currentParams.omega(j,proposedOmega);
      currentParams.rho(j,proposedRho,covariateType,varSelectType);
      proposedLogPost = logCondPostRhoOmegaj(currentParams,model,j);
      double logAcceptRatio = proposedLogPost - currentLogPost;
      double runiftemp = unifRand(rndGenerator);
      logAcceptRatio += logPdfBeta(currentRho[j],hyperParams.aRho(),hyperParams.bRho());

      if(runiftemp<exp(logAcceptRatio)){
        // Move accepted
        if(varSelectType.compare("Continuous")==0){
          currentLogPost=proposedLogPost;
        }
        nAccept++;

      }else{
        // Move rejected, reset parameters
        currentParams.omega(j,currentOmega[j]);
        currentParams.rho(j,currentRho[j],covariateType,varSelectType);

      }
    }else{
      if(currentOmega[j]==1){
        proposedRho  = truncNormalRand(rndGenerator,currentRho[j],stdDev,"B",0,1);
        currentParams.rho(j,proposedRho,covariateType,varSelectType);
        proposedLogPost = logCondPostRhoOmegaj(currentParams,model,j);
        double logAcceptRatio = proposedLogPost - currentLogPost;
        logAcceptRatio += logPdfTruncatedNormal(currentRho[j],proposedRho,stdDev,"B",0,1);
        logAcceptRatio -= logPdfTruncatedNormal(proposedRho,currentRho[j],stdDev,"B",0,1);
        propParams.rhoAddTry(j);

        double runiftemp = unifRand(rndGenerator);
        if(runiftemp<exp(logAcceptRatio)){
          // Move accepted
          if(varSelectType.compare("Continuous")==0){
            currentLogPost=proposedLogPost;
          }

          nAccept++;
          propParams.rhoAddAccept(j);
          // Also update the proposal standard deviation
          if(propParams.nTryRho(j)%propParams.rhoUpdateFreq()==0){
            stdDev += 0.1*(propParams.rhoLocalAcceptRate(j)-propParams.rhoAcceptTarget())/
              pow((double)(propParams.nTryRho(j)/propParams.rhoUpdateFreq())+2.0,0.75);
            propParams.rhoAnyUpdates(true);
            if(stdDev>propParams.rhoStdDevUpper(j)||stdDev<propParams.rhoStdDevLower(j)){
              propParams.rhoStdDevReset(j);
            }
            propParams.rhoLocalReset(j);
          }

        }else{
          // Move rejected, reset parameters
          currentParams.omega(j,currentOmega[j]);
          currentParams.rho(j,currentRho[j],covariateType,varSelectType);
          // Also update the proposal standard deviation
          if(propParams.nTryRho(j)%propParams.rhoUpdateFreq()==0){
            stdDev += 0.1*(propParams.rhoLocalAcceptRate(j)-propParams.rhoAcceptTarget())/
              pow((double)(propParams.nTryRho(j)/propParams.rhoUpdateFreq())+2.0,0.75);
            propParams.rhoAnyUpdates(true);
            if(stdDev>propParams.rhoStdDevUpper(j)||stdDev<propParams.rhoStdDevLower(j)){
              propParams.rhoStdDevReset(j);
            }
            propParams.rhoLocalReset(j);
          }
        }
      }else{
        proposedRho = betaRand(rndGenerator,hyperParams.aRho(),hyperParams.bRho());
        proposedOmega=1;
        currentParams.omega(j,proposedOmega);
        currentParams.rho(j,proposedRho,covariateType,varSelectType);
        proposedLogPost = logCondPostRhoOmegaj(currentParams,model,j);
        double logAcceptRatio = proposedLogPost - currentLogPost;
        logAcceptRatio -= logPdfBeta(proposedRho,hyperParams.aRho(),hyperParams.bRho());

        double runiftemp = unifRand(rndGenerator);
        if(runiftemp<exp(logAcceptRatio)){
          // Move accepted
          if(varSelectType.compare("Continuous")==0){
            currentLogPost=proposedLogPost;
          }
          nAccept++;

        }else{
          // Move rejected, reset parameters
          currentParams.omega(j,currentOmega[j]);
          currentParams.rho(j,currentRho[j],covariateType,varSelectType);
        }
      }
    }
  }
}

// Gibbs for update of sigmaSqY (Normal response case)
void gibbsForSigmaSqY(mcmcChain<pReMiuMParams>& chain,
                      unsigned int& nTry,unsigned int& nAccept,
                      const mcmcModel<pReMiuMParams,
                                      pReMiuMOptions,
                                      pReMiuMData>& model,
                                      pReMiuMPropParams& propParams,
                                      baseGeneratorType& rndGenerator){
  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();

  const pReMiuMData& dataset = model.dataset();

  unsigned int nSubjects = currentParams.nSubjects();
  unsigned int nFixedEffects = dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = dataset.nCategoriesY();
  Rprintf("Verify that nCategoriesY=1 in this case");

  nTry++;
  nAccept++;

  double sumVal=0.0;
  for(unsigned int i=0;i<nSubjects;i++){
    int Zi = currentParams.z(i);

    double mu = currentParams.theta(Zi,0);
    for(unsigned int j=0;j<nFixedEffects;j++){
      mu+=currentParams.beta(0,j,0,nCategoriesY)*dataset.W(i,j);
    }
    // for(unsigned int j=0;j<nFixedEffects_mix;j++){
    //   mu+=currentParams.beta_mix(Zi,j, 0, nCategoriesY)*dataset.W_mix(i,j);
    // }
    sumVal+=pow(dataset.continuousY(i)-mu,2.0);
  }

  double posteriorShape = hyperParams.shapeSigmaSqY()+(double)nSubjects/2.0;
  double posteriorScale = hyperParams.scaleSigmaSqY()+0.5*sumVal;

  randomGamma gammaRand(posteriorShape,1.0/posteriorScale);
  double sigmaSqY=1.0/gammaRand(rndGenerator);
  currentParams.sigmaSqY(sigmaSqY);


}

// Gibbs for update of nu (survival response case) using adaptive rejection sampling
void gibbsForNu(mcmcChain<pReMiuMParams>& chain,
                unsigned int& nTry,unsigned int& nAccept,
                const mcmcModel<pReMiuMParams,
                                pReMiuMOptions,
                                pReMiuMData>& model,
                                pReMiuMPropParams& propParams,
                                baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const bool weibullFixedShape=model.options().weibullFixedShape();
  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();

  nTry++;
  nAccept++;
  if (weibullFixedShape){
    double nu = ARSsampleNu2(currentParams, model, 0,logNuPostSurvival,rndGenerator);
    currentParams.nu(0,nu);
  } else {
    for (unsigned int c=0;c<=maxZ;c++){
      double nu = ARSsampleNu2(currentParams, model, c,logNuPostSurvival,rndGenerator);
      currentParams.nu(c,nu);
    }
  }
}

//RJ metropolisHastingsForL
void metropolisHastingsForL(mcmcChain<pReMiuMParams>& chain,
                            unsigned int& nTry,unsigned int& nAccept,
                            const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                            pReMiuMPropParams& propParams,
                            baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  const string outcomeType = model.dataset().outcomeType();

  // Find the number of clusters
  unsigned int maxZ = currentParams.workMaxZi();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);
  // Define a normal random number generator
  randomNormal normRand(0,1);

  double LTargetRate = propParams.LAcceptTarget(); //0.44
  unsigned int LUpdateFreq = propParams.LUpdateFreq(); //25

  double currentCondLogPost = 0.0;
  //int ratio = 1; // if ratio=1, exp(L2)=exp(L1)/4 L2= L1- log(4), otherwise L2 separate
  bool ratio_estim = model.options().estim_ratio();

  for(unsigned int c=0;c<=maxZ;c++){

    unsigned int nL=4;
    if(model.dataset().kernelType().compare("SQexponential")==0)
      nL-=1;

    if(model.options().sampleGPmean()){//AR model.options().sampleGPmean()

      if(ratio_estim){
        for (unsigned int l=0;l<nL;l++){
          if(l==0)
            currentParams.L(c,2,currentParams.L(c,0)+ log(currentParams.ratio(c)));

          if(l < 2 ){ //l<2
            currentCondLogPost = logCondPostL_covGP(currentParams,model,c); // p(L_k)p(f_k|L_k,k)
            //if(l==0 & ratio_estim != 0)
            //currentCondLogPost += logCondPostL(currentParams,model,c,0); // p(L_k)p(f_k|L_k,k) * p(Y_k|f_k, L_k,k)


            nTry++;
            propParams.LAddTry(l);
            double& stdDev = propParams.LStdDev(l);
            double LOrig = currentParams.L(c,l);
            double LProp = LOrig+stdDev*normRand(rndGenerator);
            currentParams.L(c,l,LProp);
            if(l==0)
              currentParams.L(c,2,LProp+ log(currentParams.ratio(c)));

            double propCondLogPost = 0;
            double logAcceptRatio = 0;

            propCondLogPost = logCondPostL_covGP(currentParams,model,c); // p(L_013k)p(f_k|L_013k,k)

            // if(l==0 & ratio_estim != 0)
            //  propCondLogPost += logCondPostL(currentParams,model,c,0); // p(L_013k)p(f_k|L_013k,k)

            logAcceptRatio = propCondLogPost - currentCondLogPost;

            double uni = unifRand(rndGenerator);

            if(uni<exp(logAcceptRatio)){
              nAccept++;
              propParams.LAddAccept(l);
              currentCondLogPost = propCondLogPost;
            }else{
              currentParams.L(c,l,LOrig);
              if(l==0)
                currentParams.L(c,2,LOrig+ log(currentParams.ratio(c)));
            }

            // Update the std dev of the proposal
            if(propParams.nTryL(l)%LUpdateFreq==0){
              //if(propParams.LLocalAcceptRate(l)>LTargetRate)
              //	stdDev *= std::exp(1.0/(propParams.LLocalAcceptRate(l)*LUpdateFreq));
              //else
              //	stdDev /= std::exp(1.0/(LUpdateFreq-propParams.LLocalAcceptRate(l)*LUpdateFreq));

              stdDev += 10*(propParams.LLocalAcceptRate(l)-LTargetRate)/
                pow((double)(propParams.nTryL(l)/LUpdateFreq)+2.0,0.75);

              propParams.LAnyUpdates(true);
              if(stdDev>propParams.LStdDevUpper(l)||stdDev<propParams.LStdDevLower(l)){
                propParams.LStdDevReset(l);
              }
              propParams.LLocalReset(l);
            }
          }


          if((l==2) & (1>2)){ //estimation ratio with beta distribution
            currentParams.L(c,2,currentParams.L(c,0)+ log(currentParams.ratio(c)));
            currentCondLogPost = logCondPostL(currentParams,model,c,0); // p(Y^k|L_k,f_k,k) 0 to remove L prior

            currentCondLogPost += currentParams.ratio(c)*(1-currentParams.ratio(c));
            //currentCondLogPost += logPdfBeta(currentParams.ratio(c),currentParams.hyperParams().aRatio(),currentParams.hyperParams().bRatio()); // p(ratio) prior

            nTry++;
            // propParams.LAddTry(l);
            // double& stdDev = propParams.LStdDev(l);
            // double LOrig = currentParams.L(c,l);
            // double LProp = LOrig+stdDev*normRand(rndGenerator);

            propParams.LAddTry(l);
            //double& stdDev = propParams.LStdDev(l);
            double stdDev=0.5;
            double logitRatioOrig = log(currentParams.ratio(c)/(1-currentParams.ratio(c)));
            double logitRatioProp = logitRatioOrig+stdDev*normRand(rndGenerator);
            double RatioOrig = currentParams.ratio(c);
            double RatioProp = exp(logitRatioProp)/(1+exp(logitRatioProp));

            double LOrig = currentParams.L(c,2);
            double LProp = currentParams.L(c,0)+log(RatioProp);

            currentParams.ratio(c,RatioProp);
            currentParams.L(c,2,LProp);

            double propCondLogPost = 0;
            double logAcceptRatio = 0;



            propCondLogPost = logCondPostL(currentParams,model,c,0); // p(Y^k|L_2k,f_k,k)
            propCondLogPost += currentParams.ratio(c)*(1-currentParams.ratio(c));

            logAcceptRatio = propCondLogPost - currentCondLogPost;

            double uni = unifRand(rndGenerator);

            if(uni<exp(logAcceptRatio)){
              nAccept++;
              propParams.LAddAccept(l);
              currentCondLogPost = propCondLogPost;
            }else{
              currentParams.L(c,l,LOrig);
              currentParams.ratio(c,RatioOrig);
            }

            // Update the std dev of the proposal
            if(propParams.nTryL(l)%LUpdateFreq==0){
              //if(propParams.LLocalAcceptRate(l)>LTargetRate)
              //	stdDev *= std::exp(1.0/(propParams.LLocalAcceptRate(l)*LUpdateFreq));
              //else
              //	stdDev /= std::exp(1.0/(LUpdateFreq-propParams.LLocalAcceptRate(l)*LUpdateFreq));

              stdDev += 10*(propParams.LLocalAcceptRate(l)-LTargetRate)/
                pow((double)(propParams.nTryL(l)/LUpdateFreq)+2.0,0.75);

              propParams.LAnyUpdates(true);
              if(stdDev>propParams.LStdDevUpper(l)||stdDev<propParams.LStdDevLower(l)){
                propParams.LStdDevReset(l);
              }
              propParams.LLocalReset(l);
            }
          }
        }

      }else{ //ratio_estim=FALSE and sampleGPmean
        for (unsigned int l=0;l<nL;l++){ //AR change nL-1

          if(l == 2){
            currentCondLogPost = logCondPostL(currentParams,model,c,1); // p(L_k)p(Y^k|L_k,f_k,k), 1 to add L prior
          }else{ //l<2
            currentCondLogPost = logCondPostL_covGP(currentParams,model,c); // p(L_k)p(f_k|L_k,k)
            //if(l==0 & ratio_estim != 0)
            //currentCondLogPost += logCondPostL(currentParams,model,c,0); // p(L_k)p(f_k|L_k,k) * p(Y_k|f_k, L_k,k)
          }

          nTry++;
          propParams.LAddTry(l);
          double& stdDev = propParams.LStdDev(l);
          double LOrig = currentParams.L(c,l);
          double LProp = LOrig+stdDev*normRand(rndGenerator);
          currentParams.L(c,l,LProp);

          double propCondLogPost = 0;
          double logAcceptRatio = 0;

          if(l == 2){
            propCondLogPost = logCondPostL(currentParams,model,c,1); // p(L_2)p(Y^k|L_2k,f_k,k)
          }else{
            propCondLogPost = logCondPostL_covGP(currentParams,model,c); // p(L_013k)p(f_k|L_013k,k)

            // if(l==0 & ratio_estim != 0)
            //  propCondLogPost += logCondPostL(currentParams,model,c,0); // p(L_013k)p(f_k|L_013k,k)
          }
          logAcceptRatio = propCondLogPost - currentCondLogPost;

          double uni = unifRand(rndGenerator);

          if(uni<exp(logAcceptRatio)){
            nAccept++;
            propParams.LAddAccept(l);
            currentCondLogPost = propCondLogPost;
          }else{
            currentParams.L(c,l,LOrig);
          }

          // Update the std dev of the proposal
          if(propParams.nTryL(l)%LUpdateFreq==0){
            //if(propParams.LLocalAcceptRate(l)>LTargetRate)
            //	stdDev *= std::exp(1.0/(propParams.LLocalAcceptRate(l)*LUpdateFreq));
            //else
            //	stdDev /= std::exp(1.0/(LUpdateFreq-propParams.LLocalAcceptRate(l)*LUpdateFreq));

            stdDev += 10*(propParams.LLocalAcceptRate(l)-LTargetRate)/
              pow((double)(propParams.nTryL(l)/LUpdateFreq)+2.0,0.75);

            propParams.LAnyUpdates(true);
            if(stdDev>propParams.LStdDevUpper(l)||stdDev<propParams.LStdDevLower(l)){
              propParams.LStdDevReset(l);
            }
            propParams.LLocalReset(l);
          }
        }
      }



      //AR sample meanGP
      //if(model.options().sampleGPmean()){ //AR change
      unsigned int nTimes_unique = model.dataset().nTimes_unique();
      //int threads = omp_get_max_threads();
      //vector<baseGeneratorType> rngArray(threads);

      // Seed by taking random numbers from the existing generator
      //boost::random::uniform_int_distribution<> seeder;
      //for (auto &rng : rngArray)
      //rng.seed(seeder(rndGenerator));

      //#pragma omp parallel for
      //for(unsigned int c=0;c<=maxZ;c++){
      //  int threadNum = omp_get_thread_num();


      VectorXd Fval = VectorXd::Zero(nTimes_unique);
      //Fval =  Sample_GPmean(currentParams, model.dataset(), c, rngArray[threadNum],0);
      Fval =  Sample_GPmean(currentParams, model.dataset(), c, rndGenerator,0);

      for(unsigned int j=0;j<nTimes_unique;j++){
        currentParams.meanGP(c,j, Fval(j));
      }
    }else{
      currentCondLogPost = logCondPostL(currentParams,model,c,1);

      for (unsigned int l=0;l<nL;l++){
        nTry++;
        propParams.LAddTry(l);
        double& stdDev = propParams.LStdDev(l);
        double LOrig = currentParams.L(c,l);
        double ui1= normRand(rndGenerator);
        double LProp = LOrig+stdDev*ui1;
        currentParams.L(c,l,LProp);
        double propCondLogPost = logCondPostL(currentParams,model,c,1);
        double logAcceptRatio = propCondLogPost - currentCondLogPost;
        double uii= unifRand(rndGenerator);

        if(uii<exp(logAcceptRatio)){
          nAccept++;
          propParams.LAddAccept(l);
          currentCondLogPost = propCondLogPost;
        }else{
          currentParams.L(c,l,LOrig);
        }

        // Update the std dev of the proposal
        if(propParams.nTryL(l)%LUpdateFreq==0){
          //if(propParams.LLocalAcceptRate(l)>LTargetRate)
          //	stdDev *= std::exp(1.0/(propParams.LLocalAcceptRate(l)*LUpdateFreq));
          //else
          //	stdDev /= std::exp(1.0/(LUpdateFreq-propParams.LLocalAcceptRate(l)*LUpdateFreq));
          stdDev += 10*(propParams.LLocalAcceptRate(l)-LTargetRate)/
            pow((double)(propParams.nTryL(l)/LUpdateFreq)+2.0,0.75);
          propParams.LAnyUpdates(true);
          if(stdDev>propParams.LStdDevUpper(l)||stdDev<propParams.LStdDevLower(l)){
            propParams.LStdDevReset(l);
          }
          propParams.LLocalReset(l);
        }
      }
    }
  }
}

// void gibbsForSigmaEpsilonLME_invchi2(mcmcChain<pReMiuMParams>& chain,
//                                      unsigned int& nTry,unsigned int& nAccept,
//                                      const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
//                                      pReMiuMPropParams& propParams,
//                                      baseGeneratorType& rndGenerator){
//
//
//
//
//   mcmcState<pReMiuMParams>& currentState = chain.currentState();
//   pReMiuMParams& currentParams = currentState.parameters();
//   const pReMiuMData& dataset=model.dataset();
//   const string outcomeType = dataset.outcomeType();
//
//   // Find the number of clusters
//   unsigned int nFixedEffects = dataset.nFixedEffects();
//   vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
//   vector<unsigned int> nCategoriesY = dataset.nCategoriesY();
//   unsigned int nOutcomes = dataset.nOutcomes();
//   // Define a normal random number generator
//   randomNormal normRand(0,1);
//   double S2= 0.0;
//
//   unsigned int  nSubjects = dataset.nSubjects();
//   vector<double>        y = dataset.continuousY();
//   const string kernelType = dataset.kernelType(); //AR
//   vector<double>    times = dataset.times();
//   vector<int>      tStart = dataset.tStart();
//   vector<int>       tStop = dataset.tStop();
//   int                ntot = 0;
//
//   for(unsigned int m=0;m<nOutcomes;m++){
//     for(unsigned int i=0;i<nSubjects;i++){
//
//       int nmes = (tStop[i] - tStart[i] + 1);
//       VectorXd   Yi(nmes);
//       int zi   = currentParams.z(i);
//       MatrixXd   Xi(nmes, nFixedEffects[m]);
//       ntot += nmes;
//       for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
//         Yi(j) = y[tStart[i]-1+j];
//
//         for(unsigned int b=0;b<nFixedEffects[m];b++){
//           Yi(j) -= currentParams.beta(m,b,0,nCategoriesY)*dataset.W_LME(m,tStart[i]-1+j,b);
//         }
//
//         for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
//           Yi(j) -= currentParams.beta_mix(m,zi,b, 0, nCategoriesY)*dataset.W_LME_mix(m,tStart[i]-1+j,b);
//         }
//       }
//
//       MatrixXd block=dataset.W_RE(m,tStart[i]-1, 0, tStop[i]-tStart[i]+1 , dataset.nRandomEffects(m));
//       Yi -= block*currentParams.RandomEffects(m,i);
//       S2 += Yi.transpose()*Yi;
//     }
//
//     S2 +=  currentParams.hyperParams().eps_sigma2_0()*currentParams.hyperParams().eps_vu();
//
//     int nu_post = currentParams.hyperParams().eps_vu() + ntot;
//
//     // Define a ChiSquare random number generator
//     randomChiSquare chisQRand(nu_post);
//     double epsilon;
//     double temp = chisQRand(rndGenerator);
//     epsilon = S2/temp;
//     currentParams.SigmaE(m,epsilon);
//   }
//   //currentParams.SigmaE(0.01);
// }


void gibbsForSigmaEpsilonLME(mcmcChain<pReMiuMParams>& chain,
                             unsigned int& nTry,unsigned int& nAccept,
                             const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                             pReMiuMPropParams& propParams,
                             baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  const pReMiuMData& dataset=model.dataset();
  const string outcomeType = dataset.outcomeType();

  // Find the number of clusters
  vector<unsigned int> nFixedEffects = dataset.nFixedEffects();
  vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = dataset.nCategoriesY();
  unsigned int nOutcomes = dataset.nOutcomes();

  double S2= 0.0;

  unsigned int  nSubjects = dataset.nSubjects();
  vector<double>        y = dataset.continuousY();
  const string kernelType = dataset.kernelType(); //AR
  vector<double>    times = dataset.times();
  vector<int>      tStart = dataset.tStart();
  vector<int>       tStop = dataset.tStop();

  unsigned int ind_y=0;
  unsigned int ind=0;

  for(unsigned int m=0;m<nOutcomes;m++){
    int  ntot = 0;
    for(unsigned int i=0;i<nSubjects;i++){

      int nmes = (tStop[ind] - tStart[ind] + 1);
      VectorXd   Yi(nmes);
      int zi   = currentParams.z(i);
      MatrixXd   Xi(nmes, nFixedEffects[m]);
      ntot += nmes;
      for(unsigned int j=0;j<tStop[ind]-tStart[ind]+1;j++){
        Yi(j) = y[ind_y + tStart[ind]-1+j];

        for(unsigned int b=0;b<nFixedEffects[m];b++){
          Yi(j) -= currentParams.beta(m,b,0,nCategoriesY)*dataset.W_LME(m,tStart[ind]-1+j,b);
        }

        for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
          Yi(j) -= currentParams.beta_mix(m,zi,b, 0, nCategoriesY)*dataset.W_LME_mix(m,tStart[ind]-1+j,b);
        }
      }

      MatrixXd block=dataset.W_RE(m,tStart[ind]-1, 0, tStop[ind]-tStart[ind]+1 , dataset.nRandomEffects(m));
      Yi -= block*currentParams.RandomEffects(m,i);
      S2 += Yi.transpose()*Yi;
      if(i==(nSubjects=1))
        ind_y += nSubjects;
      ind++;
    }

    S2 /=2;
    double shape_post = currentParams.hyperParams().eps_shape() + ntot/2;
    //double scale_post = 2*currentParams.hyperParams().eps_rate()/(2 + S2*currentParams.hyperParams().eps_rate());
    double scale_rate = currentParams.hyperParams().eps_rate() + S2;

    // Define a inverse gamma random number generator
    randomGamma gammaRand(shape_post, 1/scale_rate);

    double temp = gammaRand(rndGenerator);
    double epsilon = 1/temp;

    currentParams.SigmaE(m,epsilon);// variance
  }
}


//AR Sigma epsilon for LME
// void metropolisHastingsForSigmaEpsilonLME(mcmcChain<pReMiuMParams>& chain,
//                                           unsigned int& nTry,unsigned int& nAccept,
//                                           const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
//                                           pReMiuMPropParams& propParams,
//                                           baseGeneratorType& rndGenerator){
//   mcmcState<pReMiuMParams>& currentState = chain.currentState();
//   pReMiuMParams& currentParams = currentState.parameters();
//   const string outcomeType = model.dataset().outcomeType();
//
//   // Define a uniform random number generator
//   randomUniform unifRand(0,1);
//   // Define a normal random number generator
//   randomNormal normRand(0,1);
//
//   double SigmaETargetRate = propParams.SigmaEAcceptTarget(); //0.44
//   unsigned int SigmaEUpdateFreq = propParams.SigmaEUpdateFreq(); //25
//
//   double currentCondLogPost = 0.0;
//   currentCondLogPost = logCondPostSigmaE(currentParams,model);
//
//   nTry++;
//   propParams.SigmaEAddTry();
//   double& stdDev = propParams.SigmaEStdDev();
//   double SigmaEOrig = currentParams.SigmaE();
//   double ui1= normRand(rndGenerator);
//   double SigmaEProp = SigmaEOrig+stdDev*ui1;
//   currentParams.SigmaE(SigmaEProp);
//   double propCondLogPost = logCondPostSigmaE(currentParams,model);
//   double logAcceptRatio = propCondLogPost - currentCondLogPost;
//   double uii= unifRand(rndGenerator);
//   if(uii<exp(logAcceptRatio)){
//     nAccept++;
//     propParams.SigmaEAddAccept();
//     currentCondLogPost = propCondLogPost;
//   }else{
//     currentParams.SigmaE(SigmaEOrig);
//   }
//
//   // Update the std dev of the proposal
//   if(propParams.nTrySigmaE()%SigmaEUpdateFreq==0){
//     //if(propParams.LLocalAcceptRate(l)>LTargetRate)
//     //	stdDev *= std::exp(1.0/(propParams.LLocalAcceptRate(l)*LUpdateFreq));
//     //else
//     //	stdDev /= std::exp(1.0/(LUpdateFreq-propParams.LLocalAcceptRate(l)*LUpdateFreq));
//     stdDev += 10*(propParams.SigmaELocalAcceptRate()-SigmaETargetRate)/
//       pow((double)(propParams.nTrySigmaE()/SigmaEUpdateFreq)+2.0,0.75);
//     propParams.SigmaEAnyUpdates(true);
//     if(stdDev>propParams.SigmaEStdDevUpper()||stdDev<propParams.SigmaEStdDevLower()){
//       propParams.SigmaEStdDevReset();
//     }
//     propParams.SigmaELocalReset();
//   }
// }

// Gibbs update for the precision of spatial random term
void gibbsForTauCAR(mcmcChain<pReMiuMParams>& chain,
                    unsigned int& nTry,unsigned int& nAccept,
                    const mcmcModel<pReMiuMParams,
                                    pReMiuMOptions,
                                    pReMiuMData>& model,
                                    pReMiuMPropParams& propParams,
                                    baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();

  //Rprintf("TauCAR before update is %f \n.", currentParams.TauCAR());

  unsigned int nSubjects=dataset.nSubjects();

  nTry++;
  nAccept++;

  double a=hyperParams.shapeTauCAR(),b=hyperParams.rateTauCAR();

  double sumCAR1 = 0.0;
  double sumCAR2 = 0.0;
  for (unsigned int i=0; i<nSubjects; i++){
    double uCARi = currentParams.uCAR(i);
    int nNeighi = dataset.nNeighbours(i);
    sumCAR1+= uCARi*uCARi*nNeighi;
    for (int j = 0; j<nNeighi; j++){
      unsigned int nj = dataset.neighbours(i,j);
      double ucarj = currentParams.uCAR(nj-1);
      sumCAR2+=uCARi*ucarj;
    }
  }
  double sumCAR=sumCAR1-sumCAR2;

  a+=(double)(nSubjects-1)/2.0;
  b+=sumCAR/2.0;


  // Boost uses shape and scale parameterisation
  randomGamma gammaRand(a,1.0/b);
  double tau = gammaRand(rndGenerator);
  currentParams.TauCAR(tau);
  //Rprintf("TauCAR after update is %f \n .", currentParams.TauCAR());
}

// Gibbs update for spatial random term using adaptive rejection sampling
void gibbsForUCAR(mcmcChain<pReMiuMParams>& chain,
                  unsigned int& nTry,unsigned int& nAccept,
                  const mcmcModel<pReMiuMParams,
                                  pReMiuMOptions,
                                  pReMiuMData>& model,
                                  pReMiuMPropParams& propParams,
                                  baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  const pReMiuMData& dataset = model.dataset();
  unsigned int nSubjects=dataset.nSubjects();
  const string& outcomeType = model.dataset().outcomeType();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = dataset.nCategoriesY();
  Rprintf("Verify that nCategoriesY=1 in this case");
  //Rprintf("TauCAR after update of uCAR is %f \n .", currentParams.TauCAR());
  nTry++;
  nAccept++;

  vector<double> tempU;
  tempU.resize(nSubjects);
  if(outcomeType.compare("Poisson")==0){
    for (unsigned int iSub=0; iSub<nSubjects; iSub++){
      double ui=ARSsampleCAR(currentParams, model, iSub,logUiPostPoissonSpatial,rndGenerator);
      tempU[iSub]=ui;
    }
  } else if(outcomeType.compare("Normal")==0){
    for (unsigned int iSub=0; iSub<nSubjects; iSub++){
      int nNeighi = dataset.nNeighbours(iSub);
      double sigmaSqUCAR = 1/(1/currentParams.sigmaSqY()+currentParams.TauCAR()*nNeighi);
      int Zi = currentParams.z(iSub);
      double betaW = 0.0;
      for(unsigned int j=0;j<nFixedEffects;j++){
        betaW+=currentParams.beta(0,j,0,nCategoriesY)*dataset.W(iSub,j);
      }
      // for(unsigned int j=0;j<nFixedEffects_mix;j++){
      //   betaW+=currentParams.beta_mix(Zi,j,0,nCategoriesY)*dataset.W_mix(iSub,j);
      // }
      double meanUi=0.0;
      for (int j = 0; j<nNeighi; j++){
        unsigned int nj = dataset.neighbours(iSub,j);
        double ucarj = currentParams.uCAR(nj-1);
        meanUi+=ucarj;
      }
      meanUi/=nNeighi;
      double mUCAR = 1/currentParams.sigmaSqY()*(dataset.continuousY(iSub)-currentParams.theta(Zi,0)-betaW)+currentParams.TauCAR()*nNeighi*meanUi;
      mUCAR = mUCAR * sigmaSqUCAR;
      randomNormal normRand(0,1);
      tempU[iSub]=sqrt(sigmaSqUCAR)*normRand(rndGenerator)+mUCAR;
    }
  }
  double meanU=0.0;
  for (unsigned int i=0; i<nSubjects; i++){meanU+=tempU[i];}
  meanU/=nSubjects;
  for (unsigned int i=0; i<nSubjects; i++){tempU[i]-=meanU;}
  currentParams.uCAR(tempU);
  //Rprintf("uCAR1 equals %f \n", currentParams.uCAR(1));
}


/*********** BLOCK 5 p(Z|.) **********************************/

// Gibbs update for the allocation variables
void gibbsForZ(mcmcChain<pReMiuMParams>& chain,
               unsigned int& nTry,unsigned int& nAccept,
               const mcmcModel<pReMiuMParams,
                               pReMiuMOptions,
                               pReMiuMData>& model,
                               pReMiuMPropParams& propParams,
                               baseGeneratorType& rndGenerator){

  mcmcState<pReMiuMParams>& currentState = chain.currentState();
  pReMiuMParams& currentParams = currentState.parameters();
  pReMiuMHyperParams hyperParams = currentParams.hyperParams();
  const pReMiuMData& dataset = model.dataset();
  const string& outcomeType = model.dataset().outcomeType();
  const string& kernelType = model.dataset().kernelType(); //AR
  const string& covariateType = model.dataset().covariateType();
  const string& samplerType = model.options().samplerType();
  bool computeEntropy = model.options().computeEntropy();
  unsigned int nSubjects=dataset.nSubjects();
  unsigned int nPredictSubjects=dataset.nPredictSubjects();
  unsigned int maxNClusters=currentParams.maxNClusters();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY=dataset.nCategoriesY();
  unsigned int nCovariates=dataset.nCovariates();
  unsigned int nDiscreteCovs=dataset.nDiscreteCovs();
  unsigned int nContinuousCovs=dataset.nContinuousCovs();
  vector<unsigned int>nCategories=dataset.nCategories();
  const vector<vector<bool> >& missingX=dataset.missingX();
  bool includeResponse = model.options().includeResponse();
  bool responseExtraVar = model.options().responseExtraVar();
  const bool includeCAR=model.options().includeCAR();
  const string& predictType = model.options().predictType();


  nTry++;
  nAccept++;

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  vector<unsigned int> nMembers(maxNClusters,0);

  vector<double> rnd(nSubjects+nPredictSubjects,0.0);
  vector<double> u(nSubjects+nPredictSubjects,0.0);

  for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
    rnd[i] = unifRand(rndGenerator);
    u[i] = currentParams.u(i);
  }

  vector<double> testBound(maxNClusters,0.0);
  vector<double> clusterWeight(maxNClusters,0.0);

  //#pragma omp parallel for
  for(unsigned int c=0;c<maxNClusters;c++){
    if(samplerType.compare("SliceDependent")==0){
      testBound[c] = exp(currentParams.logPsi(c));
      clusterWeight[c] = 0.0;
    }else if(samplerType.compare("SliceIndependent")==0){
      testBound[c] = hyperParams.workXiSlice(c);
      clusterWeight[c] = currentParams.logPsi(c)-(double)c*log(hyperParams.rSlice())-log(1-hyperParams.rSlice());
    }else if(samplerType.compare("Truncated")==0){
      testBound[c] = 1.0;
      clusterWeight[c] = currentParams.logPsi(c);
    }
  }


  // Compute the allocation probabilities in terms of the unique vectors
  vector<vector<double> > logPXiGivenZi;
  logPXiGivenZi.resize(nSubjects+nPredictSubjects);
  if(covariateType.compare("Discrete")==0){

    //#pragma omp parallel for
    for(unsigned int i=0;i<nSubjects;i++){
      logPXiGivenZi[i].resize(maxNClusters,0);
      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          if(currentParams.z(i)==(int)c){
            logPXiGivenZi[i][c]=currentParams.workLogPXiGivenZi(i);
          }else{
            logPXiGivenZi[i][c]=0;
            for(unsigned int j=0;j<nCovariates;j++){
              int Xij = currentParams.workDiscreteX(i,j);
              logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
            }
          }
        }
      }
    }

    // For the predictive subjects we do not count missing data
    for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
      logPXiGivenZi[i].resize(maxNClusters,0);
      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          for(unsigned int j=0;j<nCovariates;j++){
            if(!missingX[i][j]){
              int Xij = currentParams.workDiscreteX(i,j);
              logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
            }
          }
        }
      }
    }

  }else if(covariateType.compare("Normal")==0){
    for(unsigned int i=0;i<nSubjects;i++){
      logPXiGivenZi[i].resize(maxNClusters,0.0);
      VectorXd xi=VectorXd::Zero(nCovariates);
      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          if(currentParams.z(i)==(int)c){
            logPXiGivenZi[i][c]=currentParams.workLogPXiGivenZi(i);
          }else{
            for(unsigned int j=0;j<nCovariates;j++){
              xi(j)=currentParams.workContinuousX(i,j);
            }
            logPXiGivenZi[i][c]=logPdfMultivarNormal(nCovariates,xi,currentParams.workMuStar(c),currentParams.workSqrtTau(c),currentParams.workLogDetTau(c));
          }
        }

      }
    }
    // For the predictive subjects we do not count missing data
    for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
      logPXiGivenZi[i].resize(maxNClusters,0.0);

      unsigned int nNotMissing=dataset.nContinuousCovariatesNotMissing(i);

      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          VectorXd workMuStar=currentParams.workMuStar(c);

          VectorXd xi=VectorXd::Zero(nNotMissing);
          VectorXd muStar=VectorXd::Zero(nNotMissing);
          MatrixXd sqrtTau=MatrixXd::Zero(nNotMissing,nNotMissing);
          double logDetTau =0.0;
          if(nNotMissing==nCovariates){
            muStar=workMuStar;
            sqrtTau=currentParams.workSqrtTau(c);
            logDetTau=currentParams.workLogDetTau(c);
            for(unsigned int j=0;j<nCovariates;j++){
              xi(j)=currentParams.workContinuousX(i,j);
            }
          }else{
            MatrixXd workSigma=currentParams.Sigma(c);
            MatrixXd Sigma=MatrixXd::Zero(nNotMissing,nNotMissing);
            MatrixXd Tau=MatrixXd::Zero(nNotMissing,nNotMissing);
            unsigned int j=0;
            for(unsigned int j0=0;j0<nCovariates;j0++){
              if(!missingX[i][j0]){
                xi(j)=currentParams.workContinuousX(i,j0);
                muStar(j)=workMuStar(j0);
                unsigned int r=0;
                for(unsigned int j1=0;j1<nCovariates;j1++){
                  if(!missingX[i][j1]){
                    Sigma(j,r)=workSigma(j0,j1);
                    r++;
                  }
                }
                j++;
              }
            }
            Tau = Sigma.inverse();
            LLT<MatrixXd> llt;
            sqrtTau = (llt.compute(Tau)).matrixU();
            logDetTau = log(Tau.determinant());
          }
          logPXiGivenZi[i][c]=logPdfMultivarNormal(nNotMissing,xi,muStar,sqrtTau,logDetTau);
        }
      }
    }

  }else if(covariateType.compare("Mixed")==0){

    for(unsigned int i=0;i<nSubjects;i++){
      VectorXd xi=VectorXd::Zero(nContinuousCovs);
      logPXiGivenZi[i].resize(maxNClusters,0);
      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          if(currentParams.z(i)==(int)c){
            logPXiGivenZi[i][c]=currentParams.workLogPXiGivenZi(i);
          }else{
            logPXiGivenZi[i][c]=0;
            for(unsigned int j=0;j<nDiscreteCovs;j++){
              int Xij = currentParams.workDiscreteX(i,j);
              logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
             }
            for(unsigned int j=0;j<nContinuousCovs;j++){
              xi(j)=currentParams.workContinuousX(i,j);
            }
            logPXiGivenZi[i][c]+=logPdfMultivarNormal(nContinuousCovs,xi,currentParams.workMuStar(c),currentParams.workSqrtTau(c),currentParams.workLogDetTau(c));
          }
        }
      }
    }

    // For the predictive subjects we do not count missing data
    for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){
      logPXiGivenZi[i].resize(maxNClusters,0);
      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          for(unsigned int j=0;j<nDiscreteCovs;j++){
            if(!missingX[i][j]){
              int Xij = currentParams.workDiscreteX(i,j);
              logPXiGivenZi[i][c]+=currentParams.workLogPhiStar(c,j,(unsigned int)Xij);
            }
          }
        }
      }
    }

    // For the predictive subjects we do not count missing data

    for(unsigned int i=nSubjects;i<nSubjects+nPredictSubjects;i++){

      unsigned int nNotMissing=dataset.nContinuousCovariatesNotMissing(i);

      for(unsigned int c=0;c<maxNClusters;c++){
        if(u[i]<testBound[c]){
          VectorXd workMuStar=currentParams.workMuStar(c);
          VectorXd xi=VectorXd::Zero(nNotMissing);
          VectorXd muStar=VectorXd::Zero(nNotMissing);
          MatrixXd sqrtTau=MatrixXd::Zero(nNotMissing,nNotMissing);
          double logDetTau =0.0;
          if(nNotMissing==nContinuousCovs){
            muStar=workMuStar;
            sqrtTau=currentParams.workSqrtTau(c);
            logDetTau=currentParams.workLogDetTau(c);
            for(unsigned int j=0;j<nContinuousCovs;j++){
              xi(j)=currentParams.workContinuousX(i,j);
            }
          }else{
            MatrixXd workSigma=currentParams.Sigma(c);
            MatrixXd Sigma=MatrixXd::Zero(nNotMissing,nNotMissing);
            MatrixXd Tau=MatrixXd::Zero(nNotMissing,nNotMissing);
            unsigned int j=0;
            for(unsigned int j0=0;j0<nContinuousCovs;j0++){
              if(!missingX[i][nDiscreteCovs+j0]){
                xi(j)=currentParams.workContinuousX(i,j0);
                muStar(j)=workMuStar(j0);
                unsigned int r=0;
                for(unsigned int j1=0;j1<nContinuousCovs;j1++){
                  if(!missingX[i][nDiscreteCovs+j1]){
                    Sigma(j,r)=workSigma(j0,j1);
                    r++;
                  }
                }
                j++;
              }
            }
            Tau = Sigma.inverse();
            LLT<MatrixXd> llt;
            sqrtTau = (llt.compute(Tau)).matrixU();
            logDetTau = log(Tau.determinant());
          }
          logPXiGivenZi[i][c]+=logPdfMultivarNormal(nNotMissing,xi,muStar,sqrtTau,logDetTau);
        }
      }
    }
  }


  double (*logPYiGivenZiWi)(const pReMiuMParams&, const pReMiuMData&,
          const vector<unsigned int>&,const int&,
          const unsigned int&)=NULL;
  if(includeResponse){
    if(!responseExtraVar){
      if(outcomeType.compare("Bernoulli")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiBernoulli;
      }else if(outcomeType.compare("Binomial")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiBinomial;
      }else if(outcomeType.compare("Poisson")==0){
        if (includeCAR){
          logPYiGivenZiWi = &logPYiGivenZiWiPoissonSpatial;
        }else{
          logPYiGivenZiWi = &logPYiGivenZiWiPoisson;
        }
      }else if(outcomeType.compare("Normal")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiNormal;
      }else if(outcomeType.compare("Categorical")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiCategorical;
      }else if(outcomeType.compare("Survival")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiSurvival;
      }else if(outcomeType.compare("Longitudinal")==0){//RJ set logPYiGivenZiWi
        if(model.options().sampleGPmean()){//AR change
          logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal_meanGP; //AR
        }else{
          logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal;
        }
      }else if(outcomeType.compare("MVN")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiMVN;
      }else if(outcomeType.compare("LME")==0){
        logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal_parametric;
      }
    }
  }

  vector<double> meanVec(nSubjects,0.0);
  if(includeResponse){
    if(outcomeType.compare("Poisson")==0){
      meanVec = dataset.logOffset();
    }
  }

  unsigned int maxZ=0;
  //RJ for Longitudinal: save cluster marginal likelihoods, update only when membership changes
  vector<double> clusterMarginal, numerator, denominator;
  clusterMarginal.resize(maxNClusters,0);
  numerator.resize(maxNClusters,0);
  denominator.resize(maxNClusters,0);
  int origZi;


  //AR
  // double * det_M0;
  //det_M0 = new double [maxNClusters];
  vector<double> det_M0(maxNClusters,0);
  vector<MatrixXd> Sigma_inv_c_ord(maxNClusters);
  vector<int> sizek(maxNClusters,0);

  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> times = dataset.times();
  int Ana=1; // 0 ROB, 1 AR, 2 both
  vector<double> grid;
  //double minui=*std::min_element(std::begin(u), std::end(u) );

  for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){

    vector<double> logPyXz(maxNClusters,0.0);
    // We calculate p(Z=c|y,X) \propto p(y,X,Z=c)
    // p(y,X,z=c) = p(y|Z=c)p(X|z=c)p(z=c)
    double maxLogPyXz = -(numeric_limits<double>::max());

    if(includeResponse&&i<nSubjects){
      // Response only included in allocation probs for fitting subjects, not predicting subjects
      if(responseExtraVar){
        // In this case the Y only go into this conditional
        // through lambda
        if(outcomeType.compare("LME")!=0 ){
          for(unsigned int c=0;c<maxNClusters;c++){
            if(u[i]<testBound[c]){
              double meanVal=meanVec[i]+currentParams.theta(c,0);
              for(unsigned int j=0;j<nFixedEffects[0];j++){
                meanVal+=currentParams.beta(0,j,0,nCategoriesY)*dataset.W(i,j);
              }
              // for(unsigned int j=0;j<nFixedEffects_mix;j++){
              //   meanVal+=currentParams.beta_mix(c,j,0,nCategoriesY)*dataset.W_mix(i,j);
              // }
              logPyXz[c]+=logPdfNormal(currentParams.lambda(i),meanVal,
                                       1/sqrt(currentParams.tauEpsilon()));
            }
          }
        }else{
          printf("responseExtraVar not possible with yModel = LME.");
        }
      }else{
        // In this case the Y go in directly
        //RJ for longitudinal, pYi = marginal(cluster with gene)/marginal(cluster without gene)
        //save original cluster label
        origZi = currentParams.z(i);

        if(i==0 && outcomeType.compare("Longitudinal")==0 && Ana>0){ //grid construction
          //int grid_size = model.dataset().nTimes_unique();
          grid=dataset.times_unique();

          // double min_grid=*std::min_element(std::begin(times), std::end(times));
          // double pas_grid=(*std::max_element(std::begin(times), std::end(times) )-*std::min_element(std::begin(times), std::end(times)))/(grid_size-1);
          //
          // for(unsigned int i2=0;i2<grid_size;i2++){
          //   grid[i2]= min_grid+i2*pas_grid;
          // }
        }

        //#pragma omp parallel for
        for(unsigned int c=0;c<maxNClusters;c++){

          double temp;
          vector<double> timesk;

          if( !model.options().sampleGPmean() && i==0 && outcomeType.compare("Longitudinal")==0 ){//for first person, calculate all marginals
            // Computation of the initial cluster-specific likelihoods clusterMarginal[c]=p(Y|c=g)
            if(Ana!=1)
              clusterMarginal[c] = logPYiGivenZiWi(currentParams,dataset,nFixedEffects,c,nSubjects);

            if(Ana==2)
              temp=clusterMarginal[c];


            if(Ana>0){
              //AR Compute number of timepoints of subjects in cluster c to compute inverse Sigma_c
              sizek[c]=0;

              // set sizes based on cluster occupation
              for(unsigned int i2=0;i2<nSubjects;i2++){
                if(currentParams.z(i2) == c){
                  sizek[c] = sizek[c] + tStop[i2] - tStart[i2] + 1;
                }
              }

              if(sizek[c]>0){
                Sigma_inv_c_ord[c].setZero(sizek[c],sizek[c]);
                timesk.resize(sizek[c]);

                int counter = 0;
                for(unsigned int i2=0;i2<nSubjects;i2++){
                  if(currentParams.z(i2) == c){
                    for(unsigned int j=0;j<tStop[i2]-tStart[i2]+1;j++){
                      timesk[counter+j] = times[tStart[i2]-1+j];
                    }
                    counter = counter + tStop[i2] - tStart[i2] + 1;
                  }
                }

                //AR compute covariance matrix
                //fout << " A "<< c << " det_M0 "<< det_M0[c]<< " clusterMarg "<< clusterMarginal[c] <<endl;

                det_M0[c] = Get_Sigma_inv_GP_cov(Sigma_inv_c_ord[c],currentParams.L(c),timesk,dataset.equalTimes(),grid,kernelType);
                //fout <<  " B "<<c << " det_M0 "<< det_M0[c]<< " clusterMarg "<< clusterMarginal[c] <<endl;
                clusterMarginal[c] = logPYiGivenZiWiLongitudinal_bis(Sigma_inv_c_ord[c],det_M0[c],currentParams,dataset,nFixedEffects,c,sizek[c],nSubjects);
                //cout << " i "<< i << " "<< c << " clusMarg "<< clusterMarginal[c] <<endl ;
              }

              if((Ana==2) & (abs(temp- clusterMarginal[c])>pow (1.0, -5.0))){
                std::cout << i << " c "<< c << " det_M0[c]" << det_M0[c]<< endl;
                std::cout << i << " c "<< c << " clusterMarginal " << clusterMarginal[c]<< " temp" << temp << endl<< endl;
              }
            }
          }//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          if(u[i]<testBound[c]){

            if(outcomeType.compare("Longitudinal")==0){

              if(model.options().sampleGPmean()){//AR change model.options().sampleGPmean()

                double temp2 = logPYiGivenZiWi(currentParams,dataset,nFixedEffects,c,i); // only if u[i]<testBound[c]
                logPyXz[c] += temp2;
                //logPyXz[c] += logPdfMultivariateNormal(currentParams.meanGP(c),currentParams.L(c), dataset.times_unique(), kernelType);

                sizek[c]=0;

                // set sizes based on cluster occupation
                for(unsigned int i2=0;i2<nSubjects;i2++){
                  if(currentParams.z(i2) == c){
                    sizek[c] +=  tStop[i2] - tStart[i2] + 1;
                  }
                }

                if(std::isinf(-logPyXz[c]))
                  cout <<" i "<<i<< " c " <<c << " u[i] "<<u[i]<< " testBound[c] "<<testBound[c]<<" sizek[c] "<<sizek[c]<<" logPyXz[c] "<<logPyXz[c]<<endl;


              }else{//RJ marginal likelihood without i

                if(origZi == c){
                  //if point currently in cluster:
                  //marginal(cluster with gene) = clusterMarginal[c]
                  numerator[c] = clusterMarginal[c]; // likelihood with i

                  //get marginal(cluster without gene)
                  if(Ana!=1){
                    denominator[c] = logPYiGivenZiWi(currentParams,dataset,nFixedEffects,c,i); // likelihood without i
                    if(Ana==2)
                      temp=denominator[c];
                  }

                  if(Ana>0){
                    denominator[c] = logPYiGivenZiWiLongitudinal_bis(Sigma_inv_c_ord[c],det_M0[c],currentParams,dataset,nFixedEffects,c,sizek[c],i); // likelihood without i
                  }

                  if((Ana==2) & (abs(temp- denominator[c])>pow (1.0, -5.0))){
                    std::cout << i << " c "<< c << " det_M0[c]" << det_M0[c]<< endl;
                    std::cout << i << " c "<< c << " denominator " << denominator[c]<< " temp" << temp << endl<<endl;
                  }
                }else{
                  //take cluster marginal for denominator
                  denominator[c] = clusterMarginal[c]; // likelihood without i

                  if(Ana!=1){
                    //move point for numerator
                    currentParams.z(i,c,covariateType);//RJ marginal likelihood with i
                    numerator[c] = logPYiGivenZiWi(currentParams,dataset,nFixedEffects,c,nSubjects); // likelihood with i
                    currentParams.z(i,origZi,covariateType);//RJ set cluster to original
                    if(Ana==2)
                      temp=numerator[c];
                  }

                  if(Ana>0){
                    numerator[c] = logPYiGivenZiWiLongitudinal_bis(Sigma_inv_c_ord[c],det_M0[c],currentParams,dataset,nFixedEffects,c,sizek[c],i);// likelihood with i
                  }

                  if((Ana==2) & (abs(temp- numerator[c])>pow(1.0, -5.0))){
                    std::cout << i << " c "<< c << " det_M0[c]" << det_M0[c]<< endl;
                    std::cout << i << " c "<< c << " numerator " << numerator[c]<< " temp" << temp << " sizeS0 "<< Sigma_inv_c_ord[c].rows()<< " sizek "<< sizek[c] << endl;
                  }
                }
                logPyXz[c]+= numerator[c] - denominator[c];
              }
            }else{
              logPyXz[c]+=logPYiGivenZiWi(currentParams,dataset,nFixedEffects,c,i);
            }
          }
        }
      }
    }

    for(unsigned int c=0;c<maxNClusters;c++){

      if(u[i]<testBound[c]){
        // Make sure prediction subjects can only be allocated to one
        // of the non-empty clusters
        if(i<nSubjects||nMembers[c]>0){
          logPyXz[c]+=clusterWeight[c];
          logPyXz[c]+=logPXiGivenZi[i][c];
        }else{
          logPyXz[c]=-(numeric_limits<double>::max());
        }
      }else{
        logPyXz[c]=-(numeric_limits<double>::max());
      }
      //AR
      if(std::isinf(logPyXz[c]))
        logPyXz[c] = -(numeric_limits<double>::max());

      if(logPyXz[c]>maxLogPyXz){
        maxLogPyXz=logPyXz[c];
      }
    }

    vector<double> pzGivenXy(maxNClusters);
    double sumVal=0;
    for(unsigned int c=0;c<maxNClusters;c++){
      double exponent = logPyXz[c] - maxLogPyXz;

      // Check for negative infinity (can only be negative)
      if(std::isinf(exponent)||std::isnan(exponent)){
        exponent=-(numeric_limits<double>::max());
      }
      pzGivenXy[c]=exp(exponent);
      sumVal+=pzGivenXy[c];
    }

    vector<double> expectedTheta(nCategoriesY);
    double entropyVal=0.0;
    vector<double> cumPzGivenXy(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      pzGivenXy[c]/=sumVal;
      if(computeEntropy){
        if(pzGivenXy[c]>0){
          entropyVal-=pzGivenXy[c]*log(pzGivenXy[c]);
        }
      }

      if(c==0){
        cumPzGivenXy[c]=pzGivenXy[c];
      }else{
        cumPzGivenXy[c]=cumPzGivenXy[c-1]+pzGivenXy[c];
      }
      if (predictType.compare("RaoBlackwell")==0){
        if(includeResponse&&i>=nSubjects){
          if(outcomeType.compare("Categorical")==0){
            for (unsigned int k=0;k<nCategoriesY;k++){
              expectedTheta[k]+=currentParams.theta(c,k)*pzGivenXy[c];
            }
          } else {
            expectedTheta[0]+=currentParams.theta(c,0)*pzGivenXy[c];
          }
        }
      }
    }

    if(includeResponse&&i>=nSubjects){
      if (predictType.compare("random")==0){
        // choose which component of the mixture we are sampling from
        double u=unifRand(rndGenerator);
        unsigned int c=0;
        while(cumPzGivenXy[c]<=u){
          c++;
        }
        // draw from the normal distribution of that sample (only for yModel=Normal)
        // Create a normal random generator
        randomNormal normRand(0,1);
        expectedTheta[0]=sqrt(currentParams.sigmaSqY())*normRand(rndGenerator)+currentParams.theta(c,0);
      }
    }

    unsigned int zi;

    if(maxNClusters==1){
      zi=0;
    }else{
      zi = 0;
      for(unsigned int c=0;c<maxNClusters;c++){
        if(rnd[i]<cumPzGivenXy[c]){
          zi=c;
          break;
        }
      }
    }
    if(zi>maxZ){
      maxZ=zi;
    }

    currentParams.z(i,zi,covariateType);

    // if(i<130){
    //   currentParams.z(i,0,covariateType);
    // }else{
    //   currentParams.z(i,1,covariateType);
    // }

    //AR Compute again sizek, yk and timeskfor z(i) and currentParams.z(i)
    if(outcomeType.compare("Longitudinal")==0 && zi!=origZi && !model.options().sampleGPmean() ){//AR change !model.options().sampleGPmean()

      //RJ update clusterMarginals if a subject has swapped clusters
      clusterMarginal[zi] = numerator[zi];
      clusterMarginal[origZi] = denominator[origZi];

      if(Ana>0){
        sizek[origZi] =sizek[origZi] -(tStop[i] - tStart[i] + 1) ;
        sizek[currentParams.z(i)] =sizek[currentParams.z(i)] + (tStop[i] - tStart[i] + 1) ;

        //AR Compute  clustermarginal again
        int unsigned c=origZi;
        Sigma_inv_c_ord[c].setZero(sizek[c],sizek[c]);

        if(sizek[c]>0){
          vector<double> timesk;
          timesk.resize(sizek[c]);

          int counter = 0;

          for(unsigned int i2=0;i2<nSubjects;i2++){
            if(currentParams.z(i2) == c){
              for(unsigned int j=0;j<tStop[i2]-tStart[i2]+1;j++){
                timesk[counter+j] = times[tStart[i2]-1+j];
              }
              counter = counter + tStop[i2] - tStart[i2] + 1;
            }
          }

          det_M0[c] = Get_Sigma_inv_GP_cov(Sigma_inv_c_ord[c],currentParams.L(c),timesk,dataset.equalTimes(),grid,kernelType);
        }else{
          det_M0[c]=0;
        }

        c=zi;

        Sigma_inv_c_ord[c].setZero(sizek[c],sizek[c]);
        vector<double> timesk;
        timesk.resize(sizek[c]);

        int counter = 0;
        for(unsigned int i2=0;i2<nSubjects;i2++){
          if(currentParams.z(i2) == c){
            for(unsigned int j=0;j<tStop[i2]-tStart[i2]+1;j++){
              timesk[counter+j] = times[tStart[i2]-1+j];
            }
            counter = counter + tStop[i2] - tStart[i2] + 1;
          }
        }
        det_M0[c] = Get_Sigma_inv_GP_cov(Sigma_inv_c_ord[c],currentParams.L(c),timesk,dataset.equalTimes(),grid,kernelType);
      }
    }

    if(computeEntropy){
      currentParams.workEntropy(i,entropyVal);
    }

    if(i<nSubjects){
      nMembers[zi]++;
    }else{
      if(includeResponse){
        for (unsigned int k=0;k<nCategoriesY;k++){
          currentParams.workPredictExpectedTheta(i-nSubjects,k,expectedTheta[k]);
        }
      }
    }
  }
  currentParams.workNXInCluster(nMembers);
  currentParams.workMaxZi(maxZ);
}


void updateMissingPReMiuMData(baseGeneratorType& rndGenerator,
                              pReMiuMParams& params,
                              const pReMiuMOptions& options,
                              pReMiuMData& dataset){

  unsigned int nSubjects = dataset.nSubjects();
  unsigned int nCovariates = dataset.nCovariates();
  unsigned int nDiscreteCovs = dataset.nDiscreteCovs();
  unsigned int nContinuousCovs = dataset.nContinuousCovs();
  vector<unsigned int> nCategories = params.nCategories();
  string covariateType = options.covariateType();

  // Define a uniform random number generator
  randomUniform unifRand(0,1);

  // For the missing variables we check which cluster the subject is allocated
  // to and then sample for the appropriate
  if(covariateType.compare("Discrete")==0){
    // We don't update the predictive subjects as their X values which
    // were missing are not used anywhere
    // change: we now compute them because we are interested in looking at their posterior predictive distributions
    for(unsigned int i=0;i<nSubjects;i++){
      int zi = params.z(i);
      for(unsigned int j=0;j<nCovariates;j++){
        if(dataset.missingX(i,j)){
          vector<double> cumPhiStar(nCategories[j]);
          // Sample uniform
          double u=unifRand(rndGenerator);
          int k=0;
          cumPhiStar[k]=exp(params.workLogPhiStar(zi,j,k));
          while(u>cumPhiStar[k]){
            k++;
            // Make phi cumulative
            cumPhiStar[k]=cumPhiStar[k-1]+exp(params.workLogPhiStar(zi,j,k));
          }

          int prevX = params.workDiscreteX(i,j);
          dataset.discreteX(i,j,k);
          params.workDiscreteX(i,j,k);
          // Now we need to recompute the workLogPXiGivenZi values based
          // on the new X
          if(prevX!=k){
            double logVal = params.workLogPXiGivenZi(i);
            double oldVal, newVal;
            oldVal = params.workLogPhiStar(zi,j,prevX);
            newVal = params.workLogPhiStar(zi,j,k);
            logVal+=(newVal-oldVal);
            params.workLogPXiGivenZi(i,logVal);
          }
        }
      }
    }
  }else if(covariateType.compare("Normal")==0){
    for(unsigned int i=0;i<nSubjects;i++){
      // Check if there is anything to do
      if(dataset.nContinuousCovariatesNotMissing(i)<nCovariates){
        int zi = params.z(i);
        VectorXd newXi=multivarNormalRand(rndGenerator,params.workMuStar(zi),params.Sigma(zi));
        for(unsigned int j=0;j<nCovariates;j++){
          if(dataset.missingX(i,j)){
            dataset.continuousX(i,j,newXi(j));
            params.workContinuousX(i,j,newXi(j));
          }else{
            newXi(j)=dataset.continuousX(i,j);
          }
        }
        double logVal = logPdfMultivarNormal(nCovariates,newXi,params.workMuStar(zi),params.workSqrtTau(zi),params.workLogDetTau(zi));
        params.workLogPXiGivenZi(i,logVal);
      }
    }

  }else if(covariateType.compare("Mixed")==0){
    // Discrete part of mixed type covariates
    // We don't update the predictive subjects as their X values which
    // were missing are not used anywhere
    for(unsigned int i=0;i<nSubjects;i++){
      int zi = params.z(i);
      for(unsigned int j=0;j<nDiscreteCovs;j++){
        if(dataset.missingX(i,j)){
          vector<double> cumPhiStar(nCategories[j]);
          // Sample uniform
          double u=unifRand(rndGenerator);
          int k=0;
          cumPhiStar[k]=exp(params.workLogPhiStar(zi,j,k));
          while(u>cumPhiStar[k]){
            k++;
            // Make phi cumulative
            cumPhiStar[k]=cumPhiStar[k-1]+exp(params.workLogPhiStar(zi,j,k));
          }

          int prevX = params.workDiscreteX(i,j);
          dataset.discreteX(i,j,k);
          params.workDiscreteX(i,j,k);
          // Now we need to recompute the workLogPXiGivenZi values based
          // on the new X
          if(prevX!=k){
            double logVal = params.workLogPXiGivenZi(i);
            double oldVal, newVal;
            oldVal = params.workLogPhiStar(zi,j,prevX);
            newVal = params.workLogPhiStar(zi,j,k);
            logVal+=(newVal-oldVal);
            params.workLogPXiGivenZi(i,logVal);
          }
        }
      }
    }

    // Normal part of mixed type covariates
    for(unsigned int i=0;i<nSubjects;i++){
      // Check if there is anything to do
      if(dataset.nContinuousCovariatesNotMissing(i)<nContinuousCovs){
        int zi = params.z(i);
        VectorXd newXi=multivarNormalRand(rndGenerator,params.workMuStar(zi),params.Sigma(zi));
        for(unsigned int j=0;j<nContinuousCovs;j++){
          if(dataset.missingX(i,j)){
            dataset.continuousX(i,j,newXi(j));
            params.workContinuousX(i,j,newXi(j));
          }else{
            newXi(j)=dataset.continuousX(i,j);
          }
        }
        double logVal = logPdfMultivarNormal(nCovariates,newXi,params.workMuStar(zi),params.workSqrtTau(zi),params.workLogDetTau(zi));
        params.workLogPXiGivenZi(i,logVal);
      }
    }
  }
}


#endif /* DIPBACPROPOSALS_H_ */
