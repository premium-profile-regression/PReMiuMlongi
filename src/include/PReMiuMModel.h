/// \file PReMiuMModel.h
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


#ifndef DIPBACMODEL_H_
#define DIPBACMODEL_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<numeric>
#include<limits>
#include<map>
#include<sys/time.h>
//#include <algorithm>

#include<boost/math/distributions/normal.hpp>
#include<boost/math/distributions/gamma.hpp>
#include<boost/math/distributions/beta.hpp>
#include<boost/math/distributions/students_t.hpp>
#include<boost/math/special_functions/gamma.hpp>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>
#include<Eigen/CholmodSupport>
#include <Eigen/Eigenvalues>

// Custom includes
#include"MCMC/model.h"
#include<MCMC/state.h>
#include<Math/random.h>
#include<Math/distribution.h>
#include<Math/mathfunctions.h>
#include<PReMiuMData.h>
#include<PReMiuMOptions.h>


using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;
using std::accumulate;
using std::numeric_limits;
using std::cout;
using std::endl;
using std::max;

using namespace Eigen;

using boost::math::normal_distribution;
using boost::math::gamma_distribution;
using boost::math::beta_distribution;
using boost::math::students_t_distribution;
using boost::math::lgamma;


class pReMiuMHyperParams{

public:
  // Default constructor
  pReMiuMHyperParams() {};

  // Default virtual destructor
  ~pReMiuMHyperParams() {};

  // Set sizes
  void setSizes(const unsigned int& nCovariates,
                const unsigned int& nDiscreteCov,
                const unsigned int& nContinuousCov,
                const unsigned int& nOutcomes,
                const string& covariateType,
                const string& outcomeType,
                const vector<unsigned int>& nRandomEffects){
    if (covariateType.compare("Discrete")==0){
      _aPhi.resize(nCovariates);
    } else if (covariateType.compare("Normal")==0){
      _mu0.setZero(nCovariates);
      _Tau0.setZero(nCovariates,nCovariates);
      _workSqrtTau0.setZero(nCovariates,nCovariates);
      _R0.setZero(nCovariates,nCovariates);
      _workInverseR0.setZero(nCovariates,nCovariates);
    }
    else if (covariateType.compare("Mixed")==0 ){
      _aPhi.resize(nDiscreteCov);
      _mu0.setZero(nContinuousCov);
      _Tau0.setZero(nContinuousCov,nContinuousCov);
      _workSqrtTau0.setZero(nContinuousCov,nContinuousCov);
      _R0.setZero(nContinuousCov,nContinuousCov);
      _workInverseR0.setZero(nContinuousCov,nContinuousCov);
    }
    if (outcomeType.compare("MVN")==0){
      _MVNmu0.setZero(nOutcomes);
      _MVNTau0.setZero(nOutcomes,nOutcomes);
      _workSqrtMVNTau0.setZero(nOutcomes,nOutcomes);
      _MVNR0.setZero(nOutcomes,nOutcomes);
      _workInverseMVNR0.setZero(nOutcomes,nOutcomes);
    }

    if (outcomeType.compare("LME")==0){
      _SigmaLME_R0.resize(nOutcomes);
      _workTauLME_R0.resize(nOutcomes);
      _workLogDetTauLME_R0.resize(nOutcomes);
      _SigmaLME_kappa0.resize(nOutcomes);
    }
  }

  // Set defaults
  void setDefaults(const pReMiuMData& dataset,
                   const pReMiuMOptions& options){

    unsigned int nSubjects = dataset.nSubjects();
    unsigned int nOutcomes = dataset.nOutcomes();
    //unsigned int nCov = dataset.nCovariates();
    //unsigned int nDiscreteCov = dataset.nDiscreteCovs();
    //unsigned int nContinuousCov = dataset.nContinuousCovs();

    // For alpha
    _shapeAlpha=2.0;
    _rateAlpha=1.0;

    if(options.covariateType().compare("Discrete")==0 || options.covariateType().compare("Mixed")==0 ){
      // For Phi
      for(unsigned int j=0;j<_aPhi.size();j++){
        _aPhi[j]=1.0;
      }
    }

    if(options.covariateType().compare("Normal")==0 || options.covariateType().compare("Mixed")==0 ){
      unsigned int nContCovs=_mu0.size();
      // Values for mu, Tau0, R0, kappa0
      // In the following it is useful to have the rows of X as
      // Eigen vectors
      vector<VectorXd> xi(nSubjects);
      for(unsigned int i=0;i<nSubjects;i++){
        xi[i].setZero(nContCovs);
        for(unsigned int j=0;j<nContCovs;j++){
          xi[i](j)=dataset.continuousX(i,j);
        }
      }

      // First compute the hyper prior parameters
      // First compute the hyper parameter for mu_c
      VectorXd mu0=VectorXd::Zero(nContCovs);
      MatrixXd Sigma0=MatrixXd::Zero(nContCovs,nContCovs);
      for(unsigned int j=0;j<nContCovs;j++){
        double meanX=0.0;
        double minX=0;
        double maxX=0;
        for(unsigned int i=0;i<nSubjects;i++){
          double tmpX = xi[i](j);
          meanX+= tmpX;
          if(tmpX<minX||i==0){
            minX=tmpX;
          }
          if(tmpX>maxX||i==0){
            maxX=tmpX;
          }
        }
        meanX/=(double)(nSubjects-1);
        double rangeX = maxX-minX;
        mu0(j)=meanX;
        Sigma0(j,j)=rangeX*rangeX;
      }
      MatrixXd Tau0 = Sigma0.inverse();
      _Tau0=Tau0;
      LLT<MatrixXd> llt;
      _workSqrtTau0=(llt.compute(Tau0)).matrixU();
      _workLogDetTau0 = log(Tau0.determinant());
      _mu0=mu0;

      // Now we compute the hyper parameters for Tau_c
      // First we compute the sample covariance
      MatrixXd R=MatrixXd::Zero(nContCovs,nContCovs);

      for(unsigned int i=0;i<nSubjects;i++){
        R += (xi[i]-mu0)*((xi[i]-mu0).transpose());
      }
      R/=(double)nSubjects;
      R*=(double)nContCovs;
      _workInverseR0=R;
      R=R.inverse();
      _R0=R;
      _workLogDetR0=log(R.determinant());
      _kappa0 = nContCovs;
      _nu0=0.01;
    }
    //RJ defaults for MVN parameters
    if(options.outcomeType().compare("MVN")==0){
      unsigned int nTimes = dataset.nTimes();
      unsigned int nOutcomes = nTimes/nSubjects;
      // Values for mu, Tau0, R0, kappa0
      // In the following it is useful to have the rows of X as
      // Eigen vectors
      vector<VectorXd> yi(nSubjects);
      for(unsigned int i=0;i<nSubjects;i++){
        yi[i].setZero(nOutcomes);
        for(unsigned int j=0;j<nOutcomes;j++){
          yi[i](j)=dataset.continuousY(i*nOutcomes+j);
        }
      }

      // First compute the hyper prior parameters
      // First compute the hyper parameter for mu_c
      VectorXd mu0=VectorXd::Zero(nOutcomes);
      MatrixXd Sigma0=MatrixXd::Zero(nOutcomes,nOutcomes);
      for(unsigned int j=0;j<nOutcomes;j++){
        double meanX=0.0;
        double minX=0;
        double maxX=0;
        for(unsigned int i=0;i<nSubjects;i++){
          double tmpX = yi[i](j);
          meanX+= tmpX;
          if(tmpX<minX||i==0){
            minX=tmpX;
          }
          if(tmpX>maxX||i==0){
            maxX=tmpX;
          }
        }
        meanX/=(double)(nSubjects-1);
        double rangeX = maxX-minX;
        mu0(j)=meanX;
        Sigma0(j,j)=rangeX*rangeX;
      }
      MatrixXd Tau0 = Sigma0.inverse();
      _MVNTau0=Tau0;
      LLT<MatrixXd> llt;
      _workSqrtMVNTau0=(llt.compute(Tau0)).matrixU();
      _workLogDetMVNTau0 = log(Tau0.determinant());
      _MVNmu0=mu0;

      // Now we compute the hyper parameters for Tau_c
      // First we compute the sample covariance
      MatrixXd R=MatrixXd::Zero(nOutcomes,nOutcomes);

      for(unsigned int i=0;i<nSubjects;i++){
        R += (yi[i]-mu0)*((yi[i]-mu0).transpose());
      }
      R/=(double)nSubjects;
      R*=(double)nOutcomes;
      _workInverseMVNR0=R;
      R=R.inverse();
      _MVNR0=R;
      _workLogDetMVNR0=log(R.determinant());
      _MVNnu0 = nOutcomes;
      _MVNkappa0=0.01;
    }


    //AR defaults for LME parameters
    if(options.outcomeType().compare("LME")==0){
      vector<unsigned int> nRandomEffects = dataset.nRandomEffects();
      // Now we compute the hyper parameters for Tau_c
      // First we compute the sample covariance
      //_workLogDetTauLME_R0.resize(nOutcomes);
      //_SigmaLME_kappa0.resize(nOutcomes);
      for(unsigned int m=0;m<nOutcomes;m++){
        _SigmaLME_R0[m] = 0.01*MatrixXd::Identity(nRandomEffects[m],nRandomEffects[m]);
        _workTauLME_R0[m]=_SigmaLME_R0[m].inverse();
        //_SigmaLME_R0[m]=MatrixXd::Identity(nRandomEffects[m],nRandomEffects[m]);
        //_workTauLME_R0[m]=_SigmaLME_R0[m].inverse();
        _workLogDetTauLME_R0[m] = -log(_SigmaLME_R0[m].determinant());
        _SigmaLME_kappa0[m]=nRandomEffects[m];
      }
    }

    _muTheta = 0.0;
    _sigmaTheta = 2.5;
    _dofTheta = 7;

    _muBeta = 0.0;
    _sigmaBeta = 2.5;
    _dofBeta = 7;

    _aRatio=1;
    _bRatio=1;

    //RJ default values for muL and sigmaL
    unsigned int nL;
    if(options.kernelType().compare("SQexponential")==0){    //AR
      nL=3;
    }else{
      nL=4;
    }
    _muL.resize(nL);
    _sigmaL.resize(nL);
    for(unsigned int i=0;i<nL;i++){
      _muL[i] = 0.0;
      _sigmaL[i] = 1.0;
    }

    _eps_sigma2_0 = 0.5;
    _eps_vu = 1.0;
    _eps_shape = 1.0;
    _eps_rate = 1.0;
    _shapeTauEpsilon = 5.0;
    _rateTauEpsilon = 0.5;

    _aRho = 0.5;
    _bRho = 0.5;
    _atomRho= 0.5;

    _shapeSigmaSqY = 2.5;
    _scaleSigmaSqY = 2.5;

    _rSlice =0.75;
    _truncationEps = 0.000001;

    _shapeTauCAR=0.001;
    _rateTauCAR=0.001;

    _shapeNu = 2.5;
    _scaleNu = 1;

  }

  double shapeAlpha() const{
    return _shapeAlpha;
  }

  void shapeAlpha(const double& s){
    _shapeAlpha = s;
  }

  double rateAlpha() const{
    return _rateAlpha;
  }

  void rateAlpha(const double& r){
    _rateAlpha = r;
  }

  vector<double> aPhi() const{
    return _aPhi;
  }

  void aPhi(const vector<double>& a){
    _aPhi = a;
  }

  double aPhi(const unsigned int& j) const{
    return _aPhi[j];
  }


  /// \brief Return the hyper parameter Tau0Mu0
  const VectorXd& mu0() const{
    return _mu0;
  }
  /// \brief Set the hyper parameter Tau0Mu0
  void mu0(const VectorXd& m0){
    _mu0 = m0;
  }
  /// \brief Return the hyper parameter Tau0
  const MatrixXd& Tau0() const{
    return _Tau0;
  }
  /// \brief Set the hyper parameter Tau0
  void Tau0(const MatrixXd& T0) {
    _Tau0 = T0;
    _workLogDetTau0 = log(T0.determinant());
    LLT<MatrixXd> llt;
    _workSqrtTau0=(llt.compute(T0)).matrixU();
  }
  /// \brief Return the hyper parameter R0
  const MatrixXd& R0() const{
    return _R0;
  }
  /// \brief Set the hyper parameter R0
  void R0(const MatrixXd& R) {
    _R0 = R;
    _workLogDetR0 = log(R.determinant());
    _workInverseR0 = R.inverse();
  }
  /// \brief Return the hyper parameter kappa0
  const unsigned int& kappa0() const{
    return _kappa0;
  }
  /// \brief Set the hyper parameter kappa0
  void kappa0(const unsigned int& k0){
    _kappa0 = k0;
  }
  /// \brief Return the hyper parameter nu0
  const double& nu0() const{
    return _nu0;
  }
  /// \brief Set the hyper parameter nu0
  void nu0(const double& n0){
    _nu0 = n0;
  }
  //RJ handling functions for MVN parameters
  const VectorXd& MVNmu0() const{
    return _MVNmu0;
  }
  void MVNmu0(const VectorXd& m0){
    _MVNmu0 = m0;
  }
  const MatrixXd& MVNTau0() const{
    return _MVNTau0;
  }
  void MVNTau0(const MatrixXd& T0) {
    _MVNTau0 = T0;
    _workLogDetMVNTau0 = log(T0.determinant());
    LLT<MatrixXd> llt;
    _workSqrtMVNTau0=(llt.compute(T0)).matrixU();
  }
  const MatrixXd& MVNR0() const{
    return _MVNR0;
  }
  void MVNR0(const MatrixXd& R) {
    _MVNR0 = R;
    _workLogDetMVNR0 = log(R.determinant());
    _workInverseMVNR0 = R.inverse();
  }
  const unsigned int& MVNnu0() const{
    return _MVNnu0;
  }
  void MVNnu0(const unsigned int& n0){
    _MVNnu0 = n0;
  }
  const double& MVNkappa0() const{
    return _MVNkappa0;
  }
  void MVNkappa0(const double& k0){
    _MVNkappa0 = k0;
  }

  const vector<MatrixXd>& SigmaLME_R0() const{
    return _SigmaLME_R0;
  }

  const MatrixXd& SigmaLME_R0(const unsigned int& m) const{
    return _SigmaLME_R0[m];
  }

  const vector<MatrixXd>& workTauLME_R0() const{
    return _workTauLME_R0;
  }

  const MatrixXd& workTauLME_R0(const unsigned int& m) const{
    return _workTauLME_R0[m];
  }

  /// \brief Set the hyper parameter R0
  void SigmaLME_R0(const unsigned int& m, const MatrixXd& R) {
    _SigmaLME_R0[m] = R;
    _workLogDetTauLME_R0[m] = - log(R.determinant());
    _workTauLME_R0[m] = R.inverse();
  }

  const vector<double>& workLogDetTauLME_R0() const{
    return _workLogDetTauLME_R0;
  }

  const double& workLogDetTauLME_R0(const unsigned int& m) const{
    return _workLogDetTauLME_R0[m];
  }

   vector<double>& workLogDetTauLME_R0() {
    return _workLogDetTauLME_R0;
  }

  double& workLogDetTauLME_R0(const unsigned int& m) {
    return _workLogDetTauLME_R0[m];
  }


  const vector<double>& SigmaLME_kappa0() const{
    return _SigmaLME_kappa0;
  }
  void SigmaLME_kappa0(const unsigned int& m, const double& k0){
    _SigmaLME_kappa0[m] = k0;
  }

  const double& SigmaLME_kappa0(const unsigned int& m) const{
    return _SigmaLME_kappa0[m];
  }

  double muTheta() const{
    return _muTheta;
  }

  void muTheta(const double& mu){
    _muTheta = mu;
  }

  double sigmaTheta() const{
    return _sigmaTheta;
  }

  void sigmaTheta(const double& sigma){
    _sigmaTheta = sigma;
  }

  unsigned int dofTheta() const{
    return _dofTheta;
  }

  void dofTheta(const unsigned int& dof){
    _dofTheta = dof;
  }

  double muBeta() const{
    return _muBeta;
  }
  void muBeta(const double& mu){
    _muBeta = mu;
  }
  double sigmaBeta() const{
    return _sigmaBeta;
  }
  void sigmaBeta(const double& sigma){
    _sigmaBeta = sigma;
  }
  unsigned int dofBeta() const{
    return _dofBeta;
  }
  void dofBeta(const unsigned int& dof){
    _dofBeta = dof;
  }

  //RJ handling functions for muL and sigmaL
  vector <double> muL() const{
    return _muL;
  }
  double muL(const unsigned int& i) const{
    return _muL[i];
  }
  void muL(const unsigned int& i,const double& mu){
    _muL[i] = mu;
  }
  vector <double> sigmaL() const{
    return _sigmaL;
  }
  double sigmaL(const unsigned int& i) const{
    return _sigmaL[i];
  }
  void sigmaL(const unsigned int& i,const double& sigma){
    _sigmaL[i] = sigma;
  }

  double eps_sigma2_0() const{
    return _eps_sigma2_0;
  }
  void eps_sigma2_0(const double& sigma){
    _eps_sigma2_0 = sigma;
  }

  double eps_vu() const{
    return _eps_vu;
  }
  void eps_vu(const double& val){
    _eps_vu = val;
  }

  double eps_shape() const{
    return _eps_shape;
  }

  void eps_shape(const double& sigma){
    _eps_shape = sigma;
  }

  double eps_rate() const{
    return _eps_rate;
  }

  void eps_rate(const double& sigma){
    _eps_rate = sigma;
  }

  double shapeTauEpsilon() const{
    return _shapeTauEpsilon;
  }

  void shapeTauEpsilon(const double& a){
    _shapeTauEpsilon = a;
  }

  double rateTauEpsilon() const{
    return _rateTauEpsilon;
  }

  void rateTauEpsilon(const double& b){
    _rateTauEpsilon = b;
  }

  double aRho() const{
    return _aRho;
  }

  void aRho(const double& a){
    _aRho = a;
  }

  double bRho() const{
    return _bRho;
  }

  void bRho(const double& b){
    _bRho = b;
  }

  double aRatio() const{
    return _aRatio;
  }

  void aRatio(const double& a){
    _aRatio = a;
  }

  double bRatio() const{
    return _bRatio;
  }

  void bRatio(const double& b){
    _bRatio = b;
  }

  double atomRho() const{
    return _atomRho;
  }

  void atomRho(const double& atom){
    _atomRho = atom;
  }

  double shapeSigmaSqY() const{
    return _shapeSigmaSqY;
  }

  void shapeSigmaSqY(const double& s){
    _shapeSigmaSqY = s;
  }

  double scaleSigmaSqY() const{
    return _scaleSigmaSqY;
  }

  void scaleSigmaSqY(const double& r){
    _scaleSigmaSqY = r;
  }

  double shapeNu() const{
    return _shapeNu;
  }

  void shapeNu(const double& s){
    _shapeNu = s;
  }

  double scaleNu() const{
    return _scaleNu;
  }

  void scaleNu(const double& r){
    _scaleNu = r;
  }

  double rSlice() const{
    return _rSlice;
  }

  void rSlice(const double& rSl){
    _rSlice=rSl;
  }

  double truncationEps() const{
    return _truncationEps;
  }

  void truncationEps(const double& eps){
    _truncationEps=eps;
  }


  const MatrixXd& workSqrtTau0() const{
    return _workSqrtTau0;
  }
  double workLogDetTau0() const{
    return _workLogDetTau0;
  }
  const MatrixXd& workInverseR0() const{
    return _workInverseR0;
  }
  double workLogDetR0() const{
    return _workLogDetR0;
  }
  //RJ more MVN parameter handling functions
  const MatrixXd& workSqrtMVNTau0() const{
    return _workSqrtMVNTau0;
  }
  double workLogDetMVNTau0() const{
    return _workLogDetMVNTau0;
  }
  const MatrixXd& workInverseMVNR0() const{
    return _workInverseMVNR0;
  }
  double workLogDetMVNR0() const{
    return _workLogDetMVNR0;
  }

  double workXiSlice(unsigned int c) const{
    return (1-_rSlice)*pow(_rSlice,(double)c);
  }

  double shapeTauCAR() const{
    return _shapeTauCAR;
  }

  void shapeTauCAR(const double& a){
    _shapeTauCAR = a;
  }

  double rateTauCAR() const{
    return _rateTauCAR;
  }

  void rateTauCAR(const double& b){
    _rateTauCAR = b;
  }

  vector<double> initAlloc() const{
    return _initAlloc;
  }

  void initAlloc(const vector<double>& c){
    _initAlloc = c;
  }

  double initAlloc(const unsigned int& j) const{
    return _initAlloc[j];
  }

  VectorXd initL(const unsigned int& i) const{
      return _initL.row(i);
  }

  MatrixXd initL() const{
    return _initL;
  }

  void initL(const MatrixXd& L0) {
    _initL = L0;
  }
  // Copy operator
  // pReMiuMHyperParams& operator=(const pReMiuMHyperParams& hyperParams){
  //   _shapeAlpha = hyperParams.shapeAlpha();
  //   _rateAlpha = hyperParams.rateAlpha();
  //   _aPhi = hyperParams.aPhi();
  //   _mu0 = hyperParams.mu0();
  //   _Tau0 = hyperParams.Tau0();
  //   _R0 = hyperParams.R0();
  //   _kappa0 = hyperParams.kappa0();
  //   _nu0 = hyperParams.nu0();
  //   //RJ set MVN parameters
  //   _MVNmu0 = hyperParams.MVNmu0();
  //   _MVNTau0 = hyperParams.MVNTau0();
  //   _MVNR0 = hyperParams.MVNR0();
  //   _MVNkappa0 = hyperParams.MVNkappa0();
  //   _MVNnu0 = hyperParams.MVNnu0();
  //
  //   //AR set SigmaELME parameters
  //   _SigmaLME_R0 = hyperParams.SigmaLME_R0();
  //   _SigmaLME_kappa0 = hyperParams.SigmaLME_kappa0();
  //
  //   _muTheta = hyperParams.muTheta();
  //   _sigmaTheta = hyperParams.sigmaTheta();
  //   _dofTheta = hyperParams.dofTheta();
  //   _muBeta = hyperParams.muBeta();
  //   _sigmaBeta = hyperParams.sigmaBeta();
  //   _dofBeta = hyperParams.dofBeta();
  //   //RJ set _muL and _sigmaL
  //   _muL = hyperParams.muL();
  //   _sigmaL = hyperParams.sigmaL();
  //   _eps_vu = hyperParams.eps_vu();
  //   _eps_sigma2_0 = hyperParams.eps_sigma2_0();
  //   _shapeTauEpsilon = hyperParams.shapeTauEpsilon();
  //   _rateTauEpsilon = hyperParams.rateTauEpsilon();
  //   _aRho = hyperParams.aRho();
  //   _bRho = hyperParams.bRho();
  //   _aRatio = hyperParams.aRatio();
  //   _bRatio = hyperParams.bRatio();
  //   _atomRho = hyperParams.atomRho();
  //   _shapeSigmaSqY = hyperParams.shapeSigmaSqY();
  //   _scaleSigmaSqY = hyperParams.scaleSigmaSqY();
  //   _shapeNu = hyperParams.shapeNu();
  //   _scaleNu = hyperParams.scaleNu();
  //   _workSqrtTau0 = hyperParams.workSqrtTau0();
  //   _workLogDetTau0 = hyperParams.workLogDetTau0();
  //   _workInverseR0 = hyperParams.workInverseR0();
  //   _workLogDetR0 = hyperParams.workLogDetR0();
  //   //RJ more MVN hyperparameters
  //   _workSqrtMVNTau0 = hyperParams.workSqrtMVNTau0();
  //   _workLogDetMVNTau0 = hyperParams.workLogDetMVNTau0();
  //   _workInverseMVNR0 = hyperParams.workInverseMVNR0();
  //   _workLogDetMVNR0 = hyperParams.workLogDetMVNR0();
  //   _rSlice = hyperParams.rSlice();
  //   _truncationEps = hyperParams.truncationEps();
  //   _shapeTauCAR = hyperParams.shapeTauCAR();
  //   _rateTauCAR = hyperParams.rateTauCAR();
  //   _initAlloc = hyperParams.initAlloc();
  //   return *this;
  // }


private:
  // Hyper parameters for prior for alpha
  // prior is alpha ~ Gamma(shape,rate)
  double _shapeAlpha;
  double _rateAlpha;

  // Hyper parameters for prior for Phi_j (discrete categorical covariates)
  // Prior is Phi_j ~ Dirichlet(a,a,...,a)
  vector<double> _aPhi;

  // Hyper parameters for prior for mu (for Normal covariates)
  // Prior is mu ~ N(mu0,inv(Tau0)) or Prior is mu ~ N(mu0,inv(Tau)/nu0) if the Normal inverse Wishart prior is chosen
  VectorXd _mu0;
  MatrixXd _Tau0;
  double _nu0;

  // Hyper parameters for prior of Tau (for Normal covariates)
  // Prior is Tau ~ Wishart(R0,kappa0) (which has mean R0*kappa0).
  MatrixXd _R0;
  unsigned int _kappa0;

  //RJ Hyperparameters for MVN
  VectorXd _MVNmu0;
  MatrixXd _MVNTau0;
  double _MVNkappa0;
  MatrixXd _MVNR0;
  unsigned int _MVNnu0;

  //AR Hyperparameters for LME
  vector<MatrixXd> _SigmaLME_R0;
  vector<MatrixXd> _workTauLME_R0;
  vector<double> _SigmaLME_kappa0;
  vector<double> _workLogDetTauLME_R0;

  // Hyper parameters for prior for theta
  // Prior is location and scale T distribution theta ~ t(mu,sigma,dof)
  // http://www.mathworks.com/help/toolbox/stats/brn2ivz-145.html
  double _muTheta;
  double _sigmaTheta;
  unsigned int _dofTheta;

  // Hyper parameters for prior for beta (for fixed effects)
  // Prior is location and scale T distribution theta ~ t(mu,sigma,dof)
  // http://www.mathworks.com/help/toolbox/stats/brn2ivz-145.html
  double _muBeta;
  double _sigmaBeta;
  unsigned int _dofBeta;
  //RJ declare _muL and _sigmaL
  vector<double> _muL;
  vector<double> _sigmaL;
  double _eps_sigma2_0;
  double _eps_vu;
  double _eps_shape;
  double _eps_rate;

  double _aRatio;
  double _bRatio;

  // Hyper parameters for prior for tauEpsilon (for extra variation)
  // Prior is tauEpsilon ~ Gamma(shape,rate)
  // we use the inverse scale parameterisation so that E[tauEpsilon] = shape/rate
  double _shapeTauEpsilon;
  double _rateTauEpsilon;

  // Hyper parameters for prior for tauEpsilon (for variable selection)
  // Prior is rho ~ Beta(a,b) with a sparsity inducing atom Bernoulli(atomRho)
  double _aRho;
  double _bRho;
  double _atomRho;

  //Hyper parameter for prior for sigma_y^2 (for normal response model)
  // Prior is sigma_y^2 ~ InvGamma(shapeSigmaSqY,scaleSigmaSqY)
  double _shapeSigmaSqY;
  double _scaleSigmaSqY;

  //Hyper parameter for prior for nu, the shape parameter of the Weibull (for survival response model)
  // Prior is nu ~ Gamma(shapeNu,scaleNu)
  double _shapeNu;
  double _scaleNu;

  //Some working variables for speeding up linear algebra
  double _workLogDetTau0;
  MatrixXd _workSqrtTau0;
  double _workLogDetR0;
  MatrixXd _workInverseR0;
  //RJ Some working variables for speeding up MVN linear algebra
  double _workLogDetMVNTau0;
  MatrixXd _workSqrtMVNTau0;
  double _workLogDetMVNR0;
  MatrixXd _workInverseMVNR0;

  //Slice sampler variables for independent slice sampler
  double _rSlice;

  //Truncated sampler variables
  double _truncationEps;

  // Hyper parameters for prior for tauCAR (for spatial CAR term)
  // Prior is tauCAR ~ Gamma(shape,rate)
  // we use the inverse scale parameterisation so that E[tauCAR] = shape/rate
  double _shapeTauCAR;
  double _rateTauCAR;

  // Initial allocations
  vector<double> _initAlloc;
  MatrixXd _initL;
};


/// \class pReMiuMParams PReMiuMModel.h "PReMiuMModel.h"
/// \brief A class for PReMiuM parameters
class pReMiuMParams{

public:
  /// \brief Default constructor
  //pReMiuMParams(){};

  /// \brief Destructor
  //~pReMiuMParams(){};

  /// \brief Function to set the sizes of the various class members
  void setSizes(const unsigned int& nSubjects,
                const unsigned int& nCovariates,
                const unsigned int& nDiscreteCov,
                const unsigned int& nContinuousCov,
                const vector<unsigned int>& nFixedEffects,
                const vector<unsigned int>& nFixedEffects_mix,
                const unsigned int& nCategoriesY,
                const unsigned int& nPredictSubjects,
                const unsigned int& nTimes,//RJ
                const unsigned int& nTimes_unique, //AR
                const unsigned int& nOutcomes, //AR
                const vector<unsigned int>& nCategories,
                const unsigned int& nClusInit,
                const string& covariateType,
                const string& outcomeType,
                const bool& weibullFixedShape,
                const string& kernelType,
                const vector<unsigned int>& nRandomEffects//AR
  ){

    unsigned int nDiscrCovs = 0;
    //unsigned int nOutcomes = nTimes/nSubjects;//RJ
    // Initially make the maximum number of clusters  the bigger or
    // the initial number of clusters and 150.
    // This is only used for initial memory allocation
    // And will ensure that at this initial time sufficient
    // space is allocated to prevent future allocations
    unsigned int maxNClusters = 100;
    if(nClusInit>100){
      maxNClusters=nClusInit;
    }
    _maxNClusters = maxNClusters;

    // Resize all the objects and set to 0
    _logPsi.resize(maxNClusters);
    for(unsigned int c=0;c<maxNClusters;c++){
      _logPsi[c]=0.0;
    }
    _v.resize(maxNClusters);
    _logPhi.resize(maxNClusters);
    _workLogPhiStar.resize(maxNClusters);
    _mu.resize(maxNClusters);
    _workMuStar.resize(maxNClusters);
    _Tau.resize(maxNClusters);
    _workSqrtTau.resize(maxNClusters);
    _workLogDetTau.resize(maxNClusters);
    _Sigma.resize(maxNClusters);
    _gamma.resize(maxNClusters);
    //RJ
    _MVNmu.resize(maxNClusters);
    _workMVNMuStar.resize(maxNClusters);
    _MVNTau.resize(maxNClusters);
    _workSqrtMVNTau.resize(maxNClusters);
    _workLogDetMVNTau.resize(maxNClusters);
    _MVNSigma.resize(maxNClusters);

    if (outcomeType.compare("LME")==0) {
      _RandomEffects.resize(nOutcomes);
      _workLogDetTauLME.resize(nOutcomes, maxNClusters);
      _workSqrtTauLME.resize(nOutcomes);
      _covRE.resize(nOutcomes);
      _SigmaE.resize(nOutcomes);
       for(unsigned int m=0;m<nOutcomes;m++){
        _SigmaE[m]=0.1;
        _RandomEffects[m].setZero(nSubjects,nRandomEffects[m]);
        _workSqrtTauLME[m].resize(maxNClusters);
        _covRE[m].resize(maxNClusters);

        for(unsigned int c=0;c<maxNClusters;c++){
          _covRE[m][c]=MatrixXd::Identity(nRandomEffects[m], nRandomEffects[m]);
          //LLT<MatrixXd> llt;
          //MatrixXd mat = 0.1*MatrixXd::Identity(nRandomEffects[m], nRandomEffects[m]);
          //_covRE[m][c] = mat;
          _workSqrtTauLME[m][c]=MatrixXd::Identity(nRandomEffects[m], nRandomEffects[m]);// = (llt.compute(mat)).matrixU();
        }
       }
    }

    for(unsigned int c=0;c<maxNClusters;c++){
      if(c==0){
        if (covariateType.compare("Discrete")==0) _logNullPhi.resize(nCovariates);
        if (covariateType.compare("Normal")==0) _nullMu.setZero(nCovariates);
        if (covariateType.compare("Mixed")==0){
          _logNullPhi.resize(nDiscreteCov);
          _nullMu.setZero(nContinuousCov);
        }
      }
      if (covariateType.compare("Discrete")==0) {
        _logPhi[c].resize(nCovariates);
        _workLogPhiStar[c].resize(nCovariates);
      } else if (covariateType.compare("Normal")==0) {
        _mu[c].setZero(nCovariates);
        _workMuStar[c].setZero(nCovariates);
        _Tau[c].setZero(nCovariates,nCovariates);
        _Sigma[c].setZero(nCovariates,nCovariates);
        _workSqrtTau[c].setZero(nCovariates,nCovariates);
      } else if (covariateType.compare("Mixed")==0) {
        _logPhi[c].resize(nDiscreteCov);
        _workLogPhiStar[c].resize(nDiscreteCov);
        _mu[c].setZero(nContinuousCov);
        _workMuStar[c].setZero(nContinuousCov);
        _Tau[c].setZero(nContinuousCov,nContinuousCov);
        _Sigma[c].setZero(nContinuousCov,nContinuousCov);
        _workSqrtTau[c].setZero(nContinuousCov,nContinuousCov);
      }
      //RJ
      if (outcomeType.compare("MVN")==0) {
        _MVNmu[c].setZero(nOutcomes);
        _workMVNMuStar[c].setZero(nOutcomes);
        _MVNTau[c].setZero(nOutcomes,nOutcomes);
        _MVNSigma[c].setZero(nOutcomes,nOutcomes);
        _workSqrtMVNTau[c].setZero(nOutcomes,nOutcomes);
      }
      _gamma[c].resize(nCovariates);
      if (covariateType.compare("Discrete")==0||covariateType.compare("Mixed")==0){
        if (covariateType.compare("Mixed")==0){
          nDiscrCovs = nDiscreteCov;
        } else {
          nDiscrCovs = nCovariates;
        }
        for(unsigned int j=0;j<nDiscrCovs;j++){
          if(c==0){
            _logNullPhi[j].resize(nCategories[j]);
          }
          _logPhi[c][j].resize(nCategories[j]);
          _workLogPhiStar[c][j].resize(nCategories[j]);
          for(unsigned int p=0;p<nCategories[j];p++){
            // We set logPhi to 0.0 so that we can call
            // normal member function above when we actually
            // initialise (allowing workLogPXiGivenZi to be
            // correctly calculated)
            _logPhi[c][j][p]=0.0;
            _workLogPhiStar[c][j][p]=0.0;
            if(c==0){
              _logNullPhi[j][p]=0.0;
            }
          }
        }
      }
      for(unsigned int j=0;j<nCovariates;j++){
        // Everything in by default
        // This allows us to write general case with switches.
        _gamma[c][j]=1;
      }
    }

    _theta.resize(maxNClusters);
    //RJ resize _L
    _L.resize(maxNClusters);
    _meanGP.resize(maxNClusters); //AR
    _ratio.resize(maxNClusters); //AR
    _pdf_meanGP.resize(maxNClusters); //AR
    for (unsigned int c=0;c<maxNClusters;c++){
      _theta[c].resize(nCategoriesY);
      if(kernelType.compare("SQexponential")==0){ //AR
        _L[c].resize(3);
      }else{
        _L[c].resize(4);
      }
      _meanGP[c].resize(nTimes_unique);//AR
    }
    _beta_mix.resize(nOutcomes);
    for (unsigned int m=0;m<nOutcomes;m++){
      _beta_mix[m].resize(maxNClusters, nCategoriesY*nFixedEffects_mix[m]);
      // for (unsigned int c=0;c<maxNClusters;c++){
      //   for (unsigned int j=0;j<nFixedEffects_mix[m];j++){
      //     for (unsigned int k=0;j<nCategoriesY;k++)
      //       _beta_mix[m](c,j*nCategoriesY+k)=0;
      //   }
      // }
    }

    _beta.resize(nOutcomes);
    for (unsigned int m=0;m<nOutcomes;m++){
      _beta[m].resize(nFixedEffects[m]*nCategoriesY);
      // for (unsigned int j=0;j<nFixedEffects[m];j++){
      //   for (unsigned int k=0;j<nCategoriesY;k++)
      //     _beta[m][j*nCategoriesY+k]=0;
      // }
    }

    _u.resize(nSubjects+nPredictSubjects,0.0);
    _lambda.resize(nSubjects);
    _z.resize(nSubjects+nPredictSubjects);
    _rho.resize(nCovariates);
    _omega.resize(nCovariates);
    _workNXInCluster.resize(maxNClusters,0);
    _workDiscreteX.resize(nSubjects+nPredictSubjects);
    _workContinuousX.resize(nSubjects+nPredictSubjects);
    _workLogPXiGivenZi.resize(nSubjects);
    _workPredictExpectedTheta.resize(nPredictSubjects);
    for (unsigned int i=0;i<nPredictSubjects;i++){
      _workPredictExpectedTheta[i].resize(nCategoriesY);
    }
    _workEntropy.resize(nSubjects+nPredictSubjects,0);
    for(unsigned int i=0;i<nSubjects+nPredictSubjects;i++){
      if (covariateType.compare("Discrete")==0||covariateType.compare("Normal")==0){
        _workDiscreteX[i].resize(nCovariates,0);
        _workContinuousX[i].resize(nCovariates,0);
      } else {
        _workDiscreteX[i].resize(nDiscreteCov,0);
        _workContinuousX[i].resize(nContinuousCov,0);
      }
      if(i<nSubjects){
        _workLogPXiGivenZi[i]=0;
      }
    }
    _uCAR.resize(nSubjects);
    if (weibullFixedShape) {
      _nu.resize(1);
    } else {
      _nu.resize(maxNClusters);
    }
  }



  /// \brief Return the number of clusters
  unsigned int maxNClusters() const{
    return _maxNClusters;
  }

  /// \brief Set the number of clusters
  void maxNClusters(const unsigned int& nClus,
                    const string& covariateType,
                    const string& outcomeType,
                    const string& kernelType,//AR
                    const unsigned int& nTimes_unique, //AR
                    const vector<unsigned int>& nRandomEffects //AR
  ){
    _maxNClusters=nClus;

    // Check if we need to do a resize of the
    // number of various vectors
    unsigned int prevNClus = _logPsi.size();
    if(nClus>prevNClus){
      unsigned int nCov=nCovariates();
      unsigned int nDiscrCovs=nDiscreteCovs();
      unsigned int nDCovs = 0;
      unsigned int nContCovs=nContinuousCovs();
      unsigned int nOcs = nOutcomes();
      vector<unsigned int> nCats=nCategories();
      unsigned int nCategoriesY = _theta[0].size();

      _logPsi.resize(nClus);
      _v.resize(nClus);
      _theta.resize(nClus);
      if (_nu.size()>1) _nu.resize(nClus);
      //RJ resize _L
      if (outcomeType.compare("Longitudinal")==0){
        _L.resize(nClus);
        _meanGP.resize(nClus);
        _ratio.resize(nClus);
        _pdf_meanGP.resize(nClus);
      }
      if (outcomeType.compare("LME")==0){
        _workLogDetTauLME.resize(nOcs, nClus);
        _workSqrtTauLME.resize(nClus);
        for (unsigned int m=0;m<nClus;m++){
          for (unsigned int c=0;c<nClus;c++){
            _workSqrtTauLME[m][c].resize(nRandomEffects[m],nRandomEffects[m]);
          }
        }

      }
      for (unsigned int c=0;c<nClus;c++){
        _theta[c].resize(nCategoriesY);
        if(outcomeType.compare("Longitudinal")==0 && kernelType.compare("SQexponential")==0){ //AR
          _L[c].resize(3);
        }else{
          _L[c].resize(4);
        }
        if(outcomeType.compare("Longitudinal")==0 && kernelType.compare("SQexponential")==0)
          _meanGP[c].resize(nTimes_unique);
      }
      _workNXInCluster.resize(nClus);
      if (covariateType.compare("Discrete")==0){
        _logPhi.resize(nClus);
        _workLogPhiStar.resize(nClus);
      } else if (covariateType.compare("Normal")==0){
        _mu.resize(nClus);
        _workMuStar.resize(nClus);
        _Tau.resize(nClus);
        _workSqrtTau.resize(nClus);
        _workLogDetTau.resize(nClus);
        _Sigma.resize(nClus);
      } else if (covariateType.compare("Mixed")==0){
        _logPhi.resize(nClus);
        _workLogPhiStar.resize(nClus);
        _mu.resize(nClus);
        _workMuStar.resize(nClus);
        _Tau.resize(nClus);
        _workSqrtTau.resize(nClus);
        _workLogDetTau.resize(nClus);
        _Sigma.resize(nClus);

      }
      if (outcomeType.compare("MVN")==0){//RJ
        _MVNmu.resize(nClus);
        _workMVNMuStar.resize(nClus);
        _MVNTau.resize(nClus);
        _workSqrtMVNTau.resize(nClus);
        _workLogDetMVNTau.resize(nClus);
        _MVNSigma.resize(nClus);
      }
      _gamma.resize(nClus);
      for(unsigned int c=prevNClus;c<nClus;c++){
        _workNXInCluster[c]=0;
        if (covariateType.compare("Discrete")==0){
          _logPhi[c].resize(nCov);
          _workLogPhiStar[c].resize(nCov);
        } else if (covariateType.compare("Normal")==0){
          _mu[c].setZero(nCov);
          _workMuStar[c].setZero(nCov);
          _Tau[c].setZero(nCov,nCov);
          _workSqrtTau[c].setZero(nCov,nCov);
          _Sigma[c].setZero(nCov,nCov);
        } else if (covariateType.compare("Mixed")==0){
          _logPhi[c].resize(nDiscrCovs);
          _workLogPhiStar[c].resize(nDiscrCovs);
          _mu[c].setZero(nContCovs);
          _workMuStar[c].setZero(nContCovs);
          _Tau[c].setZero(nContCovs,nContCovs);
          _workSqrtTau[c].setZero(nContCovs,nContCovs);
          _Sigma[c].setZero(nContCovs,nContCovs);
        }
        if (outcomeType.compare("MVN")==0){
          _MVNmu[c].setZero(nOcs);
          _workMVNMuStar[c].setZero(nOcs);
          _MVNTau[c].setZero(nOcs,nOcs);
          _workSqrtMVNTau[c].setZero(nOcs,nOcs);
          _MVNSigma[c].setZero(nOcs,nOcs);
        }
        _gamma[c].resize(nCov);
        if (covariateType.compare("Discrete")==0 || covariateType.compare("Mixed")==0){
          if (covariateType.compare("Mixed")==0){
            nDCovs = nDiscrCovs;
          } else {
            nDCovs = nCov;
          }
          for(unsigned int j=0;j<nDCovs;j++){
            _logPhi[c][j].resize(nCats[j]);
            _workLogPhiStar[c][j].resize(nCats[j]);
            for(unsigned int p=0;p<nCats[j];p++){
              // We set logPhi to 0.0 so that we can call
              // normal member function above when we actually
              // initialise (allowing workLogPXiGivenZi to be
              // correctly calculated)
              _logPhi[c][j][p]=0.0;
              _workLogPhiStar[c][j][p]=0.0;
            }
            // Everything in by default
            // This allows us to write general case with switches.
            _gamma[c][j]=1;
          }
        }
      }
    }
  }

  /// \brief Return the number of subjects
  unsigned int nSubjects() const{
    return _lambda.size();
  }
  unsigned int nOutcomes() const{
    return _MVNmu[0].size();
  }

  /// \brief Return the number of covariates
  unsigned int nCovariates() const{
    return _gamma[0].size();
  }

  /// \brief Return the number of covariates
  unsigned int nDiscreteCovs() const{
    return _logPhi[0].size();
  }

  /// \brief Return the number of covariates
  unsigned int nContinuousCovs() const{
    return _mu[0].size();
  }

  /// \brief Return the number of categories of outcome Y for Categorical outcome
  unsigned int nCategoriesY() const{
    return _theta[0].size();
  }

  /// \brief Return the number of clusters
  unsigned int nPredictSubjects() const{
    return _z.size()-_lambda.size();
  }

  /// \brief Return the vector of the number of categories
  vector<unsigned int> nCategories() const{
    vector<unsigned int> output;
    unsigned int nCovariates = _logPhi[0].size();
    for(unsigned int j=0;j<nCovariates;j++){
      output.push_back(_logPhi[0][j].size());
    }
    return output;
  }

  /// \brief Return the number of categories of covariate i
  unsigned int nCategories(const unsigned int& j) const{
    return _logPhi[0][j].size();
  }

  /// \brief Return the probabilities of allocation
  vector<double> logPsi() const{
    return _logPsi;
  }

  /// \brief Return the probability of allocation to cluster c
  double logPsi(const unsigned int& c) const{
    return _logPsi[c];
  }


  /// \brief Set the probability of allocation to cluster c
  void logPsi(const unsigned int& c,const double& logPsiVal){
    _logPsi[c]=logPsiVal;
  }

  /// \brief Set the vector of probability of allocations to clusters
  void logPsi(const vector<double>& logPsiVec){
    _logPsi=logPsiVec;
  }

  /// \brief Return the stick breaking constants
  vector<double> v() const{
    return _v;
  }

  /// \brief Return the stick breaking constant for cluster c
  double v(const unsigned int& c) const{
    return _v[c];
  }

  /// \brief Set the stick breaking constant for cluster c
  void v(const unsigned int& c,const double& vVal){
    _v[c]=vVal;
  }

  /// \brief Set the stick breaking constant for cluster c
  void v(const vector<double>& vVec){
    _v=vVec;
  }

  /// \brief Return the auxilliary variables
  vector<double> u() const{
    return _u;
  }

  /// \brief Return the auxilliary variable for individual i
  double u(const unsigned int& i) const{
    return _u[i];
  }

  /// \brief Set auxilliary variable for the individual i
  void u(const unsigned int& i,const double& uVal){
    _u[i]=uVal;
  }

  /// \brief Return the conditional covariate probabilites
  const vector<vector<vector<double> > >& logPhi() const{
    return _logPhi;
  }

  /// \brief Return the conditional probabilities for cluster c, covariate j
  const vector<double>& logPhi(const unsigned int& c,const unsigned int& j) const{
    return _logPhi[c][j];
  }

  /// \brief Set the conditional probabilities for cluster c, covariate j
  void logPhi(const unsigned int& c,const unsigned int& j,const vector<double>& logPhiVec){

    unsigned int nSbj = nSubjects();
    unsigned int nCat = nCategories(j);
    // Update the working variables
    vector<double> logPhiStarNew(nCat);
    for(unsigned int p=0;p<nCat;p++){
      logPhiStarNew[p] = log(gamma(c,j)*exp(logPhiVec[p])+(1.0-gamma(c,j))*exp(logNullPhi(j,p)));
    }

    for(unsigned int i=0;i<nSbj;i++){
      if(z(i)==(int)c){
        double logPhiStar;
        int Xij=workDiscreteX(i,j);
        logPhiStar = workLogPhiStar(c,j,Xij);
        _workLogPXiGivenZi[i]+=(logPhiStarNew[Xij]-logPhiStar);
      }
    }
    _workLogPhiStar[c][j]=logPhiStarNew;
    _logPhi[c][j]=logPhiVec;

  }

  /// \brief Return the conditional probability for cluster c, covariate j,category p
  double logPhi(const unsigned int& c,const unsigned int& j,const unsigned int& p) const{
    return _logPhi[c][j][p];
  }

  /// \brief Return the conditional covariate probabilites for null covariates
  const vector<vector<double> >& logNullPhi() const{
    return _logNullPhi;
  }

  /// \brief Return the conditional probabilities for null covariate j
  const vector<double>& logNullPhi(const unsigned int& j) const{
    return _logNullPhi[j];
  }

  /// \brief Set the conditional probabilities for cluster c, covariate j
  void logNullPhi(const unsigned int& j,const vector<double>& logNullPhiVec){

    unsigned int nClusters = maxNClusters();
    unsigned int nSbj = nSubjects();
    unsigned int nCat = nCategories(j);
    // Update the working variables
    vector<vector<double> > logPhiStarNew(nClusters);
    for(unsigned int c=0;c<nClusters;c++){
      logPhiStarNew[c].resize(nCat);
      for(unsigned int p=0;p<nCat;p++){
        logPhiStarNew[c][p] = log(gamma(c,j)*exp(logPhi(c,j,p))+(1.0-gamma(c,j))*exp(logNullPhiVec[p]));
      }
    }
    for(unsigned int i=0;i<nSbj;i++){
      unsigned int c=z(i);
      int Xij=workDiscreteX(i,j);
      double logPhiStar;
      logPhiStar = workLogPhiStar(c,j,Xij);
      _workLogPXiGivenZi[i]+=(logPhiStarNew[c][Xij]-logPhiStar);
    }
    for(unsigned int c=0;c<nClusters;c++){
      _workLogPhiStar[c][j]=logPhiStarNew[c];
    }
    _logNullPhi[j]=logNullPhiVec;

  }

  /// \brief Return the conditional probability for cluster c, covariate j,category p
  double logNullPhi(const unsigned int& j,const unsigned int& p) const{
    return _logNullPhi[j][p];
  }


  /// \brief Return the vector of normal means mu
  const vector<VectorXd>& mu() const{
    return _mu;
  }

  /// \brief Return the normal mean mu for cluster c
  const VectorXd& mu(const unsigned int& c) const{
    return _mu[c];
  }

  /// \brief Return the normal mean mu for cluster c covariate j
  double mu(const unsigned int& c, const unsigned int& j) const{
    return _mu[c](j);
  }

  /// \brief Set the normal mean for cluster c
  void mu(const unsigned int& c,const VectorXd& muVec){

    _mu[c]=muVec;

    unsigned int nSbj = nSubjects();
    unsigned int nCov = nCovariates();
    unsigned int nContCov = nContinuousCovs();
    if (nCov!=nContCov) {
      nCov = nContCov;
    }

    if(Sigma(c).trace()>0){
      // This condition should stop this being evaluated when
      // mu has been initialised but Sigma hasn't
      VectorXd xi=VectorXd::Zero(nCov);
      VectorXd muStar=VectorXd::Zero(nCov);
      for(unsigned int j=0;j<nCov;j++){
        muStar(j)=gamma(c,nDiscreteCovs()+j)*muVec(j)+(1.0-gamma(c,nDiscreteCovs()+j))*nullMu(j);
      }
      _workMuStar[c] = muStar;

      for(unsigned int i=0;i<nSbj;i++){
        if(z(i)==(int)c){
          for(unsigned int j=0;j<nCov;j++){
            xi(j)=workContinuousX(i,j);
          }
          _workLogPXiGivenZi[i]=logPdfMultivarNormal(nCov,xi,muStar,workSqrtTau(c),workLogDetTau(c));
        }
      }
    }
  }

  const vector<VectorXd>& MVNmu() const{
    return _MVNmu;
  }
  const VectorXd& MVNmu(const unsigned int& c) const{
    return _MVNmu[c];
  }
  double MVNmu(const unsigned int& c, const unsigned int& j) const{
    return _MVNmu[c](j);
  }
  void MVNmu(const unsigned int& c,const VectorXd& muVec){
    _MVNmu[c]=muVec;
  }

  /// \brief Return the normal mean mu for null covariates
  const VectorXd& nullMu() const{
    return _nullMu;
  }

  /// \brief Return the normal mean mu for null covariate j
  double nullMu(const unsigned int& j) const{
    return _nullMu(j);
  }

  /// \brief Set the normal mean for null covariates
  void nullMu(const VectorXd& nullMuVec){
    _nullMu=nullMuVec;

    unsigned int nClusters = maxNClusters();
    unsigned int nSbj = nSubjects();
    unsigned int nCov = nCovariates();
    unsigned int nContCov = nContinuousCovs();
    if (nCov!=nContCov) {
      nCov = nContCov;
    }

    // This condition should stop this being evaluated when
    // mu has been initialised but Sigma hasn't
    if(Sigma(0).trace()>0){
      VectorXd xi=VectorXd::Zero(nCov);
      vector<VectorXd> muStar(nClusters);
      for(unsigned int c=0;c<nClusters;c++){
        muStar[c].setZero(nCov);
        for(unsigned int j=0;j<nCov;j++){
          muStar[c](j)=gamma(c,nDiscreteCovs()+j)*mu(c,j)+(1-gamma(c,nDiscreteCovs()+j))*nullMuVec(j);
        }
        _workMuStar[c]=muStar[c];
      }
      for(unsigned int i=0;i<nSbj;i++){
        int c=z(i);
        for(unsigned int j=0;j<nCov;j++){
          xi(j)=workContinuousX(i,j);
          // Note in the continuous or global case, we just create
          // repeats of gamma(0,j) for each cluster
        }
        _workLogPXiGivenZi[i]=logPdfMultivarNormal(nCov,xi,muStar[c],workSqrtTau(c),workLogDetTau(c));
      }
    }
  }

  /// \brief Return the vector of precision matrices Tau
  const vector<MatrixXd>& Tau() const{
    return _Tau;
  }

  /// \brief Return the precision matrix Tau for cluster c
  const MatrixXd& Tau(const unsigned int& c) const{
    return _Tau[c];
  }

  /// \brief Set the precision matrix Tau for cluster c
  void Tau(const unsigned int& c,const MatrixXd& TauMat){

    _Tau[c]=TauMat;
    Sigma(c,TauMat.inverse());
    workLogDetTau(c,log(TauMat.determinant()));
    LLT<MatrixXd> llt;
    MatrixXd sqrtTau = (llt.compute(TauMat)).matrixU();
    workSqrtTau(c,sqrtTau);

    unsigned int nSbj = nSubjects();
    unsigned int nCov = nCovariates();
    unsigned int nContCov = nContinuousCovs();
    if (nCov!=nContCov) {
      nCov = nContCov;
    }

    VectorXd muStar=workMuStar(c);

    for(unsigned int i=0;i<nSbj;i++){
      VectorXd xi=VectorXd::Zero(nCov);
      if(z(i)==(int)c){
        for(unsigned int j=0;j<nCov;j++){
          xi(j)=workContinuousX(i,j);
        }
        _workLogPXiGivenZi[i]=logPdfMultivarNormal(nCov,xi,muStar,workSqrtTau(c),workLogDetTau(c));
      }
    }
  }

  /// \brief Return the covariance matrix Sigma element j1,j2 for cluster c
  double Tau(const unsigned int& c,const unsigned int& j1,const unsigned int& j2) const{
    return _Tau[c](j1,j2);
  }

  /// \brief Return the vector of covariance matrices Sigma
  const vector<MatrixXd>& Sigma() const{
    return _Sigma;
  }

  /// \brief Return the covariance matrix Sigma for cluster c
  const MatrixXd& Sigma(const unsigned int& c) const{
    return _Sigma[c];
  }

  /// \brief Set the covariance matrix Sigma for cluster c
  void Sigma(const unsigned int& c,const MatrixXd& SigmaMat){
    _Sigma[c]=SigmaMat;
  }

  /// \brief Return the covariance matrix Sigma element j1,j2 for cluster c
  double Sigma(const unsigned int& c,const unsigned int& j1,const unsigned int& j2) const{
    return _Sigma[c](j1,j2);
  }

  //RJ
  const vector<MatrixXd>& MVNTau() const{
    return _MVNTau;
  }
  const MatrixXd& MVNTau(const unsigned int& c) const{
    return _MVNTau[c];
  }
  void MVNTau(const unsigned int& c,const MatrixXd& TauMat){
    _MVNTau[c]=TauMat;
    MVNSigma(c,TauMat.inverse());
    workLogDetMVNTau(c,log(TauMat.determinant()));
    LLT<MatrixXd> llt;
    MatrixXd sqrtTau = (llt.compute(TauMat)).matrixU();
    workSqrtMVNTau(c,sqrtTau);
  }
  double MVNTau(const unsigned int& c,const unsigned int& j1,const unsigned int& j2) const{
    return _MVNTau[c](j1,j2);
  }
  const vector<MatrixXd>& MVNSigma() const{
    return _MVNSigma;
  }
  const MatrixXd& MVNSigma(const unsigned int& c) const{
    return _MVNSigma[c];
  }
  void MVNSigma(const unsigned int& c,const MatrixXd& SigmaMat){
    _MVNSigma[c]=SigmaMat;
  }
  double MVNSigma(const unsigned int& c,const unsigned int& j1,const unsigned int& j2) const{
    return _MVNSigma[c](j1,j2);
  }

  const MatrixXd& covRE(const unsigned int& m, const unsigned int& c) const{
    return _covRE[m][c];
  }

  const vector<vector<MatrixXd>>& covRE() const{
    return _covRE;
  }

  const vector<MatrixXd>& covRE(const unsigned int& m) const{
    return _covRE[m];
  }

  double covRE(const unsigned int& m, const unsigned int& c, const unsigned int& i, const unsigned int& j) const{
    return _covRE[m][c](i,j);
  }

  void covRE(const unsigned int& m, const unsigned int& c, const MatrixXd& cov) {
    _covRE[m][c]=cov;
    LLT<MatrixXd> llt;
    _workLogDetTauLME(m,c)=-log(cov.determinant());
    _workSqrtTauLME[m][c]=(llt.compute(cov.inverse())).matrixU();
  }


  const vector<MatrixXd>& RandomEffects() const{
    return _RandomEffects;
  }

  const MatrixXd& RandomEffects(const unsigned int& m) const{
    return _RandomEffects[m];
  }

  /// \brief Return the outcome probability for cluster c
  VectorXd RandomEffects(const unsigned int& m, const unsigned int& i) const{
    return _RandomEffects[m].row(i);
  }

  /// \brief Return the outcome probability for cluster c and category k
  double RandomEffects(const unsigned int& m, const unsigned int& i,const unsigned int& j) const{
    return _RandomEffects[m](i,j);
  }

  /// \brief Set the outcome probability for cluster c and category k
  void RandomEffects(const unsigned int& m, const unsigned int& i,const unsigned int& j,const double& thetaVal){
    _RandomEffects[m](i,j)=thetaVal;
  }


  void RandomEffects(const unsigned int& m, const unsigned int& i,const VectorXd& thetaVal){
    _RandomEffects[m].row(i)=thetaVal;
  }
  /// \brief Return the outcome probabilities
  vector<vector<double>> theta() const{
    return _theta;
  }

  /// \brief Return the outcome probability for cluster c
  vector <double> theta(const unsigned int& c) const{
    return _theta[c];
  }

  /// \brief Return the outcome probability for cluster c and category k
  double theta(const unsigned int& c,const unsigned int& k) const{
    return _theta[c][k];
  }

  /// \brief Set the outcome probability for cluster c and category k
  void theta(const unsigned int& c,const unsigned int& k,const double& thetaVal){
    _theta[c][k]=thetaVal;
  }


  //AR erre variance SigmaE for LME outcome option
  double  SigmaE(const unsigned int& m) const{
    return _SigmaE[m];
  }
  void SigmaE(const unsigned int& m, double& Val){
    _SigmaE[m]=Val;
  }

  //RJ handling functions for _L
  vector<vector <double> > L() const{
    return _L;
  }
  vector <double> L(const unsigned int& c) const{
    return _L[c];
  }
  double L(const unsigned int& c,const unsigned int& k) const{
    return _L[c][k];
  }
  void L(const unsigned int& c,const unsigned int& k,const double& LVal){
    _L[c][k]=LVal;
  }

  //AR handling functions for _meanGP
  vector<vector <double> > meanGP() const{
    return _meanGP;
  }
  vector <double> meanGP(const unsigned int& c) const{
    return _meanGP[c];
  }
  double meanGP(const unsigned int& c,const unsigned int& k) const{
    return _meanGP[c][k];
  }

  //RJ handling functions for _ratio
  vector <double> ratio() const{
    return _ratio;
  }
  double ratio(const unsigned int& c) const{
    return _ratio[c];
  }
  void ratio(const unsigned int& c,const double& LVal){
    _ratio[c]=LVal;
  }

  void meanGP(const unsigned int& c,const unsigned int& k,const double& LVal){
    _meanGP[c][k]=LVal;
  }

  vector <double> pdf_meanGP() const{
    return _pdf_meanGP;
  }

  double pdf_meanGP(const unsigned int& c) const{
    return _pdf_meanGP[c];
  }

  void pdf_meanGP(const unsigned int& c, const double& LVal){
    _pdf_meanGP[c]=LVal;
  }
  // void meanGP(const unsigned int& c,const unsigned int& k,const vector<vector <double> >& LVal){
  //  for(unsigned int j=0;j<LVal.size();j++){
  //    _meanGP[c][j]=LVal[j];
  //  }
  //}
  /// \brief Return the confounder coefficients
  // vector<vector <double>> beta() const{
  //   return _beta;
  // }

  /// \brief Return the coefficient for confounder j and category k, nY number of categories
  double beta(const unsigned int& m, const unsigned int& j, const unsigned int& k, const unsigned int& nY) const{
    return _beta[m][j*nY+k];
  }

  /// \brief Set the coefficient for confounder j and category k
  void beta(const unsigned int& m, const unsigned int& j, const unsigned int& k, const unsigned int& nY, const double& betaVal){
    _beta[m][j*nY+k]=betaVal;
  }

  // const vector<MatrixXd>& beta_mix() const{
  //   return _beta_mix;
  // }

  // cluster c, confounder j, k^th of the nY categories of Y
  double beta_mix(const unsigned int& m, const unsigned int& c, const unsigned int& j, const unsigned int& k, const unsigned int& nY) const{
    return _beta_mix[m](c,j*nY+k);
  }
  // Set the coefficient for cluster c confounder j and category k
  void beta_mix(const unsigned int& m, const unsigned int& c,const unsigned int& j, const unsigned int& k, const unsigned int& nY, const double& betaVal){
    _beta_mix[m](c,j*nY+k)=betaVal;
  }


  // void beta_mix(const unsigned int& k,const unsigned int& j,const double& betaVal){
  //   _beta_mix[k][j]=betaVal;
  // }

  /// \brief Return the hyper parameter alpha
  double alpha() const{
    return _alpha;
  }

  /// \brief Set the hyper parameter alpha
  void alpha(const double& alphaVal){
    _alpha=alphaVal;
  }

  /// \brief Return the hyper parameter dPitmanYor
  double dPitmanYor() const{
    return _dPitmanYor;
  }

  /// \brief Set the hyper parameter dPitmanYor
  void dPitmanYor(const double& dPitmanYorVal){
    _dPitmanYor=dPitmanYorVal;
  }


  /// \brief Return the mean y variables
  vector<double> lambda() const{
    return _lambda;
  }

  /// \brief Return the mean y variable of the ith subject
  double lambda(const unsigned int& i) const{
    return _lambda[i];
  }

  /// \brief Set the ith mean y variable
  void lambda(const unsigned int& i,const double& lam){
    _lambda[i]=lam;
  }

  /// \brief Return the hyper parameter epsilon
  double tauEpsilon() const{
    return _tauEpsilon;
  }

  /// \brief Set the hyper parameter epsilon
  void tauEpsilon(const double& tauVal){
    _tauEpsilon=tauVal;
  }

  /// \brief Return the allocation variables
  const vector<int>& z() const{
    return _z;
  }

  /// \brief Return the allocation variable of the ith subject
  int z(const unsigned int& i) const{
    return _z[i];
  }

  /// \brief Set the ith allocation variable to cluster c
  void z(const unsigned int& i,const int& c,const string& covariateType){
    unsigned int nCov = nCovariates();
    unsigned int nDiscreteCov = nDiscreteCovs();
    unsigned int nContinuousCov = nContinuousCovs();

    if(i<nSubjects()){
      if(covariateType.compare("Discrete")==0){
        for(unsigned int j=0;j<nCov;j++){
          unsigned int zi = z(i);
          int Xij=workDiscreteX(i,j);
          double logPhiStar = 0;
          double logPhiStarNew = 0;
          logPhiStar = workLogPhiStar(zi,j,Xij);
          logPhiStarNew = workLogPhiStar(c,j,Xij);
          _workLogPXiGivenZi[i]+=(logPhiStarNew-logPhiStar);
        }
      }else if(covariateType.compare("Normal")==0){
        if(Sigma(0).trace()>0){
          VectorXd xi=VectorXd::Zero(nCov);
          for(unsigned int j=0;j<nCov;j++){
            xi(j)=workContinuousX(i,j);
          }
          _workLogPXiGivenZi[i]=logPdfMultivarNormal(nCov,xi,workMuStar(c),workSqrtTau(c),workLogDetTau(c));
        }
      }else if(covariateType.compare("Mixed")==0){

        for(unsigned int j=0;j<nDiscreteCov;j++){
          unsigned int zi = z(i);
          int Xij=workDiscreteX(i,j);
          double logPhiStar,logPhiStarNew;
          logPhiStar = workLogPhiStar(zi,j,Xij);
          logPhiStarNew = workLogPhiStar(c,j,Xij);
          _workLogPXiGivenZi[i]+=(logPhiStarNew-logPhiStar);
        }
        if(Sigma(0).trace()>0){
          VectorXd xi=VectorXd::Zero(nContinuousCov);
          for(unsigned int j=0;j<nContinuousCov;j++){
            xi(j)=workContinuousX(i,j);
          }
          _workLogPXiGivenZi[i]+=logPdfMultivarNormal(nContinuousCov,xi,workMuStar(c),workSqrtTau(c),workLogDetTau(c));
        }
      }
    }
    // This is the only original code
    _z[i]=c;

  }

  /// \brief Return the variable selection matrix
  const vector<vector<double> >& gamma() const{
    return _gamma;
  }

  /// \brief Return the variable selection vector for cluster c
  const vector<double>& gamma(const unsigned int& c) const{
    return _gamma[c];
  }

  /// \brief Return the variable selection value for cluster c covariate j
  double gamma(const unsigned int& c, const unsigned int& j) const{
    return _gamma[c][j];
  }

  /// \brief Set the variable selection vector for all clusters (used for
  /// continuous variable selection)
  void gamma(const unsigned int& j,const double& gammaVal,const string& covariateType){

    unsigned int nClusters = maxNClusters();
    unsigned int nSbj = nSubjects();
    unsigned int nCov = nCovariates();
    unsigned int nDiscreteCov = nDiscreteCovs();
    unsigned int nContinuousCov = nContinuousCovs();

    // Need to update working vector here
    if(covariateType.compare("Discrete")==0){
      unsigned int nCat = nCategories(j);
      vector<vector<double> > logPhiStarNew(nClusters);
      for(unsigned int c=0;c<nClusters;c++){
        logPhiStarNew[c].resize(nCat);
        for(unsigned int p=0;p<nCat;p++){
          logPhiStarNew[c][p] = log(gammaVal*exp(logPhi(c,j,p))+(1.0-gammaVal)*exp(logNullPhi(j,p)));
        }
      }

      for(unsigned int i=0;i<nSbj;i++){
        unsigned int c=z(i);
        int Xij=_workDiscreteX[i][j];
        double logPhiStar;
        logPhiStar = workLogPhiStar(c,j,Xij);
        _workLogPXiGivenZi[i]+=(logPhiStarNew[c][Xij]-logPhiStar);
      }
      for(unsigned int c=0;c<nClusters;c++){
        _workLogPhiStar[c][j]=logPhiStarNew[c];
      }
    }else if(covariateType.compare("Normal")==0){
      if(Sigma(0).trace()>0){
        VectorXd xi=VectorXd::Zero(nCov);
        vector<VectorXd> muStar(nClusters);
        for(unsigned int c=0;c<nClusters;c++){
          muStar[c] = workMuStar(c);
          muStar[c](j)=gammaVal*mu(c,j)+(1-gammaVal)*nullMu(j);
        }
        for(unsigned int i=0;i<nSbj;i++){
          int c=z(i);
          for(unsigned int jj=0;jj<nCov;jj++){
            xi(jj)=workContinuousX(i,jj);
          }
          _workLogPXiGivenZi[i]=logPdfMultivarNormal(nCov,xi,muStar[c],workSqrtTau(c),workLogDetTau(c));

        }
        for(unsigned c=0;c<nClusters;c++){
          _workMuStar[c]=muStar[c];
        }
      }
    } else if (covariateType.compare("Mixed")==0){
      if (j < nDiscreteCov){
        unsigned int nCat = nCategories(j);
        vector<vector<double> > logPhiStarNew(nClusters);
        for(unsigned int c=0;c<nClusters;c++){
          logPhiStarNew[c].resize(nCat);
          for(unsigned int p=0;p<nCat;p++){
            logPhiStarNew[c][p] = log(gammaVal*exp(logPhi(c,j,p))+(1.0-gammaVal)*exp(logNullPhi(j,p)));
          }
        }

        for(unsigned int i=0;i<nSbj;i++){
          unsigned int c=z(i);
          int Xij=_workDiscreteX[i][j];
          double logPhiStar;
          logPhiStar = workLogPhiStar(c,j,Xij);
          _workLogPXiGivenZi[i]+=(logPhiStarNew[c][Xij]-logPhiStar);

        }
        for(unsigned int c=0;c<nClusters;c++){
          _workLogPhiStar[c][j]=logPhiStarNew[c];
        }
      } else {
        if(Sigma(0).trace()>0){
          VectorXd xi=VectorXd::Zero(nContinuousCov);
          vector<VectorXd> muStar(nClusters);
          for(unsigned int c=0;c<nClusters;c++){
            muStar[c] = workMuStar(c);
            muStar[c](j-nDiscreteCov)=gammaVal*mu(c,j-nDiscreteCov)+(1-gammaVal)*nullMu(j-nDiscreteCov);
          }
          for(unsigned int i=0;i<nSbj;i++){
            int c=z(i);
            for(unsigned int jj=0;jj<nContinuousCov;jj++){
              xi(jj)=workContinuousX(i,jj);
            }
            _workLogPXiGivenZi[i]=logPdfMultivarNormal(nContinuousCov,xi,muStar[c],workSqrtTau(c),workLogDetTau(c));
          }
          for(unsigned c=0;c<nClusters;c++){
            _workMuStar[c]=muStar[c];
          }
        }
      }
    }
    for(unsigned int c=0;c<nClusters;c++){
      _gamma[c][j]=gammaVal;
    }
  }

  /// \brief Set the variable selection vector for cluster c (used for
  /// binary cluster specific variable selection)
  void gamma(const unsigned int& c, const unsigned int& j, const double& gammaVal,const string& covariateType){

    unsigned int nSbj = nSubjects();
    unsigned int nCov = nCovariates();
    unsigned int nDiscrCov = nDiscreteCovs();
    unsigned int nContCov = nContinuousCovs();

    // Need to update working vector here
    if(covariateType.compare("Discrete")==0){
      unsigned int nCat = nCategories(j);
      vector<double> logPhiStarNew(nCat);
      for(unsigned int p=0;p<nCat;p++){
        logPhiStarNew[p] = log(gammaVal*exp(logPhi(c,j,p))+(1.0-gammaVal)*exp(logNullPhi(j,p)));
      }
      for(unsigned int i=0;i<nSbj;i++){
        if(z(i)==(int)c){
          int Xij=workDiscreteX(i,j);
          double logPhiStar;
          logPhiStar = workLogPhiStar(c,j,Xij);
          _workLogPXiGivenZi[i]+=(logPhiStarNew[Xij]-logPhiStar);
        }
      }
      _workLogPhiStar[c][j] = logPhiStarNew;

    }else if(covariateType.compare("Normal")==0){
      if(Sigma(0).trace()>0){
        VectorXd xi=VectorXd::Zero(nCov);
        VectorXd muStar = workMuStar(c);
        muStar(j)=gammaVal*mu(c,j)+(1-gammaVal)*nullMu(j);
        _workMuStar[c]=muStar;
        for(unsigned int i=0;i<nSbj;i++){
          if(z(i)==(int)c){
            for(unsigned int jj=0;jj<nCov;jj++){
              xi(jj)=workContinuousX(i,jj);
            }
            _workLogPXiGivenZi[i]=logPdfMultivarNormal(nCov,xi,muStar,workSqrtTau(c),workLogDetTau(c));

          }
        }
      }
    }else if(covariateType.compare("Mixed")==0){
      if (j<nDiscrCov){
        unsigned int nCat = nCategories(j);
        vector<double> logPhiStarNew(nCat);
        for(unsigned int p=0;p<nCat;p++){
          logPhiStarNew[p] = log(gammaVal*exp(logPhi(c,j,p))+(1.0-gammaVal)*exp(logNullPhi(j,p)));
        }
        for(unsigned int i=0;i<nSbj;i++){
          if(z(i)==(int)c){
            int Xij=workDiscreteX(i,j);
            double logPhiStar;
            logPhiStar = workLogPhiStar(c,j,Xij);
            _workLogPXiGivenZi[i]+=(logPhiStarNew[Xij]-logPhiStar);
           }
        }
        _workLogPhiStar[c][j] = logPhiStarNew;
      } else {
        if(Sigma(0).trace()>0){
          VectorXd xi=VectorXd::Zero(nContCov);
          VectorXd muStar = workMuStar(c);
          muStar(j-nDiscrCov)=gammaVal*mu(c,j-nDiscrCov)+(1-gammaVal)*nullMu(j-nDiscrCov);
          _workMuStar[c]=muStar;
          for(unsigned int i=0;i<nSbj;i++){
            if(z(i)==(int)c){
              for(unsigned int jj=0;jj<nContCov;jj++){
                xi(jj)=workContinuousX(i,jj);
              }
              _workLogPXiGivenZi[i]=logPdfMultivarNormal(nContCov,xi,muStar,workSqrtTau(c),workLogDetTau(c));
            }
          }
        }
      }
    }
    _gamma[c][j]=gammaVal;
  }

  /// \brief Get the variable selection rho vector - in the continuous
  /// case rho here is equivalent to zeta in the paper (and we copy
  /// this value straight to
  const vector<double>& rho() const{
    return _rho;
  }

  /// \brief Get the variable selection rho vector for cluster c (used for
  /// binary cluster specific variable selection)
  double rho(const unsigned int& j) const{
    return _rho[j];
  }

  /// \brief Set the variable selection rho vector for cluster c (used for
  /// binary cluster specific variable selection)
  void rho(const unsigned int& j,const double& rhoVal,const string& covariateType,const string& varSelectType){
    _rho[j]=rhoVal;
    // For continuous variable selection propagate this through to gamma
    if(varSelectType.compare("Continuous")==0){
      gamma(j,rhoVal,covariateType);
    }
  }

  /// \brief Get the variable selection rho vector for cluster c (used for
  /// binary cluster specific variable selection)
  const vector<unsigned int> omega() const{
    return _omega;
  }

  /// \brief Get the variable selection rho vector for cluster c (used for
  /// binary cluster specific variable selection)
  unsigned int omega(const unsigned int& j) const{
    return _omega[j];
  }

  /// \brief Set the variable selection rho vector for cluster c (used for
  /// binary cluster specific variable selection)
  void omega(const unsigned int& j,const unsigned int& omegaVal){
    _omega[j]=omegaVal;
  }


  /// \brief Return the parameter sigmaSqY
  double sigmaSqY() const{
    return _sigmaSqY;
  }

  /// \brief Set the parameter sigmaSqY
  void sigmaSqY(const double& sigmaSqYVal){
    _sigmaSqY=sigmaSqYVal;
  }

  /// \brief Return the vector _nu
  vector<double> nu() const{
    return _nu;
  }

  /// \brief Return the parameter _nu for cluster c
  double nu(const unsigned int& c ) const{
    return _nu[c];
  }

  /// \brief Set the vector nu
  void nu(vector<double>& nuVal){
    _nu=nuVal;
  }

  /// \brief Set the parameter nu for cluster c
  void nu(const unsigned int& c, double& nuVal){
    _nu[c]=nuVal;
  }


  /// \brief Return the vector _uCAR
  vector<double> uCAR() const{
    return _uCAR;
  }

  /// \brief Return the value of_uCAR for subject i
  double uCAR(const unsigned int& i ) const{
    return _uCAR[i];
  }

  /// \brief Set the CAR spatial random vector
  void uCAR(vector<double>& u){
    _uCAR=u;
  }

  /// \brief Set the CAR spatial random term for subject i
  void uCAR(const unsigned int& i, double& u){
    _uCAR[i]=u;
  }

  /// \brief Return the precision of spatial random term
  double TauCAR() const{
    return _TauCAR;
  }

  /// \brief Set the precision for the spatial random term
  void TauCAR(double& t){
    _TauCAR=t;
  }

  const pReMiuMHyperParams& hyperParams() const{
    return _hyperParams;
  }

  pReMiuMHyperParams& hyperParams(){
    return _hyperParams;
  }

  void hyperParams(const pReMiuMHyperParams& hyPar){
    _hyperParams = hyPar;
  }

  const vector<unsigned int>& workNXInCluster() const{
    return _workNXInCluster;
  }

  const unsigned int& workNXInCluster(unsigned int c) const{
    return _workNXInCluster[c];
  }

  void workNXInCluster(const vector<unsigned int>& nXInCluster){
    for(unsigned int c=0;c<_workNXInCluster.size();c++){
      if(c<nXInCluster.size()){
        _workNXInCluster[c]=nXInCluster[c];
      }else{
        _workNXInCluster[c]=0;
      }
    }
  }

  void workNXInCluster(const unsigned int& c,const unsigned int& n){
    _workNXInCluster[c]=n;
  }

  const unsigned int& workMaxZi() const{
    return _workMaxZi;
  }

  void workMaxZi(const unsigned int& maxZ){
    _workMaxZi=maxZ;
  }

  const double& workMinUi() const{
    return _workMinUi;
  }

  void workMinUi(const double& minUi){
    _workMinUi=minUi;
  }

  const vector<vector<int> >& workDiscreteX() const {
    return _workDiscreteX;
  }

  const vector<int>& workDiscreteX(const unsigned int& i) const{
    return _workDiscreteX[i];
  }

  int workDiscreteX(const unsigned int& i, const unsigned int& j) const {
    return _workDiscreteX[i][j];
  }

  void workDiscreteX(const vector<vector<int> >& X){
    _workDiscreteX=X;
  }

  void workDiscreteX(const unsigned int& i,const unsigned int& j,const int& x){
    _workDiscreteX[i][j]=x;
  }

  const vector<vector<double> >& workContinuousX() const {
    return _workContinuousX;
  }

  double workContinuousX(const unsigned int& i, const unsigned int& j) const {
    return _workContinuousX[i][j];
  }


  void workContinuousX(const vector<vector<double> >& X){
    _workContinuousX=X;
  }

  void workContinuousX(const unsigned int& i,const unsigned int& j,const double& x){
    _workContinuousX[i][j]=x;
  }

  /// \brief Return the conditional probabilities
  const vector<double>& workLogPXiGivenZi() const {
    return _workLogPXiGivenZi;
  }

  double workLogPXiGivenZi(const unsigned int& i) const{
    return _workLogPXiGivenZi[i];
  }

  void workLogPXiGivenZi(const unsigned int& i, const double& newVal){
    _workLogPXiGivenZi[i] = newVal;
  }

  const vector<vector<vector<double> > >& workLogPhiStar() const{
    return _workLogPhiStar;
  }

  double workLogPhiStar(const unsigned int& c,const unsigned int& j, const unsigned int& p) const{
    return _workLogPhiStar[c][j][p];
  }

  void workLogPhiStar(const unsigned int& c,const unsigned int& j, const unsigned int& p, const double& logPhiStar){
    _workLogPhiStar[c][j][p]=logPhiStar;
  }

  const VectorXd& workMuStar(const unsigned int& c) const{
    return _workMuStar[c];
  }

  const vector<VectorXd>& workMuStar() const{
    return _workMuStar;
  }

  void workMuStar(const unsigned int& c,const VectorXd& muStar){
    _workMuStar[c]=muStar;
  }

  //RJ
  const VectorXd& workMVNMuStar(const unsigned int& c) const{
    return _workMVNMuStar[c];
  }
  const vector<VectorXd>& workMVNMuStar() const{
    return _workMVNMuStar;
  }
  void workMVNMuStar(const unsigned int& c,const VectorXd& muStar){
    _workMVNMuStar[c]=muStar;
  }

  double workPredictExpectedTheta(const unsigned int& j,const unsigned int& k) const{
    return _workPredictExpectedTheta[j][k];
  }

  const vector<vector<double> >& workPredictExpectedTheta() const{
    return _workPredictExpectedTheta;
  }

  void workPredictExpectedTheta(const unsigned int& j,const unsigned int& k,const double& expectedVal){
    _workPredictExpectedTheta[j][k]=expectedVal;
  }

  double workEntropy(const unsigned int& i) const{
    return _workEntropy[i];
  }

  const vector<double>& workEntropy() const{
    return _workEntropy;
  }

  void workEntropy(const unsigned int& i,const double& entropyVal){
    _workEntropy[i]=entropyVal;
  }

  unsigned int workNClusInit() const{
    return _workNClusInit;
  }

  void workNClusInit(const unsigned int& nClusInit){
    _workNClusInit=nClusInit;
  }

  const vector<MatrixXd>& workSqrtTau() const{
    return _workSqrtTau;
  }


  const MatrixXd& workSqrtTau(const unsigned int& c) const{
    return _workSqrtTau[c];
  }

  void workSqrtTau(const unsigned int& c, const MatrixXd& sqrtTau){
    _workSqrtTau[c] = sqrtTau;
  }

  const MatrixXd& workSqrtTauLME(const unsigned int& m, const unsigned int& c) const{
    return _workSqrtTauLME[m][c];
  }

  void workSqrtTauLME(const unsigned int& m, const unsigned int& c, const MatrixXd& sqrtTau){
    _workSqrtTauLME[m][c] = sqrtTau;
  }

  const vector<double>& workLogDetTau() const{
    return _workLogDetTau;
  }

  double workLogDetTau(const unsigned int& c) const{
    return _workLogDetTau[c];
  }

  void workLogDetTau(const unsigned int& c, const double& logDetTau){
    _workLogDetTau[c] = logDetTau;
  }

  //RJ
  const vector<MatrixXd>& workSqrtMVNTau() const{
    return _workSqrtMVNTau;
  }
  const MatrixXd& workSqrtMVNTau(const unsigned int& c) const{
    return _workSqrtMVNTau[c];
  }
  void workSqrtMVNTau(const unsigned int& c, const MatrixXd& sqrtTau){
    _workSqrtMVNTau[c] = sqrtTau;
  }
  const vector<double>& workLogDetMVNTau() const{
    return _workLogDetMVNTau;
  }
  double workLogDetMVNTau(const unsigned int& c) const{
    return _workLogDetMVNTau[c];
  }
  void workLogDetMVNTau(const unsigned int& c, const double& logDetTau){
    _workLogDetMVNTau[c] = logDetTau;
  }

  const MatrixXd& workLogDetTauLME() const{
    return _workLogDetTauLME;
  }

  const double workLogDetTauLME(const unsigned int& m, const unsigned int& c) const{
    return _workLogDetTauLME(m,c);
  }

  double workLogDetTauLME(const unsigned int& m, const unsigned int& c){
    return _workLogDetTauLME(m,c);
  }

  void workLogDetTauLME(const unsigned int& m, const unsigned int& c, const double& logDetTau){
    _workLogDetTauLME(m,c) = logDetTau;
  }



  void switchLabels(const unsigned int& c1,const unsigned int& c2,
                    const string& covariateType,
                    const string& outcomeType,
                    const string& varSelectType,
                    const vector<unsigned int>& nFixedEffects,
                    const vector<unsigned int>& nFixedEffects_mix,
                    const unsigned int& nCategoriesY,
                    const unsigned int& nOutcomes){

    //Covariate parameters including working parameters
    if(covariateType.compare("Discrete")==0){
      _logPhi[c1].swap(_logPhi[c2]);
      _workLogPhiStar[c1].swap(_workLogPhiStar[c2]);
    }else if(covariateType.compare("Normal")==0){
      VectorXd muTmp = _mu[c1];
      _mu[c1]=_mu[c2];
      _mu[c2]=muTmp;
      VectorXd workMuStarTmp = _workMuStar[c1];
      _workMuStar[c1]=_workMuStar[c2];
      _workMuStar[c2]=workMuStarTmp;
      MatrixXd SigmaTmp = _Sigma[c1];
      _Sigma[c1]=_Sigma[c2];
      _Sigma[c2]=SigmaTmp;
      MatrixXd TauTmp = _Tau[c1];
      _Tau[c1]=_Tau[c2];
      _Tau[c2]=TauTmp;
      MatrixXd sqrtTauTmp = _workSqrtTau[c1];
      _workSqrtTau[c1]=_workSqrtTau[c2];
      _workSqrtTau[c2]=sqrtTauTmp;
      double logDetTauTmp = _workLogDetTau[c1];
      _workLogDetTau[c1]=_workLogDetTau[c2];
      _workLogDetTau[c2]=logDetTauTmp;
    }else if(covariateType.compare("Mixed")==0){
      _logPhi[c1].swap(_logPhi[c2]);
      _workLogPhiStar[c1].swap(_workLogPhiStar[c2]);
      VectorXd muTmp = _mu[c1];
      _mu[c1]=_mu[c2];
      _mu[c2]=muTmp;
      VectorXd workMuStarTmp = _workMuStar[c1];
      _workMuStar[c1]=_workMuStar[c2];
      _workMuStar[c2]=workMuStarTmp;
      MatrixXd SigmaTmp = _Sigma[c1];
      _Sigma[c1]=_Sigma[c2];
      _Sigma[c2]=SigmaTmp;
      MatrixXd TauTmp = _Tau[c1];
      _Tau[c1]=_Tau[c2];
      _Tau[c2]=TauTmp;
      MatrixXd sqrtTauTmp = _workSqrtTau[c1];
      _workSqrtTau[c1]=_workSqrtTau[c2];
      _workSqrtTau[c2]=sqrtTauTmp;
      double logDetTauTmp = _workLogDetTau[c1];
      _workLogDetTau[c1]=_workLogDetTau[c2];
      _workLogDetTau[c2]=logDetTauTmp;
    }

    //Variable selection parameters
    if(varSelectType.compare("BinaryCluster")==0){
      _gamma[c1].swap(_gamma[c2]);
    }
    //Response parameters
    vector<double> thetaTmp = _theta[c1];
    _theta[c1]=_theta[c2];
    _theta[c2]=thetaTmp;

    //RJ switch outcome parameters
    if(outcomeType.compare("MVN")==0){
      VectorXd muTmp = _MVNmu[c1];
      _MVNmu[c1]=_MVNmu[c2];
      _MVNmu[c2]=muTmp;
      MatrixXd SigmaTmp = _MVNSigma[c1];
      _MVNSigma[c1]=_MVNSigma[c2];
      _MVNSigma[c2]=SigmaTmp;
      MatrixXd TauTmp = _MVNTau[c1];
      _MVNTau[c1]=_MVNTau[c2];
      _MVNTau[c2]=TauTmp;
      MatrixXd sqrtTauTmp = _workSqrtMVNTau[c1];
      _workSqrtMVNTau[c1]=_workSqrtMVNTau[c2];
      _workSqrtMVNTau[c2]=sqrtTauTmp;
      double logDetTauTmp = _workLogDetMVNTau[c1];
      _workLogDetMVNTau[c1]=_workLogDetMVNTau[c2];
      _workLogDetMVNTau[c2]=logDetTauTmp;
    }

    if(outcomeType.compare("LME")==0){
      for(unsigned int m=0;m<nOutcomes;m++){
        // MatrixXd CovRETmp = _covRE[m][c1];
        // _covRE[m][c1]=_covRE[m][c2];
        // _covRE[m][c2]=CovRETmp;
        double btemp=0;
        for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
          btemp=_beta_mix[m](c1,b);//[m](c,j*nY+k)
          _beta_mix[m](c1,b)=_beta_mix[m](c2,b);
          _beta_mix[m](c2,b)=btemp;
        }

        // MatrixXd sqrtTauTmp = _workSqrtTauLME[m][c1];
        // _workSqrtTauLME[m][c1]=_workSqrtTauLME[m][c2];
        // _workSqrtTauLME[m][c2]=sqrtTauTmp;


        // double logDetTauTmp = _workLogDetTauLME(m,c1);
        // _workLogDetTauLME(m,c1)=_workLogDetTauLME(m,c2);
        // _workLogDetTauLME(m,c2)=logDetTauTmp;
      }
    }

    if(outcomeType.compare("Longitudinal")==0){
      vector<double> LTmp = _L[c1];
      _L[c1]=_L[c2];
      _L[c2]=LTmp;

      vector<double> LTmp2 = _meanGP[c1];
      _meanGP[c1]=_meanGP[c2];
      _meanGP[c2]=LTmp2;

      double LTmp3 = _pdf_meanGP[c1];
      _pdf_meanGP[c1]=_pdf_meanGP[c2];
      _pdf_meanGP[c2]=LTmp3;

    }

    //Allocation parameters (including counts)
    for(unsigned int i=0;i<nSubjects()+nPredictSubjects();i++){
      if(_z[i]==(int)c1){
        _z[i]=c2;
      }else if(_z[i]==(int)c2){
        _z[i]=c1;
      }
    }
    double workNXInClusterTmp = _workNXInCluster[c1];
    _workNXInCluster[c1]=_workNXInCluster[c2];
    _workNXInCluster[c2]=workNXInClusterTmp;
  }

  /// \brief Copy operator
  /*
  pReMiuMParams& operator=(const pReMiuMParams& params){
    _maxNClusters=params.maxNClusters();
    _logPsi = params.logPsi();
    _u = params.u();
    _v = params.v();
    _logPhi = params.logPhi();
    _logNullPhi = params.logNullPhi();
    _mu = params.mu();
    _nullMu = params.nullMu();
    _Tau = params.Tau();
    _Sigma = params.Sigma();
    _theta = params.theta();
    _beta = params.beta();
    _beta_mix = params.beta_mix();
    //RJ set _L and MVN parameters
    _L = params.L();
    _SigmaE = params.SigmaE();
    _ratio = params.ratio();
    _meanGP = params.meanGP(); //AR
    _pdf_meanGP = params.pdf_meanGP(); //AR
    _MVNmu = params.MVNmu();
    _MVNTau = params.MVNTau();
    _MVNSigma = params.MVNSigma();
    _alpha = params.alpha();
    _dPitmanYor = params.dPitmanYor();
    _lambda = params.lambda();
    _tauEpsilon = params.tauEpsilon();
    _z = params.z();
    _gamma = params.gamma();
    _rho = params.rho();
    _omega = params.omega();
    _sigmaSqY = params.sigmaSqY();
    _nu = params.nu();
    _hyperParams = params.hyperParams();
    _workNXInCluster=params.workNXInCluster();
    _workMaxZi=params.workMaxZi();
    _workMinUi=params.workMinUi();
    _workDiscreteX=params.workDiscreteX();
    _workContinuousX=params.workContinuousX();
    _workLogPXiGivenZi = params.workLogPXiGivenZi();
    _workLogPhiStar = params.workLogPhiStar();
    _workMuStar = params.workMuStar();
    _workPredictExpectedTheta = params.workPredictExpectedTheta();
    _workEntropy = params.workEntropy();
    _workNClusInit = params.workNClusInit();
    _workLogDetTau = params.workLogDetTau();
    _workSqrtTau = params.workSqrtTau();
    _uCAR = params.uCAR();
    _TauCAR = params.TauCAR();
    _workMVNMuStar = params.workMVNMuStar();
    _workLogDetMVNTau = params.workLogDetMVNTau();
    _workSqrtMVNTau = params.workSqrtMVNTau();
    _RandomEffects = params.RandomEffects();
    _workLogDetTauLME = params.workLogDetTauLME();
    _covRE = params.covRE();
    return *this;
  }
*/


private:
  /// \brief The current maximum number of clusters
  unsigned int _maxNClusters;

  /// \brief Vector of probabilities of assignment to clusters
  vector<double> _logPsi;

  /// \brief Vector of constants for stick breaking generation of logPsi
  vector<double> _v;

  /// \brief Vector of auxilliary variables for slice sampler
  vector<double> _u;

  /// \brief An array containing conditional covariate probabilities
  vector<vector<vector<double> > > _logPhi;

  /// \brief A matrix containing conditional covariate probabilities
  /// for the "null" covariates
  vector<vector<double> > _logNullPhi;

  /// \brief A vector of Eigen dynamic vectors containing covariate means
  /// for the case of Normal covariates
  vector<VectorXd> _mu;
  vector<VectorXd> _MVNmu;

  /// \brief An Eigen dynamic vector containing covariate means
  /// for the case of Normal covariates for the null covariates
  VectorXd _nullMu;

  /// \brief A vector of Eigen dynamic vectors containing covariate covariance
  /// matrices for the case of Normal covariates
  vector<MatrixXd> _Sigma;
  vector<MatrixXd> _MVNSigma;

  vector<MatrixXd> _RandomEffects;
  vector<vector<MatrixXd>> _covRE;

  /// \brief A vector of Eigen dynamic vectors containing covariate precision
  /// matrices for the case of Normal covariates
  vector<MatrixXd> _Tau;
  vector<MatrixXd> _MVNTau;

  /// \brief A vector of outcome probabilities for each cluster
  vector< vector <double> > _theta;

  /// \brief A vector of coefficients for confounding covariates
  vector<vector <double>> _beta; //m ; nCategoriesY*nFixedEffects
  vector<MatrixXd> _beta_mix; //m outcomes, MatrixXd _beta_mix[c].resize(nCategoriesY*nFixedEffects_mix)

  //RJ declare _L
  vector< vector<double> > _L;
  vector<double> _SigmaE;
  vector<double>  _ratio;
  //AR declare _meanGP
  vector< vector<double> > _meanGP;

  //AR declare _pdf_meanGP
  vector<double>  _pdf_meanGP;

  /// \brief The hyper parameter for dirichlet model
  double _alpha;

  /// \brief The discount hyper parameter for the pitman-yor process prior
  double _dPitmanYor;

  /// \brief The mean Y values (if required)
  vector<double> _lambda;

  /// \brief The precision for the extra deviation (if required)
  double _tauEpsilon;

  /// \brief A vector of allocation variables
  vector<int> _z;

  /// \brief A matrix containing the binary variable selection switches
  /// NOTE: in the continuous or global case, we just create
  /// repeats of gamma(0,j) for each cluster to make computation faster
  /// due to less function calls. In our notation, in the continuous case
  /// gamma is equal (deterministically) to rho.
  vector<vector<double> > _gamma;

  /// \brief Prior parameters for variable selection (binary and continuous)
  /// NOTE: in the continous case, rho is equivalent to zeta in the paper
  /// and we juect copy this value into gamma.
  vector<double> _rho;

  /// \brief Prior parameters for rho
  vector<unsigned int> _omega;

  /// \brief Prior variance for Y model when response is Normal
  double _sigmaSqY;

  /// \brief Prior shape parameter for Y model when response is weibull (survival)
  vector<double> _nu;

  /// NOTE: hyper parameters
  /// \brief The hyper parameter object
  pReMiuMHyperParams _hyperParams;

  /// \brief Integer determining the maximum cluster label with some members
  /// in
  unsigned int _workMaxZi;

  /// \brief Double determining the minimum ui
  double _workMinUi;

  /// \brief A vector containing the number of subjects in each cluster
  vector<unsigned int> _workNXInCluster;

  /// \brief A matrix containing a copy of X
  vector<vector<int> > _workDiscreteX;

  /// \brief A matrix containing a copy of X
  vector<vector<double> > _workContinuousX;

  /// \brief A matrix containing P(Xi|Zi)
  vector<double>  _workLogPXiGivenZi;

  /// \brief An array of phi star for variable selection
  vector<vector<vector<double> > > _workLogPhiStar;

  /// \brief An array of mu star for variable selection
  vector<VectorXd> _workMuStar;
  vector<VectorXd> _workMVNMuStar;

  /// \brief A vector of expected theta values for the prediction subjects
  vector<vector<double> > _workPredictExpectedTheta;

  /// \brief A vector of entropy values for each of the subjects
  vector<double> _workEntropy;

  /// \brief The actual number of cluster the individuals are individually
  /// allocated to. Note this is just needed for writing the log file.
  unsigned int _workNClusInit;

  /// \brief Working vector of matrices containing matrix square root of Tau
  vector<MatrixXd> _workSqrtTau;
  vector<MatrixXd> _workSqrtMVNTau;
  vector<vector<MatrixXd>> _workSqrtTauLME;
  /// \brief Working vector containing the log determinants of Tau
  vector<double> _workLogDetTau;
  vector<double> _workLogDetMVNTau;
  MatrixXd _workLogDetTauLME;
  /// \brief vector of CAR spatial random term U
  vector<double> _uCAR;

  /// \brief double for the precision of the CAR random term _uCAR
  double _TauCAR;
};


double logPYiGivenZiWiBernoulli(const pReMiuMParams& params, const pReMiuMData& dataset,
                                const vector<unsigned int>& nFixedEffects,const int& zi,
                                const unsigned int& i){

  double lambda;
  unsigned int nCategoriesY = params.nCategoriesY();
  lambda=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    lambda+=params.beta(0, j,0, nCategoriesY)*dataset.W(i,j);
  }

  double p=1.0/(1.0+exp(-lambda));
  return logPdfBernoulli(dataset.discreteY(i),p);
}

double logPYiGivenZiWiBernoulliExtraVar(const pReMiuMParams& params,
                                        const pReMiuMData& dataset,
                                        const vector<unsigned int>& nFixedEffects,const int& zi,
                                        const unsigned int& i){

  double lambda=params.lambda(i);

  double p=1.0/(1.0+exp(-lambda));
  return logPdfBernoulli(dataset.discreteY(i),p);
}

double logPYiGivenZiWiBinomial(const pReMiuMParams& params, const pReMiuMData& dataset,
                               const vector<unsigned int>& nFixedEffects,const int& zi,
                               const unsigned int& i){

  double lambda;
  unsigned int nCategoriesY = params.nCategoriesY();
  lambda=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    lambda+=params.beta(0, j,0, nCategoriesY)*dataset.W(i,j);
  }

  double p=1.0/(1.0+exp(-lambda));
  return logPdfBinomial(dataset.discreteY(i),dataset.nTrials(i),p);
}

double logPYiGivenZiWiBinomialExtraVar(const pReMiuMParams& params,
                                       const pReMiuMData& dataset,
                                       const vector<unsigned int>& nFixedEffects,const int& zi,
                                       const unsigned int& i){

  double lambda=params.lambda(i);

  double p=1.0/(1.0+exp(-lambda));
  return logPdfBinomial(dataset.discreteY(i),dataset.nTrials(i),p);
}

double logPYiGivenZiWiPoisson(const pReMiuMParams& params, const pReMiuMData& dataset,
                              const vector<unsigned int>& nFixedEffects,const int& zi,
                              const unsigned int& i){

  double lambda;
  unsigned int nCategoriesY = params.nCategoriesY();
  lambda=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    lambda+=params.beta(0, j,0, nCategoriesY)*dataset.W(i,j);
  }
  lambda+=dataset.logOffset(i);

  double mu =exp(lambda);
  return logPdfPoisson(dataset.discreteY(i),mu);
}

double logPYiGivenZiWiPoissonExtraVar(const pReMiuMParams& params,
                                      const pReMiuMData& dataset,
                                      const vector<unsigned int>& nFixedEffects,const int& zi,
                                      const unsigned int& i){

  double lambda=params.lambda(i);

  double mu=exp(lambda);
  return logPdfPoisson(dataset.discreteY(i),mu);
}

double logPYiGivenZiWiPoissonSpatial(const pReMiuMParams& params, const pReMiuMData& dataset,
                                     const vector<unsigned int>& nFixedEffects,const int& zi,
                                     const unsigned int& i){

  double lambda;
  unsigned int nCategoriesY = params.nCategoriesY();
  lambda=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    lambda+=params.beta(0, j,0, nCategoriesY)*dataset.W(i,j);
  }
  lambda+=dataset.logOffset(i);
  lambda+=params.uCAR(i);

  double mu =exp(lambda);
  return logPdfPoisson(dataset.discreteY(i),mu);
}

double logPYiGivenZiWiLongitudinal(const pReMiuMParams& params, const pReMiuMData& dataset,
                                   const vector<unsigned int>& nFixedEffects, const int& c,
                                   const unsigned int& ii){

  //vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nSubjects=dataset.nSubjects();
  const string kernelType = dataset.kernelType(); //AR
  vector<double> y = dataset.continuousY();
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> timesk;
  unsigned int nCategoriesY = params.nCategoriesY();

  VectorXd yk;
  VectorXd meanVec;
  MatrixXd precMat;
  double logDetMat =0.0;
  int sizek = 0;
  double dmvnorm = 0.0;
  int counter = 0;
  // set sizes based on cluster occupation
  for(unsigned int i=0;i<nSubjects;i++){
    if(params.z(i) == c && ii!=i){
      sizek = sizek + tStop[i] - tStart[i] + 1;
    }
  }

  timesk.resize(sizek);
  yk.resize(sizek);
  meanVec.resize(sizek);
  for(unsigned int i=0;i<nSubjects;i++){
    if(params.z(i) == c && ii!=i){
      for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
        yk(counter+j) = y[tStart[i]-1+j];
        timesk[counter+j] = times[tStart[i]-1+j];
        meanVec(counter+j) = 0.0;
        for(unsigned int b=0;b<nFixedEffects[0];b++){
          yk(counter+j)-=params.beta(0, b,0, nCategoriesY)*dataset.W(i,b);
        }
        // for(unsigned int b=0;b<nFixedEffects_mix;b++){
        //   yk(counter+j)-=params.beta_mix(c,b,0,nCategoriesY)*dataset.W_mix(i,b);
        // }
      }
      //RJ tidy up vector copying
      //timesk.insert(timesk.begin()+counter,times.begin()+tStart[i]-1,times.begin()+tStop[i]);
      counter = counter + tStop[i] - tStart[i] + 1;
    }
  }

  precMat.setZero(sizek,sizek);

  if(sizek > 0){
    int dimBlock = dataset.equalTimes();
    if(dimBlock<0){
      }else{

      // sort yk and times with corresponding ordering in idx
      std::vector<int> idx(timesk.size());
      int x=0;
      iota(idx.begin(), idx.end(), x++);
      stable_sort(idx.begin(), idx.end(), // sort indexes based on comparing values in times
                  [&](int i1,int i2) { return (timesk[i1] < timesk[i2]); });

      stable_sort(timesk.begin(), timesk.end());

      VectorXd yk_sorted = yk;
      for(unsigned int j=0;j<sizek;j++)
        yk_sorted(j)=yk(idx[j]);
      yk=yk_sorted;

      logDetMat=Get_Sigma_inv_GP_cov(precMat,params.L(c),timesk,dataset.equalTimes(),dataset.times_unique(),kernelType);
      dmvnorm = -0.5*yk.transpose()*precMat*yk - 0.5*sizek*log(2.0*pi<double>()) - 0.5*logDetMat;
    }
  }
  return dmvnorm;

  if(std::isnan(dmvnorm) || isinf(dmvnorm)){

    std::cout << " dmvnorm  "<< dmvnorm<< " logDetMat "<<logDetMat << endl<<
      " zi "<<  params.z(ii) << " c "<< c <<endl<<
        " prod "<< yk.transpose()*precMat*yk <<
          " sizek "<<sizek<<endl;
    dmvnorm=-(std::numeric_limits<double>::max());
  }
}

//AR Compute logPYiGivenZiWiLongitudinal from covariance matrix and det
double logPYiGivenZiWiLongitudinal_bis(const MatrixXd& Sigma_inv_ord, const double& logDetSigma, const pReMiuMParams& params, const pReMiuMData& dataset,
                                       const vector<unsigned int>& nFixedEffects, const unsigned int& c, unsigned int size_k,
                                       const  unsigned int& ii// subject to remove
){

  //vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nSubjects=dataset.nSubjects();
  vector<double> y = dataset.continuousY();
  const string kernelType = dataset.kernelType(); //AR
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> timesk ;
  unsigned int nCategoriesY = params.nCategoriesY();
  //int nTimes_unique = dataset.nTimes_unique();

  VectorXd yk;
  //VectorXd meanVec;
  MatrixXd precMat;
  double logDetPrecMat =0.0;

  int sizek = 0;
  double dmvnorm = 0.0;
  int counter = 0;
  unsigned int ind_permut;

  sizek=size_k;

  // set sizes based on cluster occupation
  if(ii < nSubjects){
    if(params.z(ii) == c){ // AR remove ii
      sizek = sizek - (tStop[ii] - tStart[ii] + 1); // +1 if add=1 params.z(ii) != c, -1 if add=0 params.z(ii) = c
    }else{                 // AR add ii
      sizek = sizek + tStop[ii] - tStart[ii] + 1; // +1 if add=1 params.z(ii) != c, -1 if add=0 params.z(ii) = c
    }
  }

  if(sizek > 0){
    yk.resize(sizek);
    timesk.resize(sizek);

    for(unsigned int i=0;i<nSubjects;i++){
      if(params.z(i) == c && ii!=i){
        for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
          yk(counter+j) = y[tStart[i]-1+j];
          timesk[counter+j] = times[tStart[i]-1+j];

          for(unsigned int b=0;b<nFixedEffects[0];b++){
            yk(counter+j)-=params.beta(0, b,0, nCategoriesY)*dataset.W(i,b);
          }
          // for(unsigned int b=0;b<nFixedEffects_mix;b++){
          //   yk(counter+j)-=params.beta_mix(c,b, 0, nCategoriesY)*dataset.W_mix(i,b);
          // }
        }
        counter = counter + tStop[i] - tStart[i] + 1;
      }else if((params.z(i) == c) & (ii==i)){
        ind_permut=counter; //index (in times and y) of first visit of subject to permut in M0
      }
    }

    if(ii>=nSubjects){
      precMat = Sigma_inv_ord; //ascending order of times
      logDetPrecMat = logDetSigma;

      // sort yk with corresponding ordering
      std::vector<int> idx(timesk.size());
      int x=0;
      iota(idx.begin(), idx.end(), x++);
      stable_sort(idx.begin(), idx.end(), // sort indexes based on comparing values in times
                  [&](int i1,int i2) { return (timesk[i1] < timesk[i2]); });

      stable_sort(timesk.begin(), timesk.end());

      VectorXd yk_sorted = yk;
      for(unsigned int j=0;j<sizek;j++)
        yk_sorted(j)=yk(idx[j]);
      yk=yk_sorted;

    }else{

      if(params.z(ii) != c ){ // Add one subject

        std::vector<double> times_check=timesk; // ordered with indiv to add at the end
        stable_sort(times_check.begin(), times_check.begin()+counter);

        for(unsigned int j=0;j<tStop[ii]-tStart[ii]+1;j++){
          yk(counter+j) = y[tStart[ii]-1+j];
          timesk[counter+j] = times[tStart[ii]-1+j];
          times_check[counter+j] = times[tStart[ii]-1+j];//check

          for(unsigned int b=0;b<nFixedEffects[0];b++){
            yk(counter+j)-=params.beta(0, b,0, nCategoriesY)*dataset.W(ii,b);
          }
          // for(unsigned int b=0;b<nFixedEffects_mix;b++){
          //   yk(counter+j)-=params.beta_mix(c,b,0,nCategoriesY)*dataset.W_mix(ii,b);
          // }
        }

        precMat.setZero(sizek,sizek);
        GP_cov(precMat,params.L(c),times_check,dataset.equalTimes(),kernelType,1);

        if(Sigma_inv_ord.rows()>0 && sizek>100){//sizek>nTimes_unique

          std::vector<int> idx_back(sizek);
          int x=0;
          iota(idx_back.begin(), idx_back.end(), x++);
          stable_sort(idx_back.begin(), idx_back.end(), // sort indexes based on comparing values in times
                      [&](int i1,int i2) { return (times_check[i1] < times_check[i2]); });

          logDetPrecMat=Inverse_woodbury(Sigma_inv_ord,logDetSigma, precMat, idx_back);
          if(std::isnan(logDetPrecMat) || isinf(logDetPrecMat) || logDetPrecMat== -(std::numeric_limits<double>::max())){
            //fout << "indiv: "<< ii << " z(i) "<<params.z(ii)<< " c "<< c<<endl;
          }
          //fout << "Inverse_woodbury add "<<logDetSigma << endl;
          // Sigma2 and Sigma_inv_ord2 including measurement errors
          // ordered, with subject to remove at the end

        }else{
          LLT<MatrixXd> lltOfA(precMat); // compute the Cholesky decomposition of A
          MatrixXd L = lltOfA.matrixL();
          logDetPrecMat=  2*L.diagonal().array().log().sum();
          precMat = L.inverse().transpose()*L.inverse();
        }

        // sort yk with corresponding ordering
        std::vector<int> idx(timesk.size());
        int x=0;
        iota(idx.begin(), idx.end(), x++);
        stable_sort(idx.begin(), idx.end(), // sort indexes based on comparing values in times
                    [&](int i1,int i2) { return (timesk[i1] < timesk[i2]); });

        sort(timesk.begin(), timesk.end());

        VectorXd yk_sorted = yk;
        for(unsigned int j=0;j<sizek;j++)
          yk_sorted(j)=yk(idx[j]);
        yk=yk_sorted;

      }else{ // Remove one subject
        precMat.setZero(sizek,sizek);

        if(Sigma_inv_ord.rows()>0 && sizek>100){ //sizek>nTimes_unique

          vector<double> timesk2; //timesk (t_i1, ///, t_in) + indiv to remove at the end
          timesk2.resize(Sigma_inv_ord.rows()); // sizek+tStop[i]-tStart[i]+1 + indiv in end
          for(unsigned int j=0;j<timesk.size();j++)
            timesk2[j]= timesk[j];

          // Remove one subject
          for(unsigned int j=0;j<tStop[ii]-tStart[ii]+1;j++){
            timesk2[counter+j] = times[tStart[ii]-1+j];
          }

          vector<double> timesk3; //times_init (timkes+indiv) ordered
          timesk3.resize(sizek + tStop[ii]-tStart[ii]+1); // sizek+tStop[i]-tStart[i]+1 + indiv in middle
          counter=0;
          for(unsigned int i=0;i<nSubjects;i++){
            if(params.z(i) == c ){
              for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
                timesk3[counter+j] = times[tStart[i]-1+j];
              }
              counter = counter + tStop[i] - tStart[i] + 1;
            }
          }

          vector<double> times_init= timesk3;

          // sort times with corresponding ordering
          std::vector<int> idx(timesk3.size()); //times_init[idx[j]] = timesk3
          int x=0;
          iota(idx.begin(), idx.end(), x++);
          stable_sort(idx.begin(), idx.end(), // sort indexes based on comparing values in times
                      [&](int i1,int i2) { return (timesk3[i1] < timesk3[i2]); });
          stable_sort(timesk3.begin(), timesk3.end());

          std::vector<double> times_permut(timesk3.size()); //timesk ordered with indiv to remove at the end
          std::vector<double> idx_i; //indexes of timepoints of subject i in timesk3
          int ai=0;
          int iii;
          for(int i=0;i<timesk3.size();i++){
            if(idx[i+idx_i.size()] >= ind_permut && idx[i+idx_i.size()] < ind_permut+ tStop[ii]-tStart[ii]+1)
              idx_i.push_back(i+idx_i.size());
            if(i< timesk3.size()-(tStop[ii]-tStart[ii]+1)){
              iii=i+idx_i.size();
            }else {
              iii=idx_i[ai];
              ai++;
            }
            times_permut[i]=timesk3[iii];
          }

          MatrixXd Sigma_inv_ord2=Sigma_inv_ord;
          Permut_cov_sorted(Sigma_inv_ord2, ind_permut, tStop[ii]-tStart[ii]+1, idx); // timesk3=times_init[idx]
          MatrixXd Sigma_ord2;
          Sigma_ord2.setZero(Sigma_inv_ord.rows(),Sigma_inv_ord.rows());
          GP_cov(Sigma_ord2,params.L(c),times_permut,dataset.equalTimes(),kernelType,1);
          logDetPrecMat=Inverse_woodbury(Sigma_ord2,logDetSigma, precMat, Sigma_inv_ord2);
          if(std::isnan(logDetPrecMat) || isinf(logDetPrecMat) || logDetPrecMat== -(std::numeric_limits<double>::max())){
            //fout << "indiv: "<< ii << " z(i) "<<params.z(ii)<< " c "<< c<<endl;
          }
          // Sigma_ord2 and Sigma_inv_ord2 including measurement errors
          // ordered, with subject to remove at the end

          std::vector<int> idx2(timesk.size());
          x=0;
          iota(idx2.begin(), idx2.end(), x++);
          stable_sort(idx2.begin(), idx2.end(), // sort indexes based on comparing values in times
                      [&](int i1,int i2) { return (timesk[i1] < timesk[i2]); });

          VectorXd yk_sorted = yk;
          for(unsigned int j=0;j<sizek;j++)
            yk_sorted(j)=yk(idx2[j]);
          yk=yk_sorted;

        }else{
          GP_cov(precMat,params.L(c),timesk,dataset.equalTimes(),kernelType,1);
          LLT<MatrixXd> lltOfA(precMat); // compute the Cholesky decomposition of A
          MatrixXd L = lltOfA.matrixL();
          logDetPrecMat=  2*L.diagonal().array().log().sum();
          precMat = L.inverse().transpose()*L.inverse();
        }
      }
    }
    dmvnorm =  -0.5*yk.transpose()*precMat*yk - 0.5*sizek*log(2.0*pi<double>()) - 0.5*logDetPrecMat ;
  }else{
    //std::cout << "no one remaining, empty class\n";
    dmvnorm=0;
  }

  if(std::isnan(dmvnorm) || isinf(dmvnorm)){

    std::cout << " dmvnorm  "<< dmvnorm<< " logdetprecMat "<<logDetPrecMat << endl<<
      " zi "<<  params.z(ii) << " c "<< c <<endl<<
        " prod "<< yk.transpose()*precMat*yk <<
          " sizek "<<sizek<<endl;
    dmvnorm=-(std::numeric_limits<double>::max());
  }

  return dmvnorm;
}

//AR Compute logPYiGivenZiWiLongitudinal from sampled meanGP f(Yk|k, fk, B, sigma, w, t)
double logPYiGivenZiWiLongitudinal_meanGP(const pReMiuMParams& params, const pReMiuMData& dataset,
                                          const vector<unsigned int>& nFixedEffects, const int& c,
                                          const  unsigned int& ii){

  //vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  vector<double> y = dataset.continuousY();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  unsigned int nCategoriesY = params.nCategoriesY();
  //AR for meanGP
  vector<double> times_corr = dataset.times_corr();
  VectorXd yi;
  MatrixXd Vi_inv;
  double dmvnorm = 0.0;
  unsigned int ni = 0;
  int unsigned l = 2; // SQexponential or SQuadratic


  for(unsigned int i=0; i<dataset.nSubjects(); i++){
    if((params.z(i) == c &&  ii == dataset.nSubjects()) || ( i == ii ) ){

      ni =  (tStop[i] - tStart[i] + 1);
      yi.resize(ni);
      Vi_inv.setZero(ni,ni);

      for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
        yi(j) = y[tStart[i]-1+j];
        yi(j) -= params.meanGP(c, times_corr[tStart[i]-1+j]);

        for(unsigned int b=0;b<nFixedEffects[0];b++){
          yi(j)-=params.beta(0, b,0, nCategoriesY)*dataset.W(i,b);
        }
        // for(unsigned int b=0;b<nFixedEffects_mix;b++){
        //   yi(j)-=params.beta_mix(c,b,0,nCategoriesY)*dataset.W_mix(i,b);
        // }
      }
      Vi_inv = MatrixXd::Identity(ni, ni)*exp(-params.L(c,l));
      //double eL0 = exp(L[0]);
      //double eL1 = exp(L[1])*2.0;
      //double eL2 = exp(L[2]);

      double out =  -0.5*yi.transpose()*Vi_inv*yi - 0.5*ni*log(2.0*pi<double>())
        - 0.5*ni*params.L(c,l);//*abs(params.L(c,l)); //log(exp(params.L(c,2)*ni)) ;

      dmvnorm += out ;
      if((dmvnorm==0) & (ii <dataset.nSubjects())){
        cout <<  " i "<< ii  << " c "<< c <<" dmvnorm "<<dmvnorm<< " ni "<< ni << " prod "<<yi.transpose()*Vi_inv*yi
             << " L0 "<< params.L(c,0) << " L1 "<< params.L(c,1) << " L2 "<< params.L(c,2)
             << " yi "<<yi<<endl;

      }
    }
  }
  return dmvnorm;
}

//RJ
double logPYiGivenZiWiMVN(const pReMiuMParams& params, const pReMiuMData& dataset,
                          const vector<unsigned int>& nFixedEffects,const int& zi,
                          const unsigned int& i){

  //vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nSubjects = dataset.nSubjects();
  unsigned int nTimes = dataset.nTimes();
  unsigned int nOutcomes = nTimes/nSubjects;
  unsigned int nCategoriesY = params.nCategoriesY();
  VectorXd yi = VectorXd::Zero(nOutcomes);
  VectorXd mu = params.MVNmu(zi);
  MatrixXd workSqrtMVNTau = params.workSqrtMVNTau(zi);
  double workLogDetMVNTau = params.workLogDetMVNTau(zi);
  for(unsigned int j=0;j<nOutcomes;j++){
    yi(j)=dataset.continuousY(i*nOutcomes+j);
    for(unsigned int k=0;k<nFixedEffects[0];k++){
      mu(j)+=params.beta(0, k,0, nCategoriesY)*dataset.W(i,k);
    }
    // for(unsigned int k=0;k<nFixedEffects_mix;k++){
    //   mu(j)+=params.beta_mix(zi,k,0,nCategoriesY)*dataset.W_mix(i,k);
    // }
  }

  return logPdfMultivarNormal(nOutcomes,yi,mu,workSqrtMVNTau,workLogDetMVNTau);
}


//AR Parametric Longitudinal model - Latent class model
double logPYiGivenZiWiLongitudinal_parametric(const pReMiuMParams& params, const pReMiuMData& dataset,
                                   const vector<unsigned int>& nFixedEffects, const int& c,
                                   const unsigned int& ii){

  vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  vector<double> y = dataset.continuousY();
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  unsigned int nOutcomes = dataset.nOutcomes();
  unsigned int nCategoriesY = params.nCategoriesY();
  //AR for meanGP
  VectorXd yi;
  MatrixXd Vi_inv;
  double dmvnorm = 0.0;
  unsigned int ni_m = 0;

  unsigned int sum_ind_m = 0;
  unsigned int curr_i=0;
  for(unsigned int m=0; m<nOutcomes; m++){
    for(unsigned int i=0; i<dataset.nSubjects(); i++){

      ni_m =  (tStop[curr_i] - tStart[curr_i] + 1);

      if((params.z(i) == c &&  ii == dataset.nSubjects()) || ( i == ii ) ){

        yi.resize(ni_m);
        for(unsigned int j=0;j<tStop[curr_i]-tStart[curr_i]+1;j++){
          yi(j) = y[sum_ind_m + tStart[curr_i]-1+j];

          for(unsigned int b=0;b<nFixedEffects[m];b++){
            yi(j)-=params.beta(m,b, 0, nCategoriesY)*dataset.W_LME(m, tStart[curr_i]-1+j,b);
          }
          for(unsigned int b=0;b<nFixedEffects_mix[m];b++){
            yi(j)-=params.beta_mix(m, c,b,0,nCategoriesY)*dataset.W_LME_mix(m, tStart[curr_i]-1+j,b);

          }
        }

        //Random effects of subject i
        bool cond_RE=true;
        if(!cond_RE){
          MatrixXd block=dataset.W_RE(m, tStart[curr_i]-1, 0, ni_m, dataset.nRandomEffects(m));

          MatrixXd Vi=block*params.covRE(m,0)*block.transpose()+MatrixXd::Identity(ni_m, ni_m) * params.SigmaE(m);

          LLT<MatrixXd> lltOfA(Vi);
          MatrixXd L = lltOfA.matrixL();
          double logDetPrecMat=  2*L.diagonal().array().log().sum();
          MatrixXd precMat = L.inverse().transpose()*L.inverse(); //Vi^-1

          dmvnorm +=  -0.5*yi.transpose()*precMat*yi - 0.5*ni_m*log(2.0*pi<double>()) - 0.5*logDetPrecMat ;
          // yi -= block*params.RandomEffects(i);

        }else{
          MatrixXd block=dataset.W_RE(m,tStart[curr_i]-1, 0, ni_m, dataset.nRandomEffects(m));
          yi -= block*params.RandomEffects(m,i);

          MatrixXd Vi=MatrixXd::Identity(ni_m, ni_m) * params.SigmaE(m);
          dmvnorm += -0.5*yi.transpose()*Vi.inverse()*yi - 0.5*ni_m*log(2.0*pi<double>()) - 0.5*ni_m*log(params.SigmaE(m));
        }




         // If random effects integrated out:
         // MatrixXd V = block *params.covRE(c)* block.transpose() + MatrixXd::Identity(ni, ni) * exp(params.SigmaE());
         // LLT<MatrixXd> lltOfA(V); // compute the Cholesky decomposition of A
         // MatrixXd L = lltOfA.matrixL();
         // double logDetPrecMat=  2*log(L.determinant());
         // MatrixXd Vinv = L.inverse().transpose()*L.inverse();
         //
         // dmvnorm +=  -0.5*yi.transpose()*Vinv*yi - 0.5*ni*log(2.0*pi<double>())
         //   - 0.5*logDetPrecMat;
       }
       curr_i ++;
       sum_ind_m += ni_m;
     }
   }
  return dmvnorm;
}

double logPYiGivenZiWiNormal(const pReMiuMParams& params, const pReMiuMData& dataset,
                             const vector<unsigned int>& nFixedEffects,const int& zi,
                             const unsigned int& i){

  double mu;
  //vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();

  mu=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    mu+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j); //marker 0, effect j, Ycategory 0 nCategoriesY=1
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   mu+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
  // }

  return logPdfNormal(dataset.continuousY(i),mu,sqrt(params.sigmaSqY()));
}

double logPYiGivenZiWiNormalSpatial(const pReMiuMParams& params, const pReMiuMData& dataset,
                                    const vector<unsigned int>& nFixedEffects,const int& zi,
                                    const unsigned int& i){

  double mu;
  //vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();

  mu=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    mu+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   mu+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
  // }
  mu+=params.uCAR(i);

  return logPdfNormal(dataset.continuousY(i),mu,sqrt(params.sigmaSqY()));
}

double logPYiGivenZiWiCategorical(const pReMiuMParams& params, const pReMiuMData& dataset,
                                  const vector<unsigned int>& nFixedEffects,
                                  const int& zi,
                                  const unsigned int& i){

  vector<double> lambda;
  //vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();

  lambda.resize(nCategoriesY);

  double lambdaSum = 1.0;

  double value = 0.0;
  for (unsigned int k=0;k<nCategoriesY;k++){
    for (unsigned int j=0;j<nFixedEffects[0];j++){
      value+=params.beta(0, j, k, nCategoriesY)*dataset.W(i,j);
    }
    // for (unsigned int j=0;j<nFixedEffects_mix;j++){
    //   value+=params.beta_mix(zi,j,k,nCategoriesY)*dataset.W_mix(i,j);
    // }
    lambda[k] = exp(value + params.theta(zi,k));
    lambdaSum += exp(value + params.theta(zi,k));
    value = 0.0;
  }
  vector<double> p;
  p.resize(nCategoriesY+1);
  p[0]=1.0/lambdaSum;
  for (unsigned int k=0;k<nCategoriesY;k++){
    p[k+1]=lambda[k]/lambdaSum;
  }
  return logPdfMultinomialSizeOne(dataset.discreteY(i),p);
}


double logPYiGivenZiWiSurvival(const pReMiuMParams& params, const pReMiuMData& dataset,
                               const vector<unsigned int>& nFixedEffects,const int& zi,
                               const unsigned int& i){

  unsigned int weibullFixedShape=params.nu().size();
  //vector<unsigned int> nFixedEffects_mix = dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();

  double lambda = params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    lambda+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   lambda+=params.beta_mix(zi, j, 0,nCategoriesY)*dataset.W_mix(i,j);
  // }
  double nu=0;

  if (weibullFixedShape==1) {
    nu=params.nu(0);
  } else {
    nu=params.nu(zi);
  }
  return logPdfWeibullCensored(dataset.continuousY(i), nu, exp(lambda), dataset.censoring(i));
}


vector<double> pReMiuMLogPost(const pReMiuMParams& params,
                              const mcmcModel<pReMiuMParams,
                                              pReMiuMOptions,
                                               pReMiuMData>& model){
  const pReMiuMData& dataset = model.dataset();
  const string outcomeType = model.dataset().outcomeType();
  const string kernelType = model.options().kernelType(); //AR
  const string covariateType = model.dataset().covariateType();
  const bool includeResponse = model.options().includeResponse();
  const bool responseExtraVar = model.options().responseExtraVar();
  const string varSelectType = model.options().varSelectType();
  double fixedAlpha = model.options().fixedAlpha();
  unsigned int nSubjects=dataset.nSubjects();
  //unsigned int nTimes=dataset.nTimes();
  unsigned int nOutcomes=dataset.nOutcomes();
  unsigned int maxNClusters=params.maxNClusters();
  unsigned int nCovariates=dataset.nCovariates();
  unsigned int nDiscreteCov=dataset.nDiscreteCovs();
  unsigned int nContinuousCov=dataset.nContinuousCovs();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  vector<unsigned int> nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY=dataset.nCategoriesY();
  vector<unsigned int> nCategories = dataset.nCategories();
  const pReMiuMHyperParams& hyperParams = params.hyperParams();
  const bool includeCAR=model.options().includeCAR();
  const bool weibullFixedShape=model.options().weibullFixedShape();
  const bool useNormInvWishPrior=model.options().useNormInvWishPrior();

  // (Augmented) Log Likelihood first
  // We want p(y,X|z,params,W) = p(y|z=c,W,params)p(X|z=c,params)
  // if(outcomeType.compare("Longitudinal")==0)
  //   nOutcomes=nTimes/nSubjects;

  double logLikelihood=0.0;
  //std::fstream foutL("compare.txt", std::ios::in | std::ios::out | std::ios::app);

  // Add in contribution from X
  for(unsigned int i=0;i<nSubjects;i++){
    //P(xi|zi,params)
    logLikelihood+=params.workLogPXiGivenZi(i);
  }

  // Add in contribution from Y
  vector<double> extraVarPriorVal(nSubjects,0.0);
  vector<double> extraVarPriorMean(nSubjects,0.0);

  if(includeResponse){

    double (*logPYiGivenZiWi)(const pReMiuMParams&,const pReMiuMData&,
            const vector<unsigned int>&,const int&,
            const unsigned int&)=NULL;

    if(outcomeType.compare("Bernoulli")==0){
      if(!responseExtraVar){
        logPYiGivenZiWi = &logPYiGivenZiWiBernoulli;
      }else{
        logPYiGivenZiWi = &logPYiGivenZiWiBernoulliExtraVar;
        for(unsigned int i=0;i<nSubjects;i++){
          extraVarPriorVal[i]=params.lambda(i);
          int zi=params.z(i);
          extraVarPriorMean[i]=params.theta(zi,0);
          for(unsigned int j=0;j<nFixedEffects[0];j++){
            extraVarPriorMean[i]+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
          }
          // for(unsigned int j=0;j<nFixedEffects_mix;j++){
          //   extraVarPriorMean[i]+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
          // }
        }
      }
    }else if(outcomeType.compare("Binomial")==0){
      if(!responseExtraVar){
        logPYiGivenZiWi = &logPYiGivenZiWiBinomial;
      }else{
        logPYiGivenZiWi = &logPYiGivenZiWiBinomialExtraVar;
        for(unsigned int i=0;i<nSubjects;i++){
          extraVarPriorVal[i]=params.lambda(i);
          int zi=params.z(i);
          extraVarPriorMean[i]=params.theta(zi,0);
          for(unsigned int j=0;j<nFixedEffects[0];j++){
            extraVarPriorMean[i]+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
          }
          // for(unsigned int j=0;j<nFixedEffects_mix;j++){
          //   extraVarPriorMean[i]+=params.beta_mix(zi,j, 0,nCategoriesY)*dataset.W_mix(i,j);
          // }
        }
      }
    }else if(outcomeType.compare("Poisson")==0){
      if(!responseExtraVar){
        if (includeCAR){
          logPYiGivenZiWi = &logPYiGivenZiWiPoissonSpatial;
        }else{
          logPYiGivenZiWi = &logPYiGivenZiWiPoisson;
        }
      }else{
        logPYiGivenZiWi = &logPYiGivenZiWiPoissonExtraVar;
        for(unsigned int i=0;i<nSubjects;i++){
          extraVarPriorVal[i]=params.lambda(i);
          int zi=params.z(i);
          extraVarPriorMean[i]=params.theta(zi,0);
          for(unsigned int j=0;j<nFixedEffects[0];j++){
            extraVarPriorMean[i]+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
          }
          // for(unsigned int j=0;j<nFixedEffects_mix;j++){
          //   extraVarPriorMean[i]+=params.beta_mix(zi,j, 0,nCategoriesY)*dataset.W_mix(i,j);
          // }
          extraVarPriorMean[i]+=dataset.logOffset(i);
        }
      }
    }else if(outcomeType.compare("Categorical")==0){
      logPYiGivenZiWi = &logPYiGivenZiWiCategorical;
    }else if(outcomeType.compare("Normal")==0){
      if (includeCAR){
        logPYiGivenZiWi = &logPYiGivenZiWiNormalSpatial;
      }else{
        logPYiGivenZiWi = &logPYiGivenZiWiNormal;
      }
    }else if(outcomeType.compare("Survival")==0){
      logPYiGivenZiWi = &logPYiGivenZiWiSurvival;
    }else if(outcomeType.compare("MVN")==0){
      logPYiGivenZiWi = &logPYiGivenZiWiMVN;
    }else if(outcomeType.compare("Longitudinal")==0){
      if(model.options().sampleGPmean()){
        logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal_meanGP;
      }else{
        logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal;
      }
      for(unsigned int c=0;c<maxNClusters;c++){
        logLikelihood+=logPYiGivenZiWi(params,dataset,nFixedEffects,c,nSubjects);
      }
    }    else if(outcomeType.compare("LME")==0){
      logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal_parametric;
    }

    //RJ logLikelihood is a sum only if Yi are conditionally independent
    if(outcomeType.compare("Longitudinal")!=0){
      for(unsigned int i=0;i<nSubjects;i++){
        int zi = params.z(i);
         double temp = logPYiGivenZiWi(params,dataset,nFixedEffects,zi,i);
         logLikelihood+= temp;
      }
    }

  // Now need to add in the prior (i.e. p(z,params) in above notation)
  // Note we integrate out u and all parameters in Theta with no members
  double logPrior=0.0;

  // Prior for z
  for(unsigned int i=0;i<nSubjects;i++){
    int zi = params.z(i);
    logPrior+=params.logPsi(zi);
  }


  // Prior for V (we only need to include these up to maxNCluster, but we do need
  // to include all V, whether or not a cluster is empty, as the V themselves
  //don't correspond to a cluster
  for(unsigned int c=0;c<maxNClusters;c++){
    logPrior+=logPdfBeta(params.v(c),1.0-params.dPitmanYor(),params.alpha()+params.dPitmanYor()*(c+1));
  }

  if(isinf(logPrior)){
    //for(unsigned int c=0;c<maxNClusters;c++){
    //  foutL <<c << " v " << params.v(c) << " pitmanyor " <<params.dPitmanYor() << " alpha "<< params.alpha()<<endl;
    //}
  }
  // Prior for alpha
  if(fixedAlpha<=-1){
    logPrior+=logPdfGamma(params.alpha(),hyperParams.shapeAlpha(),hyperParams.rateAlpha());
  }

  // Prior for phi
  if(covariateType.compare("Discrete")==0){
    // If covariate type is discrete
    for(unsigned int c=0;c<maxNClusters;c++){
      if(params.workNXInCluster(c)>0){
        for(unsigned int j=0;j<nCovariates;j++){
          // We use a Dirichlet prior
          vector<double> dirichParams(nCategories[j]);
          for(unsigned int k=0;k<nCategories[j];k++){
            dirichParams[k]=hyperParams.aPhi(j);
          }
          logPrior+=logPdfDirichlet(params.logPhi(c,j),dirichParams,true);
        }
      }
    }
  }else if(covariateType.compare("Normal")==0){
    // If covariate type is Normal
    // Add in the prior for mu_c and Sigma_c for each c
    for(unsigned int c=0;c<maxNClusters;c++){
      if(params.workNXInCluster(c)>0){
        if (useNormInvWishPrior){
          logPrior+=logPdfMultivarNormal(nCovariates,params.mu(c),hyperParams.mu0(),hyperParams.nu0()*params.Tau(c),nCovariates*hyperParams.nu0()+params.workLogDetTau(c));
        }else{
          logPrior+=logPdfMultivarNormal(nCovariates,params.mu(c),hyperParams.mu0(),hyperParams.workSqrtTau0(),hyperParams.workLogDetTau0());
        }
        logPrior+= logPdfWishart(nCovariates, params.Tau(c), params.workLogDetTau(c), hyperParams.workInverseR0(), hyperParams.workLogDetR0(), (double)hyperParams.kappa0());
      }
    }

  }else if(covariateType.compare("Mixed")==0){
    // If covariate type is discrete
    for(unsigned int c=0;c<maxNClusters;c++){
      if(params.workNXInCluster(c)>0){
        for(unsigned int j=0;j<nDiscreteCov;j++){
          // We use a Dirichlet prior
          vector<double> dirichParams(nCategories[j]);
          for(unsigned int k=0;k<nCategories[j];k++){
            dirichParams[k]=hyperParams.aPhi(j);
          }
          logPrior+=logPdfDirichlet(params.logPhi(c,j),dirichParams,true);
        }

        if (useNormInvWishPrior){
          logPrior+=logPdfMultivarNormal(nContinuousCov,params.mu(c),hyperParams.mu0(),hyperParams.nu0()*params.Tau(c),nContinuousCov*hyperParams.nu0()+params.workLogDetTau(c));
        }else{
          logPrior+=logPdfMultivarNormal(nContinuousCov,params.mu(c),hyperParams.mu0(),hyperParams.workSqrtTau0(),hyperParams.workLogDetTau0());
        }
        logPrior+=logPdfWishart(nContinuousCov,params.Tau(c),params.workLogDetTau(c),hyperParams.workInverseR0(),hyperParams.workLogDetR0(),(double)hyperParams.kappa0());
      }
    }
  }

  // Prior for variable selection parameters
  if(varSelectType.compare("None")!=0){
    if(varSelectType.compare("BinaryCluster")==0){
      for(unsigned int j=0;j<nCovariates;j++){
        if(params.rho(j)>0){
          for(unsigned int c=0;c<maxNClusters;c++){
            if(params.workNXInCluster(c)>0){
              logPrior+=params.gamma(c,j)*log(params.rho(j))+(1.0-params.gamma(c,j))*log(1.0-params.rho(j));
            }
          }
        }
      }
    }

    // We can add in the prior for rho and omega
    logPrior+=log(hyperParams.atomRho());
    for(unsigned int j=0;j<nCovariates;j++){
      if(params.omega(j)==1){
        logPrior+=logPdfBeta(params.rho(j),hyperParams.aRho(),hyperParams.bRho());
      }
    }

  }

  if(includeResponse){

    if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
      // Prior for theta
      // We use a location/scale t distribution
      // http://www.mathworks.com/help/toolbox/stats/brn2ivz-145.html
      // as per Molitor et al. 2008 (from Gelman et al. 2008)
      // This is different from Papathomas who uses a normal
      for(unsigned int c=0;c<maxNClusters;c++){
        if(params.workNXInCluster(c)>0){
          for (unsigned int k=0;k<nCategoriesY;k++){
            logPrior+=logPdfLocationScaleT(params.theta(c,k),hyperParams.muTheta(),
                                           hyperParams.sigmaTheta(),hyperParams.dofTheta());
          }
        }
      }
    }

    // Prior for beta
    // There were no fixed effects in the Molitor paper but to be consistent with
    // theta we use the same distribution (based on same reasoning).
    // Note we should pre-process variables to standardise the fixed effects to
    // have 0 mean and standard dev of 0.5
    for (unsigned int m=0;m<nOutcomes;m++){
      for(unsigned int j=0;j<nFixedEffects[m];j++){
        for (unsigned int k=0;k<nCategoriesY;k++){
          if(outcomeType.compare("LME")!=0){
            logPrior+=logPdfLocationScaleT(params.beta(0, j, k, nCategoriesY),hyperParams.muBeta(),
                                           hyperParams.sigmaBeta(),hyperParams.dofBeta());
          }else{
            logPrior+=logPdfNormal(params.beta(m, j, k, nCategoriesY),hyperParams.muBeta(),
                                   hyperParams.sigmaBeta());
          }
        }
      }

      for(unsigned int j=0;j<nFixedEffects_mix[m];j++){
        for(unsigned int c=0;c<maxNClusters;c++){
          for (unsigned int k=0;k<nCategoriesY;k++){
            //if(outcomeType.compare("LME")!=0){
            logPrior+=logPdfLocationScaleT(params.beta_mix(m, c, j, k, nCategoriesY),hyperParams.muBeta(),
                                           hyperParams.sigmaBeta(),hyperParams.dofBeta());
            //}else{
            // logPrior+=logPdfNormal(params.beta_mix(c, j+k*nFixedEffects_mix),hyperParams.muBeta(),
            //                        hyperParams.sigmaBeta());
            //}
          }
        }
      }
    }




    // Take account priors for epsilon and tauEpsilon if there is extra variation
    if(responseExtraVar){
      double sigma = 1.0/sqrt(params.tauEpsilon());
      for(unsigned int i=0;i<nSubjects;i++){
        logPrior+=logPdfNormal(extraVarPriorVal[i],extraVarPriorMean[i],sigma);
      }
      logPrior+=logPdfGamma(params.tauEpsilon(),hyperParams.shapeTauEpsilon(),
                            hyperParams.rateTauEpsilon());
    }

    if(outcomeType.compare("Normal")==0){
      double tau = 1.0/params.sigmaSqY();
      logPrior+=logPdfGamma(tau,hyperParams.shapeSigmaSqY(),hyperParams.scaleSigmaSqY());
    }

    if(outcomeType.compare("Survival")==0){
      if (weibullFixedShape){
        logPrior+=logPdfGamma(params.nu(0),hyperParams.shapeNu(),hyperParams.scaleNu());
      } else {
        for (unsigned int c=0;c<maxNClusters;c++) {
          logPrior+=logPdfGamma(params.nu(c),hyperParams.shapeNu(),hyperParams.scaleNu());
        }
      }
    }

    //RJ add up prior for L
    if(outcomeType.compare("Longitudinal")==0){
      for(unsigned int c=0;c<maxNClusters;c++){
        unsigned int nL;
        if(kernelType.compare("SQexponential")==0){    //AR
          nL=3;
        }else{
          nL=4;
        }

        for(unsigned int l=0;l<nL;l++){
          if(hyperParams.sigmaL(l)>0)
            logPrior+= logPdfNormal(params.L(c,l),hyperParams.muL(l),
                                    hyperParams.sigmaL(l));
        }

        if(model.options().sampleGPmean()){//AR
          double a =logPdfMultivariateNormal(params.meanGP(c),params.L(c), dataset.times_unique(), kernelType);
          logPrior+= a;
          //if(logPrior>pow(10,10))
          //  foutL << c <<" logPrior7_f  "<< logPrior << " a "<< a <<endl;

          if(a>pow(10,10)){
            std::cout << c <<" writeoutput "<<endl;
            std::cout << " p(f|L) "<< a << " logPrior "<< logPrior <<endl;
            std::cout << "L "<<params.L(c,0) << " "<<params.L(c,1) << " "<<params.L(c,2) << " "<<endl;
            std::cout << " meanGP "<<endl;
            for(unsigned int l=0;l<dataset.times_unique().size();l++)
              std::cout << params.meanGP(c,l)<<" ";
            std::cout << endl;
          }
        }
      }
    }

    if(outcomeType.compare("MVN")==0){
      for(unsigned int c=0;c<maxNClusters;c++){
        //logPrior+=logPdfMultivarNormal(nOutcomes,params.MVNmu(c),hyperParams.MVNmu0(),hyperParams.MVNkappa0()*params.MVNTau(c),nOutcomes*hyperParams.MVNkappa0()+params.workLogDetMVNTau(c));
        //logPrior+=logPdfWishart(nOutcomes,params.MVNTau(c),params.workLogDetMVNTau(c),hyperParams.workInverseMVNR0(),hyperParams.workLogDetMVNR0(),(double)hyperParams.MVNnu0());
        logPrior += logPdfNIW(nOutcomes, hyperParams.workInverseMVNR0(), hyperParams.workLogDetMVNR0(), params.MVNTau(c), -params.workLogDetMVNTau(c), hyperParams.MVNkappa0(), hyperParams.MVNnu0(), params.MVNmu(c),hyperParams.MVNmu0());
      }
    }

     if(outcomeType.compare("LME")==0){
       vector<unsigned int> nRandomEffects=dataset.nRandomEffects();
       for(unsigned int m=0;m<nOutcomes;m++){
        for(unsigned int c=0;c<maxNClusters;c++){
    //     //logPrior+=logPdfMultivarNormal(nOutcomes,params.MVNmu(c),hyperParams.MVNmu0(),hyperParams.MVNkappa0()*params.MVNTau(c),nOutcomes*hyperParams.MVNkappa0()+params.workLogDetMVNTau(c));
          logPrior+= logPdfInverseWishart(nRandomEffects[m], params.covRE(m,0), -params.workLogDetTauLME(m,0), hyperParams.SigmaLME_R0(m),
                                 -hyperParams.workLogDetTauLME_R0(m),(double)hyperParams.SigmaLME_kappa0(m));
         // double logPdfInverseWishart(const unsigned int& dimA, const MatrixXd& A, const double& logDetA, const MatrixXd& covR,
        // const double& logDetR, const double& kappa){
        // MatrixXd work = covR*A.inverse();
        // return 0.5*(kappa*logDetR-(kappa+(double)dimA+1)*logDetA-work.trace())
        //   - (0.5*kappa*(double)dimA*log(2.0)+logMultivarGammaFn(kappa/2.0,dimA));


         //MatrixXd Tau = params.covRE(c).inverse();
        //logPrior+=logPdfWishart(dataset.nRandomEffects(), Tau, params.workLogDetMVNTau(c)
        }
        logPrior += logPdfGamma(1/params.SigmaE(m),hyperParams.eps_shape(),
                                hyperParams.eps_rate()); //Gamma(shape,rate=1/scale)

         for(unsigned int i=0;i<nSubjects;i++){
           //unsigned int zi = params.z(i);
           VectorXd RE=params.RandomEffects(m,i);
           VectorXd mu;
           mu.setZero(nRandomEffects[m]);
           logPrior +=  logPdfMultivarNormal(RE.size(),RE,mu,params.workSqrtTauLME(m,0),params.workLogDetTauLME(m,0));
         }
       }

       //logPrior += logPdfNormal(params.SigmaE(),hyperParams.muSigmaE(),
      //                         hyperParams.sigmaSigmaE());
       //logPrior += logPdfScaleInvChiSquare(params.SigmaE(),hyperParams.eps_vu(),
      //                          hyperParams.eps_sigma2_0());


    }

    // Prior for TauCAR and UCAR
    if (includeCAR){
      logPrior+=logPdfGamma(params.TauCAR(),hyperParams.shapeTauCAR(),hyperParams.rateTauCAR());
      logPrior+=logPdfIntrinsicCAR(params.uCAR(), dataset.neighbours() , params.TauCAR());
    }
  }

   vector<double> outVec(3);
   outVec[0]=logLikelihood+logPrior;
   outVec[1]=logLikelihood;
   outVec[2]=logPrior;

   std::cout << " logLikelihood "<<logLikelihood <<  " logPrior "<<logPrior<<std::endl;
  }
   //return outVec;
   vector<double> outVec={0,1};
   return outVec;
}


// Log conditional posterior for phi (only used in continuous variable selection)
double logCondPostPhicj(const pReMiuMParams& params,
                        const mcmcModel<pReMiuMParams,
                                        pReMiuMOptions,
                                        pReMiuMData>& model,
                                        const unsigned int& c,
                                        const unsigned int& j){


  const pReMiuMData& dataset = model.dataset();
  unsigned int nSubjects=dataset.nSubjects();
  vector<unsigned int> nCategories = dataset.nCategories();
  const pReMiuMHyperParams& hyperParams = params.hyperParams();

  double out=0.0;

  // Add in contribution from likelihood
  for(unsigned int i=0;i<nSubjects;i++){
    //P(xi|zi,params)
    if(params.z(i)==(int)c){
      out+=params.workLogPXiGivenZi(i);
    }
  }

  // Add in contribution from prior
  vector<double> dirichParams(nCategories[j]);
  for(unsigned int k=0;k<nCategories[j];k++){
    dirichParams[k]=hyperParams.aPhi(j);
  }
  out+=logPdfDirichlet(params.logPhi(c,j),dirichParams,true);

  return out;
}

// Log conditional posterior for rho and omega (only used in variable selection)
double logCondPostRhoOmegaj(const pReMiuMParams& params,
                            const mcmcModel<pReMiuMParams,
                                            pReMiuMOptions,
                                            pReMiuMData>& model,
                                            const unsigned int& j){


  const pReMiuMData& dataset = model.dataset();
  unsigned int nSubjects=dataset.nSubjects();
  unsigned int maxNClusters=params.maxNClusters();
  unsigned int nCovariates=dataset.nCovariates();
  string varSelectType = model.options().varSelectType();

  const pReMiuMHyperParams& hyperParams = params.hyperParams();

  double out=0.0;

  // Add in contribution from likelihood (only in the case of continuous switches)
  if(varSelectType.compare("Continuous")==0){
    for(unsigned int i=0;i<nSubjects;i++){
      //P(xi|zi,params)
      out+=params.workLogPXiGivenZi(i);
    }
  }else{

    // Add in contribution from prior
    if(params.omega(j)==1){
      for(unsigned int c=0;c<maxNClusters;c++){
        out+=params.gamma(c,j)*log(params.rho(j))+(1.0-params.gamma(c,j))*log(1.0-params.rho(j));
      }
    }else{
      for(unsigned int c=0;c<maxNClusters;c++){
        if(params.gamma(c,j)>0.5){
          out = -(numeric_limits<double>::max());
          return out;
        }
      }
    }

  }
  // We can add in the prior for rho and omega
  // We keep the loop here because it saves evaluations in the continuous case
  for(unsigned int j1=0;j1<nCovariates;j1++){
    out+=log(hyperParams.atomRho());
    if (params.omega(j1)==1){
      out+=logPdfBeta(params.rho(j1),hyperParams.aRho(),hyperParams.bRho());
    }
  }
  return out;
}

double logCondPostThetaBeta(const pReMiuMParams& params,
                            const mcmcModel<pReMiuMParams,
                                            pReMiuMOptions,
                                            pReMiuMData>& model){
  const pReMiuMData& dataset = model.dataset();
  const string outcomeType = model.dataset().outcomeType();
  const bool responseExtraVar = model.options().responseExtraVar();
  unsigned int nSubjects=dataset.nSubjects();
  unsigned int maxNClusters=params.maxNClusters();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY=dataset.nCategoriesY();
  const pReMiuMHyperParams& hyperParams = params.hyperParams();
  const bool includeCAR=model.options().includeCAR();

  double out=0.0;

  // Add in contribution from Y
  vector<double> extraVarPriorVal(nSubjects,0.0);
  vector<double> extraVarPriorMean(nSubjects,0.0);
  double (*logPYiGivenZiWi)(const pReMiuMParams&,const pReMiuMData&,
          const vector<unsigned int>&,const int&,
          const unsigned int&) = NULL;

  if(outcomeType.compare("Bernoulli")==0){
    if(!responseExtraVar){
      logPYiGivenZiWi = &logPYiGivenZiWiBernoulli;
    }else{
      logPYiGivenZiWi = &logPYiGivenZiWiBernoulliExtraVar;
      for(unsigned int i=0;i<nSubjects;i++){
        extraVarPriorVal[i]=params.lambda(i);
        int zi=params.z(i);
        extraVarPriorMean[i]=params.theta(zi,0);
        for(unsigned int j=0;j<nFixedEffects[0];j++){
          extraVarPriorMean[i]+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
        }
        // for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //   extraVarPriorMean[i]+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
        // }
      }
    }
  }else if(outcomeType.compare("Binomial")==0){
    if(!responseExtraVar){
      logPYiGivenZiWi = &logPYiGivenZiWiBinomial;
    }else{
      logPYiGivenZiWi = &logPYiGivenZiWiBinomialExtraVar;
      for(unsigned int i=0;i<nSubjects;i++){
        extraVarPriorVal[i]=params.lambda(i);
        int zi=params.z(i);
        extraVarPriorMean[i]=params.theta(zi,0);
        for(unsigned int j=0;j<nFixedEffects[0];j++){
          extraVarPriorMean[i]+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
        }
        // for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //   extraVarPriorMean[i]+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
        // }
      }
    }
  }else if(outcomeType.compare("Poisson")==0){
    if(!responseExtraVar){
      if (includeCAR){
        logPYiGivenZiWi = &logPYiGivenZiWiPoissonSpatial;
      }else{
        logPYiGivenZiWi = &logPYiGivenZiWiPoisson;
      }
    }else{
      logPYiGivenZiWi = &logPYiGivenZiWiPoissonExtraVar;
      for(unsigned int i=0;i<nSubjects;i++){
        extraVarPriorVal[i]=params.lambda(i);
        int zi=params.z(i);
        extraVarPriorMean[i]=params.theta(zi,0);
        for(unsigned int j=0;j<nFixedEffects[0];j++){
          extraVarPriorMean[i]+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
        }
        // for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //   extraVarPriorMean[i]+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
        // }
        extraVarPriorMean[i]+=dataset.logOffset(i);
      }
    }
  }else if(outcomeType.compare("Categorical")==0){
    logPYiGivenZiWi = &logPYiGivenZiWiCategorical;
  }else if(outcomeType.compare("Normal")==0){
    if (includeCAR){
      logPYiGivenZiWi = &logPYiGivenZiWiNormalSpatial;
    }else{
      logPYiGivenZiWi = &logPYiGivenZiWiNormal;
    }
  }else if(outcomeType.compare("Survival")==0){
    logPYiGivenZiWi = &logPYiGivenZiWiSurvival;
  }else if(outcomeType.compare("MVN")==0){
    logPYiGivenZiWi = &logPYiGivenZiWiMVN;
  }else if(outcomeType.compare("Longitudinal")==0){
    if( model.options().sampleGPmean()){//AR
      logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal_meanGP; //AR
    }else{
      logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal;
    }
    for(unsigned int c=0;c<maxNClusters;c++){
      out+=logPYiGivenZiWi(params,dataset,nFixedEffects,c,nSubjects);
    }
  }else if(outcomeType.compare("LME")==0){
    logPYiGivenZiWi = &logPYiGivenZiWiLongitudinal_parametric;
    for(unsigned int c=0;c<maxNClusters;c++){
      out+=logPYiGivenZiWi(params,dataset,nFixedEffects,c,nSubjects);
    }
  }
  //RJ logLikelihood is a sum only if Yi are conditionally independent
  if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
    for(unsigned int i=0;i<nSubjects;i++){
      int zi = params.z(i);
      out+=logPYiGivenZiWi(params,dataset,nFixedEffects,zi,i);
    }
  }

  // Prior for theta
  // We use a location/scale t distribution
  // http://www.mathworks.com/help/toolbox/stats/brn2ivz-145.html
  // as per Molitor et al. 2008 (from Gelman et al. 2008)
  // This is different from Papathomas who uses a normal
  if(outcomeType.compare("Longitudinal")!=0 && outcomeType.compare("LME")!=0){
    for(unsigned int c=0;c<maxNClusters;c++){
      for (unsigned int k=0;k<nCategoriesY;k++){
        double temppp = logPdfLocationScaleT(params.theta(c,k),hyperParams.muTheta(),
                                             hyperParams.sigmaTheta(),hyperParams.dofTheta());
        out+=temppp;
      }
    }
  }

  if(outcomeType.compare("LME")!=0){
    for(unsigned int j=0;j<nFixedEffects[0];j++){
      for (unsigned int k=0;k<nCategoriesY;k++){
        out+=logPdfLocationScaleT(params.beta(0, j, k, nCategoriesY),hyperParams.muBeta(),
                                  hyperParams.sigmaBeta(),hyperParams.dofBeta());
        //out+= logPdfNormal(params.beta(j,k), hyperParams.muBeta(),
        //                   hyperParams.sigmaBeta());

      }
    }

    // for(unsigned int j=0;j<nFixedEffects_mix;j++){
    //   for(unsigned int c=0;c<maxNClusters;c++){
    //     for (unsigned int k=0;k<nCategoriesY;k++){
    //       out+=logPdfLocationScaleT(params.beta_mix(0,c,j,k,nCategoriesY),hyperParams.muBeta(),
    //                                 hyperParams.sigmaBeta(),hyperParams.dofBeta());
    //       //out+= logPdfNormal(params.beta_mix(c,j+k*nFixedEffects_mix),hyperParams.muBeta(),
    //       //                   hyperParams.sigmaBeta());
    //     }
    //   }
    // }
  }else{
    printf("should not use this function when outcome = LME");
  }
  // Prior for beta
  // There were no fixed effects in the Molitor paper but to be consistent with
  // theta we use the same distribution (based on same reasoning).
  // Note we should pre-process variables to standardise the fixed effects to
  // have 0 mean and standard dev of 0.5


  if(responseExtraVar){
    for(unsigned int i=0;i<nSubjects;i++){
      out+=logPdfNormal(extraVarPriorVal[i],extraVarPriorMean[i],1/sqrt(params.tauEpsilon()));
    }
  }
  return out;
}

double logCondPostLambdaiBernoulli(const pReMiuMParams& params,
                                   const mcmcModel<pReMiuMParams,
                                                   pReMiuMOptions,
                                                   pReMiuMData>& model,
                                                   const unsigned int& i){

  const pReMiuMData& dataset = model.dataset();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();
  //Rprintf("Verify that nCategoriesY=1 in this case");

  int zi = params.z(i);
  double meanVal = params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    meanVal+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   meanVal+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
  // }
  return logPYiGivenZiWiBernoulliExtraVar(params,dataset,nFixedEffects,zi,i)
    + logPdfNormal(params.lambda(i),meanVal,1.0/sqrt(params.tauEpsilon()));

}

double logCondPostLambdaiBinomial(const pReMiuMParams& params,
                                  const mcmcModel<pReMiuMParams,
                                                  pReMiuMOptions,
                                                  pReMiuMData>& model,
                                                  const unsigned int& i){

  const pReMiuMData& dataset = model.dataset();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();
  //Rprintf("Verify that nCategoriesY=1 in this case");

  int zi = params.z(i);
  double meanVal = params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    meanVal+=params.beta(0, j, 0, nCategoriesY)*dataset.W(i,j);
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   meanVal+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
  // }
  return logPYiGivenZiWiBinomialExtraVar(params,dataset,nFixedEffects,zi,i)
    + logPdfNormal(params.lambda(i),meanVal,1.0/sqrt(params.tauEpsilon()));

}



  //AR logCondPostSigmaE function definition: p(sigma)p(Y^k|sigma_k,k)
  //addprior = 0 if prior should not be computed
  double logCondPostSigmaE(const pReMiuMParams& params,
                      const mcmcModel<pReMiuMParams,
                                      pReMiuMOptions,
                                      pReMiuMData>& model){

    const pReMiuMData& dataset = model.dataset();
    const pReMiuMHyperParams& hyperParams = params.hyperParams();
    vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
    unsigned int nSubjects=dataset.nSubjects();
    unsigned int nOutcomes=dataset.nOutcomes();
    unsigned int maxNClusters=params.maxNClusters();

    double logPrior = 0.0;
    double out = 0;

    //logPrior+=logPdfScaleInvChiSquare(params.SigmaE(),hyperParams.eps_vu(),
    //                      hyperParams.eps_sigma2_0());

    for(unsigned int m=0;m<nOutcomes;m++){
      logPrior += logPdfGamma(1/params.SigmaE(m),hyperParams.eps_shape(),
                              hyperParams.eps_rate()); //Gamma(shape,rate=1/scale)
    }


    // double sigma = 1.0/sqrt(params.tauEpsilon());
    // for(unsigned int i=0;i<nSubjects;i++){
    //   logPrior+=logPdfNormal(extraVarPriorVal[i],extraVarPriorMean[i],sigma);
    // }
    // logPrior+=logPdfGamma(params.tauEpsilon(),hyperParams.shapeTauEpsilon(),
    //                       hyperParams.rateTauEpsilon());

    for(unsigned int c=0;c<maxNClusters;c++){
      out +=  logPYiGivenZiWiLongitudinal_parametric(params,dataset,nFixedEffects,c,nSubjects);
    }

    out += logPrior;
    return out;
  }

//RJ logCondPostL function definition: p(sigma)p(Y^k|sigma_k,k)
//addprior = 0 if prior should not be computed
double logCondPostL(const pReMiuMParams& params,
                    const mcmcModel<pReMiuMParams,
                                    pReMiuMOptions,
                                    pReMiuMData>& model,
                                    const unsigned int& c,
                                    const int addprior){

  const pReMiuMData& dataset = model.dataset();
  const pReMiuMHyperParams& hyperParams = params.hyperParams();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  const string kernelType = model.options().kernelType(); //AR

  double logPrior = 0.0;
  unsigned int nL;
  if(kernelType.compare("SQexponential")==0){    //AR
    nL=3;
  }else{
    nL=4;
  }

  if(model.options().sampleGPmean()){//AR model.options().sampleGPmean()
    unsigned int l = 2;
    if(addprior==1)
      logPrior+=logPdfNormal(params.L(c,l),hyperParams.muL(l),
                             hyperParams.sigmaL(l));
  }else{
    for(unsigned int l=0;l<nL;l++){
      logPrior+=logPdfNormal(params.L(c,l),hyperParams.muL(l),
                             hyperParams.sigmaL(l));
    }
  }

  //AR
  double out = 0;
  if(model.options().sampleGPmean()){//AR model.options().sampleGPmean()
    out = logPrior +  logPYiGivenZiWiLongitudinal_meanGP(params,dataset,nFixedEffects,c,dataset.nSubjects());
    //double a= logPYiGivenZiWiLongitudinal_meanGP(params,dataset,nFixedEffects,c,dataset.nSubjects());
    //double b=logPYiGivenZiWiLongitudinal(params,dataset,nFixedEffects,c,dataset.nSubjects());
    //foutL << " lik Y|f "<< a<< " lik Y|L "<< b<< endl;
  }else{
    out = logPrior +  logPYiGivenZiWiLongitudinal(params,dataset,nFixedEffects,c,dataset.nSubjects());
  }

  return out;

}

//AR logPdfPostMultivariateNormal P(fk|L_k, k, yk)
double logPdfPostMultivariateNormal(const pReMiuMParams& params, const pReMiuMData& dataset, const unsigned int& c){
  //const vector<double>& muGP, const vector<double>& L, const vector<double>& times, const string& kernelType){

  vector<double> times_unique = dataset.times_unique();
  unsigned int nTimes_unique=dataset.nTimes_unique();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();

  VectorXd postM(nTimes_unique);

  unsigned int nSubjects=dataset.nSubjects();
  vector<double> y = dataset.continuousY();
  const string kernelType = dataset.kernelType(); //AR
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> timesk ;
  VectorXd yk;

  int sizek = 0;
  int counter = 0;

  for(unsigned int i=0;i<nSubjects;i++){
    if(params.z(i) == c){
      sizek +=   (tStop[i] - tStart[i] + 1);
    }
  }

  yk.resize(sizek);
  timesk.resize(sizek);

  for(unsigned int i=0;i<nSubjects;i++){
    if(params.z(i) == c){
      for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
        yk(counter+j) = y[tStart[i]-1+j];
        timesk[counter+j] = times[tStart[i]-1+j];

        for(unsigned int b=0;b<nFixedEffects;b++){
          yk(counter+j)-=params.beta(0,b,0,nCategoriesY)*dataset.W(i,b);
        }
        // for(unsigned int b=0;b<nFixedEffects_mix;b++){
        //   yk(counter+j)-=params.beta_mix(0,c,b,0,nCategoriesY)*dataset.W_mix(i,b);
        // }
      }
      counter = counter + tStop[i] - tStart[i] + 1;
    }
  }

  MatrixXd priorCor(sizek, nTimes_unique); // take transpose
  priorCor.setZero(sizek, nTimes_unique);
  MatrixXd invC(sizek, sizek);

  GP_cov_star(priorCor,params.L(c),timesk,times_unique, kernelType);
  GP_cov(invC, params.L(c), timesk, 1, kernelType,1);

  postM = priorCor.transpose() * (invC.inverse() * yk);

  //Computation posterior covariance matrix of GP: postV
  MatrixXd postV(nTimes_unique, nTimes_unique);
  MatrixXd priorCor_star(nTimes_unique, nTimes_unique);
  GP_cov(priorCor_star, params.L(c), times_unique,  1, kernelType,0);
  postV = priorCor_star - priorCor.transpose() * invC.inverse() * priorCor;

  LLT<MatrixXd> lltOfA(postV); // compute the Cholesky decomposition of A
  MatrixXd L = lltOfA.matrixL();
  double logDetV=  2*L.diagonal().array().log().sum();
  MatrixXd postV_inv= L.inverse().transpose()*L.inverse();

  double dmvnorm = 0;
  VectorXd GPmean(nTimes_unique);
  for(int i = 0; i < nTimes_unique; ++i)
    GPmean(i)=params.meanGP(c)[i];

  VectorXd diff = GPmean - postM;
  dmvnorm = -0.5*diff.transpose()*postV_inv*diff - 0.5*nTimes_unique*log(2.0*pi<double>()) - 0.5*logDetV;

  return dmvnorm;
}

//AR logCondPostL function definition: p(L_12)p(f^k|L_12k,k,yk)
double logCondPostL_covGP(const pReMiuMParams& params,
                          const mcmcModel<pReMiuMParams,
                                          pReMiuMOptions,
                                          pReMiuMData>& model,
                                          const unsigned int& c){

  const pReMiuMData& dataset = model.dataset();
  const pReMiuMHyperParams& hyperParams = params.hyperParams();
  double logPrior = 0.0;
  bool ratio_estim = model.options().estim_ratio();

  unsigned int nL;
  if(dataset.kernelType().compare("SQexponential")==0){    //AR
    nL=3;
  }else{
    nL=4;
  }

  for(unsigned int l=0;l<nL;l++){
    if(! ratio_estim || l != 2)
      logPrior+=logPdfNormal(params.L(c,l),hyperParams.muL(l),
                             hyperParams.sigmaL(l));
  }
  double out = logPrior + logPdfMultivariateNormal(params.meanGP(c),params.L(c), dataset.times_unique(), dataset.kernelType());
  //double out = logPrior +  logPYiGivenZiWiLongitudinal(params,dataset,dataset.nFixedEffects(),c,dataset.nSubjects());

  return out;
}


double logCondPostLambdaiPoisson(const pReMiuMParams& params,
                                 const mcmcModel<pReMiuMParams,
                                                 pReMiuMOptions,
                                                 pReMiuMData>& model,
                                                 const unsigned int& i){

  const pReMiuMData& dataset = model.dataset();
  vector<unsigned int> nFixedEffects=dataset.nFixedEffects();
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();
  //Rprintf("Verify that nCategoriesY=1 in this case");

  int zi = params.z(i);
  double meanVal = params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects[0];j++){
    meanVal+=params.beta(0,j,0,nCategoriesY)*dataset.W(i,j);
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   meanVal+=params.beta_mix(zi,j,0,nCategoriesY )*dataset.W_mix(i,j);
  // }
  meanVal+=dataset.logOffset(i);

  return logPYiGivenZiWiPoissonExtraVar(params,dataset,nFixedEffects,zi,i)
    + logPdfNormal(params.lambda(i),meanVal,1.0/sqrt(params.tauEpsilon()));

}


// Evaluation log of conditional density of U_i given U_{-i} and Y for adaptive rejection
void logUiPostPoissonSpatial(const pReMiuMParams& params,
                             const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                             const unsigned int& iSub,
                             const double& x,
                             double* Pt_y1, double* Pt_y2){

  const pReMiuMData& dataset=model.dataset();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();
  //Rprintf("Verify that nCategoriesY=1 in this case");

  double y1;
  double y2;
  int Yi=dataset.discreteY(iSub);
  int zi=params.z(iSub);
  double meanVal=params.theta(zi,0);
  for(unsigned int j=0;j<nFixedEffects;j++){
    meanVal+=params.beta(0,j,0,nCategoriesY)*dataset.W(iSub,j);
  }
  // for(unsigned int j=0;j<nFixedEffects_mix;j++){
  //   meanVal+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(iSub,j);
  // }
  int nNeighi=dataset.nNeighbours(iSub);
  // mean of Ui is mean of Uj where j are the neighbours of i
  double meanUi=0.0;
  for (int j = 0; j<nNeighi; j++){
    unsigned int nj = dataset.neighbours(iSub,j);
    double ucarj = params.uCAR(nj-1);
    meanUi+=ucarj;
  }
  meanUi/=nNeighi;
  y1=Yi*x-exp(dataset.logOffset(iSub)+meanVal+x)-0.5*params.TauCAR()*nNeighi*(x-meanUi)*(x-meanUi);
  // derivative of y1
  y2=Yi-exp(dataset.logOffset(iSub)+meanVal+x)-params.TauCAR()*nNeighi*(x-meanUi);
  *Pt_y1=y1;
  *Pt_y2=y2;
}


// Evaluation log of conditional density of \nu (shape parameter of Weibull for Survival response)
// given \gamma for adaptive rejection sampling
void logNuPostSurvival(const pReMiuMParams& params,
                       const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
                       const unsigned int& cluster,
                       const double& x,
                       double* Pt_y1, double* Pt_y2){

  const pReMiuMData& dataset=model.dataset();
  const pReMiuMHyperParams& hyperParams = params.hyperParams();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nSubjects=dataset.nSubjects();
  const vector<unsigned int> censoring = dataset.censoring();
  vector<double> y = dataset.continuousY();
  const bool weibullFixedShape=model.options().weibullFixedShape();
  unsigned int nCategoriesY = params.nCategoriesY();

  double y1;
  double y2;
  // compute d
  double dCensored = 0;

  if (weibullFixedShape){
    for (unsigned int i=0;i<nSubjects;i++){
      dCensored += censoring[i];
    }
  } else {
    for (unsigned int i=0;i<nSubjects;i++){
      int zi=params.z(i);
      if (zi==(int)cluster) dCensored += censoring[i];
    }
  }
  // sum of yi^nu
  double yNu = 0;
  double yNulogy = 0;
  for (unsigned int i=0;i<nSubjects;i++){
    int zi=params.z(i);
    if (weibullFixedShape){
      double lambda = params.theta(zi,0);
      for(unsigned int j=0;j<nFixedEffects;j++){
        lambda+=params.beta(0,j,0,nCategoriesY)*dataset.W(i,j);
      }
      // for(unsigned int j=0;j<nFixedEffects_mix;j++){
      //   lambda+=params.beta_mix(zi,j, 0, nCategoriesY)*dataset.W_mix(i,j);
      // }
      yNulogy += pow(y[i],x) * log(y[i]) * exp(lambda);
      yNu += pow(y[i],x) * exp(lambda);
    } else {
      if (zi==(int)cluster) {
        double lambda = params.theta(zi,0);
        for(unsigned int j=0;j<nFixedEffects;j++){
          lambda+=params.beta(0,j,0, nCategoriesY)*dataset.W(i,j);
        }
        // for(unsigned int j=0;j<nFixedEffects_mix;j++){
        //   lambda+=params.beta_mix(zi,j,0,nCategoriesY)*dataset.W_mix(i,j);
        // }
        yNulogy += pow(y[i],x) * log(y[i]) * exp(lambda);
        yNu += pow(y[i],x) * exp(lambda);
      }
    }
  }
  // sum of di*log yi
  double dlogY = 0;
  for (unsigned int i=0;i<nSubjects;i++){
    if (weibullFixedShape){
      dlogY += log(y[i]) * censoring[i];
    } else {
      int zi=params.z(i);
      if (zi==(int)cluster) dlogY += log(y[i]) * censoring[i];
    }
  }
  y1=dCensored * log(x) - yNu + x * dlogY + (hyperParams.shapeNu()-1) * log(x) - hyperParams.scaleNu() * x;
  // derivative of y1
  y2=dCensored / x - yNulogy + dlogY + (hyperParams.shapeNu()-1) / x - hyperParams.scaleNu();

  *Pt_y1=y1;
  *Pt_y2=y2;
}



//void EvalHXHPrimaXPoissonSpatial(const pReMiuMParams& params,
//                                 const mcmcModel<pReMiuMParams,pReMiuMOptions,pReMiuMData>& model,
//                                 const unsigned int& iSub,
//                                 const double& x,
//                                 double* Pt_y1, double* Pt_y2){
//
//    const pReMiuMData& dataset=model.dataset();
//    unsigned int nFixedEffects=dataset.nFixedEffects();
//
//    double y1;
//    double y2;
//    y1=-0.5*x*x;
//    y2=-x;
//    *Pt_y1=y1;
//    *Pt_y2=y2;
//}

//AR to sample the GP mean for all unique timepoints
VectorXd Sample_GPmean(pReMiuMParams& params, const pReMiuMData& dataset,
                       const unsigned int c, baseGeneratorType& rndGenerator,
                       const unsigned int init){
// If init == 1 or sizek == 0, f drawn from prior
  vector<double> times_unique = dataset.times_unique();
  unsigned int nTimes_unique=dataset.nTimes_unique();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  VectorXd GPmean = VectorXd::Zero(nTimes_unique);

  // Define a normal random number generator
  randomNormal normRand(0,1);

  VectorXd postM(nTimes_unique);

  unsigned int nSubjects=dataset.nSubjects();
  vector<double> y = dataset.continuousY();
  const string kernelType = dataset.kernelType(); //AR
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> timesk ;
  unsigned int nCategoriesY = params.nCategoriesY();

  VectorXd yk;
  MatrixXd GPSigma;

  int sizek = 0;
  int counter = 0;

  // set sizes based on cluster occupation
  if(init == 0){
    for(unsigned int i=0;i<nSubjects;i++){
      if(params.z(i) == c){
        sizek +=   (tStop[i] - tStart[i] + 1);
      }
    }
  }

  if(sizek > 0 && init == 0){
    yk.resize(sizek);
    timesk.resize(sizek);

    for(unsigned int i=0;i<nSubjects;i++){
      if(params.z(i) == c){
        for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
          yk(counter+j) = y[tStart[i]-1+j];
          timesk[counter+j] = times[tStart[i]-1+j];

          for(unsigned int b=0;b<nFixedEffects;b++){
            yk(counter+j)-=params.beta(0,b,0,nCategoriesY)*dataset.W(i,b);
          }
          // for(unsigned int b=0;b<nFixedEffects_mix;b++){
          //   yk(counter+j)-=params.beta_mix(0,c,b, 0, nCategoriesY)*dataset.W_mix(i,b);
          // }
        }
        counter = counter + tStop[i] - tStart[i] + 1;
      }
    }

    // sort yk and times with corresponding ordering in idx
    std::vector<int> idx(timesk.size());
    int x=0;
    iota(idx.begin(), idx.end(), x++);
    stable_sort(idx.begin(), idx.end(), // sort indexes based on comparing values in times
                [&](int i1,int i2) { return (timesk[i1] < timesk[i2]); });

    stable_sort(timesk.begin(), timesk.end());

    VectorXd yk_sorted = yk;
    for(unsigned int j=0;j<sizek;j++)
      yk_sorted(j)=yk(idx[j]);
    yk=yk_sorted;


    //Computation posterior mean of GP: postM
    MatrixXd priorCor(sizek, nTimes_unique); // K_{}k,*}
    priorCor.setZero(sizek, nTimes_unique);
    MatrixXd invC(sizek, sizek);

    GP_cov_star(priorCor,params.L(c),timesk,times_unique, kernelType);
    // cout << " ---- " <<endl<< " check sans le GP_cov: " << log(invC.determinant())<<endl;
    //double logDetPrecMat= Get_Sigma_inv_GP_cov(invC,params.L(c),timesk,dataset.equalTimes(),times_unique,kernelType);
    //cout << " check sans le GP_cov logDetPrecMat"<<logDetPrecMat<< " vs " << log(invC.determinant())<<endl;

    invC.setZero(sizek, sizek);
    double logDetPrecMat=0;
    GP_cov(invC, params.L(c), timesk, dataset.equalTimes(), kernelType,1);
    //cout << " check sans le GP_cov 2: " << log(invC.determinant())<<endl;
     logDetPrecMat= Get_Sigma_inv_GP_cov(invC,params.L(c),timesk,dataset.equalTimes(),times_unique,kernelType);
    //cout << " check sans le GP_cov logDetPrecMat 2: "<<logDetPrecMat<< " vs " << log(invC.determinant())<<endl
    //    <<" ---- " <<endl;


    if(std::isnan(logDetPrecMat) || isinf(logDetPrecMat)){
      for (unsigned int j=0;j<GPmean.size();j++)
        GPmean(j)=0;

    }else{

      postM = priorCor.transpose() * (invC * yk);


      // foutC << " yk "<< yk<<endl<< " times "<< endl;
      // for(unsigned int j=0;j<yk.size();j++)
      //   foutC << timesk[j]<<" " ;
      //
      //   foutC << endl<<" postM "<< postM<<endl<< " times_unique "<<endl;
      //   for(unsigned int j=0;j<times_unique.size();j++)
      //     foutC << times_unique[j]<<" " ;
      //     foutC << endl;

      //Computation posterior covariance matrix of GP: postV
      MatrixXd postV(nTimes_unique, nTimes_unique);
      MatrixXd priorCor_star(nTimes_unique, nTimes_unique);
      GP_cov(priorCor_star, params.L(c), times_unique,  1, kernelType,0);
      postV = priorCor_star - priorCor.transpose() * invC * priorCor;

      VectorXd random_vector(nTimes_unique);
      for(int i = 0; i < nTimes_unique; ++i){
        random_vector(i) = normRand(rndGenerator);
      }


      //Eigen decomposition of postV
      EigenSolver<MatrixXd> es(postV, true);
      VectorXd eigenvalues(es.eigenvalues().size()) ;

      //bool sdp=true;
      for(int i = 0; i < nTimes_unique; ++i){
        eigenvalues(i) = es.eigenvalues()(i).real() ;

        if(es.eigenvalues()(i).real()<0){
          eigenvalues(i) =0;
          //sdp=false;
        }

        if(es.eigenvalues()(i).imag() != 0 || es.eigenvectors()(i).imag()!=0 ){
          std::cout << c <<" c post complex eigenvalue " <<es.eigenvalues()(i).imag()<< endl;
          std::cout << c <<" c post complex eigenvector " <<es.eigenvectors()(i).imag()<< endl;
        }
      }

      // if(!sdp){
      //   postV = postV +0.001*MatrixXd::Identity(nTimes_unique, nTimes_unique);
      //   EigenSolver<MatrixXd> es2(postV, true);
      //
      //   bool sdp=true;
      //   for(int i = 0; i < nTimes_unique; ++i){
      //     if(es2.eigenvalues()(i).real()<0)
      //       sdp=false;
      //     if(es2.eigenvalues()(i).imag() != 0 || es.eigenvectors()(i).imag()!=0 ){
      //       fout << c <<" c post complex eigenvalue bis " <<es2.eigenvalues()(i).imag()<< endl;
      //       fout << c <<" c post complex eigenvector bis " <<es2.eigenvectors()(i).imag()<< endl;
      //     }
      //   }
      //   GPSigma= es2.eigenvectors().real() * es2.eigenvalues().real().cwiseSqrt().asDiagonal();
      //}else{
      GPSigma= es.eigenvectors().real() * eigenvalues.cwiseSqrt().asDiagonal();
      //}


      //LLT<MatrixXd> llt;
      //MatrixXd B = (llt.compute(covMat)).matrixL();
      //V = meanVec + B*V;


      GPmean = postM  +  GPSigma* random_vector;

      if(std::isnan(GPmean(1))){
        std::cout << c <<" GPmean " << endl <<GPmean(0)<< endl;
        std::cout << c <<" L " << params.L(c,0)<<" L " << params.L(c,1) <<" L " << params.L(c,2)<<endl;
        std::cout << c <<" postM " << endl <<postM<< endl;
        for(int i = 0; i < nTimes_unique; ++i){
          std::cout << c <<" eigenvalue" <<es.eigenvalues()(i)<< endl<< endl;
        }
        for(int i = 0; i < nTimes_unique; ++i){
          std::cout << c <<" eigenvalue2" <<es.eigenvalues()(i)<< endl<< endl;
        }
      }
    }
  }else{ // init == 1 || sizek==0
    MatrixXd priorCor(nTimes_unique, nTimes_unique);
    GP_cov(priorCor, params.L(c), times_unique, 1, kernelType,0);
    EigenSolver<MatrixXd> es(priorCor, true);

    //bool sdp=true;
    VectorXd eigenvalues(es.eigenvalues().size());
    for(int i = 0; i < nTimes_unique; ++i){
      eigenvalues(i)=es.eigenvalues()(i).real();
      if(es.eigenvalues()(i).real()<0){
        //sdp=false;
        eigenvalues(i)=0;
      }
      if(es.eigenvalues()(i).imag() != 0 || es.eigenvectors()(i).imag()!=0 ){
        std::cout <<"eigenvalue"<<endl;
        cout <<"eigenvalue"<<endl;
        //fout << c <<" c prior eigenvalue" <<es.eigenvalues()(i).imag()<< endl;
        //fout << c <<" c prior eigenvector" <<es.eigenvectors()(i).imag()<< endl;
      }
    }


    // if(!sdp){
    //   priorCor = priorCor +0.001*MatrixXd::Identity(nTimes_unique, nTimes_unique);
    //   EigenSolver<MatrixXd> es2(priorCor, true);
    //
    //   for(int i = 0; i < nTimes_unique; ++i){
    //     if(es2.eigenvalues()(i).imag() != 0 || es.eigenvectors()(i).imag()!=0 ){
    //       fout <<"eigenvalue"<<endl;
    //       //fout << c <<" c post eigenvalue bis " <<es2.eigenvalues()(i).imag()<< endl;
    //       //fout << c <<" c post eigenvector bis " <<es2.eigenvectors()(i).imag()<< endl;
    //     }
    //   }
    //   GPSigma= es2.eigenvectors().real() * es2.eigenvalues().real().cwiseSqrt().asDiagonal();
    // }else{

    GPSigma= es.eigenvectors().real() * eigenvalues.cwiseSqrt().asDiagonal();
    //}

    VectorXd random_vector(nTimes_unique);
    for(int i = 0; i < nTimes_unique; ++i)
      random_vector(i) = normRand(rndGenerator);

    GPmean =  GPSigma * random_vector;

  }
  return(GPmean);
}



//AR logPdfPostMultivariateNormal P(f_k|L_k, k, yk)
double logPdfPostMultivariateNormal(pReMiuMParams& params, const pReMiuMData& dataset, const unsigned int c){
    //const vector<double>& muGP, const vector<double>& L, const vector<double>& times, const string& kernelType){

  vector<double> times_unique = dataset.times_unique();
  unsigned int nTimes_unique=dataset.nTimes_unique();
  unsigned int nFixedEffects=dataset.nFixedEffects(0);
  //unsigned int nFixedEffects_mix=dataset.nFixedEffects_mix();
  unsigned int nCategoriesY = params.nCategoriesY();

  VectorXd postM(nTimes_unique);

  unsigned int nSubjects=dataset.nSubjects();
  vector<double> y = dataset.continuousY();
  const string kernelType = dataset.kernelType(); //AR
  vector<double> times = dataset.times();
  vector<int> tStart = dataset.tStart();
  vector<int> tStop = dataset.tStop();
  vector<double> timesk ;
  VectorXd yk;

  int sizek = 0;
  int counter = 0;

  for(unsigned int i=0;i<nSubjects;i++){
    if(params.z(i) == c){
      sizek +=   (tStop[i] - tStart[i] + 1);
    }
  }

  yk.resize(sizek);
  timesk.resize(sizek);

  for(unsigned int i=0;i<nSubjects;i++){
    if(params.z(i) == c){
      for(unsigned int j=0;j<tStop[i]-tStart[i]+1;j++){
        yk(counter+j) = y[tStart[i]-1+j];
        timesk[counter+j] = times[tStart[i]-1+j];

        for(unsigned int b=0;b<nFixedEffects;b++){
          yk(counter+j)-=params.beta(0,b,0,nCategoriesY)*dataset.W(i,b);
        }
        // for(unsigned int b=0;b<nFixedEffects_mix;b++){
        //   yk(counter+j)-=params.beta_mix(0,c, b,0,nCategoriesY)*dataset.W_mix(i,b);
        // }
      }
      counter = counter + tStop[i] - tStart[i] + 1;
    }
  }

  MatrixXd priorCor(sizek, nTimes_unique); // take transpose
  priorCor.setZero(sizek, nTimes_unique);
  MatrixXd invC(sizek, sizek);

  GP_cov_star(priorCor,params.L(c),timesk,times_unique, kernelType);
  GP_cov(invC, params.L(c), timesk, 1, kernelType,1);

  postM = priorCor.transpose() * (invC.inverse() * yk);

  //Computation posterior covariance matrix of GP: postV
  MatrixXd postV(nTimes_unique, nTimes_unique);
  MatrixXd priorCor_star(nTimes_unique, nTimes_unique);
  GP_cov(priorCor_star, params.L(c), times_unique,  1, kernelType,0);
  postV = priorCor_star - priorCor.transpose() * invC.inverse() * priorCor;

  LLT<MatrixXd> lltOfA(postV); // compute the Cholesky decomposition of A
  MatrixXd L = lltOfA.matrixL();
  double logDetV=  2*L.diagonal().array().log().sum();
  MatrixXd postV_inv =L.inverse().transpose()*L.inverse();

  double dmvnorm = 0;
  VectorXd GPmean(nTimes_unique);
  for(int i = 0; i < nTimes_unique; ++i)
    GPmean(i)=params.meanGP(c)[i];

  VectorXd diff = GPmean - postM;
  dmvnorm = -0.5*diff.transpose()*postV_inv*diff - 0.5*nTimes_unique*log(2.0*pi<double>()) - 0.5*logDetV;

  return dmvnorm;
}


#endif /* DIPBACMODEL_H_ */

