/// \file PReMiuMData.h
/// \author David Hastie
/// \brief Header file for data specification for PReMiuMpp

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

#ifndef DIPBACDATA_H_
#define DIPBACDATA_H_

// Standard includes
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<stdexcept>


#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>
#include<Eigen/CholmodSupport>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using std::vector;
using std::ifstream;
using std::string;

/// \class pReMiuMData PReMiuMData.h "PReMiuMModel.h"
/// \brief A class for PReMiuMpp Data
class pReMiuMData{

	public:
		/// \brief Default constructor //RJ add _nTimes
		pReMiuMData():  _nSubjects(0), _nOutcomes(0), _nTimes(0),_nTimes_unique(0), _nCovariates(0), _nCategoriesY(0), _nPredictSubjects(0) {};

		/// \brief Default destructor
		~pReMiuMData(){};

		/// \brief Return the number of subjects
		unsigned int nSubjects() const{
			return _nSubjects;
		}

		/// \brief Return the number of subjects
		unsigned int& nSubjects(){
			return _nSubjects;
		}

		/// \brief Set the number of subjects
		void nSubjects(const unsigned int& nSubj){
			_nSubjects = nSubj;
		}



		/// \brief Return the number of outcomes
		unsigned int nOutcomes() const{
		  return _nOutcomes;
		}

		/// \brief Return the number of outcomes
		unsigned int& nOutcomes(){
		  return _nOutcomes;
		}

		/// \brief Set the number of outcomes
		void nOutcomes(const unsigned int& n){
		  _nOutcomes = n;
		}


		unsigned int size() const{
			return _nSubjects;
		}

		void nTimes_m(const unsigned int& m, const unsigned int& nm) {
		  _nTimes_m[m] = nm;
		}

		vector<unsigned int> nTimes_m()const{
		  return _nTimes_m;
		}
		vector<unsigned int>& nTimes_m(){
		  return _nTimes_m;
		}

		/// \brief Set the output vector
		void nTimes_m(const vector<unsigned int>& yVec){
		  _nTimes_m.clear();
		  _nTimes_m.resize(yVec.size());
		  _nTimes_m.insert(_nTimes_m.begin(),yVec.begin(),yVec.end());
		}

		/// \brief Return the output value for the ith subject
		unsigned int nTimes_m(const unsigned int& i) const{
		  if(i>_nOutcomes){
		    throw std::range_error("y subscript i out of range");
		  }
		  return _nTimes_m[i];
		}


		//RJ handling functions for _nTimes
		unsigned int nTimes() const{
			return _nTimes;
		}
		unsigned int& nTimes(){
			return _nTimes;
		}
		void nTimes(const unsigned int& nT){
			_nTimes = nT;
		}
		unsigned int sizeT() const{
			return _nTimes;
		}

		//AR handling functions for _nTimes_unique
		unsigned int nTimes_unique() const{
		  return _nTimes_unique;
		}
		unsigned int& nTimes_unique(){
		  return _nTimes_unique;
		}
		void nTimes_unique(const unsigned int& nT){
		  _nTimes_unique = nT;
		}
		unsigned int sizeT_unique() const{
		  return _nTimes_unique;
		}

		/// \brief Return the number of covariates
		unsigned int nCovariates() const{
			return _nCovariates;
		}

		/// \brief Return the number of covariates
		unsigned int& nCovariates(){
			return _nCovariates;
		}

		/// \brief Set the number of covariates
		void nCovariates(const unsigned int& nCov){
			_nCovariates = nCov;
		}

		/// \brief Return the number of discrete covariates
		unsigned int nDiscreteCovs() const{
			return _nDiscreteCovs;
		}

		/// \brief Return the number of discrete covariates
		unsigned int& nDiscreteCovs(){
			return _nDiscreteCovs;
		}

		/// \brief Set the number of discrete covariates
		void nDiscreteCovs(const unsigned int& nDiscrCovs){
			_nDiscreteCovs = nDiscrCovs;
		}
		/// \brief Return the number of continuous covariates
		unsigned int nContinuousCovs() const{
			return _nContinuousCovs;
		}

		/// \brief Return the number of continuous covariates
		unsigned int& nContinuousCovs(){
			return _nContinuousCovs;
		}


		/// \brief Set the number of continuous covariates
		void nContinuousCovs(const unsigned int& nContCovs){
			_nContinuousCovs = nContCovs;
		}

		/// \brief Return the number of fixed effectss
		unsigned int& nFixedEffects_mix(const unsigned int& m){
		  return _nFixedEffects_mix[m];
		}

		unsigned int nFixedEffects_mix(const unsigned int& m) const{
			return _nFixedEffects_mix[m];
		}

		void nFixedEffects_mix(const unsigned int& m, const unsigned int& nConf){
		  _nFixedEffects_mix[m] = nConf;
		}


		vector<unsigned int>& nFixedEffects_mix() {
		  return _nFixedEffects_mix;
		}

		vector<unsigned int> nFixedEffects_mix() const{
		  return _nFixedEffects_mix;
		}

		unsigned int& nFixedEffects(const unsigned int& m){
			return _nFixedEffects[m];
		}

		unsigned int nFixedEffects(const unsigned int& m)const{
		  return _nFixedEffects[m];
		}

		void nFixedEffects(const unsigned int& m, const unsigned int& nConf){
			_nFixedEffects[m] = nConf;
		}

		vector<unsigned int>& nFixedEffects() {
		  return _nFixedEffects;
		}

		vector<unsigned int> nFixedEffects() const{
		  return _nFixedEffects;
		}

		unsigned int& nRandomEffects(const unsigned int& m){
		  return _nRandomEffects[m];
		}

		unsigned int nRandomEffects(const unsigned int& m)const{
		  return _nRandomEffects[m];
		}

		void nRandomEffects(const unsigned int& m, const unsigned int& nConf){
		  _nRandomEffects[m] = nConf;
		}

		vector<unsigned int>& nRandomEffects() {
		  return _nRandomEffects;
		}

		vector<unsigned int> nRandomEffects() const{
		  return _nRandomEffects;
		}
		//RJ handling functions for equalTimes
		unsigned int equalTimes() const{
			return _equalTimes;
		}
		unsigned int& equalTimes(){
			return _equalTimes;
		}
		void equalTimes(const unsigned int& equalTimes){
			_equalTimes = equalTimes;
		}

		/// \brief Return the number of categories
		unsigned int nCategoriesY() const{
			return _nCategoriesY;
		}

		/// \brief Return the number of categories
		unsigned int& nCategoriesY(){
			return _nCategoriesY;
		}

		/// \brief Set the number of categories
		void nCategoriesY(const unsigned int& nCat){
			_nCategoriesY = nCat;
		}

		/// \brief Return the number of subjects
		unsigned int nPredictSubjects() const{
			return _nPredictSubjects;
		}

		/// \brief Return the number of subjects
		unsigned int& nPredictSubjects(){
			return _nPredictSubjects;
		}

		/// \brief Set the number of subjects
		void nPredictSubjects(const unsigned int& nPredSubj){
			_nPredictSubjects = nPredSubj;
		}



		/// \brief Return the vector of the number of categories
		vector<unsigned int> nCategories() const{
			return _nCategories;
		}

		/// \brief Return the vector of the number of categories
		vector<unsigned int>& nCategories(){
			return _nCategories;
		}


		/// \brief Set the vector of the number of categories
		void nCategories(const vector<unsigned int>& nCats){
			_nCategories.clear();
			_nCategories.resize(nCats.size());
			_nCategories.insert(_nCategories.begin(),nCats.begin(),nCats.end());
		}

		/// \brief Return the number of categories for covariate j
		unsigned int nCategories(const unsigned int& j) const{
			if(j>_nCovariates){
				throw std::range_error("nCategories subscript j out of range");
			}
			return _nCategories[j];
		}

		/// \brief Return the vector of the covariate names
		vector<string> covariateNames() const{
			return _covariateNames;
		}

		/// \brief Return the vector of the covariate names
		vector<string>& covariateNames(){
			return _covariateNames;
		}


		/// \brief Set the vector of the covariate names
		void covariateNames(const vector<string>& covNames){
			_covariateNames.clear();
			_covariateNames.resize(covNames.size());
			_covariateNames.insert(_covariateNames.begin(),covNames.begin(),covNames.end());
		}

		/// \brief Return name for covariate j
		const string& covariateNames(const unsigned int& j) const{
			return _covariateNames[j];
		}

		vector<vector<string>> fixedEffectNames() const {
		  return _fixedEffectNames;
		}

		/// \brief Return the vector of the fixed effects names
		vector<vector<string>>& fixedEffectNames() {
			return _fixedEffectNames;
		}

		/// \brief Set the vector of the fixed effects names
		void fixedEffectNames(const unsigned int& m, const vector<string>& fixEffNames){
		  _fixedEffectNames[m].clear();
		  _fixedEffectNames[m].resize(fixEffNames.size());
		  _fixedEffectNames[m].insert(_fixedEffectNames[m].begin(),fixEffNames.begin(),fixEffNames.end());
		}

		/// \brief Return name for fixed effects j
		const string& fixedEffectNames(const unsigned int& m, const unsigned int& j) const{
		  return _fixedEffectNames[m][j];
		}

		vector<string>& fixedEffectNames(const unsigned int& m){
		  return _fixedEffectNames[m];
		}

		/// \brief Return the vector of the fixed effects names
	  vector<vector<string>> fixedEffectNames_mix() const {
		  return _fixedEffectNames_mix;
		}

		/// \brief Return the vector of the fixed effects names
		vector<vector<string>>& fixedEffectNames_mix()  {
		  return _fixedEffectNames_mix;
		}

		/// \brief Set the vector of the fixed effects names
		void fixedEffectNames_mix(const unsigned int& m, const vector<string>& fixEffNames){
			_fixedEffectNames_mix[m].clear();
			_fixedEffectNames_mix[m].resize(fixEffNames.size());
			_fixedEffectNames_mix[m].insert(_fixedEffectNames_mix[m].begin(),fixEffNames.begin(),fixEffNames.end());
		}

		/// \brief Return name for fixed effects j
		const string& fixedEffectNames_mix(const unsigned int& m, const unsigned int& j) const{
			return _fixedEffectNames_mix[m][j];
		}

		vector<string> fixedEffectNames_mix(const unsigned int& m) const{
		  return _fixedEffectNames_mix[m];
		}

		vector<string>& fixedEffectNames_mix(const unsigned int& m){
		  return _fixedEffectNames_mix[m];
		}

		vector<vector<string>> randomEffectNames() const{
		  return _randomEffectNames;
		}
		vector<vector<string>>& randomEffectNames() {
		  return _randomEffectNames;
		}

		 void randomEffectNames(const unsigned int& m, const vector<string>& fixEffNames){
		   _randomEffectNames[m].clear();
		   _randomEffectNames[m].resize(fixEffNames.size());
		   _randomEffectNames[m].insert(_randomEffectNames[m].begin(),fixEffNames.begin(),fixEffNames.end());
		 }

		string randomEffectNames(const unsigned int& m, const unsigned int& j) const{
		   return _randomEffectNames[m][j];
		 }
		string& randomEffectNames(const unsigned int& m, const unsigned int& j) {
		  return _randomEffectNames[m][j];
		}

		vector<string>& randomEffectNames(const unsigned int& m) {
		  return _randomEffectNames[m];
		}

		vector<string> randomEffectNames(const unsigned int& m) const{
		  return _randomEffectNames[m];
		}
		/// \brief Return the outcome model type
		const string& outcomeType() const{
			return _outcomeType;
		}

		/// \brief Set the outcome model type
		void outcomeType(const string& outType){
			_outcomeType=outType;
		}

		/// \brief Return the kernel type //AR
		const string& kernelType() const{
		  return _kernelType;
		}

		/// \brief Set the kernel type
		void kernelType(const string& outType){
		  _kernelType=outType;
		}

		/// \brief Return the outcome model type
		const string& covariateType() const{
			return _covariateType;
		}

		/// \brief Set the outcome model type
		void covariateType(const string& covType){
			_covariateType=covType;
		}

		/// \brief Return the output vector
		const vector<unsigned int>& discreteY() const{
			return _discreteY;
		}

		/// \brief Return the output vector
		vector<unsigned int>& discreteY() {
			return _discreteY;
		}

		/// \brief Set the output vector
		void discreteY(const vector<unsigned int>& yVec){
			_discreteY.clear();
			_discreteY.resize(yVec.size());
			_discreteY.insert(_discreteY.begin(),yVec.begin(),yVec.end());
		}

		/// \brief Return the output value for the ith subject
		unsigned int discreteY(const unsigned int& i) const{
			if(i>_nSubjects){
				throw std::range_error("y subscript i out of range");
			}
			return _discreteY[i];
		}

		/// \brief Return the output vector
		const vector<double>& continuousY() const{
			return _continuousY;
		}
		/// \brief Return the output vector
		vector<double>& continuousY() {
			return _continuousY;
		}
		/// \brief Set the output vector
		void continuousY(const vector<double>& yVec){
			_continuousY.clear();
			_continuousY.resize(yVec.size());
			_continuousY.insert(_continuousY.begin(),yVec.begin(),yVec.end());
		}
		/// \brief Return the output value for the ith subject
		double continuousY(const unsigned int& i) const{
			if(i>_nTimes){//RJ was i>_nSubjects; now i>_nTimes
				throw std::range_error("y subscript i out of range");
			}
			return _continuousY[i];
		}

		//RJ handling functions for _times
		vector<double>& times() {
			return _times;
		}
		const vector<double>& times() const{
			return _times;
		}
		void times(const vector<double>& timepoints){
			_times.clear();
			_times.resize(timepoints.size());
			_times.insert(_times.begin(),timepoints.begin(),timepoints.end());
		}
		double times(const unsigned int& i)const{
			if(i>_nTimes){
				throw std::range_error("t subscript i out of range.\n");
			}
			return _times[i];
		}

		//AR handling functions for _times_unique
		vector<double>& times_unique() {
		  return _times_unique;
		}
		const vector<double>& times_unique() const{
		  return _times_unique;
		}
		void times_unique(const vector<double>& timepoints){
		  _times_unique.clear();
		  _times_unique.resize(timepoints.size());
		  _times_unique.insert(_times_unique.begin(),timepoints.begin(),timepoints.end());
		}
		double times_unique(const unsigned int& i)const{
		  if(i>_nTimes_unique){
		    throw std::range_error("t subscript i out of range.\n");
		  }
		  return _times_unique[i];
		}

		//AR handling functions for _times_corr
		vector<double>& times_corr() {
		  return _times_corr;
		}
		const vector<double>& times_corr() const{
		  return _times_corr;
		}
		void times_corr(const vector<double>& timepoints){
		  _times_corr.clear();
		  _times_corr.resize(timepoints.size());
		  _times_corr.insert(_times_corr.begin(),timepoints.begin(),timepoints.end());
		}
		double times_corr(const unsigned int& i)const{
		  if(i>_nTimes){
		    throw std::range_error("t subscript i out of range.\n");
		  }
		  return _times_corr[i];
		}
		vector<int>& tStart() {
			return _tStart;
		}
		const vector<int>& tStart() const{
			return _tStart;
		}
		void tStart(const vector<int>& indices){
			_tStart.clear();
			_tStart.resize(indices.size());
			_tStart.insert(_tStart.begin(),indices.begin(),indices.end());
		}
		int tStart(const unsigned int& i)const{
			if(i>_nSubjects){
				throw std::range_error("index subscript i out of range.\n");
			}
			return _tStart[i];
		}
		vector<int>& tStop() {
			return _tStop;
		}
		const vector<int>& tStop() const{
			return _tStop;
		}
		void tStop(const vector<int>& indices){
			_tStop.clear();
			_tStop.resize(indices.size());
			_tStop.insert(_tStop.begin(),indices.begin(),indices.end());
		}
		int tStop(const unsigned int& i)const{
			if(i>_nSubjects){
				throw std::range_error("index subscript i out of range.\n");
			}
			return _tStop[i];
		}

		/// \brief Return the covariate matrix
		const vector<vector<int> >& discreteX() const{
			return _discreteX;
		}

		/// \brief Return the covariate matrix
		vector<vector<int> >& discreteX(){
			return _discreteX;
		}

		/// \brief Return the jth covariate for subject i
		int discreteX(const unsigned int& i,const unsigned int& j) const{
			return _discreteX[i][j];
		}

		/// \brief Set the jth covariate for subject i
		void discreteX(const unsigned int& i,const unsigned int& j,const int& x){
			_discreteX[i][j]=x;
		}

		/// \brief Return the covariate matrix
		const vector<vector<double> >& continuousX() const{
			return _continuousX;
		}

		/// \brief Return the covariate matrix
		vector<vector<double> >& continuousX(){
			return _continuousX;
		}

		/// \brief Return the jth covariate for subject i
		double continuousX(const unsigned int& i,const unsigned int& j) const{
			return _continuousX[i][j];
		}

		/// \brief Set the jth covariate for subject i
		void continuousX(const unsigned int& i,const unsigned int& j,const double& x){
			_continuousX[i][j]=x;
		}

		/// \brief Return the missing covariate matrix
		const vector<vector<bool> >& missingX() const{
			return _missingX;
		}

		/// \brief Return the missing covariate matrix
		vector<vector<bool> >& missingX(){
			return _missingX;
		}

		/// \brief Return the jth covariate for subject i
		bool missingX(const unsigned int& i,const unsigned int& j) const{
			return _missingX[i][j];
		}

		/// \brief Return the number of covariates not missing for subject i
		unsigned int nContinuousCovariatesNotMissing(const unsigned int& i) const{
			return _nContinuousCovariatesNotMissing[i];
		}

		/// \brief Return the number of covariates not missing for each subject
		vector<unsigned int>& nContinuousCovariatesNotMissing(){
			return _nContinuousCovariatesNotMissing;
		}

		/// \brief Return the fixed effects matrix
		const vector<vector<double> >& W() const{
			return _W;
		}

		/// \brief Return the fixed effects matrix
		vector<vector<double> >& W(){
			return _W;
		}

		/// \brief Return the jth covariate for subject i
		double W(const unsigned int& i,const unsigned int& j) const{
			return _W[i][j];
		}

		const vector<vector<double> >& W_mix() const{
		  return _W_mix;
		}

		/// \brief Return the fixed effects matrix
		vector<vector<double> >& W_mix(){
		  return _W_mix;
		}

		/// \brief Return the jth covariate for subject i
		double W_mix(const unsigned int& i,const unsigned int& j) const{
		  return _W_mix[i][j];
		}

		/// \brief Return the jth covariate for subject i and marker m
		vector<MatrixXd>& W_RE(){
		  return _W_RE;
		}
		const vector<MatrixXd>& W_RE()const {
		  return _W_RE;
		}

		/// \brief Return the jth covariate for subject i and marker m
		double W_RE(const unsigned int& m, const unsigned int& i,const unsigned int& j) const{
		  return _W_RE[m](i,j);
		}

		const MatrixXd& W_RE(const unsigned int& m) const{
		  return _W_RE[m];
		}

		/// \brief Return the fixed effects matrix
		MatrixXd& W_RE(const unsigned int& m){
		  return _W_RE[m];
		}

//
// 		void W_RE(const unsigned int& i,const unsigned int& j, const double& temp){
// 		   _W_RE(i,j)=temp;
// 		}
//
// 		// const Matrix2d& W_RE() const{
// 		//   return _W_RE;
// 		// }
//
// 		MatrixXd W_RE() const{
// 		  return _W_RE;
// 		}
//
 		VectorXd W_RE(const unsigned int& m, const unsigned int& i) const{
 		  return _W_RE[m].row(i);
 		}

		// MatrixXd W_RE(const unsigned int& i) {
		//
		//   MatrixXd r(_tStop[i]-_tStart[i], _W_RE.cols());
		//   for(unsigned int j=0;j<_tStop[i]-_tStart[i]+1;j++){
		//     r.row(j) << _W_RE.row(_tStop[i]+j);
		//   }
		//   return r;
		// }


		MatrixXd W_RE(const unsigned int& m,const unsigned int& i,const unsigned int& j,const unsigned int& n,const unsigned int& p) const{
		  return _W_RE[m].block(i, j, n, p);
		}


		/// \brief Return the fixed effects matrix for LME
		const MatrixXd& W_LME(const unsigned int& m) const{
		  return _W_LME[m];
		}

		/// \brief Return the fixed effects matrix for LME
		MatrixXd& W_LME(const unsigned int& m){
		  return _W_LME[m];
		}

		/// \brief Return the jth covariate for subject i for LME
		double W_LME(const unsigned int& m, const unsigned int& i,const unsigned int& j) const{
		  return _W_LME[m](i,j);
		}

		vector<MatrixXd>& W_LME(){
		  return _W_LME;
		}
		const vector<MatrixXd>& W_LME()const {
		  return _W_LME;
		}

		const MatrixXd& W_LME_mix(const unsigned int& m) const{
		  return _W_LME_mix[m];
		}

		/// \brief Return the fixed effects matrix for LME
		MatrixXd& W_LME_mix(const unsigned int& m){
		  return _W_LME_mix[m];
		}

		/// \brief Return the jth covariate for subject i for LME
		double W_LME_mix(const unsigned int& m, const unsigned int& i,const unsigned int& j) const{
		  return _W_LME_mix[m](i,j);
		}

		vector<MatrixXd>& W_LME_mix(){
		  return _W_LME_mix;
		}
		const vector<MatrixXd>& W_LME_mix()const {
		  return _W_LME_mix;
		}

		/// \brief Return the fixed effects matrix


		/// \brief Return the logOffset vector
		const vector<double>& logOffset() const{
			return _logOffset;
		}

		/// \brief Return the logOffset vector
		vector<double>& logOffset(){
			return _logOffset;
		}

		/// \brief Return the logOffset for subject i
		double logOffset(const unsigned int& i) const{
			return _logOffset[i];
		}

		/// \brief Return the n vector
		const vector<unsigned int>& nTrials() const{
			return _nTrials;
		}

		/// \brief Return the n vector
		vector<unsigned int>& nTrials(){
			return _nTrials;
		}

		/// \brief Return n for subject i
		unsigned int nTrials(const unsigned int& i) const{
			return _nTrials[i];
		}


		/// \brief Return the n vector for Survival data
		const vector<unsigned int>& censoring() const{
			return _censoring;
		}

		/// \brief Return the n vector for Survival data
		vector<unsigned int>& censoring(){
			return _censoring;
		}

		/// \brief Return n for subject i for Survival data
		unsigned int censoring(const unsigned int& i) const{
			return _censoring[i];
		}
		/// \brief Return the vector of lists of neighbours
		const vector<vector<unsigned int> >& neighbours() const{
			return _neighbours;
		}

		/// \brief Return the vector of lists of neighbours
		vector<vector<unsigned int> >& neighbours() {
			return _neighbours;
		}

		/// \brief Return the list of neighbours for subject i
		vector<unsigned int> neighbours(const unsigned int& i) const{
			return _neighbours[i];
		}

		/// \brief Return the j-th neighbour for subject i
		unsigned int neighbours(const unsigned int& i, const unsigned& j) const{
			return _neighbours[i][j];
		}

		/// \brief Set the vector of lists of neighbours
		void neighbours(const vector<vector<unsigned int> >&  neighvec){
			_neighbours = neighvec;
		}

		/// \brief Set the list of neighbours for subject i
		void neighbours(const vector<unsigned int>&  neighvec, const unsigned int& i){
			_neighbours[i] = neighvec;
		}

		/// \brief Return the number of neighbours for all subjects
		const vector<unsigned int>& nNeighbours() const{
			return _nNeighbours;
		}

		/// \brief Return the number of neighbours for all subjects
		vector<unsigned int>& nNeighbours(){
			return _nNeighbours;
		}

		/// \brief Return the number of neighbours for subject i
		unsigned int nNeighbours(const unsigned int& i) const{
			return _nNeighbours[i];
		}

		/// \brief Set the number of neighbours for all subjects
		void nNeighbours(const vector<unsigned int>& nNeigh){
			_nNeighbours=nNeigh;
		}

		/// \brief Set the number of neighbours for subject i
		void nNeighbours(const unsigned int& nNeigh, const unsigned int& i ){
			_nNeighbours[i]=nNeigh;
		}

		/// \brief Set includeCAR
		void includeCAR(const bool& incl ){
			_includeCAR=incl;
		}

		/// \brief return includeCAR
		bool& includeCAR(){
			return _includeCAR;
		}

		/// \brief return includeCAR
		bool includeCAR() const{
			return _includeCAR;
		}

	private:
		/// \brief The number of subjects
		unsigned int _nSubjects;
	  unsigned int _nOutcomes;
		//RJ declare _nTimes
		unsigned int _nTimes;
		//AR declare _nTimes_unique
		unsigned int _nTimes_unique;
		/// \brief The number of covariates
		unsigned int _nCovariates;

		/// \brief The number of discrete covariates
		unsigned int _nDiscreteCovs;

		/// \brief The number of continuous covariates
		unsigned int _nContinuousCovs;

		/// \brief The number of fixed effects covariates
		vector<unsigned int> _nFixedEffects;
		vector<unsigned int> _nFixedEffects_mix;
		vector<unsigned int> _nRandomEffects;

		//RJ indicator: 0, unequal times; 1, equal times
		unsigned int _equalTimes;

		/// \brief The number of categories for outcome discreteY when outcome is categorical
		unsigned int _nCategoriesY;

		/// \brief The number of subjects we are making predictions for
		unsigned int _nPredictSubjects;

		/// \brief A vector of the number of categories for each covariate
		vector<unsigned int> _nCategories;

		vector<unsigned int> _nTimes_m; //If yModel=LME

		/// \brief A string describing the model for y
		string _outcomeType;

		/// \brief A string describing the kernel //AR
		string _kernelType;

		/// \brief A string describing the model for X
		string _covariateType;


		/// \brief A vector of the output variables
		vector<unsigned int> _discreteY;

		/// \brief A vector of the output variables
		vector<double> _continuousY;
		//RJ declare _times, _tStart, _tStop
		vector<double> _times;
		vector<double> _times_corr; //AR
		vector<double> _times_unique; //AR
		vector<int> _tStart;
		vector<int> _tStop;


		/// \brief A matrix (vector of vectors) of the covariate data
		/// \note this is a signed int because missing values are typically stored
		/// as negative values
		vector<vector<int> > _discreteX;

		/// \brief A matrix (vector of vectors) of the covariate data
		/// \note this is a signed int because missing values are typically stored
		/// as negative values
		vector<vector<double> > _continuousX;

		/// \brief A vector of covariate names
		vector<string> _covariateNames;

		/// \brief A matrix (vector of vectors) of where there are missing
		/// covariate values
		vector<vector<bool> > _missingX;

		/// \brief A matrix of the number of non missing covariates for each subject
		vector<unsigned int> _nContinuousCovariatesNotMissing;

		/// \brief A matrix of the fixed effects covariates
		/// \note This may need to changed to be signed or double
		vector<vector<double> > _W;
		vector<vector<double> > _W_mix;
		//MatrixXd _W_RE;
		//MatrixXd _W_LME;
		//MatrixXd _W_LME_mix;
		vector<MatrixXd> _W_RE;
		vector<MatrixXd> _W_LME;
		vector<MatrixXd> _W_LME_mix; //m x [i x j]

		/// \brief A vector of fixed effects names
		vector<vector<string>> _fixedEffectNames;
		vector<vector<string>> _fixedEffectNames_mix;
		vector<vector<string>> _randomEffectNames;

		/// \brief A vector of logOffsets (only used in the Poisson model)
		vector<double> _logOffset;

		/// \brief A vector of n for each individual (only used in the Binomial model)
		vector<unsigned int> _nTrials;

		/// \brief A vector of n for each individual (only used in the survival model)
		vector<unsigned int> _censoring;

		/// \brief A vector of vector of neighbours for each subject
		vector<vector<unsigned int> > _neighbours;

		/// \brief A containing the number of neighbours for each subject
		vector<unsigned int> _nNeighbours;

		/// \brief Is the CAR term is included
		bool _includeCAR;
};


#endif //DIPBACDATA_H_
