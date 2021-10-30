/// \file distribution.h
/// \author David Hastie and Silvia Liverani
/// \date 19 Mar 2012
/// \brief Header file to define distributions

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


#ifndef DISTRIBUTION_H_
#define DISTRIBUTION_H_

#include<cmath>
#include<string>
#include <vector>
#include<boost/math/distributions/normal.hpp>
#include<boost/math/special_functions/gamma.hpp>
#include<boost/math/constants/constants.hpp>

#include<Eigen/Core>
#include<Eigen/Cholesky>
#include<Eigen/LU>
#include <Eigen/Dense>

#include<Math/mathfunctions.h>

using boost::math::lgamma;
using std::string;
using std::vector;

using namespace Eigen;
using namespace boost::math::constants;

double logPdfBernoulli(const unsigned int& x,const double& p){
	return x*log(p)+(1-x)*log(1.0-p);
}

double logPdfBeta(const double& x,const double& a,const double& b){
	return lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(x)+(b-1)*log(1-x);
}

double logPdfBinomial(const unsigned int& x, const unsigned int& n, const double& p){
	return (double)x*log(p)+(double)(n-x)*log(1.0-p)
        		  + lgamma((double)(n+1)) - lgamma((double)(x+1)) - lgamma((double)(n-x+1));
}

double logPdfDirichlet(const vector<double>& x, const vector<double>& alpha,const bool& xOnLogScale){

	double out = 0.0;
	unsigned int K = x.size();
	double sumAlpha = 0.0;
	for(unsigned int k=0;k<K;k++){
		double logxk=x[k];
		if(!xOnLogScale){
			logxk=log(logxk);
		}

		out+=(alpha[k]-1)*logxk-lgamma(alpha[k]);
		sumAlpha+=alpha[k];
	}
	out+=lgamma(sumAlpha);

	return out;
}

double logPdfExponential(const double& x,const double& lambda){
	return log(lambda)-lambda*x;
}

double logPdfGamma(const double& x,const double& shape,const double& rate){
	return shape*log(rate)+(shape-1)*log(x)-rate*x-lgamma(shape);
}

double logPdfNormal(const double& x, const double& mu, const double& sigma){
	return -0.5*log(2*pi<double>())-log(sigma)-0.5*(x-mu)*(x-mu)/(sigma*sigma);
}

double logPdfScaleInvChiSquare(const double& x, const int& nu, const double& sigma2_0){
  return nu/2*(log(sigma2_0) + log(nu) - log(2)) - lgamma((double)nu/2) - (nu*sigma2_0/(2*x)) - (1+nu/2)*log(x);
}

double logPdfPoisson(const unsigned int &x, const double& mu){
	return (double)x*log(mu) -mu - lgamma((double)x+1);
}

double logPdfStudentsT(const double &x, const unsigned int& dof){
	return -0.5*log((double)dof*pi<double>())-(((double)dof+1)/2.0)*log(1.0+x*x/(double)dof)
				+lgamma((double)(dof+1)/2.0) - lgamma((double)dof/2.0);
}

double logPdfLocationScaleT(const double &x, const double& mu,
										const double& sigma, const unsigned int& dof){
	// http://www.mathworks.com/help/toolbox/stats/brn2ivz-145.html
	// location param=mu, scale param=sigma, shape param=dof
	return logPdfStudentsT((x-mu)/sigma,dof);

}

double logPdfTruncatedNormal(const double& x,const double& mean,const double& stdDev,
								const string& distType, const double& lower,const double& upper){

	boost::math::normal_distribution<double> normDist(mean,stdDev);
	double pLower,pUpper;

	if(distType.compare("U")==0){
		pLower=0.0000000001;
		pUpper=cdf(normDist,upper);
	}else if(distType.compare("L")==0){
		pLower=cdf(normDist,lower);
		pUpper=1.0-0.0000000001;
	}else{
		pLower=cdf(normDist,lower);
		pUpper=cdf(normDist,upper);
	}

	return log(pdf(normDist,x))-log(pUpper-pLower);

}

double logPdfMultinomialSizeOne(const unsigned int& x, const vector<double>& p){
       return log(p[x]);
}

double logPdfWeibullCensored(const double& x, const double& shape, const double& scale, const unsigned int& censoring){

	double out = - scale * pow(x,shape);
	if (censoring == 1) out += log(scale) + log(shape) + (shape-1) * log(x) ;
	return out;
}

double logPdfMultivarNormal(const unsigned int& sizeX, const VectorXd& x,const VectorXd& meanVec,const MatrixXd& sqrtPrecMat,const double& logDetPrecMat){
	// If S is the (upper triangular matrix) square root of the precision P, then S'S=P

	VectorXd work = sqrtPrecMat*(x-meanVec);
	double tmp = work.squaredNorm();
	return -0.5*((double)sizeX*log(2*pi<double>())-logDetPrecMat+tmp);
}

//RJ logPdfNIW
double logPdfNIW(const unsigned int& d, const MatrixXd& Lambda0, const double& logDetInvLambda0, const MatrixXd& invSigma, const double& logDetSigma, const double& kappa0, const int& nu0, const VectorXd& mu,const VectorXd& mu0){

	double Z = 0.5*(double)nu0*(double)d*log(2.0) + logMultivarGammaFn((double)nu0/2.0,d) + (double)d/2.0*log(2.0*pi<double>()/kappa0) + (double)nu0/2.0*logDetInvLambda0;
	double logPmuSigma = - (((double)nu0+(double)d)/2.0+1.0)*logDetSigma - 0.5*(Lambda0*invSigma).trace() - 0.5*kappa0*(mu-mu0).transpose()*invSigma*(mu-mu0) - Z;
	//cout << " logPmuSigma "<< logPmuSigma << " logDetSigma "<<logDetSigma<< " trace "<<(Lambda0*invSigma).trace()<< " kappa "<<kappa0*(mu-mu0).transpose()*invSigma*(mu-mu0)<< " Z "<< Z <<endl;
	return logPmuSigma;
}

//AR logPdfMultivariateNormal p(fk|Lk,k)
double logPdfMultivariateNormal(const vector<double>& muGP, const vector<double> L, const vector<double>& times, const string& kernelType){

  VectorXd y(muGP.size());
  //std::fstream foutC("compare.txt", std::ios::in | std::ios::out | std::ios::app);

  for(unsigned int i=0;i<muGP.size();i++)
    y(i)=muGP[i];


  MatrixXd precMat(muGP.size(),muGP.size());
  precMat.setZero(muGP.size(),muGP.size());

  GP_cov(precMat,L,times,0,kernelType,0);
  double logDetPrecMat=0 ;//=log(precMat.determinant());
  LLT<MatrixXd> lltOfA(precMat); // compute the Cholesky decomposition of A
  MatrixXd Llt = lltOfA.matrixL();

  for(unsigned int i=0;i<muGP.size();i++)
    logDetPrecMat += 2*log(Llt(i,i));
   // logDetPrecMat=  2*log(Llt.determinant());

   MatrixXd precMat1=precMat.inverse();
   //fout << "det L_inv 2 "<< det << " mat "<<log(Mat.determinant()) <<endl;
   MatrixXd Mat2 = Llt.inverse().transpose()*Llt.inverse();



  //double dmvnorm = -0.5*y.transpose()*precMat*y - 0.5*y.size()*log(2.0*pi<double>()) - 0.5*logDetPrecMat;
  double bb=y.transpose()*precMat1*y;
  double dmvnorm = -0.5*bb - 0.5*y.size()*log(2.0*pi<double>()) - 0.5*logDetPrecMat;

  if(bb < 0){
    //LLT<MatrixXd> llt;
    //llt.compute(precMat);
    //MatrixXd Sigma= llt.solve(MatrixXd::Identity(precMat.rows(), precMat.rows())).transpose()*Llt.solve(MatrixXd::Identity(precMat.rows(), precMat.rows()));

    //MatrixXd Sigma = precMat+0.01*MatrixXd::Identity(precMat.rows(), precMat.rows());
    //MatrixXd precMat2=PartialPivLU(precMat);
    //double bb2=y.transpose()*precMat2*y;
    //dmvnorm = -0.5*bb2 - 0.5*y.size()*log(2.0*pi<double>()) - 0.5*logDetPrecMat;
    //foutC << "bb "<< bb<<" dmvnorm "<<  dmvnorm<< " sizek " <<y.size()<<
     // " L "<< L[0] << " "<< L[1] << " "<< L[2] <<endl<< " y "<< y.transpose()<<endl<< " det "<< precMat.determinant() <<" vs "<< 2*log(Llt.determinant()) ;

//    double a=-(times[1]-times[0])*(times[1]-times[0])/(2*exp(L[1]));
//    double b = exp(L[0])*exp(a);
//    foutC << endl<<" times[1] "<< std::setprecision(9)<<times[1] << " times[0] "<< std::setprecision(9)<< times[0]<< " b "<< b  << " vs "<<precMat<<endl<<endl;

    dmvnorm = -pow(10,10);
  }
  return dmvnorm;
}

double logPdfWishart(const unsigned int& dimA, const MatrixXd& A, const double& logDetA, const MatrixXd& inverseR, const double& logDetR, const double& kappa){
	MatrixXd work = inverseR*A;
	return -0.5*(kappa*logDetR-(kappa-(double)dimA-1)*logDetA+work.trace())
			  - (0.5*kappa*(double)dimA*log(2.0)+logMultivarGammaFn(kappa/2.0,dimA));
}

double logPdfInverseWishart(const unsigned int& dimA, const MatrixXd& A, const double& logDetA, const MatrixXd& covR, const double& logDetR, const double& kappa){
  MatrixXd work = covR*A.inverse();
  return 0.5*(kappa*logDetR-(kappa+(double)dimA+1)*logDetA-work.trace())
    - (0.5*kappa*(double)dimA*log(2.0)+logMultivarGammaFn(kappa/2.0,dimA));
}

double logPdfIntrinsicCAR(const vector<double>& x, const vector<vector<unsigned int> >& Neighbours, const double& precision){
	// compute xtPx with P the precision Matrix of the intrinsic CAR
	double sumCAR1 = 0.0;
	double sumCAR2 = 0.0;
	int n=x.size();
	for (int i=0; i<n; i++){
		double xi = x[i];
		int nNeighi = (Neighbours[i]).size();
		sumCAR1+= xi*xi*nNeighi;
		for (int j = 0; j<nNeighi; j++){
			unsigned int nj = Neighbours[i][j];
			double xj = x[nj-1];
			sumCAR2+=xi*xj;
	        }
	}
	double xtPx=sumCAR1-sumCAR2;

	return 0.5*(n-1)*log(precision)-0.5*precision*xtPx;
}

#endif /*DISTRIBUTION_H_*/
