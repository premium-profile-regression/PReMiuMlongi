/// \file mathfunctions.h
/// \author David Hastie
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


#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include<cmath>
#include <Eigen/Eigen>
#include<iostream>
#include<fstream>
#include <Eigen/Dense>
#include<boost/math/special_functions/gamma.hpp>
#include<string>
#include <algorithm>
#include <iterator>


using std::ifstream;
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::string;
using namespace boost::math::constants;
using namespace Eigen;

using boost::math::lgamma;

double logMultivarGammaFn(const double x,const unsigned int p){

  double out;
  out = 0.25*(double)(p*(p-1))*log(pi<double>());
  for(unsigned int i=1;i<=p;i++){
    out += lgamma(x+(1.0-(double)i)/2.0);
  }
  return out;
}

double logit(const double& lambda){
  return 1.0/(1.0+exp(-lambda));
}

//RJ matrix (fast) inversion function
void invert(MatrixXd& invSigma,MatrixXd& Sigma,const unsigned int dimBlock, double noise){
  if(dimBlock>0){
    int dimSigma = Sigma.rows();
    int nBlocks = dimSigma/dimBlock;
    MatrixXd C;
    C = Sigma.block(0, 0, dimBlock, dimBlock);
    for(int i=0;i<dimBlock;i++)
      C(i,i) = C(i,i) - noise;
    double invNoise = 1.0/noise;
    MatrixXd E;
    E = nBlocks*C;
    for(int i=0;i<dimBlock;i++)
      E(i,i) = E(i,i) + noise;
    E = - invNoise * C * E.inverse();
    invSigma = E.replicate(nBlocks,nBlocks);
    for(int i=0;i<dimSigma;i++)
      invSigma(i,i) = invSigma(i,i) + invNoise;
  }else{
    invSigma = Sigma.inverse();
  }
}
//RJ GP covariance function
void GP_cov(MatrixXd& Mat,std::vector<double> L,std::vector<double> times,const unsigned int dimBlock, const string& kernel, const unsigned int error){
  double a;
  int i,j;
  int nTimes = times.size();
  Mat.setZero(Mat.rows(),Mat.cols());

  if(kernel.compare("SQexponential")==0){

    double eL0 = exp(L[0]);
    double eL1 = exp(L[1])*2.0;
    double eL2 = exp(L[2]);
    if(error==0)
      eL2 = 0.001;


    // if(dimBlock<0){
    //   int nBlocks = nTimes/dimBlock;
    //   MatrixXd unit;
    //   unit.resize(dimBlock,dimBlock);
    //   unit.fill(0.0);
    //
    //   for(i=1;i<dimBlock;i++){
    //     for(j=0;j<i;j++){
    //       a=-(times[i]-times[j])*(times[i]-times[j])/eL1;
    //       unit(i,j)=eL0*std::exp(a);
    //     }
    //   }
    //   unit = unit + unit.transpose() + eL0*MatrixXd::Identity(dimBlock, dimBlock);
    //   Mat = unit.replicate(nBlocks,nBlocks);
    //   for(i=0;i<nTimes;i++) Mat(i,i) = Mat(i,i) + eL2;
    // }else{

    for(i=0 ;i<nTimes;i++){
      for(j=0;j<i;j++){
        a=-(times[i]-times[j])*(times[i]-times[j])/eL1;
        Mat(i,j)=eL0*std::exp(a);
      }
    }
    Mat = Mat + Mat.transpose();
    for(int i=0;i<nTimes;i++)
      Mat(i,i) = eL0+eL2;


  }else{
    //[sigma_b^2 + sigma_v^2(t-l)(t-l)]^2
    double eL0 = exp(L[0]); //sigma_b^2
    double eL1 = exp(L[1]); //sigma_v^2
    double eL2 = exp(L[2]); // sigma_e^2
    if(error==0)
      eL2=0;
    double eL3 = exp(L[3]); // l

    for(i=0 ;i<nTimes;i++){
      for(j=0;j<i;j++){
        a=eL0+eL1*(times[i]-eL3)*(times[j]-eL3);
        Mat(i,j)=a*a;
      }
    }
    Mat = Mat + Mat.transpose();

    for(int i=0;i<nTimes;i++){
      a=eL0+eL1*(times[i]-eL3)*(times[i]-eL3);
      Mat(i,i)=a*a + eL2;
    }
  }

}

//AR covariance matrix for Linear Mixed Model
// void ME_cov(MatrixXd& Mat,std::vector<double> L,std::vector<double> times){
//   ME_cov(Sigma,params,timesk); // Sigma ordered;
//
//
// }

//AR to obtain Sigma_permut^{-1} from Sigma^{-1} = Mat
// in: Mat a1 a2 a3 b1 b2 b3 c1 c2 c3
// out: Mat a1 a2 a3 c1 c2 c3 b1 b2 b3
void Permut_cov(MatrixXd& Mat,  const int start_permut, const int length_permut){
  int i,j,ii,jj,ai,aj;
  int nTimes = Mat.rows();
  MatrixXd Mat_permut;
  Mat_permut.setZero(nTimes,nTimes);

  ai=0;
  for(i=0;i<nTimes;i++){

    if(i< start_permut){
      ii=i;
    }else if(i<nTimes-length_permut){
      ii=i+length_permut;
    }else{
      ii=start_permut+ai;
      ai++;
    }

    aj=0;
    for(j=0;j<nTimes;j++){
      if(j< start_permut){
        jj=j;
      }else if(j<nTimes-length_permut){
        jj=j+length_permut;
      }else{
        jj=start_permut+aj;
        aj++;
      }
      Mat_permut(i,j)=Mat(ii,jj);
    }
  }
  Mat=Mat_permut;
}

//AR Permut_cov with M0 (Mat) sorted.
// in: Mat a1 b1 c1 a2 b2 c2 a3 b3 c3
// out: Mat a1 c1 a2 c2 a3 c3 b1 b2 b3
void Permut_cov_sorted(MatrixXd& Mat_sorted,  const int start_permut, const int length_permut, const std::vector<int> & idx){
  int i,j,ii,jj,ai,aj;
  int nTimes = Mat_sorted.rows();

  MatrixXd Mat_permut;
  Mat_permut.setZero(nTimes,nTimes);
  std::vector<double> idx_i;
  ai=0;
  for(i=0;i<nTimes;i++){
    if(idx[i+idx_i.size()] >= start_permut && idx[i+idx_i.size()] < (start_permut+ length_permut))
      idx_i.push_back(i+idx_i.size());
    if(i< (nTimes-length_permut)){
      ii=i+idx_i.size();
    }else {
      ii=idx_i[ai];
      ai++;
    }

    aj=0;
    std::vector<double> idx_j;
    for(j=0;j<nTimes;j++){
      if(idx[j+idx_j.size()] >= start_permut && idx[j+idx_j.size()] < (start_permut+ length_permut))
        idx_j.push_back(j+idx_j.size());
      if(j< (nTimes-length_permut)){
        jj=j+idx_j.size();
      }else {
        jj=idx_j[aj];
        aj++;
      }
      Mat_permut(i,j)=Mat_sorted(ii,jj);
    }
  }
  Mat_sorted=Mat_permut;
  //cout << endl<<endl;
}

void Permut_cov_sorted(MatrixXd& Mat,   const std::vector<int> & idx){
  //Permut according to idx

  int i,j;
  int nTimes = Mat.rows();
  MatrixXd Mat_permut;
  Mat_permut.setZero(nTimes,nTimes);
  std::vector<double> idx_i;

  for(i=0;i<nTimes;i++){
    for(j=0;j<nTimes;j++){
      Mat_permut(i,j)=Mat(idx[i],idx[j]);
    }
  }
  Mat=Mat_permut;
}

//AR
double Get_Sigma_inv_GP_cov(MatrixXd& Mat, std::vector<double> L, std::vector<double> &times,
                            const unsigned int dimBlock, const std::vector<double>& grid,
                            const string& kernel){

  double a,eL0,eL1,eL2,eL3;
  int i,j;
  int nTimes = times.size();
  double det=0;
  std::fstream fout("file_output.txt", std::ios::in | std::ios::out | std::ios::app);

  Mat.setZero(nTimes,nTimes);

  // ascending order of times
  // std::vector<int> idx(times.size());
  // int x=0;
  // iota(idx.begin(), idx.end(), x++);
  // stable_sort(idx.begin(), idx.end(), // sort indexes based on comparing values in times
  //      [&](int i1,int i2) { return (times[i1] < times[i2]); }
  //      );
  sort(times.begin(), times.end());

  if(kernel.compare("SQexponential")==0){
    eL0 = exp(L[0]);
    eL1 = exp(L[1])*2.0;
    eL2 = exp(L[2]);

    for(i=1;i<nTimes;i++){
      for(j=0;j<i;j++){
        a=-(times[i]-times[j])*(times[i]-times[j])/eL1;
        Mat(i,j)=eL0*std::exp(a);
      }
    }
    Mat = Mat + Mat.transpose();
    for(int i=0;i<nTimes;i++)
      Mat(i,i) = Mat(i,i) + eL0+eL2;

  }else{

    //[sigma_b^2 + sigma_v^2(t-l)(t-l)]^2
    eL0 = exp(L[0]); //sigma_b^2
    eL1 = exp(L[1]); //sigma_v^2
    eL2 = exp(L[2]);
    eL3 = exp(L[3]); // l

    for(i=1;i<nTimes;i++){
      for(j=0;j<i;j++){
        a=eL0+eL1*(times[i]-eL3)*(times[j]-eL3);
        Mat(i,j)=a*a;
      }
    }
    Mat = Mat + Mat.transpose();

    for(int i=0;i<nTimes;i++){
      a=eL0+eL1*(times[i]-eL3)*(times[i]-eL3);
      Mat(i,i)=a*a+eL2;
    }
  }

  int trick=0;
  if(*std::max_element(std::begin(grid), std::end(grid))==0 || grid.size() >= nTimes || nTimes < 100){
    LLT<MatrixXd> lltOfA(Mat); // compute the Cholesky decomposition of A
    MatrixXd L = lltOfA.matrixL();
    det=  2*L.diagonal().array().log().sum();
    Mat = L.inverse().transpose()*L.inverse();
    if(isnan(det))
      trick=1;
  }else{
    trick=1;
  }


  if(trick==1 && nTimes >= grid.size()){

    MatrixXd Ktu;
    MatrixXd Kuu;
    Ktu.setZero(nTimes,grid.size());
    Kuu.setZero(grid.size(),grid.size());

    for(i=0;i<nTimes;i++){
      for(j=0;j<grid.size();j++){
        if(kernel.compare("SQexponential")==0){
          a=-(times[i]-grid[j])*(times[i]-grid[j])/eL1;
          Ktu(i,j)=eL0*std::exp(a);
        }else{
          a=eL0+eL1*(times[i]-eL3)*(grid[j]-eL3);
          Ktu(i,j)=a*a;
        }
      }
    }

    for(i=1;i<grid.size();i++){
      for(j=0;j<i;j++){
        if(kernel.compare("SQexponential")==0){
          a=-(grid[i]-grid[j])*(grid[i]-grid[j])/eL1;
          Kuu(i,j)=eL0*std::exp(a);
        }else{
          a=eL0+eL1*(grid[i]-eL3)*(grid[j]-eL3);
          Kuu(i,j)=a*a;
        }
      }
    }
    Kuu = Kuu + Kuu.transpose();

    for(int i=0;i<grid.size();i++){
      if(kernel.compare("SQexponential")==0){
        Kuu(i,i) = Kuu(i,i) + eL0;//+0.00;
      }else{
        for(int i=0;i<grid.size();i++){
          a=eL0+eL1*(grid[i]-eL3)*(grid[i]-eL3);
          Kuu(i,i)=a*a ;//+ eL2;
        }
      }
    }

    if(Kuu.determinant()<0){
      for(int i=0;i<grid.size();i++)
        Kuu(i,i) = Kuu(i,i)+ 0.01;
    }

    LLT<MatrixXd> lltOfA(Kuu); // compute the Cholesky decomposition of A
    MatrixXd Lm = lltOfA.matrixL();
    MatrixXd L_inv = Lm.inverse();
    MatrixXd Kuu_inv = L_inv.transpose()*L_inv;

    //Mat = Mat - eL2*MatrixXd::Identity(nTimes,nTimes);
    MatrixXd Qtt=Ktu*Kuu_inv*Ktu.transpose();
    MatrixXd Lambda=(Mat.diagonal()-Qtt.diagonal()).asDiagonal();
    //Lambda = Lambda + eL2*MatrixXd::Identity(nTimes,nTimes);
    MatrixXd Lambda_inv = Lambda.inverse();

    double logdet_lambda=0;

    LLT<MatrixXd> lltOfAL(Lambda); // compute the Cholesky decomposition of A
    MatrixXd L_Lambda = lltOfAL.matrixL();
    logdet_lambda=  2*L_Lambda.diagonal().array().log().sum();

    if(isnan(logdet_lambda)){
      logdet_lambda=0;
      PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu(Lambda);
      auto& LU = lu.matrixLU();
      double c = lu.permutationP().determinant(); // -1 or 1
      for (unsigned i = 0; i < LU.rows(); ++i) {
        const auto& lii = LU(i,i);
        if (lii < double(0)) c *= -1;
        logdet_lambda += log(abs(lii));
      }
      logdet_lambda += log(c);

      if(isnan(logdet_lambda)){
        cout <<  " logdet_lambda=nan, diag Lambda "<< LU.diagonal()<<endl;
      }
    }
    //for(int i=0;i<nTimes;i++)
    //  logdet_lambda += log(Lambda(i,i));

    MatrixXd Aut=L_inv* Ktu.transpose();
    MatrixXd T=(MatrixXd::Identity(grid.size(), grid.size())+
      Aut*Lambda_inv*Aut.transpose());

    //MatrixXd Mat2=Mat;
    Mat = Lambda_inv-Lambda_inv*Aut.transpose()*T.inverse()*Aut*Lambda_inv;
    MatrixXd Prod1(grid.size(), grid.size());

    Prod1=MatrixXd::Identity(grid.size(), grid.size())+Aut*Lambda_inv*Aut.transpose();

    LLT<MatrixXd> lltOfA2(Prod1); // compute the Cholesky decomposition of A
    MatrixXd P = lltOfA2.matrixL();
    double logDetP=  2*P.diagonal().array().log().sum();

    if(isnan(logDetP)){
      logDetP=0;
      PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu(Prod1);
      auto& LU = lu.matrixLU();
      double c = lu.permutationP().determinant(); // -1 or 1
      for (unsigned i = 0; i < LU.rows(); ++i) {
        const auto& lii = LU(i,i);
        if (lii < double(0)) c *= -1;
        logDetP += log(abs(lii));
      }
      logDetP += log(c);

      if(isnan(logDetP)){
        cout << " Prod1.det "<< Prod1.determinant() << " Mat.det "<< Mat.determinant() << " c "<<c<<" logDetP=nan, diag Lambda "<< LU.diagonal().transpose()<<endl<<endl;
      }
    }

    det=logDetP+logdet_lambda;

    if(std::isnan(det) || isinf(det)){
      fout << "det Get_Sigma_inv_GP_cov "<< det << " logdet_lambda "<<logdet_lambda<<
        " + logDetP:  "<< logDetP<<" grid.size() "<< grid.size()
                       << " eL0 " << eL0<< " eL1 " << eL1<<  " eL2 " << eL2<<endl <<
        " logdet Lambda "<< log(Lambda.determinant()) <<
            " Prod1.logdet "<< log(Prod1.determinant()) <<endl;

      det=-(std::numeric_limits<double>::max());
      // Add an offset on time vector as equal times can lead to instability

      //std::vector<int>::iterator it;
      //it=std::unique_copy (times.begin(), it, times.begin(), myfunction);
//
//       std::vector<double> times_unique = times;
//       auto last = std::unique(times_unique.begin(), times_unique.end());
//       times_unique.erase(last, times_unique.end());
//
//       cout << " ~~~~~~~~~~~~" <<endl<< "det Get_Sigma_inv_GP_cov "<< det << " logdet_lambda "<<logdet_lambda<<
//         " + logDetP:  "<< logDetP<<" grid.size() "<< grid.size()
//                        << " eL0 " << eL0<< " eL1 " << eL1<<  " eL2 " << eL2<<endl <<
//         " logdet Lambda "<< log(Lambda.determinant()) <<
//             " Prod1.logdet "<< log(Prod1.determinant()) <<endl;
//
//       if(times_unique.size() != times.size()){
//
//         std::vector<double> times_offset = times;
//
//         int j=0;
//         for(i=0;i<times_unique.size();i++){
//           int offset=0;
//             while(times[j]==times_unique[i]){
//               times_offset[j] += 0.001*offset;
//               j++;
//               offset++;
//             }
//         }
//
//         if(kernel.compare("SQexponential")==0){
//           Mat.setZero(nTimes,nTimes);
//
//           for(i=1;i<nTimes;i++){
//             for(j=0;j<i;j++){
//               a=-(times_offset[i]-times_offset[j])*(times_offset[i]-times_offset[j])/eL1;
//               Mat(i,j)=eL0*std::exp(a);
//             }
//           }
//           Mat = Mat + Mat.transpose();
//           for(int i=0;i<nTimes;i++)
//             Mat(i,i) = Mat(i,i) + eL0 + eL2;
//
//         }else{
//
//           for(i=1;i<nTimes;i++){
//             for(j=0;j<i;j++){
//               a=eL0+eL1*(times_offset[i]-eL3)*(times_offset[j]-eL3);
//               Mat(i,j)=a*a;
//             }
//           }
//           Mat = Mat + Mat.transpose();
//
//           for(int i=0;i<nTimes;i++){
//             a=eL0+eL1*(times_offset[i]-eL3)*(times_offset[i]-eL3);
//             Mat(i,i)=a*a+eL2;
//           }
//         }
//
//         cout << " 2 - logdet mat " << Mat.determinant()<< " nTimes "<< nTimes<<" time "<< times_offset[2]<<endl  <<endl ;
//
//
//         MatrixXd Ktu;
//         MatrixXd Kuu;
//         Ktu.setZero(nTimes,grid.size());
//         Kuu.setZero(grid.size(),grid.size());
//
//         for(i=0;i<nTimes;i++){
//           for(j=0;j<grid.size();j++){
//             if(kernel.compare("SQexponential")==0){
//               a=-(times_offset[i]-grid[j])*(times_offset[i]-grid[j])/eL1;
//               Ktu(i,j)=eL0*std::exp(a);
//             }else{
//               a=eL0+eL1*(times_offset[i]-eL3)*(grid[j]-eL3);
//               Ktu(i,j)=a*a;
//             }
//           }
//         }
//
//         for(i=1;i<grid.size();i++){
//           for(j=0;j<i;j++){
//             if(kernel.compare("SQexponential")==0){
//               a=-(grid[i]-grid[j])*(grid[i]-grid[j])/eL1;
//               Kuu(i,j)=eL0*std::exp(a);
//             }else{
//               a=eL0+eL1*(grid[i]-eL3)*(grid[j]-eL3);
//               Kuu(i,j)=a*a;
//             }
//           }
//         }
//         Kuu = Kuu + Kuu.transpose();
//
//         for(int i=0;i<grid.size();i++){
//           if(kernel.compare("SQexponential")==0){
//             Kuu(i,i) = Kuu(i,i) + eL0;//+0.00;
//           }else{
//             for(int i=0;i<grid.size();i++){
//               a=eL0+eL1*(grid[i]-eL3)*(grid[i]-eL3);
//               Kuu(i,i)=a*a +0.0001;//+ eL2;
//             }
//           }
//         }
//
//         if(Kuu.determinant()<0){
//           for(int i=0;i<grid.size();i++)
//             Kuu(i,i) = Kuu(i,i)+ 0.01;
//         }
//
//         LLT<MatrixXd> lltOfA(Kuu); // compute the Cholesky decomposition of A
//         MatrixXd Lm = lltOfA.matrixL();
//         MatrixXd L_inv = Lm.inverse();
//         MatrixXd Kuu_inv = L_inv.transpose()*L_inv;
//
//         //Mat = Mat - eL2*MatrixXd::Identity(nTimes,nTimes);
//         MatrixXd Qtt=Ktu*Kuu_inv*Ktu.transpose();
//         //Mat = Mat - eL2*MatrixXd::Identity(nTimes,nTimes);
//         MatrixXd Lambda=(Mat.diagonal()-Qtt.diagonal()).asDiagonal();
//         //Lambda = Lambda + eL2*MatrixXd::Identity(nTimes,nTimes);
//         MatrixXd Lambda_inv = Lambda.inverse();
//
//         double logdet_lambda=0;
//
//         LLT<MatrixXd> lltOfAL(Lambda); // compute the Cholesky decomposition of A
//         MatrixXd L_Lambda = lltOfAL.matrixL();
//         logdet_lambda=  2*L_Lambda.diagonal().array().log().sum();
//
//
//         //for(int i=0;i<nTimes;i++)
//         //  logdet_lambda += log(Lambda(i,i));
//
//         double ld = 0;
//         PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu(Lambda);
//         auto& LU = lu.matrixLU();
//         double c = lu.permutationP().determinant(); // -1 or 1
//         for (unsigned i = 0; i < LU.rows(); ++i) {
//           const auto& lii = LU(i,i);
//           if (lii < double(0)) c *= -1;
//           ld += log(abs(lii));
//         }
//         ld += log(c);
//
//
//
//         MatrixXd Aut=L_inv* Ktu.transpose();
//         MatrixXd T=(MatrixXd::Identity(grid.size(), grid.size())+
//           Aut*Lambda_inv*Aut.transpose());
//
//         //MatrixXd Mat2=Mat;
//         Mat = Lambda_inv-Lambda_inv*Aut.transpose()*T.inverse()*Aut*Lambda_inv;
//         MatrixXd Prod1(grid.size(), grid.size());
//
//         Prod1=MatrixXd::Identity(grid.size(), grid.size())+Aut*Lambda_inv*Aut.transpose();
//
//         LLT<MatrixXd> lltOfA2(Prod1); // compute the Cholesky decomposition of A
//         MatrixXd P = lltOfA2.matrixL();
//         double logDetP=  2*P.diagonal().array().log().sum();
//
//
//         double ld2 = 0;
//         PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu2(Prod1);
//         auto& LU2 = lu2.matrixLU();
//         double c2 = lu2.permutationP().determinant(); // -1 or 1
//         for (unsigned i = 0; i < LU2.rows(); ++i) {
//           const auto& lii2 = LU2(i,i);
//           if (lii2 < double(0)) c2 *= -1;
//           ld2 += log(abs(lii2));
//         }
//         ld2 += log(c2);
//
//
//         det=logDetP+logdet_lambda;
//
//
//         cout <<"bis Get_Sigma_inv_GP_cov "<< det << " logdet_lambda "<<logdet_lambda<<
//           " + logDetP:  "<< logDetP<<" grid.size() "<< grid.size()
//                          << " eL0 " << eL0<< " eL1 " << eL1<<  " eL2 " << eL2<<endl <<
//           " logdet Lambda "<< log(Lambda.determinant()) <<
//             " ld "<< ld <<
//               " Prod1.logdet "<< log(Prod1.determinant()) <<
//                 " ld2 "<< ld2 <<endl <<" ~~~~~~~~~~"<<endl<<endl;
//
        }



      // <<
      //     " Lambda diag "<<endl<<Lambda.diagonal().transpose()<<endl<<
      //       " Qtt diag "<<endl<<Qtt.diagonal().transpose()<<endl<<
      //         " Mat diag "<<endl<<Mat2.diagonal().transpose()<<endl<<
      //           " Lambda "<<Lambda<< endl<<endl<<
      //             " Aut "<<Aut<< endl<<endl;
    }

    //
    //     if(nTimes ==5){
    //       fout <<  " times "<<endl;
    //       for(unsigned int j=0; j<times.size(); j++)
    //         fout << times[j]<< " ";
    //       fout <<  " grid "<<endl;
    //       for(unsigned int j=0; j<grid.size(); j++)
    //         fout << grid[j]<< " ";
    //       fout << endl << " params.L(c) "<<endl;
    //       for(unsigned int j=0; j<3; j++)
    //         fout << L[j]<< " ";
    //       fout << endl<< " Mat inv " << log(Mat.determinant()) <<endl<< Mat.inverse() << endl;
    //
    //       Mat =Lambda.inverse()-Lambda.inverse()*Aut.transpose()*(MatrixXd::Identity(grid.size(), grid.size())+Aut*Lambda.inverse()*Aut.transpose()).inverse()*Aut*Lambda.inverse();
    //       det=log((MatrixXd::Identity(grid.size(), grid.size())+Aut*Lambda.inverse()*Aut.transpose()).determinant()*Lambda.determinant());
    //
    //       fout <<" precMat "<< det<<endl<< Mat << endl << " Ktu "<< endl<<Ktu  << endl<<" Kuu "<< endl<<Kuu  << endl << " Kut "<<endl<< Ktu.transpose() <<endl << " Ktu*Kuu-1 "<<endl<<Ktu* Kuu.inverse() << " Kuu*Kuu-1 "<<endl<<Kuu* Kuu.inverse() <<endl<<endl;
    //     }

  return det;
}

//AR Inverse function for block matrix with Woodbury identity
// double Inverse_woodbury(const MatrixXd& M0, const double& log_det_M0, MatrixXd& Mat,std::vector<double> times);
//
// double Inverse_woodbury(const MatrixXd& M0, const double& log_det_M0, MatrixXd& Mat){
//   std::vector<double> times;
//   return Inverse_woodbury( M0,  log_det_M0,  Mat,  times);
// }

double Inverse_woodbury(const MatrixXd& M0_inv, const double& log_det_M0, MatrixXd& Mat, const std::vector<int> idx){
  // Function to add one subject
  double log_DetPrecMat=0.0;
  MatrixXd kno;
  MatrixXd Knew;
  std::fstream fout("file_output.txt", std::ios::in | std::ios::out | std::ios::app);

  int i,j;
  int nTimes = Mat.rows();

  i=M0_inv.rows();

  if(i<nTimes){ // add one subject

    Knew.setZero(nTimes-i,nTimes-i);
    kno.setZero(nTimes-i,i);

    for(int i2=i;i2<nTimes;i2++){
      for(j=0;j<i2;j++){
        if(j<i){
          kno(i2-i,j)=Mat(i2,j);
        }else{
          Knew(i2-i,j-i)=Mat(i2,j);
        }
      }
    }

    Knew = Knew + Knew.transpose();
    for(int i2=0;i2<nTimes-i;i2++)
      Knew(i2,i2) = Mat(i+i2,i+i2);

    MatrixXd A(nTimes-i,nTimes-i), A_inv(nTimes-i,nTimes-i);
    MatrixXd B(i,nTimes-i);
    A.setZero(nTimes-i,nTimes-i);
    B.setZero(i,nTimes-i);
    B=M0_inv*kno.transpose();
    A=Knew-kno*B;

    if(A.determinant()<=0){
      for(int i=0;i<A.cols();i++)
        A(i,i) = A(i,i)+ 0.01;
    }

    LLT<MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
    MatrixXd La = lltOfA.matrixL();
    double logdetA=  2*La.diagonal().array().log().sum();
    A_inv = La.inverse().transpose()*La.inverse();
    log_DetPrecMat=log_det_M0+logdetA;

    Mat.setZero(nTimes,nTimes);
    Mat.topRows(i)<<M0_inv+B*A_inv*kno*M0_inv, -B*A_inv;
    Mat.bottomRows(nTimes-i)<<-A_inv*B.transpose(), A_inv;

    //Permut back
    Permut_cov_sorted(Mat, idx);

    if(std::isnan(log_DetPrecMat) || isinf(log_DetPrecMat)){
      logdetA=0;
      PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu(A);
      auto& LU = lu.matrixLU();
      double c = lu.permutationP().determinant(); // -1 or 1
      for (unsigned i = 0; i < LU.rows(); ++i) {
        const auto& lii = LU(i,i);
        if (lii < double(0)) c *= -1;
        logdetA += log(abs(lii));
      }
      logdetA += log(c);
      log_DetPrecMat=log_det_M0+logdetA;

      if(std::isnan(log_DetPrecMat) || isinf(log_DetPrecMat)){
        //fout << "in inverse_Woodbury to add: log_DetPrecMat="<<log_DetPrecMat << " = log_det_M0 "<< log_det_M0 <<" +logdetA "<< logdetA<<  " log(detA)"<< log(A.determinant())<<" ";
        log_DetPrecMat=-(std::numeric_limits<double>::max());
      }
    }

  }else{ // remove one subject i>nTimes
    cout << " problem inverse_Woodbury !!!!!" <<endl;
  }

  return(log_DetPrecMat);
}

double Inverse_woodbury(const MatrixXd& M0, const double& log_det_M0, MatrixXd& Mat, MatrixXd& M0_inv){
  // Function to remove one subject
  std::fstream fout("file_output.txt", std::ios::in | std::ios::out | std::ios::app);

  double log_DetPrecMat=0.0;
  MatrixXd kno;
  MatrixXd Knew;

  int i,j;
  int nTimes = Mat.rows();

  i=M0.rows();
  if(i<nTimes){ cout << " problem Inverse_woodbury, should be removing one subject"<<endl;}else{ // remove one subject i>nTimes

    Knew.setZero(i-nTimes,i-nTimes);
    kno.setZero(i-nTimes,nTimes);

    for(int i2=nTimes;i2<i;i2++){
      for(j=0;j<i2;j++){
        if(j<nTimes){
          kno(i2-nTimes,j)=M0(i2,j);
        }else{
          Knew(i2-nTimes,j-nTimes)=M0(i2,j);
        }
      }
    }

    Knew = Knew + Knew.transpose();


    for(int i2=0;i2<i-nTimes;i2++)
      Knew(i2,i2) = M0(nTimes+i2,nTimes+i2);



    MatrixXd U(i,(i-nTimes));
    U << MatrixXd::Zero(nTimes,i-nTimes), MatrixXd::Identity(i-nTimes,i-nTimes);
    MatrixXd V(i-nTimes,i);
    V << -kno, MatrixXd::Zero(i-nTimes,i-nTimes);

    MatrixXd Mnew_inv=M0_inv; // To improve

    MatrixXd T=(MatrixXd::Identity(i-nTimes,i-nTimes)+V*Mnew_inv*U);
    MatrixXd M1_inv(i,i);
    M1_inv.setZero(i,i);
    M1_inv=Mnew_inv - Mnew_inv*U*T.inverse()*V*Mnew_inv;

    MatrixXd U1=-V.transpose();
    MatrixXd V1=-U.transpose();

    MatrixXd M2(i,i);
    M2.setZero(i,i);
    MatrixXd T1=(MatrixXd::Identity(i-nTimes,i-nTimes)+V1*M1_inv*U1);
    M2=M1_inv - M1_inv*U1*T1.inverse()*V1*M1_inv;

    Mat.setZero(nTimes,nTimes);
    for(int i2=0;i2<nTimes;i2++){
      for(j=0;j<i2;j++){
        Mat(i2,j)=M2(i2,j);
      }
    }
    Mat = Mat + Mat.transpose();
    for(int i2=0;i2<nTimes;i2++)
      Mat(i2,i2) = M2(i2,i2);

    MatrixXd A(i-nTimes,i-nTimes);
    MatrixXd B(nTimes,i-nTimes);

    A.setZero(i-nTimes,i-nTimes);
    B.setZero(nTimes,i-nTimes);
    B=Mat*kno.transpose();
    A=Knew-kno*B;

    if(A.determinant()<=0){
      for(int i=0;i<A.cols();i++)
        A(i,i) = A(i,i)+ 0.01;
    }

    LLT<MatrixXd> lltOfA(A); // compute the Cholesky decomposition of A
    MatrixXd La = lltOfA.matrixL();
    double logdetA=  2*La.diagonal().array().log().sum();
    log_DetPrecMat=log_det_M0-logdetA;

    if(std::isnan(logdetA) || isinf(logdetA)){
      logdetA=0;
      PartialPivLU<Matrix<double,Dynamic,Dynamic>> lu(A);
      auto& LU = lu.matrixLU();
      double c = lu.permutationP().determinant(); // -1 or 1
      for (unsigned i = 0; i < LU.rows(); ++i) {
        const auto& lii = LU(i,i);
        if (lii < double(0)) c *= -1;
        logdetA += log(abs(lii));
      }
      logdetA += log(c);
      log_DetPrecMat=log_det_M0+logdetA;

      if(std::isnan(log_DetPrecMat) || isinf(log_DetPrecMat)){
        //fout << "in inverse_Woodbury to remove: log_DetPrecMat="<<log_DetPrecMat << " = log_det_M0 "<< log_det_M0 <<" +logdetA "<< logdetA<< " log(detA)"<< log(A.determinant()) <<" ";
        log_DetPrecMat=-(std::numeric_limits<double>::max());
      }
    }
  }
  return(log_DetPrecMat);
}


//RJ GP covariance star function
void GP_cov_star(MatrixXd& Mat,std::vector<double> L,std::vector<double> times,std::vector<double> timestar, const string& kernel){

  double a;
  int i,j;
  int nTimes = times.size();
  int nStar = timestar.size();
  Mat.setZero(Mat.rows(),Mat.cols());

  for(i=0;i<nTimes;i++){
    for(j=0;j<nStar;j++){
      if(kernel.compare("SQexponential")==0){
        a=-1.0/2.0/exp(L[1])*(times[i]-timestar[j])*(times[i]-timestar[j]);
        Mat(i,j)=exp(L[0])*std::exp(a);
      }else{
        a=exp(L[0])+exp(L[1])*(times[i]-exp(L[3]))*(timestar[j]-exp(L[3]));
        Mat(i,j)=a*a ;//+ exp(L[3]);
      }
    }
  }
}
bool myfunction (int i, int j) {
  return (i==j);
}

#endif /*MATHFUNCTIONS_H_*/
