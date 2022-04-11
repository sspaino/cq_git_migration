/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#pragma once

#include <itersolver.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/eig.hpp>
#include <cqlinalg.hpp>
#include <cerr.hpp>

namespace ChronusQ {


  template <typename _F>
  bool GPLHR<_F>::runMicro() {


    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;


    if( isRoot ) {
      std::cout << "\n\n\n";
      std::cout << "  * GPLHR Settings\n"
                << "    * Sigma                        = " << std::real(sigma) << "\n"
                << "    * Min                          = " << std::real(hardLim) << "\n"
                << "    * Kyrlov-Arnoldi Parameter (M) = " << this->m
                << std::endl;


      std::cout << "\n\n" << std::endl;
    }
    
    const size_t N = this->N_;
    const size_t MSS = this->mSS_;
    const size_t nR = this->nRoots_;


    const size_t MSS2 = MSS * MSS;
    const size_t NMSS = N * MSS;
    const size_t NNR  = N * nR;
    const size_t M_NR = this->m * nR;


    // FIXME Reassign sigma (TEMP)
    if( std::imag(sigma) > 1e-10 )
      CErr("GPLHR Only Supports Real Sigma for the Time Being",std::cout);
    double sigmaD   = std::real(sigma);
    double hardLimD = std::real(hardLim);
    //sigmaD = hardLimD;


    // Initialize pointers
    _F *V = nullptr , *W = nullptr , *S = nullptr , *P = nullptr ;
    _F *AV = nullptr, *AW = nullptr, *AS = nullptr, *AP = nullptr;
    _F *Q = nullptr , *Qp = nullptr, *Q1 = nullptr, *Q2 = nullptr, 
       *Q3_p = nullptr;


    _F *BETA = nullptr; dcomplex *ALPHA = nullptr;
    _F *RMAT = nullptr, *PHI = nullptr, *PSI = nullptr;
    _F *VSR  = nullptr, *VSL = nullptr, *MA  = nullptr, *MB = nullptr;
    _F *VSCR = nullptr, *VSCR2 = nullptr;

    if( isRoot ) {

      // V  = [V W S1...Sm P]
      V = this->memManager_.template malloc<_F>(NMSS);
      W = V + NNR;
      S = W + NNR;
      P = S + this->m * NNR;

      // AV = [AV AW AS1...ASm AP]
      AV = this->memManager_.template malloc<_F>(NMSS);
      AW = AV + NNR;
      AS = AW + NNR;
      AP = AS + this->m * NNR;


      // Q = [Q Q1 Q2 Q3]
      Q = this->memManager_.template malloc<_F>(NMSS);
      Qp = Q + NNR;

      Q1   = Qp;
      Q2   = Q1 + NNR;
      Q3_p = Q2 + this->m * NNR;





      // EVAL(I) = ALPHA(I) / BETA(I)
      ALPHA = this->memManager_.template malloc<dcomplex>(MSS);
      BETA  = this->memManager_.template malloc<_F>(MSS);
      //double *RITZ  = this->memManager_.template malloc<double>(MSS);
      //double *nRITZ = this->memManager_.template malloc<double>(nR);

      RMAT = this->memManager_.template malloc<_F>(MSS2); 
      PHI  = this->memManager_.template malloc<_F>(MSS2); 
      PSI  = this->memManager_.template malloc<_F>(MSS2); 

      VSR = this->memManager_.template malloc<_F>(MSS2);
      VSL = this->memManager_.template malloc<_F>(MSS2);
      MA  = this->memManager_.template malloc<_F>(MSS2);
      MB  = this->memManager_.template malloc<_F>(MSS2);

    } // ROOT only






    // Need a NNR scratch on all processes
    VSCR  = this->memManager_.template malloc<_F>(NNR); 
    VSCR2 = this->memManager_.template malloc<_F>(NNR); 








    if( isRoot ) {

      // Initailize V as Guess, if no Guess set, init to identity
      if( Guess ) std::copy_n(Guess,NNR,V);
      else {
        std::fill_n(V, NNR, 0.);
        for(auto i = 0ul; i < nR; i++) V[i * (N+1)] = 1.0;
      }


      // V <- QR(V)
      QR(N,nR,V,N,this->memManager_);

    } // ROOT only





    // Sync processes
    MPI_Barrier(this->comm_);


    _F * VSend  = isRoot ? V  : nullptr;
    _F * AVRecv = isRoot ? AV : nullptr;


    // AV <- A * V
    this->linearTrans_(nR,VSend,AVRecv);

    // Sync processes
    MPI_Barrier(this->comm_);





    double nrmA = 0.;

    if( isRoot ) {


      // Q <- AV - sig * V 
      std::copy_n(AV, NNR, Q);
      blas::axpy( NNR, -sigma, V, 1, Q, 1 );

      // Q <- QR(Q)
      QR(N,nR,Q,N,this->memManager_);
        



      // PHI <- Q**H * AV
      // PSI <- Q**H * V
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nR,nR,N,_F(1.),Q,N,AV,N,_F(0.),PHI,nR);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nR,nR,N,_F(1.),Q,N,V,N,_F(0.),PSI,nR);


      // VSR, VSL, ALPHA, BETA <- ORDQZ(PHI,PSI,sigma)
      OrdQZ2('V','V',nR,PHI,nR,PSI,nR,ALPHA,BETA,hardLimD,sigmaD,
             VSL,nR,VSR,nR,this->memManager_);


      // Update eigs with ALPHA/BETA
      for (auto i = 0ul; i < nR; i++)
        this->eigVal_[i] = ALPHA[i] / BETA[i];

      // Print warinings if needed
      if( std::is_same<double, _F>::value ) {

        for( auto i = 0ul; i < nR; i++ )
          if( std::abs(std::imag(this->eigVal_[i])) > 1e-10 )
            std::cout << "  *** WARNING: INTERMEDIATE COMPLEX EIGENVALUE "
                      << "INCURRED ***" << std::endl;

      }

      /*
      // Adaptive sigma
      double limDiff = std::numeric_limits<double>::infinity();
      double omega0 = 0.;
      double midPt = 0.;
      bool belowTh = false;
      double maxW = -std::numeric_limits<double>::infinity();
      */


      // Get MA and MB ('Q-free' Schur form)
      getTriU(nR,PHI,nR,PSI,nR,MA,nR,MB,nR);


      // Right and Left Schur Vectors
      // VR  <- VR  * VSR
      // VL  <- VL  * VSL
      // AVR <- AVR * VSR
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),V,N,VSR,nR,_F(0.),VSCR,N);
      std::copy_n(VSCR,NNR,V);
      
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),Q,N,VSL,nR,_F(0.),VSCR,N);
      std::copy_n(VSCR,NNR,Q);
      
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),AV,N,VSR,nR,_F(0.),VSCR,N);
      std::copy_n(VSCR,NNR,AV);

      // Form initial residuals in W
      // W = AV * MB - V * MA
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.) ,AV,N,MB,nR,_F(0.),W,N);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(-1.),V ,N,MA,nR,_F(1.),W,N);


      // Get Residual Norms
      nrmA =  lapack::lange(lapack::Norm::Fro,N,nR,AV,N);
      getResidualNorms(N,nR,W,RelRes,this->eigVal_,nrmA);




      std::cout << "  * Initial Matrix Norm Estimate = " << nrmA << std::endl;

      std::cout << "\n\n  * Initial Eigenvalues\n\n";
      std::cout << std::setprecision(10) << std::scientific;
      for(auto iRt = 0ul; iRt < nR; iRt++) {

        std::cout << "    IRt = " << std::setw(5) << std::left << iRt;

        std::cout << std::right;
        if( std::is_same<double, _F>::value ) std::cout << "EigVal = "; 
        else                                  std::cout << "Re(EigVal) = ";


        std::cout << std::setw(20) << std::real(this->eigVal_[iRt]);

        if( std::is_same<dcomplex,_F>::value ) {
          std::cout << "    Im(EigVal) = ";
          std::cout << std::setw(20) << std::imag(this->eigVal_[iRt]);
        }


        std::cout << "    RelResNorm = ";
        std::cout << std::setw(20) << RelRes[iRt];


        std::cout << std::endl;
      }


      std::cout << "\n\n\n" << std::endl;

      std::cout << "  * Starting GPLHR Iterations\n" << std::endl;



    } // ROOT only

    size_t iter = 0;

    // ****************************
    // ** Begin GPLHR iterations **
    // ****************************
    for( iter = 0; iter < this->maxMicroIter_; iter++ ) {

      ProgramTimer::tick("Diagonalize Iter");

      if( isRoot ) {

        std::cout << "    GPLHRIter " << std::setw(5) << iter+1;
    
        // V, RMAT <- QR(V)
        QR(N,nR,V,N,RMAT,nR,this->memManager_);

        // AV <- X : [X * RMAT = AV]
        blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
          N,nR,_F(1.),RMAT,nR,AV,N);

        // Q <- QR(Q)
        QR(N,nR,Q,N,this->memManager_);


        // W = (I - V * V**H) * T * (I - V * V**H) * W 
        newSMatrix(N,nR,V,N,V,N,W,N,RMAT,nR);

        // W <- QR(W)
        QR(N,nR,W,N,this->memManager_);

      } // ROOT only


      // Sync processes
      MPI_Barrier(this->comm_);

      // AW <- A * W
      auto LTst = tick();

      _F * WSend  = isRoot ? W  : nullptr;
      _F * AWRecv = isRoot ? AW : nullptr;
      this->linearTrans_(nR,WSend,AWRecv);


      double LTdur = tock(LTst);

      // Sync processes
      MPI_Barrier(this->comm_);
    

      // Form S-blocks
        
      // S(0) = W, AS(0) = AW
      for( auto k = 1; k <= this->m; k++ ) {

        _F *SSend = nullptr, *ASRecv = nullptr;

        if( isRoot ) {

          // S(k-1) / AS(k-1)
          _F *Sprev  = W  + (k-1) * NNR;
          _F *ASprev = AW + (k-1) * NNR;

          // S(k) / AS(k)
          _F *Scur  = W  + k * NNR;
          _F *AScur = AW + k * NNR;

          SSend  = Scur;
          ASRecv = AScur;


          // S(k) = AS(k-1) * MB - S(k-1) * MA
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.) ,ASprev,N,MB,nR,_F(0.),Scur,N);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(-1.),Sprev ,N,MA,nR,_F(1.),Scur,N);

          // S(k) = (I - V * V**H) * T * (I - V * V**H) * S(k-1) 
          newSMatrix(N,nR,V,N,V,N,Scur,N,RMAT,nR);


          // Project out previous S's
          // S(k) = (I - S(l) * S(l)**H) * S(k)  for l = [0,k)
          for(auto l = 0; l < k; l++) 
            halfProj(N,nR,W + l * NNR,N,Scur,N,RMAT,nR);
       

          // S(k) = QR(S(k))
          QR(N,nR,Scur,N,this->memManager_);

        } else { // ROOT only

          SSend  = nullptr;
          ASRecv = nullptr;

        }



        // Sync processes
        MPI_Barrier(this->comm_);

        // AS(k) = A * S(k) 
        LTst = tick();

        this->linearTrans_(nR,SSend,ASRecv);

        LTdur += tock(LTst);

        // Sync processes
        MPI_Barrier(this->comm_);

      }




      const size_t nQp = iter ? 2 + this->m : 1 + this->m;
      const size_t nQ  = nQp + 1;

      if( isRoot ) {

        // Conjugate direction orthogonalization
        if (iter) {
          
          // P = (I - V * V**H) * P
          // P = (I - W * W**H) * P
          // P = (I - S * S**H) * P
          halfProj2(N,     nR,V,N,AV,N,P,N,AP,N,RMAT,nR  ); // Project out V
          halfProj2(N,     nR,W,N,AW,N,P,N,AP,N,RMAT,nR  ); // Project out W
          halfProj2(N,M_NR,nR,S,N,AS,N,P,N,AP,N,RMAT,M_NR); // Project out S

          // P, RMAT <- QR(P)
          QR(N,nR,P,N,RMAT,nR,this->memManager_);

          // AP <- X : [X * RMAT = AP]
          blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
            N,nR,_F(1.),RMAT,nR,AP,N);

        }

        _F * Q3 = iter ? Q3_p : nullptr;

        // Q' = (A - sigma * I) [W, S, [P]]
        std::copy_n(AW,nQp * NNR, Qp);
        blas::axpy( nQp * NNR, -sigma, W, 1, Qp, 1);




        // Q1 = (I - Q * Q**H) * Q1
        // Q1 = QR(Q1)
        halfProj(N,nR,Q,N,Q1,N,RMAT,nR);
        QR(N,nR,Q1,N,this->memManager_);

        // Q2 = (I - Q  * Q**H ) * Q2
        // Q2 = (I - Q1 * Q1**H) * Q2
        // Q2 = QR(Q2)
        halfProj(N,nR,M_NR,Q ,N,Q2,N,RMAT,nR);
        halfProj(N,nR,M_NR,Q1,N,Q2,N,RMAT,nR);
        QR(N,M_NR,Q2,N,this->memManager_);


        if( Q3 ) {
          // Q3 = (I - Q  * Q**H ) * Q3
          // Q3 = (I - Q1 * Q1**H) * Q3
          // Q3 = (I - Q2 * Q2**H) * Q3
          // Q3 = QR(Q3)
          halfProj(N,     nR,Q ,N,Q3,N,RMAT,nR);
          halfProj(N,     nR,Q1,N,Q3,N,RMAT,nR);
          halfProj(N,M_NR,nR,Q2,N,Q3,N,RMAT,M_NR);
          QR(N,nR,Q3,N,this->memManager_);
        }



        // PHI <- Q**H * AV
        // PSI <- Q**H * V
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nQ*nR, nQ*nR, N, _F(1.), Q,N, AV,N, _F(0.), PHI, nQ*nR); 
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nQ*nR, nQ*nR, N, _F(1.), Q,N, V ,N, _F(0.), PSI, nQ*nR); 

        // VSR, VSL, ALPHA, BETA <- ORDQZ(PHI,PSI,sigma)
        OrdQZ2('V','V',nQ*nR,PHI,nQ*nR,PSI,nQ*nR,ALPHA,BETA,hardLimD,sigmaD,
          VSL,nQ*nR,VSR,nQ*nR,this->memManager_);


        // VSRt is thick-restart
        _F * VSRt = VSR + nR * (nQ * nR);

        _F * VSR_V = VSR;
        _F * VSR_W = VSR_V + nR;
        _F * VSR_S = VSR_W + nR;
        _F * VSR_P = VSR_S + M_NR;

        _F * VSL_V = VSL;
        _F * VSL_W = VSL_V + nR;
        _F * VSL_S = VSL_W + nR;
        _F * VSL_P = VSL_S + M_NR;

        _F * VSRt_V = VSRt;
        _F * VSRt_W = VSRt_V + nR;
        _F * VSRt_S = VSRt_W + nR;
        _F * VSRt_P = VSRt_S + M_NR;

        // VSCR = V * VSR_V + W * VSR_W + S * VSR_S + P * VSR_P
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),V,N,VSR_V,nQ*nR,_F(0.),VSCR,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),W,N,VSR_W,nQ*nR,_F(1.),VSCR,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,M_NR,_F(1.),S,N,VSR_S,nQ*nR,_F(1.),VSCR,N);
        if( iter )
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),P,N,VSR_P,nQ*nR,_F(1.),VSCR,N);


        // VSCR2 = V * VSRt_V + W * VSRt_W + S * VSRt_S + P * VSRt_P
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),V,N,VSRt_V,nQ*nR,_F(0.),VSCR2,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),W,N,VSRt_W,nQ*nR,_F(1.),VSCR2,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,M_NR,_F(1.),S,N,VSRt_S,nQ*nR,_F(1.),VSCR2,N);
        if( iter )                          
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),P,N,VSRt_P,nQ*nR,_F(1.),VSCR2,N);

        // V = VSCR
        // P = VSCR2
        std::copy_n(VSCR ,NNR,V);
        std::copy_n(VSCR2,NNR,P);


        // VSCR = AV * VSR_V + AW * VSR_W + AS * VSR_S + AP * VSR_P
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),AV,N,VSR_V,nQ*nR,_F(0.),VSCR,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),AW,N,VSR_W,nQ*nR,_F(1.),VSCR,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,M_NR,_F(1.),AS,N,VSR_S,nQ*nR,_F(1.),VSCR,N);
        if( iter )
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),AP,N,VSR_P,nQ*nR,_F(1.),VSCR,N);

        // VSCR2 = AV * VSRt_V + AW * VSRt_W + AS * VSRt_S + AP * VSRt_P
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),AV,N,VSRt_V,nQ*nR,_F(0.),VSCR2,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),AW,N,VSRt_W,nQ*nR,_F(1.),VSCR2,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,M_NR,_F(1.),AS,N,VSRt_S,nQ*nR,_F(1.),VSCR2,N);
        if( iter )                           
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),AP,N,VSRt_P,nQ*nR,_F(1.),VSCR2,N);


        // AV = VSCR
        // AP = VSCR2
        std::copy_n(VSCR ,NNR,AV);
        std::copy_n(VSCR2,NNR,AP);

        // VSCR = Q * VSR_V + Q1 * VSR_W + Q2 * VSR_S + Q3 * VSR_P
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),Q ,N,VSL_V,nQ*nR,_F(0.),VSCR,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR  ,_F(1.),Q1,N,VSL_W,nQ*nR,_F(1.),VSCR,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,M_NR,_F(1.),Q2,N,VSL_S,nQ*nR,_F(1.),VSCR,N);
        if( Q3 )
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),Q3,N,VSL_P,nQ*nR,_F(1.),VSCR,N);

        // Q = VSCR
        std::copy_n(VSCR,NNR,Q);

      } // ROOT only

      // Refresh AV
      if( (iter+1) % 100 == 0 ) {

        if( isRoot ) 
          std::cout << "  * Refreshing AV at iteration " << iter+1 
            << std::endl;

        // Sync processes
        MPI_Barrier(this->comm_);

        LTst = tick();


        this->linearTrans_(nR,VSend,AVRecv);

        LTdur += tock(LTst);

        // Sync processes
        MPI_Barrier(this->comm_);
      }

      if( isRoot ) {

        // Update MA, MB
        getTriU(nR,PHI,nQ*nR,PSI,nQ*nR,MA,nR,MB,nR);
        
        // Update eigenvalues
        for(auto i = 0; i < nR; i++)
          this->eigVal_[i] = ALPHA[i] / BETA[i];
        

        // W = AV * MB - V * MA
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.) ,AV,N,MB,nR,_F(0.),W,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(-1.),V ,N,MA,nR,_F(1.),W,N);
        
        /*
        // Adaptive sigma
        limDiff = std::numeric_limits<double>::infinity();
        maxW = -std::numeric_limits<double>::infinity();

        
        for(auto i = 0ul; i < nQ*nR; i++)
          RITZ[i] = std::real(ALPHA[i] / BETA[i]);

        std::sort(RITZ,RITZ + nQ*nR);
        //prettyPrintSmart(std::cout,"RITZ SORT",RITZ,MSS,1,MSS);

        for( auto i = 0ul; i < MSS; i++ ) {    
          if( RITZ[i] > hardLimD) {
            std::copy(RITZ + i, RITZ + i + nR, nRITZ);
            omega0 = RITZ[i-1];
            maxW = nRITZ[nR-1];
            break;
          }
        }
        
        belowTh = true;
        

        
        midPt = (maxW + omega0) / 2.;
        if (midPt > hardLimD and belowTh) {
          sigmaD = midPt;
        } else {
          sigmaD = hardLimD;
        }
        std::cout << "midPt: " << midPt << std::endl;
        std::cout << "New Sigma: " << sigmaD << std::endl;

        */


        // Get Residual norms
        getResidualNorms(N,nR,W,RelRes,this->eigVal_,nrmA);
        

        // Check convergence
        isConverged = checkConv(nR,RelRes);

      } // ROOT only


      // Bcast converged
      MPIBCast(isConverged,0,this->comm_);





      ProgramTimer::tock("Diagonalize Iter");


      if( isRoot ) {

        auto GPLHRdur = ProgramTimer::getDurationTotal<CQSecond>(
          "Diagonalize Iter").count();

        double perLT = LTdur * 100 / GPLHRdur;

        std::cout << "  DURATION = " << std::setprecision(8) << GPLHRdur 
          << " s  ( " << perLT << " % LT )" << std::endl;

        std::cout << std::setprecision(10) << std::scientific;
        for(auto iRt = 0ul; iRt < nR; iRt++) {

          std::cout << "      IRt = " << std::setw(5) << std::left << iRt;

          std::cout << std::right;
          if( std::is_same<double, _F>::value ) std::cout << "EigVal = "; 
          else                                  std::cout << "Re(EigVal) = ";


          std::cout << std::setw(20) << std::real(this->eigVal_[iRt]);

          if( std::is_same<dcomplex,_F>::value ) {
            std::cout << "    Im(EigVal) = ";
            std::cout << std::setw(20) << std::imag(this->eigVal_[iRt]);
          }


          std::cout << "    RelResNorm = ";
          std::cout << std::setw(20) << RelRes[iRt];


          std::cout << std::endl;
        }

        std::cout << "\n" << std::endl;

      } // ROOR only



      if( isConverged ) break; // Break loop on convergence


      

    } // end for


    MPI_Barrier(this->comm_); // Sync processes

    


    if( isRoot ) {

      // Reconstruct Eigen vectors
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nR,nR,N,_F(1.),V,N,AV,N,_F(0.),PSI,nR);
      GeneralEigenSymm('N','V',nR,PSI,nR,ALPHA,VSL,nR,VSR,nR);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nR,_F(1.),V,N,VSR,nR,_F(0.),this->VR_,N);

    }


    // Free Scratch space
    
    if(V)      this->memManager_.free(V);
    if(AV)     this->memManager_.free(AV);
    if(Q)      this->memManager_.free(Q);
    if(ALPHA)  this->memManager_.free(ALPHA);
    if(BETA)   this->memManager_.free(BETA);
    if(RMAT)   this->memManager_.free(RMAT);
    if(PHI)    this->memManager_.free(PHI);
    if(PSI)    this->memManager_.free(PSI);
    if(VSL)    this->memManager_.free(VSL);
    if(VSR)    this->memManager_.free(VSR);
    if(MA)     this->memManager_.free(MA);
    if(MB)     this->memManager_.free(MB);
    if(VSCR)   this->memManager_.free(VSCR);
    if(VSCR2)  this->memManager_.free(VSCR2);

    if( isRoot ) {

      std::cout << "  * ";
      if( isConverged )
        std::cout << "GPLHR Converged in " << iter+1 << " Iterations" 
                  << std::endl;
      else
        std::cout << "GPLHR Failed to Converged in " << iter+1 << " Iterations" 
                  << std::endl;

    }


    return isConverged;

  }



  /*
    Calculate relative norm for each residual vector 
    Assuming B = Identity.
    Assuming eigenvalues are finite.
  */
  template <typename _F>
  void GPLHR<_F>::getResidualNorms(size_t N, size_t nR, _F *WMAT, double *RelRes, dcomplex *LAMBDA, double nrmA) { 

    ROOT_ONLY(this->comm_);


    for (auto i = 0; i < nR; i++) {

      double nrmI = blas::nrm2(N,WMAT + i*N,1);
      RelRes[i] = nrmI / (nrmA + std::abs(LAMBDA[i]));

    }

  }

  template <typename _F>
  bool GPLHR<_F>::checkConv(size_t nR, double *RelRes) { 

    const bool conv = std::none_of(RelRes,RelRes + nR,
        [&]( double x ) -> bool { 
          return x > this->convCrit_ or std::isnan(x);
        });

    return conv;
  }



  template <typename _F>
  void GPLHR<_F>::getTriU(size_t N, _F *TRIUA, size_t LDTRIUA, _F *TRIUB, 
    size_t LDTRIUB, _F *MA, size_t LDMA, _F *MB, size_t LDMB){

    ROOT_ONLY(this->comm_);
#if 0

    // G(TRIUB) = TRIUA * CA + TRIUB * CB
    for(auto i = 0ul; i < N; i++){ 
      blas::scal(N,CB[i*(N+1)],                TRIUB + i*ldm,1);
      blas::axpy( N,CA[i*(N+1)],TRIUA + i*ldm,1,TRIUB + i*ldm,1);
    }

    // TRIUA = inv(TRIUB) * TRIUA
    blas::trsm(blas::Layout::ColMajor,blas::Side::Left,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
      N,N,_F(1.),TRIUB,ldm,TRIUA,ldm);

    //for(auto i = 0ul; i < N; i++){ 

    //  std::copy_n(TRIUA + i*ldm, N, MA + i*N);
    //  std::copy_n(TRIUA + i*ldm, N, MB + i*N);

    //  blas::scal(N, CB[i*(N+1)],MA + i*N,1);
    //  blas::scal(N,-CA[i*(N+1)],MB + i*N,1);

    //  MB[ i*(N + 1) ] += 1.;

    //}

#else

    // Initialize CA and CB as identity FIXME: memory
    _F *CA    = this->memManager_.template malloc<_F>(N*N);
    _F *CB    = this->memManager_.template malloc<_F>(N*N);
    _F *Ident = this->memManager_.template malloc<_F>(N*N);
    _F *G     = this->memManager_.template malloc<_F>(N*N);
    std::fill_n(CA,N*N,_F(0.)); 
    std::fill_n(CB,N*N,_F(0.)); 
    std::fill_n(Ident,N*N,_F(0.)); 
    for (auto i = 0; i < N; i++) {
      CA[i + i*N]    = _F(1.0);
      CB[i + i*N]    = _F(1.0);
      Ident[i + i*N] = _F(1.0);
    }

    // Create CA and CB
    for (auto i = 0; i < N; i++) {
      if (std::abs(TRIUA[i + i*LDTRIUA]) >= std::abs(TRIUB[i + i*LDTRIUB])) {
        CA[i + i*N] = (_F(1.) - TRIUB[i + i*LDTRIUB]) / TRIUA[i + i*LDTRIUA];
      } else {
        CA[i + i*N] = _F(0.);
        CB[i + i*N] = _F(1.) / TRIUB[i + i*LDTRIUA];
      }
    }

    // Form G
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),TRIUA,LDTRIUA,CA,N,_F(0.),G,N);  
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),TRIUB,LDTRIUB,CB,N,_F(1.),G,N);
    
    // Form MA
    blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
      N,N,_F(1.),G,N,CB,N);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),CB,N,TRIUA,LDTRIUA,_F(0.),MA,LDMA);

    // Form MB   
    blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
      N,N,_F(1.),G,N,CA,N);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),CA,N,TRIUA,LDTRIUA,_F(0.),MB,LDMB);
    blas::axpy(N*N,_F(-1.),MB,1,Ident,1);
    std::copy_n(Ident, N*N, MB);


    // Free SCR Mem
    this->memManager_.free(CA);
    this->memManager_.free(CB);
    this->memManager_.free(G);
#endif


  }








  /**
   *  \brief Projects out set of vectors, V, from another set of vectors, S
   *
   *  S = (I - V * V**H) * S
   */
  template <typename _F>
  void GPLHR<_F>::halfProj(size_t N, size_t nV, size_t nS, _F *V, size_t LDV,
    _F *S, size_t LDS, _F *SCR, size_t LDSCR) {

    if( LDSCR < nV ) CErr("nV MUST be >= LDSCR");

    ROOT_ONLY(this->comm_);

    // SCR = V**H * S
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nV,nS,N,_F(1.) ,V,LDV,S  ,LDS  ,_F(0.) ,SCR,LDSCR);

    // S = S - V * SCR
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nS,nV, _F(-1.),V,LDV,SCR,LDSCR,_F(1.),S,  LDS  );

  }

  /**
   *  \brief Projects out set of vectors, V, from another set of vectors, S
   *  and updates the linear transformed AS with the same projection on AV
   *
   *  S  = (I - V * V**H) * S
   *  AS = AS - V * V**H  * S
   */
  template <typename _F>
  void GPLHR<_F>::halfProj2(size_t N, size_t nV, size_t nS, _F *V, size_t LDV,
    _F* AV, size_t LDAV, _F *S, size_t LDS, _F* AS, size_t LDAS, _F *SCR, 
    size_t LDSCR) {

    if( LDSCR < nV ) CErr("nV MUST be >= LDSCR");

    ROOT_ONLY(this->comm_);

    // SCR = V**H * S
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nV,nS,N,_F(1.) ,V,LDV,S  ,LDS  ,_F(0.) ,SCR,LDSCR);

    // S = S - V * SCR
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nS,nV, _F(-1.),V,LDV,SCR,LDSCR,_F(1.),S,  LDS  );

    // AS = AS - AV * SCR
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nS,nV, _F(-1.),AV,LDAV,SCR,LDSCR,_F(1.),AS,LDAS  );

  }



  /**
   *  \brief Forms a new S matrix in the Krylov-Arnoldi space
   *
   *  S = (I - V * V**H) * T * (I - U * U**H ) * S
   *
   *  where T is the preconditioner
   */
  template <typename _F>
  void GPLHR<_F>::newSMatrix(size_t N, size_t nR, _F *V, size_t LDV, _F *Q, 
    size_t LDQ, _F *S, size_t LDS, _F *SCR, size_t LDSCR) {


    ROOT_ONLY(this->comm_);

    halfProj(N,nR,V,LDV,S,LDS,SCR,LDSCR);
    this->preCondWShift_(nR,sigma,S,S);
    halfProj(N,nR,Q,LDQ,S,LDS,SCR,LDSCR);

  }



  template <typename _F>
  void GPLHR<_F>::restart() { }

}; // namespace ChronusQ

