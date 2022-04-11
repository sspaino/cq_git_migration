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

namespace ChronusQ {


  template <typename _F>
  void GMRES<_F>::runBatch(size_t nRHS, size_t nShift, _F* RHS, _F *shifts, 
    _F* SOL, double *RHSNorm ) {

    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;

    auto topGMRES = tick();

    // Do the standard stuff...
    IterLinearSolver<_F>::runBatch(nRHS,nShift,RHS,shifts,SOL,RHSNorm);


    // Zero out the scratch allocations
    if( isRoot ) {

      std::fill_n(J_, 2 * this->mSS_ * nRHS * nShift         , 0.);
      std::fill_n(U_, this->N_ * this->mSS_ * nRHS * nShift  , 0.);  
      std::fill_n(R_, this->mSS_ * this->mSS_ * nRHS * nShift, 0.);
      std::fill_n(W_, (this->mSS_ + 1) * nRHS * nShift       , 0.);

    }






    // Construct the initial Householder reflectors in place
      
    if( isRoot ) {

      // Copy over residuals to HHR
      std::copy_n(this->RES_, nRHS * nShift * this->N_, HHR_);

      for(auto iDo = 0ul; iDo < nShift * nRHS; iDo++) {

        _F * curHHR = HHR_ + iDo*this->N_;

        _F beta = curHHR[0];
        if(std::abs(curHHR[0]) < std::numeric_limits<double>::epsilon())
          beta = this->resNorm_.back()[iDo];
        else
          beta *= this->resNorm_.back()[iDo] / std::abs(curHHR[0]);

        curHHR[0] += beta;

        // Normalize the HHR column
        double norm = Normalize(this->N_,curHHR,1);

        // Copy the HHR to the first column of U
        std::copy_n(curHHR,this->N_, U_ + iDo*this->N_*this->mSS_);

        // Apply HHR projection to residual
        W_[iDo * (this->mSS_ + 1)] = -beta;

        //std::cout << "W0 " << W_[iDo * (this->mSS_ + 1)] << "\n";

      }

    }



    // Update tracking / counting
    std::vector<size_t> mDim( nRHS * nShift, 0 );
    std::vector<bool>   solConv( nRHS * nShift, false );

    // AX Scratch
    _F * VContract  = nullptr;
    _F * AVContract = nullptr;
    if( nRHS * nShift > 1 and isRoot ) {

      VContract = 
        this->memManager_.template malloc<_F>(this->N_*nRHS*nShift);
      AVContract = 
        this->memManager_.template malloc<_F>(this->N_*nRHS*nShift);

    }


    // Start the micro iterations
    size_t maxMicroIter = this->mSS_;

    if( isRoot ) std::cout << "    * Starting GMRES iterations\n\n";
    size_t iMicro;


    for(iMicro = 0; iMicro < maxMicroIter; iMicro++) {

      ProgramTimer::tick("Lin Solve Iter");

      auto topMicro = tick();

      for(auto iDo = 0; iDo < nRHS * nShift; iDo++)
        if( not solConv[iDo] ) mDim[iDo]++;

      if( isConverged ) break;

      size_t nConv(0), nNotConv(0);
      if( isRoot ) {
        nNotConv = std::count(solConv.begin(),solConv.end(),false);
        nConv = solConv.size() - nNotConv;
      }


      if( isRoot ) // Only root process
      for(auto iDo = 0ul; iDo < nRHS * nShift; iDo++) {

        if( solConv[iDo] ) continue;

        _F * curV   = this->V_ + iDo*this->N_;
        _F * curHHR = HHR_     + iDo*this->N_;

        std::copy_n(curHHR,this->N_,curV);
        blas::scal(this->N_, -2. * SmartConj(curHHR[iMicro]), curV, 1);
        curV[iMicro] += 1.; 

        if( iMicro > 0 )
        for(int k = iMicro-1; k >= 0; k--) {

          _F * curU = U_ + (k+iDo*this->mSS_)*this->N_;

          _F inner = blas::dot(this->N_, curU, 1, curV, 1);
          blas::axpy(this->N_, -2.*inner, curU, 1, curV, 1);

        }

        // Explicitly normalize V column
        Normalize(this->N_,curV,1);

      }

      bool CopyBuffer = isRoot and VContract and (nConv != 0);

      size_t nContract = nShift * nRHS;
      _F * VSend  = this->V_;
      _F * AVRecv = this->AV_;
      if( CopyBuffer ) {

        VSend  = VContract;
        AVRecv = AVContract;

        nContract = 0;
        for(auto iDo = 0; iDo < nRHS * nShift; iDo++)
        if( not solConv[iDo] ) {

          std::copy_n(this->V_ + iDo*this->N_, this->N_, 
            VContract + nContract*this->N_);
          nContract++;

        }


      }

      if( MPISize(this->comm_) > 1 ) MPIBCast(nContract,0,this->comm_);

      // Form (A-sB)X product and precondition product
      MPI_Barrier(this->comm_);
      auto topLT = tick();

      this->linearTrans_(nContract, VSend, AVRecv);

      double durLT = tock(topLT);
      MPI_Barrier(this->comm_);

      if( CopyBuffer ) {

        size_t iContract = 0;
        for(auto iDo = 0; iDo < nRHS * nShift; iDo++)
        if( not solConv[iDo] ) {

          std::copy_n(AVRecv + iContract*this->N_, this->N_, 
            this->AV_ + iDo*this->N_);
          iContract++;

        }

      }

      if( isRoot )
      for(auto iS = 0; iS < nShift; iS++) {
        _F *curV  = this->V_  + iS*nRHS*this->N_;
        _F *curAV = this->AV_ + iS*nRHS*this->N_;
        this->shiftVec_     (nRHS, -shifts[iS], curV , curAV);
        this->preCondWShift_(nRHS, shifts[iS], curAV, curAV);
      }


      // Copy AV -> V
      if( isRoot ) 
        std::copy_n(this->AV_,nRHS * nShift * this->N_,this->V_);

      
      if( isRoot ) this->resNorm_.emplace_back(this->resNorm_.back());

      if( isRoot )
      for(auto iDo = 0ul; iDo < nRHS * nShift; iDo++) {

        if( solConv[iDo] ) continue;

        _F * curV   = this->V_ + iDo*this->N_;
        _F * curHHR = HHR_     + iDo*this->N_;
        _F * curJ   = J_       + iDo*2*this->mSS_;
        _F * curW   = W_       + iDo*(this->mSS_ + 1);
        _F * curR   = R_       + iDo*this->mSS_*this->mSS_;

        for(auto k = 0; k <= iMicro; k++) {

          _F * curU = U_ + (k+iDo*this->mSS_)*this->N_;

          _F inner = blas::dot(this->N_, curU, 1, curV, 1);
          blas::axpy(this->N_, -2.*inner, curU, 1, curV, 1);

        }


        _F * curU = U_ + (iMicro + 1 + iDo*this->mSS_)*this->N_;

        // Determine next projector
        if( iMicro < maxMicroIter - 1 ) {

          std::copy_n(curV,this->N_,curHHR);
          std::fill_n(curHHR,iMicro+1,0.);

          _F alpha = blas::nrm2(this->N_,curHHR,1);
          if( std::abs(alpha) > 1e-10 ) {

            if (std::abs(curV[iMicro+1]) > 0)
            alpha *= curV[iMicro+1] / std::abs(curV[iMicro+1]);

            curHHR[iMicro+1] += alpha;
            double norm = Normalize(this->N_,curHHR,1);
            
            std::copy_n(curHHR,this->N_,curU);

            std::fill_n(curV + iMicro + 2, this->N_ - iMicro - 2, 0.);
            curV[iMicro+1] = -alpha;

          }

        }


        // Apply Given's rotation
        if( iMicro > 0 ) 
        for(auto j = 0ul; j < iMicro; j++) {

          auto tmp = curV[j];
          
          curV[j] =
            SmartConj(curJ[2*j    ]) * tmp + 
            SmartConj(curJ[2*j + 1]) * curV[j+1];

          curV[j + 1] = curJ[2*j] * curV[j+1] - curJ[2*j + 1] * tmp;

        }

        if( iMicro < maxMicroIter - 1 ) {

          double rho = blas::nrm2(2,curV + iMicro,1);

          curJ[2*iMicro]   = curV[iMicro]   / rho;
          curJ[2*iMicro+1] = curV[iMicro+1] / rho;

          curW[iMicro+1]   = -curJ[2*iMicro+1]         * curW[iMicro];
          curW[iMicro]     = SmartConj(curJ[2*iMicro]) * curW[iMicro];

          curV[iMicro]     = rho;
          curV[iMicro + 1] = 0.;

        }

        std::copy_n(curV,this->mSS_,curR + iMicro*this->mSS_);

        double normR = std::abs(curW[iMicro+1]);
        this->resNorm_.back()[iDo] = normR;


      } // loop over vectors


      double durMicro = tock(topMicro);



      if( isRoot ) {
        std::cout << "      GMRESIter " << std::setw(5) << iMicro + 1;
        if( nRHS * nShift > 1 ) {

          std::cout << "  :  ";
          std::cout << "  IHAVE   = " << nConv ;
          std::cout << "  NUPDATE = " << nNotConv;
          std::cout << "  DURATION = " << std::scientific 
            << durMicro << " s ( " << std::fixed 
            << durLT * 100. / durMicro << "% LT )\n";

        }
      }


      if( isRoot )
      for(auto iDo = 0; iDo < nShift*nRHS; iDo++) {

        double normR = this->resNorm_.back()[iDo];
        if(nRHS * nShift > 1) 
          std::cout << "        iDo = " << std::setw(6) << iDo << "  ";

        std::cout << std::scientific << std::setprecision(8);
        std::cout << "  ResNorm = " << normR;
        std::cout << "  RelResNorm = " << normR / RHSNorm[iDo % nRHS]; 
        std::cout << "\n";

      }
      if( isRoot and (nRHS * nShift > 1) ) std::cout << "\n";


      // Evaluate convergence
      if( isRoot ) {

        auto OldSolConv(solConv);
        size_t nNewConv(0);
        for(auto iDo = 0; iDo < nRHS*nShift; iDo++) {
          solConv[iDo] = 
            (this->resNorm_.back()[iDo] / RHSNorm[iDo % nRHS]) < 
            this->convCrit_;

          if( nRHS * nShift > 1 and solConv[iDo] and not OldSolConv[iDo] ) {

            nNewConv++;
            std::cout << "          "
              <<  "*** CEASING TO UPDATE IDO = " << std::setw(4) << iDo 
              << " ***\n";

          }
        }

        if(nNewConv) std::cout << "\n";



        isConverged = std::all_of(solConv.begin(),solConv.end(),
            [&](bool x){ return x; });

      }

      // Broadcast the convergence result to all the mpi processes
      if(MPISize(this->comm_) > 1) MPIBCast(isConverged,0,this->comm_);

      ProgramTimer::tock("Lin Solve Iter");

    } // Micro iterations


    double durGMRES = tock(topGMRES);

    // Cleanup memory
    if( VContract )  this->memManager_.free(VContract);
    if( AVContract ) this->memManager_.free(AVContract);

    if( isRoot )
      std::cout << "\n    * GMRES Converged in " << iMicro  
        << " Iterations (" << durGMRES << " s) \n\n";


    // Reconstruct the solution for the batch
          
    if( isRoot )
    for(auto iDo = 0; iDo < nShift*nRHS; iDo++) {

      _F *curR   = R_   + iDo * this->mSS_ * this->mSS_;
      _F *curW   = W_   + iDo * (this->mSS_ + 1);
      _F *curU   = U_   + iDo * this->mSS_ * this->N_;
      _F *curHHR = HHR_ + iDo * this->N_;

      int nMicro = mDim[iDo];

      // Linear Solve
      int64_t* IPIV = this->memManager_.template malloc<int64_t>(nMicro);
      lapack::gesv(nMicro,1,curR,this->mSS_,IPIV,curW,this->mSS_+1);
      this->memManager_.free(IPIV);

      std::copy_n(curU + (nMicro-1)*this->N_, this->N_, curHHR);

      _F fact = curW[nMicro-1] * SmartConj(curU[(nMicro-1)*(this->N_+1)]);
      blas::scal(this->N_,-2.*fact,curHHR,1);
      curHHR[nMicro-1] += curW[nMicro-1];

      if( nMicro > 1 )
      for(int k = nMicro - 2; k >= 0; k--) {

        curHHR[k] += curW[k];
        _F inner = blas::dot(this->N_,curU + k*this->N_,1,curHHR,1);
        blas::axpy(this->N_,-2.*inner,curU + k*this->N_,1,curHHR,1);

      }

      // FIXME: assumes 0 guess
      _F *curSOL = SOL + iDo * this->N_;
      std::copy_n(curHHR,this->N_,curSOL);

    }


  }




};

