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
#ifndef __INCLUDED_DAVIDSON_HPP__
#define __INCLUDED_DAVIDSON_HPP__

#include <itersolver.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/ortho.hpp>
#include <cqlinalg/eig.hpp>
#include <cerr.hpp>

// #define DEBUG_DAVIDSON


namespace ChronusQ {


  template <typename _F>
  bool Davidson<_F>::runMicro() {
    
    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;
    
    if( isRoot ) {
      std::cout  << "\n\n";
      std::cout  << "  * Davidson Settings:\n";
      std::cout  << "    * Right Eigenvector is Requested. \n";
      
      if(this->DoLeftEigVec) {
        std::cout  << "    * Left Eigenvector is Requested.\n";
        CErr("Do Left Eig Vec is not implemented yet");
      }
      
      if (this->EnergySpecific) {          
        std::cout<< "    * Use Energy Specific:           " << this->EnergySpecific << "\n"
                 << "      * Number of Low  Energy Roots = " << this->nHighERoots    << "\n"
                 << "      * Number of High Energy Roots = " << this-> nLowERoots    << "\n";
        if(this->adaptiveERef) {
          std::cout<< "      * Use Ground State Energy in each iteration as Reference \n";
        } else {
          std::cout<< "      * Energy Referene  = " << this->EnergyRef    << "\n";
    }
        CErr("Energy specific is not implemented yet");
      }          
      std::cout << "\n\n" << std::endl;
    }
    
    std::cout << std::setprecision(10) << std::scientific;
    
    const size_t N   = this->N_;
    const size_t MSS = this->mSS_;
    const size_t nR  = this->nRoots_;
    
    const size_t MSS2 = MSS * MSS;
    const size_t NMSS = N * MSS;
    const size_t NNR  = N * nR;
    
    const size_t nG  = this->nGuess_;
    const size_t NNG = N * nG;  
    
    // Initialize pointers
    _F *VR  = nullptr,  *XR  = nullptr, *XRPrev = nullptr; 
    _F *VL  = nullptr,  *XL  = nullptr; // for left search space    
    _F *AVR = nullptr,  *Ovlp = nullptr; // overlap is using right eigenvectors  
    _F *R   = nullptr,  *S  = nullptr; // Scratch space for residue and perturbed vector
    _F *SubA = nullptr;
    _F *SCR = nullptr;

    dcomplex *Eig = nullptr, *EPrev = nullptr;
    
    if( isRoot ) {
       
      VR     = this->memManager_.template malloc<_F>(NMSS); 
      AVR    = this->memManager_.template malloc<_F>(NMSS); 
      XR     = this->memManager_.template malloc<_F>(MSS2);
      XRPrev = this->memManager_.template malloc<_F>(MSS2);
      SubA   = this->memManager_.template malloc<_F>(MSS2);
      Eig    = this->memManager_.template malloc<dcomplex>(MSS);
      EPrev  = this->memManager_.template malloc<dcomplex>(MSS);
      Ovlp   = this->memManager_.template malloc<_F>(MSS2); 
      R      = this->memManager_.template malloc<_F>(NNG); 
      S      = this->memManager_.template malloc<_F>(NNG);
      SCR    = this->memManager_.template malloc<_F>(MSS2);
      
      if(this->DoLeftEigVec) {
        VL   = this->memManager_.template malloc<_F>(NMSS); 
        XL   = this->memManager_.template malloc<_F>(MSS2);
      }

    } // Root Only

    // generate guess
    if( isRoot ) {
      // Initailize VR as Guess, if no Guess set, 
      if( Guess ) std::copy_n(Guess,NNG,VR);
      else {
        std::cout << "  * use unit vector guess" << std::endl;
        std::fill_n(VR, NNG, 0.);
        for(auto i = 0ul; i < nG; i++) VR[i * (N+1)] = 1.0;
      } // right vector guess
      
      // left vector guess 
      if(this->DoLeftEigVec) { 
        // VL and VR should be biothogonalized    
        CErr("Do Left Eig Vec is not implemented yet");
      }
    } // Root Only
    
    // Sync processes
    MPI_Barrier(this->comm_);
    
    if( isRoot ) {
      std::cout << "\n\n  * Starting Davidson Iterations" << std::endl;
    } // Root Only
    
    // variables during iteration
    std::vector<bool> SiConv(nG); // state_i converged?
    size_t iter   = 0;
    size_t nDo    = nG;
    size_t nExam  = nG;
    size_t nVPrev = 0;     // number of vectors at previous iteration 
    size_t nVCur  = nG;    // number of vectors at current iteration 
    char   JOBVL  = this->DoLeftEigVec ? 'V': 'N';
    double VecNear = 0.1;  // criterion to determine if two vector are similar by their overlap 
    double VecConv = this->convCrit_;    
    double EConv   = VecConv * 0.01;   

    // ****************************
    // ** Begin Davidson iterations **
    // ****************************
    for( iter = 0; iter < this->maxMicroIter_; iter++) {
      
      auto DavidsonSt = tick();
      
      if( isRoot ) {
        std::cout << "\n    DavidsonIter " << std::setw(5) << iter+1  
                  << ": Number of new vectors = " << std::setw(5) << nDo << std::endl;
      } // Root Only
      
      
      // Sync processes
      MPI_Barrier(this->comm_);

      auto LTst = tick();
      
      // AVR <- A * VR
      _F * VRSend  = isRoot ? VR  + N*nVPrev: nullptr;
      _F * AVRRecv = isRoot ? AVR + N*nVPrev: nullptr;
      this->linearTrans_(nDo,VRSend,AVRRecv);
      
      double LTdur = tock(LTst);
      
      // Sync processes
      MPI_Barrier(this->comm_);
       
      if( isRoot ) {
        
        // Construct submatrix of A 
        if(this->DoLeftEigVec) { 
        //SubA <- VL * AVR 
          CErr("Do Left Eig Vec is not implemented yet");
        } else {
        // SubA <- VR_\dagger * AVR 
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nVCur,nVCur,N,_F(1.),VR,N,AVR,N,_F(0.),SubA,nVCur); 
        }
    
#ifdef DEBUG_DAVIDSON
        prettyPrintSmart(std::cout,"HH Davidson SubMatix ",SubA,nVCur,nVCur,nVCur);
#endif 

        // Diagonalize SubA 
        GeneralEigenSymm(JOBVL,'V',nVCur,SubA,nVCur,Eig,XL,nVCur,XR,nVCur); 

        // swap high energy roots for energy specific
        if(this->EnergySpecific) {
          CErr("Energy Specific is not implemented yet");
        }   

        // print eigenvalues at current iteration
        std::cout << "      - Eigenvalues at the current iteration:" << std::endl;
        for(auto i = 0; i < nExam; i++) {
          std::cout << "        Root " << std::setw(5) << std::right << i << ":"
                    << std::right << std::setw(20) << std::real(Eig[i]);
              
          if( std::is_same<dcomplex,_F>::value ) {
            std::cout << " + " << std::setw(20) << std::imag(Eig[i]) << " i";
          }
              
          std::cout << std::endl;
        }
          
        // Exam Eigenvalues and eigenvectors and do mapping if iter > 0
        std::fill_n(SiConv.begin(),nExam,false);
        if( iter > 0) {
            
          // overlap = (VR XR)_old ^\dagger * (VR XR)_new
          std::fill_n(Ovlp, nVCur*nVPrev, _F(0.)); 
              for(auto i = 0ul; i < nVCur; i++) Ovlp[i + i*nVPrev] = _F(1.); 
              
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nExam,nVCur,nVPrev,_F(1.),XRPrev,nVPrev,Ovlp,nVPrev,_F(0.),SCR,nExam);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nExam,nExam,nVCur ,_F(1.),SCR,nExam,XR,nVCur,_F(0.),Ovlp,nExam);
          
          // mapping old vector to new vectors based on overlap
          std::vector<int> StMap(nExam);
          std::fill_n(StMap.begin(),nExam,-1);
          std::copy_n(Ovlp,nExam*nExam,SCR);
              
          // i -> new state, j -> old state
          int i = 0, j = 0;
          for(i = 0; i < nExam; i++) {
            auto iO = Ovlp + i*nExam;
            j = std::distance(iO, std::max_element(iO, iO+nExam, 
                  [&] (_F A, _F B) {return std::norm(A) < std::norm(B); 
                                  }));
            if( std::abs(iO[j]) > VecNear) {
              StMap[i] = j;
              for (auto k = i+1; k < nExam; k++) Ovlp[j + k*nExam] = _F(0.);
            }
          }
        
          std::copy_n(SCR,nExam*nExam,Ovlp);
          
#ifdef DEBUG_DAVIDSON
          prettyPrintSmart(std::cout,"HH Davidson States Overlap",Ovlp,nExam,nExam,nExam);
#endif 
      
          std::cout << "\n      - Comparison to the previous iteration: " << std::endl;  
          // exam eigenvectors
          // R and S as scratch space to hold full vector old and new repectively
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nExam,nVPrev,_F(1.),VR,N,XRPrev,nVPrev,_F(0.),R,N);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nExam,nVCur,_F(1.),VR,N,XR,nVCur,_F(0.),S,N);
              
          _F phase;
          _F maxDel;  // maximum differece in vectors
          for(i = 0; i < nExam; i++) {
            j = StMap[i];
            if(j < 0) {
              std::cout << "          New root" << std::setw(5) << i << " is brand new" << std::endl;
              continue; 
            } else if (j!=i) { 
              std::cout << "          New root" << std::setw(5) << i
                        << " was old root" << std::setw(5) << j << std::endl; 
            }

            phase = Ovlp[j + i*nExam];
            phase /= std::abs(phase); 
            
            auto iNew = S + i*N;
            auto jOld = R + j*N;

            MatAdd ('N','N',N,1,_F(1.),iNew,N,-phase,jOld,N,jOld,N);
            maxDel = *std::max_element(jOld, jOld+N, 
                        [&] (_F A, _F B) { return std::norm(A) < std::norm(B); });
 
            SiConv[i] = std::abs(maxDel) < VecConv; 
            if(SiConv[i]) std::cout << "        Root "   << std::setw(5) << std::right << i
                                    << " has converged"  << std::endl;
            else std::cout << "        Root " <<  std::right << std::setw(5) << i
                           << " has not converged, maximum delta is"<<  std::right  
                           << std::setw(20) << std::abs(maxDel) << std::endl;
          }

          // exam eigenvalues
          dcomplex EDiff;
          for(i = 0; i < nExam; i++) {
            j = StMap[i];
            EDiff = Eig[i] - EPrev[j];
            SiConv[i] = SiConv[i] and (std::abs(EDiff) < EConv);  
          } 
        }   // Exam eigenvales and eigenvectors
          
        isConverged = std::all_of(SiConv.begin(), SiConv.begin()+nExam, [&] (bool i) { return i; });
        if(isConverged or nVCur >= MSS) {
        
          double DavidsonDur = tock(DavidsonSt);
          double perLT = LTdur * 100 / DavidsonDur;
        
          std::cout << "\n      - DURATION = " << std::setprecision(8) << DavidsonDur 
            << " s  ( " << perLT << " % LT )" << std::endl;
          break;
        }

        // Form residue vectors only for unconverged vectors
        std::vector<int> unConvS(nExam);
        nDo = 0ul;
        std::fill_n(unConvS.begin(),nExam,-1);
        std::fill_n(SCR,nVCur*nExam,_F(0.));
        for (auto i = 0ul; i < nExam; i++) {
          if(not SiConv[i]) {
            unConvS[nDo] = i;
            std::copy_n(XR + i*nVCur,nVCur,SCR + nDo*nVCur);
            nDo++;
          }
        }

        if(nVCur+nDo > MSS)  nDo = MSS - nVCur;
        
        // R <- AVR * XR
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nDo,nVCur,_F(1.),AVR,N,SCR,nVCur,_F(0.),R,N);

        // S as scratch space, <- eig_i * (VR * XR_i)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nDo,nVCur,_F(1.),VR,N,SCR,nVCur,_F(0.),S,N);
            
        for (auto i = 0ul; i < nDo; i++) 
          blas::scal(N,this->dcomplexTo_F(Eig[unConvS[i]]),S + i*N,1);
            
        // Compute the residue norm and generate perturbbed vectors
        MatAdd ('N','N',N,nDo,_F(1.),R,N,_F(-1.),S,N,S,N);
        std::cout << "\n      - Residues of non-converged roots: " << std::endl;  
        
        if(EigForT) std::fill_n(EigForT,nG,dcomplex(0.));
        std::fill_n(RelRes, nG, 0.);
        for (auto i = 0ul; i < nDo; i++) {
          auto j = unConvS[i];
          RelRes[i] = blas::nrm2(N,S + i*N,1); 
          std::cout << "        Root " << std::setw(5) << std::right << j+1 
                    << " 2nd order lowering " << std::right << std::setw(20) << RelRes[i]*RelRes[i] 
                    << " norm " << std::right << std::setw(20) << RelRes[i] << std::endl;
              
          if(EigForT) this->EigForT[i] = Eig[j];
        }
                
        this->preCondNoShift_(nDo,S,S); 
            
        // Append the new vectors to VR and orthogoalize against existing ones 
        // Also update the dimensions and save the XRPrev
        std::copy_n(S,nDo*N,VR+nVCur*N);
        std::copy_n(XR,nVCur*nVCur,XRPrev);
        std::copy_n(Eig,nVCur,EPrev);
        nVPrev = nVCur;
        if(this->DoLeftEigVec) {
          CErr("Do Left Eig Vec is not implemented yet");
        } else {
          // disable printing from GramSchmidt
          
#ifndef DEBUG_DAVIDSON
          std::cout.setstate(std::ios_base::failbit);
#endif

          nVCur = GramSchmidt(N,nVCur,nDo,VR,N,this->memManager_);   
          
#ifndef DEBUG_DAVIDSON
          std::cout.clear();
#endif
          
          nDo   = nVCur - nVPrev;
        }
         
        if(iter+1 == this->whenSc) { 
          nExam = nR; 
          nDo   = std::min(nDo, nExam);
          nVCur = nVPrev + nDo;
        }
        
        double DavidsonDur = tock(DavidsonSt);
        double perLT = LTdur * 100 / DavidsonDur;
    
        std::cout << "\n      - DURATION = " << std::setprecision(8) << DavidsonDur 
          << " s  ( " << perLT << " % LT )" << std::endl;
    
        if(nDo == 0) {
          isConverged = true;
          break;
        }
      } // Root Only 
    
    } // Davidson iteration    

    if( isRoot ) {

      // move data before exit runMicro      
      std::copy_n(Eig,nR,this->eigVal_);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nR,nVCur,_F(1.),VR,N,XR,nVCur,_F(0.),this->VR_,N);
      //size_t nVSave = isConverged ? nR: nG;  
      //std::copy_n(Eig,nVSave,this->eigVal_);
      //blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVSave,nVCur,_F(1.),VR,N,XR,nVCur,_F(0.),this->VR_,N);
      
      std::cout << "\n  * ";
      if( isConverged and nDo !=0)
        std::cout << "Davidson Converged in " << iter+1 << " Iterations" 
                  << std::endl;
      else if (nDo ==0) {
        double maxResDel = *std::max_element(RelRes, RelRes+nR);
        std::cout << "Davidson expansion finished, and wavefunction converged "
          << "below " << std::setw(20) << std::right << maxResDel << std::endl;
      } else
        std::cout << "Davidson Failed to Converged in " << iter+1 << " Iterations" 
                  << std::endl;
    } // Root Only 

    // Free Scratch space
      
    if(VR)       this->memManager_.template free(VR);    
    if(AVR)      this->memManager_.template free(AVR);   
    if(XR)       this->memManager_.template free(XR);    
    if(XRPrev)   this->memManager_.template free(XRPrev);
    if(SubA)     this->memManager_.template free(SubA);
    if(Eig)      this->memManager_.template free(Eig);
    if(EPrev)    this->memManager_.template free(EPrev); 
    if(Ovlp)     this->memManager_.template free(Ovlp);  
    if(R)        this->memManager_.template free(R);     
    if(S)        this->memManager_.template free(S);     
    if(SCR)      this->memManager_.template free(SCR);   
    if(VL)       this->memManager_.template free(VL);    
    if(XL)       this->memManager_.template free(XL);    
    
    return isConverged;
  
  } // Davidson::runMicro

  template <typename _F>
  void Davidson<_F>::restart() {
    // copy full vectors as new guess  
    std::cout << "\n  * Restarting Davidson..." << std::endl;
    if(Guess) this->memManager_.template free(Guess);
    
    // restart with only nRoots_ of guess
    // as now it's more close to the solution
    this->kG      = 1;
    this->whenSc  = 1;
    this->nGuess_ = this->nRoots_;
    const size_t NNS = this->N_ * this->nGuess_;
    Guess = this->memManager_.template malloc<_F>(NNS);
    std::copy_n(this->VR_, NNS, this->Guess);
  } // Davidson::restart
  
}; // namespace ChronusQ

#endif
