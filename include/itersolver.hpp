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

#include <chronusq_sys.hpp>
#include <memmanager.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/solve.hpp>
#include <cerr.hpp>

#include <util/mpi.hpp>
#include <util/math.hpp>

#include <iostream>

namespace ChronusQ {


  template <typename _F>
  class IterSolver {

  protected:

    typedef std::function< void(size_t,_F*,_F*) >    LinearTrans_t;
    typedef std::function< void(size_t,_F,_F*,_F*) > Shift_t;

    CQMemManager &memManager_;
    MPI_Comm     comm_;

    size_t N_;       ///< Problem dimension
    size_t mSS_;     ///< Max subspace dimension
    size_t maxMacroIter_; ///< Maximum number of Macro iterations
    size_t maxMicroIter_; ///< Maximum number of Micro iterations


    double convCrit_; ///< Convergence criteria

    LinearTrans_t linearTrans_;    ///< AX Product
    LinearTrans_t preCondNoShift_; ///< Unshifted preconditioner

    Shift_t shiftVec_;       ///< (A - sB)X given AX
    Shift_t preCondWShift_;  ///< Shifted preconditioner

    inline void defaultShiftVec_(size_t nVec, _F shift, _F *V, _F *AV) {

      blas::axpy(nVec*N_,shift,V,1,AV,1);

    };

    inline void defaultPreCondShift_(size_t nVec, _F shift, _F *V, _F *AV) {

      preCondNoShift_(nVec,V,AV);
      if(std::abs(std::abs(shift)) > 1e-15) shiftVec_(nVec,-1./shift,V,AV);

    };

  public:

    // Constructor (unshifted preconditioner)
    IterSolver(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans, 
      const LinearTrans_t &preNoShift, 
      const Shift_t &shiftVec = Shift_t()) : 
        comm_(c), memManager_(mem), N_(N), mSS_(MSS), maxMacroIter_(MAXMACRO),
        maxMicroIter_(MAXMICRO), convCrit_(conv), 
        linearTrans_(linearTrans), shiftVec_(shiftVec), 
        preCondNoShift_(preNoShift) { 
    
      using namespace std::placeholders;
      if( not shiftVec_ ) {
        shiftVec_ = 
          std::bind(&IterSolver<_F>::defaultShiftVec_,this,_1,_2,_3,_4); 
      }

      preCondWShift_ = 
        std::bind(&IterSolver<_F>::defaultPreCondShift_,this,_1,_2,_3,_4); 

    }

    // Constructor (shifted preconditioner)
    IterSolver(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans, 
      const Shift_t &preShift, 
      const Shift_t &shiftVec = Shift_t()) : 
        comm_(c), memManager_(mem), N_(N), mSS_(MSS), maxMacroIter_(MAXMACRO),
        maxMicroIter_(MAXMICRO), convCrit_(conv), 
        linearTrans_(linearTrans), shiftVec_(shiftVec), 
        preCondWShift_(preShift) { 

      if( not shiftVec_ ) {
        using namespace std::placeholders;
        shiftVec_ = 
          std::bind(&IterSolver<_F>::defaultShiftVec_,this,_1,_2,_3,_4); 
      }

    }



  };




  template <typename _F>
  class IterLinearSolver : public IterSolver<_F> {

  protected:

    size_t nRHS_;
    _F*    RHS_ = nullptr;
    _F*    SOL_ = nullptr;
    _F*    V_   = nullptr;
    _F*    AV_  = nullptr;
    _F*    RES_ = nullptr;

    std::vector<_F> shifts_;

    std::vector<double>              rhsNorm_;
    std::vector<std::vector<double>> resNorm_;

  public:


    using LinearTrans_t = typename IterSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterSolver<_F>::Shift_t;

    size_t shiftBS = 1;
    size_t rhsBS   = 1;

    IterLinearSolver(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS, 
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv, 
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift, 
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,preNoShift,
                    shiftVec){ } 

    IterLinearSolver(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS, 
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv, 
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift, 
      const Shift_t &shiftVec = Shift_t()) : 
     IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,preShift,
                    shiftVec){ } 


    ~IterLinearSolver() {

      if(RHS_) this->memManager_.free(RHS_);
      if(SOL_) this->memManager_.free(SOL_);
      if(V_)   this->memManager_.free(V_);
      if(AV_)  this->memManager_.free(AV_);
      if(RES_) this->memManager_.free(RES_);

    };




    template <typename T>
    void setRHS(size_t nRHS, T*RHS, size_t LDRHS);

    template <typename T>
    void setShifts(size_t nShift, T* shifts); 


    inline void getSol(_F *SOL) {

      std::copy_n(SOL_, nRHS_ * shifts_.size() * this->N_, SOL);

    };




    virtual void alloc() {

      // NO MPI
      ROOT_ONLY(this->comm_);

      size_t nRHSnSN       = nRHS_ * shifts_.size() * this->N_;
      size_t nRHSnSN_batch = rhsBS * shiftBS * this->N_;

      SOL_ = this->memManager_.template malloc<_F>(nRHSnSN);
      V_   = this->memManager_.template malloc<_F>(nRHSnSN_batch*this->mSS_);
      AV_  = this->memManager_.template malloc<_F>(nRHSnSN_batch*this->mSS_);
      RES_  = this->memManager_.template malloc<_F>(nRHSnSN_batch);

    }



    void run();

    virtual void runBatch(size_t nRHS, size_t nShift, _F* RHS, _F *shifts, 
      _F* SOL, double *RHSNorm);

  };




  template <typename _F>
  class IterDiagonalizer : public IterSolver<_F> {


  protected:

    size_t nRoots_;
    size_t nGuess_;
    _F     *VR_   = nullptr;
    _F     *VL_   = nullptr;
    _F     *AVR_  = nullptr;
    _F     *AVL_  = nullptr;
    _F     *RESR_ = nullptr;
    _F     *RESL_ = nullptr;


    dcomplex * eigVal_ = nullptr;

    void RayleighRitz();



  public:

    using LinearTrans_t = typename IterSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterSolver<_F>::Shift_t;

    IterDiagonalizer(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS, 
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv, 
      size_t nR, 
      size_t nG, 
      const LinearTrans_t &linearTrans, 
      const LinearTrans_t &preNoShift, 
      const Shift_t &shiftVec = Shift_t()) :
      IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,
                     preNoShift,shiftVec), 
        nRoots_(nR), nGuess_(nG){ } 

    IterDiagonalizer(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS, 
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv, 
      size_t nR, 
      size_t nG, 
      const LinearTrans_t &linearTrans, 
      const Shift_t &preShift, 
      const Shift_t &shiftVec = Shift_t()) : 
     IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,
                    preShift,shiftVec),
       nRoots_(nR), nGuess_(nG){ } 


    ~IterDiagonalizer() {

      if(VR_)   this->memManager_.free(VR_);
      if(VL_)   this->memManager_.free(VL_);
      if(AVR_)  this->memManager_.free(AVR_);
      if(AVL_)  this->memManager_.free(AVL_);
      if(RESR_) this->memManager_.free(RESR_);
      if(RESL_) this->memManager_.free(RESL_);

    };




    virtual void alloc() {

      // NO MPI
      ROOT_ONLY(this->comm_);

      size_t NNR = this->N_ * this->mSS_;

      eigVal_ = this->memManager_.template malloc<dcomplex>(this->nRoots_);
      VR_     = this->memManager_.template malloc<_F>(NNR);
      //AVR_    = this->memManager_.template malloc<_F>(NNR);
      //RESR_   = this->memManager_.template malloc<_F>(NNR);

    }


    dcomplex* eigVal() const { return eigVal_; }
    _F*       VR()     const { return VR_;     }


    virtual void setGuess(size_t nGuess, 
        std::function<void(size_t,_F*,size_t)>) = 0;

    void run();

    virtual bool runMicro() = 0;
    virtual void restart() = 0;

  };






  

  // GPLHR
  template <typename _F>
  class GPLHR : public IterDiagonalizer<_F> {

    double *RelRes = nullptr;
    _F *Guess = nullptr;


  protected:


    bool useAdaptiveSigma() const { return false; }


    void getResidualNorms(size_t N, size_t nR, _F *WMAT, double *RelRes, 
        dcomplex *LAMBDA, double nrmA);

    void getTriU(size_t N, _F *TRIUA, size_t LDTRIUA, _F *TRIUB, 
      size_t LDTRIUB, _F *MA, size_t LDMA, _F *MB, size_t LDMB);
  
    bool checkConv(size_t nR, double *RelRes);




    void halfProj(size_t N, size_t nV, size_t nS, _F *V, size_t LDV, _F *S, 
        size_t LDS, _F *smSCR, size_t LDSCR);

    /**
     *  \brief Overload of halfProj where nV == nS
     */
    inline void halfProj(size_t N, size_t nV, _F *V, size_t LDV, _F *S, 
        size_t LDS, _F *smSCR, size_t LDSCR) {

      halfProj(N,nV,nV,V,LDV,S,LDS,smSCR,LDSCR);

    }

    void halfProj2(size_t N, size_t nV, size_t nS, _F *V, size_t LDV, _F *AV, 
        size_t LDAV, _F *S, size_t LDS, _F *AS, size_t LDAS, _F *smSCR, 
        size_t LDSCR);

    /**
     *  \brief Overload of halfProj2 where nV == nS
     */
    inline void halfProj2(size_t N, size_t nV, _F *V, size_t LDV, _F *AV, 
        size_t LDAV, _F *S, size_t LDS, _F *AS, size_t LDAS, _F *smSCR, 
        size_t LDSCR) {

      halfProj2(N,nV,nV,V,LDV,AV,LDAV,S,LDS,AS,LDAS,smSCR,LDSCR);

    }


    void newSMatrix(size_t N, size_t nR, _F *V, size_t LDV, _F *Q, size_t LDQ,
        _F *S, size_t LDS, _F *smSCR, size_t LDSCR);


  public:


    size_t m = 1;
    _F sigma = 0.;
    _F hardLim = -std::numeric_limits<double>::infinity();

    using LinearTrans_t = typename IterDiagonalizer<_F>::LinearTrans_t;
    using Shift_t       = typename IterDiagonalizer<_F>::Shift_t;

    GPLHR(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N,
      const size_t MAXITER,
      double conv, 
      size_t nR, 
      const LinearTrans_t &linearTrans, 
      const LinearTrans_t &preNoShift, 
      const Shift_t &shiftVec = Shift_t()) :
      IterDiagonalizer<_F>(c,mem,N,2 + m*nR,1,MAXITER,conv,nR,nR,
          linearTrans,preNoShift,shiftVec){ } 

    GPLHR(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N,
      const size_t MAXITER,
      double conv, 
      size_t nR, 
      const LinearTrans_t &linearTrans, 
      const Shift_t &preShift, 
      const Shift_t &shiftVec = Shift_t()) :
      IterDiagonalizer<_F>(c,mem,N,2 + m*nR,1,MAXITER,conv,nR,nR,
          linearTrans,preShift,shiftVec){ } 

    ~GPLHR() { 
    
      if( RelRes ) this->memManager_.free(RelRes);
      if( Guess  ) this->memManager_.free(Guess);

    }

    void setM(size_t _m) {

      m = _m;
      this->mSS_ = (3 + m)*this->nRoots_;

    }

    void alloc() {

      IterDiagonalizer<_F>::alloc();

      // NO MPI
      ROOT_ONLY(this->comm_);


      // Allocate GPLHR specific Memory
        
      this->RelRes = this->memManager_.template malloc<double>(this->nRoots_);

    }

    bool runMicro();

    void restart();

    void setGuess(size_t nGuess, std::function<void(size_t,_F*,size_t)> func) {

      if( nGuess != this->nRoots_ ) 
        CErr("GPLHR Requires nGuess = nRoots",std::cout);

      // NO MPI
      ROOT_ONLY(this->comm_);

      Guess = this->memManager_.template malloc<_F>(nGuess * this->N_);

      func(nGuess, Guess, this->N_);

    }

  }; // GPLHR





  // Davidson
  template <typename _F>
  class Davidson : public IterDiagonalizer<_F> {

    bool DoLeftEigVec   = false;
   
    // Energy specific options 
    bool EnergySpecific = false;
    size_t nHighERoots  = 0; 
    size_t  nLowERoots  = 0; 
    double EnergyRef    = 0.; 
    bool adaptiveERef   = false; // use loweset eigenvalues at current iteration

    double   *RelRes  = nullptr;
    _F       *Guess   = nullptr;
    dcomplex *EigForT = nullptr; // Eigenvalues can be used for preconditioner 

    // creating double(dcomplex) operator to accomodate the
    // complex eigenvalue to real arithmetric 
    _F dcomplexTo_F(dcomplex &a) { return * reinterpret_cast<_F*>(&a);}

  protected:
    
  public:

    size_t m = 50;
    size_t whenSc = 2;
    size_t kG = 3;

    using LinearTrans_t = typename IterDiagonalizer<_F>::LinearTrans_t;
    using Shift_t       = typename IterDiagonalizer<_F>::Shift_t;

    Davidson(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N,
      const size_t MAXMACROITER,
      const size_t MAXMICROITER,
      double conv, 
      size_t nR, 
      const LinearTrans_t &linearTrans, 
      const LinearTrans_t &preNoShift, 
      const Shift_t &shiftVec = Shift_t()):
      IterDiagonalizer<_F>(c,mem,N,m*nR,MAXMACROITER,MAXMICROITER,conv,nR,nR*kG,
          linearTrans,preNoShift,shiftVec){ } 
    
    ~Davidson() { 
    
      if( RelRes  ) this->memManager_.free(RelRes);
      if( Guess   ) this->memManager_.free(Guess);

    }

    void setM(size_t _m) {
      m = _m;
      this->mSS_ = m*this->nRoots_;
    }
    
    void setkG(size_t _kG) {
      kG = _kG;
      this->nGuess_ = kG*this->nRoots_;
    }

    void setWhenSc(size_t _WhenSc) { whenSc = _WhenSc;}
    
    void setEigForT(dcomplex * _Eig) {
      
      if( this->memManager_.getSize(_Eig) < this->nGuess_ ) 
        CErr("Davison EigForT requires a memory block with size at least nGuess ",std::cout);
        
      EigForT = _Eig;
    } 

    void alloc() {

      IterDiagonalizer<_F>::alloc();

      // NO MPI
      ROOT_ONLY(this->comm_);

      // Allocate Davidson specific Memory
      this->RelRes  = this->memManager_.template malloc<double>(this->nGuess_);
    }

    bool runMicro();

    void restart();
    
    void useEnergySpecific(size_t nHER, double energyThd, double ERef) {
      assert (this->nRoots_ >= nHER);
      assert (energyThd > 0);
      
      this->EnergySpecific = true;
      this->nHighERoots    = nHER;
      this->nLowERoots     = this->nRoots_ - nHER;
      this->EnergyRef      = ERef;
    }

    void useEnergySpecific(size_t nHER, double energyThd) {
      this->adaptiveERef = true;
      this->useEnergySpecific(nHER, energyThd, 0.);
    }

    void doLeftEigenvector() {this->DoLeftEigVec = true; }
    
    void setGuess(size_t nGuess, std::function<void(size_t,_F*,size_t)> func) {

      if( nGuess != this->nGuess_ ) 
        CErr("Davison Requires nGuess = nGuess_",std::cout);

      // NO MPI
      ROOT_ONLY(this->comm_);

      Guess = this->memManager_.template malloc<_F>(nGuess * this->N_);

      func(nGuess, Guess, this->N_);

    }


  }; // Davidson 


  template <typename _F>
  class GMRES : public IterLinearSolver<_F> {


    _F * HHR_ = nullptr;
    _F * U_   = nullptr;
    _F * W_   = nullptr;
    _F * J_   = nullptr;
    _F * R_   = nullptr;

  public:

    using LinearTrans_t = typename IterLinearSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterLinearSolver<_F>::Shift_t;


    GMRES(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS,
      double conv, 
      const LinearTrans_t &linearTrans, 
      const LinearTrans_t &preNoShift, 
      const Shift_t &shiftVec = Shift_t()) :
     IterLinearSolver<_F>(c,mem,N,MSS,1,MSS,conv,linearTrans,preNoShift,
         shiftVec){ } 

    GMRES(
      MPI_Comm c, 
      CQMemManager &mem, 
      const size_t N, 
      const size_t MSS, 
      double conv, 
      const LinearTrans_t &linearTrans, 
      const Shift_t &preShift, 
      const Shift_t &shiftVec = Shift_t()) : 
     IterLinearSolver<_F>(c,mem,N,MSS,1,MSS,conv,linearTrans,preShift,
         shiftVec){ } 


    ~GMRES() {

      if(HHR_)  this->memManager_.free(HHR_);
      if(U_  )  this->memManager_.free(U_);
      if(W_  )  this->memManager_.free(W_);
      if(J_  )  this->memManager_.free(J_);
      if(R_  )  this->memManager_.free(R_);

    }


    void alloc() {

      // Standard linear solver allocations
      IterLinearSolver<_F>::alloc();

      // No MPI for GMRES
      ROOT_ONLY(this->comm_);
      size_t nBatch = this->rhsBS * this->shiftBS;
      size_t MSSnBatch = this->mSS_ * nBatch;

      HHR_ = this->memManager_.template malloc<_F>(nBatch * this->N_);

      J_   = this->memManager_.template malloc<_F>(2 * MSSnBatch);
      U_   = this->memManager_.template malloc<_F>(MSSnBatch * this->N_);
      R_   = this->memManager_.template malloc<_F>(MSSnBatch * this->mSS_);
      W_   = this->memManager_.template malloc<_F>(MSSnBatch + nBatch);

    }

    void runBatch(size_t nRHS, size_t nShift, _F* RHS, _F *shifts, 
      _F* SOL, double *RHSNorm );

  };

  

};

