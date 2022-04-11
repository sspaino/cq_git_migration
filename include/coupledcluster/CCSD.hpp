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

#include <singleslater.hpp>
#include <particleintegrals.hpp>
#include <wavefunction.hpp>
#include <chronusq_sys.hpp>
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <util/files.hpp>
#include <util/math.hpp>
#include <util/threads.hpp>
#include <integrals.hpp>
#include <cerr.hpp>
#include <coupledcluster.hpp>
#include <matrix.hpp>
#include <cqlinalg.hpp>
#include <typeinfo> 
#include <stdlib.h>
#include "DiskDIIS.hpp"


namespace ChronusQ{


  // Safe always
  void startTAThreads() {
    madness::ThreadPool::begin(GetLAThreads() - 1);
    SetLAThreads(1);
  };

  // Will throw if ThreadPool not initialized
  void endTAThreads() {
    SetLAThreads(madness::ThreadPool::size() + 1);
    if ( madness::ThreadPool::size() != 0 )
      madness::ThreadPool::end();
  };

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::run() {
    if(std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(this->ref_.aoints.TPI)){
      runConventional();
    }
    else {
      CErr("Only GHF/X2C-CCSD is supported!");
    }
  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initIntermediates() {
      size_t NO = this->ref_.nO;
      size_t NV = this->ref_.nV;  
    if (not this->A_1.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,orange});
    tmp.fill(0.0);
    this->A_1("p,q") = tmp("p,q");
    }
    if (not this->A_2.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{vrange,vrange});
    tmp.fill(0.0);
    this->A_2("p,q") = tmp("p,q");
    }
    if (not this->A_3.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,vrange});
    tmp.fill(0.0);
    this->A_3("p,q") = tmp("p,q");
    }
    if (not this->A_4.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,orange,vrange,orange});
    tmp.fill(0.0);
    this->A_4("p,q,r,s") = tmp("p,q,r,s");
    }
  
    if (not this->B_1.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,orange});
    tmp.fill(0.0);
    this->B_1("p,q") = tmp("p,q");
    }
    if (not this->B_2.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{vrange,vrange});
    tmp.fill(0.0);
    this->B_2("p,q") = tmp("p,q");
    }
    if (not this->B_3.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,vrange,orange,orange});
    tmp.fill(0.0);
    this->B_3("p,q,r,s") = tmp("p,q,r,s");
    }
    if (not this->B_5.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,orange,orange,orange});
    tmp.fill(0.0);
    this->B_5("p,q,r,s") = tmp("p,q,r,s");
    }
    if (not this->B_6.is_initialized()){
    TArray tmp(TA::get_default_world(),TA::TiledRange{orange,vrange,vrange,orange});
    tmp.fill(0.0);
    this->B_6("p,q,r,s") = tmp("p,q,r,s");
    }
  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initAmplitudes() {
      size_t NO = this->ref_.nO;
      size_t NV = this->ref_.nV;  
      if (not this->T1_.is_initialized()){
        TArray tmp(TA::get_default_world(),TA::TiledRange{vrange,orange});
        tmp.fill(0.0);
        this->T1_("p,q") = tmp("p,q");
      }
      if (not this->T2_.is_initialized()){
        TArray tmp(TA::get_default_world(),TA::TiledRange{vrange,vrange,orange,orange});
        tmp.fill(0.0);
        this->T2_("p,q,r,s") = tmp("p,q,r,s");
      }    
  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initFock() {
      size_t nO = this->ref_.nO;
      size_t nV = this->ref_.nV;  
      std::vector<std::string> fockTypes{"oo","vo","vv","ov"};
      for(const auto& focktype:fockTypes){
  
        std::vector<int> trange;
        std::vector<int> offset;
  
        for(const auto& ftype:focktype){
          if(ftype == 'o'){
            trange.push_back(nO);
            offset.push_back(0);
          }
          else{
            trange.push_back(nV);
            offset.push_back(nO);
          }
        }
  
        TArray tmp(TA::get_default_world(),TA::TiledRange({{0,trange[0]},{0,trange[1]}}));
  
        for(auto it = std::begin(tmp); it != std::end(tmp); ++it) {
          // Construct a tile
        typename TA::Array<MatsT,2>::value_type tile(
        tmp.trange().make_tile_range(it.ordinal()));
        const auto& lobound = tile.range().lobound();
        const auto& upbound = tile.range().upbound();
        
        std::size_t i[] = {0,0};  
        for(i[0] = lobound[0]; i[0] != upbound[0]; ++i[0])
          for(i[1] = lobound[1]; i[1] != upbound[1]; ++i[1])
            tile[i] = 0.0 ;
  //            tile[i] = ref_.fockMO[0](i[0] + offset[0], i[1] + offset[1]);              
        *it = tile;
        }
  
        this->fockMatrix_ta.insert(std::make_pair(focktype,tmp));
      }
  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::transformInts(){
    std::size_t nMO;
    std::size_t nO;
    std::size_t nV;      
    if (this->ref_.nC == 2){      
      nMO = this->ref_.nOrbital();
      nO = this->ref_.nO;
      nV = this->ref_.nV;
    }    
    else {
      CErr("Only GHF/X2C-CCSD is supported!");
    }  

    std::vector<std::string> fockTypes{"oo","vo","vv","ov"};
    for(const auto& focktype:fockTypes){

      std::vector<TA::TiledRange1> fockRangeTypes;

      for(const auto& ftype:focktype){
        if(ftype == 'o'){
          fockRangeTypes.push_back(orange);
        }
        else{
          fockRangeTypes.push_back(vrange);
        }
      }

      TArray tmp(TA::get_default_world(),TA::TiledRange{fockRangeTypes[0],fockRangeTypes[1]});

      for(auto it = std::begin(tmp); it != std::end(tmp); ++it) {
        // Construct a tile
      typename TA::Array<MatsT,2>::value_type tile(
      tmp.trange().make_tile_range(it.ordinal()));
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();
      
      std::size_t i[] = {0,0};  
      for(i[0] = lobound[0]; i[0] != upbound[0]; ++i[0])
        for(i[1] = lobound[1]; i[1] != upbound[1]; ++i[1])
          tile[i] = 0.0 ;
//            tile[i] = ref_.fockMO[0](i[0] + offset[0], i[1] + offset[1]);              
      *it = tile;
      }

      this->fockMatrix_ta.insert(std::make_pair(focktype,tmp));
    }

    // Start parallel without TA
    endTAThreads();

    auto aoeri = std::make_shared<InCore4indexTPI<IntsT>>(
              std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(this->ref_.aoints.TPI)
               ->template spatialToSpinBlock<IntsT>());
    auto moeri = std::make_shared<InCore4indexTPI<MatsT>>(
        aoeri->transform('N', this->ref_.mo[0].pointer(), nMO, nMO));
    Integrals<MatsT> moints;
    moints.TPI = moeri;    

    // Go back to TA threads
    startTAThreads();

    if (std::shared_ptr<InCore4indexTPI<dcomplex>> eri = 
      std::dynamic_pointer_cast<InCore4indexTPI<dcomplex>>(moints.TPI)){
    std::vector<std::string> mointsTypes{"abij","iabj","aibj","ijkl","abcd","abci","aijk"};
    std::set<char> hole{'i','j','k','l'};
    std::set<char> particle{'a','b','c','d'};  
    for(const auto& motype:mointsTypes){

      std::vector<TA::TiledRange1> MORangeTypes;
      std::vector<unsigned> offset;
      
      for (const auto& ch:motype){
       if (hole.count(ch)){
           MORangeTypes.push_back(orange);
           offset.push_back(0);
       }
       else if(particle.count(ch)){
           MORangeTypes.push_back(vrange);
           offset.push_back(nO);
       }
      }

               
      TArray tmp(TA::get_default_world(),TA::TiledRange{MORangeTypes[0],MORangeTypes[1],MORangeTypes[2],MORangeTypes[3]});

      for(auto it = std::begin(tmp); it != std::end(tmp); ++it) {
        // Construct a tile
        typename TA::Array<MatsT,4>::value_type tile(tmp.trange().make_tile_range(it.ordinal()));
      
        const auto& lobound = tile.range().lobound();
        const auto& upbound = tile.range().upbound();
      
        std::size_t i[] = {0,0,0,0};  
        for(i[0] = lobound[0]; i[0] != upbound[0]; ++i[0])
          for(i[1] = lobound[1]; i[1] != upbound[1]; ++i[1])
            for(i[2] = lobound[2]; i[2] != upbound[2]; ++i[2])
              for(i[3] = lobound[3]; i[3] != upbound[3]; ++i[3])
                tile[i] = (*eri)(i[0] + offset[0],i[2] + offset[2],i[1] + offset[1],i[3] + offset[3]) 
                - (*eri)(i[0] + offset[0], i[3] + offset[3],i[1] + offset[1], i[2] + offset[2]);
      
        *it = tile;
      }

      this->antiSymMoints.insert(std::make_pair(motype,tmp));

    }
    } else { CErr("Only GHF/X2C-CCSD is supported!");}

  }  

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initRanges(){
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV; 

    size_t NOblocks = NO % blksize == 0 ? NO / blksize : NO / blksize + 1;
    size_t NVblocks = NV % blksize == 0 ? NV / blksize : NV / blksize + 1; 

    std::vector<std::size_t> vblk, oblk; 

    for (auto i = 0 ; i < NOblocks; i++)
      oblk.push_back(blksize * i);
    oblk.push_back(NO);    
    for (auto i = 0 ; i < NVblocks; i++)
      vblk.push_back(blksize * i);
    vblk.push_back(NV); 

    TA::TiledRange1 orangetmp(oblk.begin(), oblk.end());
    TA::TiledRange1 vrangetmp(vblk.begin(), vblk.end());

    orange = orangetmp;
    vrange = vrangetmp;
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::doDIIS(size_t t1offset, size_t NT1blks, size_t NT2blks, TArray T1_old, TArray T2_old, std::shared_ptr<DiskDIIS<MatsT>> diis ){
      size_t t2offset = 0;
      for(size_t i = 0; i < NT2blks; i++){
        diis->setVector(this->T2_.find(i).get().data(),this->T2_.find(i).get().size(),t2offset);
        t2offset += this->T2_.find(i).get().size();
      }      
      
      size_t t1_offset = t1offset;
      for(size_t i = 0; i < NT1blks; i++){
        diis->setVector(this->T1_.find(i).get().data(),this->T1_.find(i).get().size(),t1_offset);
        t1_offset += this->T1_.find(i).get().size();
      }

      for(size_t i = 0; i < NT2blks; i++){
        blas::axpy(this->T2_.find(i).get().size(),-1.0,this->T2_.find(i).get().data(),1,T2_old.find(i).get().data(),1);
        blas::scal(this->T2_.find(i).get().size(),-1.0,T2_old.find(i).get().data(),1); 
      }      
      

      for(size_t i = 0; i < NT1blks; i++){
        blas::axpy(this->T1_.find(i).get().size(),-1.0,this->T1_.find(i).get().data(),1,T1_old.find(i).get().data(),1);
        blas::scal(this->T1_.find(i).get().size(),-1.0,T1_old.find(i).get().data(),1);
      }

      t2offset = 0;
      for(size_t i = 0; i < NT2blks; i++){
        diis->setErrorVector(T2_old.find(i).get().data(),T2_old.find(i).get().size(),t2offset);
        t2offset += T2_old.find(i).get().size();
      }      
      
      t1_offset = t1offset;
      for(size_t i = 0; i < NT1blks; i++){
        diis->setErrorVector(T1_old.find(i).get().data(),T1_old.find(i).get().size(),t1_offset);
        t1_offset += T1_old.find(i).get().size();
      }

      diis->extrapolate();

      t2offset = 0;
      for(size_t i = 0; i < NT2blks; i++){
        diis->getVector(this->T2_.find(i).get().data(),this->T2_.find(i).get().size(),t2offset);
        t2offset += this->T2_.find(i).get().size();
      }      
    
      t1_offset = t1offset;
      for(size_t i = 0; i < NT1blks; i++){
        diis->getVector(this->T1_.find(i).get().data(),this->T1_.find(i).get().size(),t1_offset);
        t1_offset += this->T1_.find(i).get().size();
      }    

  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::runConventional(){

    // Has to be *before* TA::initialize
    startTAThreads();

    int argc = 1;
    char **argv = NULL;
    TA::World &world = TA::initialize(argc,argv);

    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV; 
    std::cout << BannerTop << '\n';
    std::cout << "Coupled Cluster (CC) Settings:\n\n";

    std::cout << std::setw(45) << std::left << "  Type:"
              << "Singles and Doubles" << '\n';

    std::cout << std::setw(45) << std::left << "  Occupied Orbitals:"
              << this->ref_.nO << '\n';
    std::cout << std::setw(45) << std::left << "  Virtual Orbitals:"
              << this->ref_.nV << '\n';

    std::cout << std::setprecision(6) << std::scientific;

    std::cout << std::setw(45) << std::left << "  Energy Convergence Tolerance:"
              << this->ccSettings.eConv << '\n';
    std::cout << std::setw(45) << std::left << "  Amplitude Convergence Tolerance:"
              << this->ccSettings.tConv << '\n';

    std::cout << std::setw(45) << std::left << "  Direct Inversion of Iterative Subspace:";

    if ( this->useDIIS ) {
      std::cout << "On\n";
      std::cout << std::left << "    * DIIS will track up to " << this->nDIIS
                << " previous iterations\n";
    }
    else
      std::cout << "Off\n";

    std::cout << std::endl;   

    //Calculate the estimated usage of memory
    // amplitudes: OV + O^2V^2,OV + 2O^2V^2 if diis is used
    // A intermediates: O^2 + V^2 + OV + O^3V
    // B intermediates: O^2 + V^2 + O^3V + O^4 + O^2V^2
    double mem = 3*NO*NV + 2*NO*NO + 2*NV*NV + 2*NO*NO*NV*NV + 2*NO*NO*NO*NV + NO*NO*NO*NO;
    if (this->useDIIS){
      mem += NO*NO*NV*NV;
    }
    mem = mem * sizeof(MatsT)/1e9;
    std::cout << "  " << std::fixed << mem << " GB"
              << " will be used by the CC calculation" << std::endl;

    std::cout << BannerMid << '\n' << std::endl;


    int diis_dim = NO*NO*NV*NV + NO*NV;;
    std::shared_ptr<DiskDIIS<MatsT> > diis (new DiskDIIS<MatsT> (diis_dim, this->nDIIS, this->ref_.savFile, this->ref_.memManager));
    
    size_t NOblocks = NO % blksize == 0 ? NO / blksize : NO / blksize + 1;
    size_t NVblocks = NV % blksize == 0 ? NV / blksize : NV / blksize + 1; 

    size_t NT1blks = NOblocks * NVblocks;
    size_t NT2blks = NT1blks * NT1blks;
    size_t t1offset = NO * NO * NV * NV;

    initRanges();
    transformInts();
    initIntermediates();
    initAmplitudes();

    TArray T1_old(TA::get_default_world(),TA::TiledRange{vrange,orange});
    TArray T2_old(TA::get_default_world(),TA::TiledRange{vrange,vrange,orange,orange});
    TArray T2_old_copy(TA::get_default_world(),TA::TiledRange{vrange,vrange,orange,orange});
    double tnorm = 0.0;
    double tnorm_old = 0.0;

    std::cout << std::setw(18) << std::left <<  "  CC Iterations";
    std::cout << std::setw(22) << std::right << "Corr. Energy (Eh)";
    std::cout << std::setw(19) << std::right << "\u0394Ec (Eh)";
    std::cout << std::setw(19) << std::right << "|\u0394T|";
    std::cout << '\n';
    std::cout << std::setw(18) << std::left <<  "  -------------";
    std::cout << std::setw(22) << std::right << "-----------------";
    std::cout << std::setw(18) << std::right << "--------";
    std::cout << std::setw(18) << std::right << "----";
    std::cout << '\n' << std::endl;



    for (auto iter = 0; iter < this->ccSettings.maxiter; iter++){

      T1_old("p,q") = this->T1_("p,q");
      if(this->useDIIS){
        T2_old("p,q,r,s") = this->T2_("p,q,r,s");
      }
      T2_old_copy("p,q,r,s") = this->T2_("p,q,r,s");

      formA_1();
      formA_2();
      formA_3();
      formA_4();

      formB_1();    
      formB_2();
      formB_3();
      formB_5();
      formB_6();     
            
      updateT1(T1_old,T2_old_copy);
      //T2_old_copy will be modified in updateT2()!
      updateT2(T1_old,T2_old_copy);
      
      if (this->useDIIS){
      // diis
        doDIIS(t1offset, NT1blks, NT2blks, T1_old, T2_old, diis);
      }

      double Eold = this->CorrE;

      getCorrEnergy();     
      tnorm = sqrt(TA::norm2(this->T1_) + TA::norm2(this->T2_));
      double dE = this->CorrE - Eold;
      double dT = tnorm - tnorm_old;

      std::cout << std::setprecision(12) << std::fixed;
      std::cout << "  Iteration "  << std::setw(6) << std::left << iter;
      std::cout << std::setw(22) << std::right << std::fixed << this->CorrE;
      std::cout << std::setw(18) << std::right << std::fixed << dE;
      std::cout << std::setw(18) << std::right << std::fixed << std::abs(dT);
      std::cout << std::endl;

      if (std::abs(dE) < this->ccSettings.eConv and std::abs(dT) < this->ccSettings.tConv) { 

        std::cout << "\n  CC Completed: CorrE is " << std::setprecision(12) << this->CorrE << "Eh" << std::endl;
        
        if (this->ref_.savFile.exists()) {
        this->ref_.savFile.safeWriteData("/CC/CORRELATION_ENERGY",&this->CorrE, {1});
        this->ref_.savFile.safeWriteData("/CC/T_NORM",&tnorm, {1});
        }
        
        std::cout << BannerMid << std::endl;
        
        break;
        }
      else { 
        Eold = this->CorrE;
        tnorm_old = tnorm;
        }

      if(iter == ccSettings.maxiter-1){
        CErr(std::string("CC iterations didn't converge in ") + std::to_string(ccSettings.maxiter) + " steps." );
      }
    }

    endTAThreads();
    TA::finalize();

  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::getCorrEnergy() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    MatsT CorrEnergy = dot(conj(this->antiSymMoints["abij"]("c,d,k,l")),this->T2_("c,d,k,l"));
    MatsT CorrE1 = conj(this->antiSymMoints["abij"]("c,d,k,l")).dot(this->T1_("c,k") * this->T1_("d,l"));
    CorrE1 = CorrE1 * 2.0;
    this->CorrE = 0.25 * std::real(CorrEnergy + CorrE1);
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::updateT1(const TArray T1_old, const TArray T2_old) {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
  
    this->T1_("a,i") =   fockMatrix_ta["vo"]("a,i");
    this->T1_("a,i") +=    this->antiSymMoints["iabj"]("k,a,c,i") * T1_old("c,k");  
    this->T1_("a,i") +=  0.5  *   conj(this->antiSymMoints["abci"]("d,c,a,k")) * T2_old("c,d,k,i");  
    this->T1_("a,i") += A_1("k,i") * T1_old("a,k");   
    this->T1_("a,i") += A_2("a,c") * T1_old("c,i");
    this->T1_("a,i") += A_3("k,c") * T2_old("c,a,k,i");
    this->T1_("a,i") += A_4("k,l,c,i") * T2_old("c,a,k,l");
  
    TArray oneoverD_ai(TA::get_default_world(),TA::TiledRange{vrange,orange});
    for(auto it = std::begin(oneoverD_ai); it != std::end(oneoverD_ai); ++it) {
      // Construct a tile
    typename TA::Array<MatsT,2>::value_type tile(
    oneoverD_ai.trange().make_tile_range(it.ordinal()));
    const auto& lobound = tile.range().lobound();
    const auto& upbound = tile.range().upbound();
    std::size_t i[] = {0,0};  
    for(i[0] = lobound[0]; i[0] != upbound[0]; ++i[0])
      for(i[1] = lobound[1]; i[1] != upbound[1]; ++i[1])
        tile[i] = 1.0 / (this->ref_.eps1[i[1]] - this->ref_.eps1[i[0]+NO]);
    *it = tile;
    }
  
    this->T1_("a,i") = this->T1_("a,i") * oneoverD_ai("a,i");
  
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_1() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
  
    A_1("k,i") =   conj(this->antiSymMoints["aijk"]("c,i,k,l")) * this->T1_("c,l"); 
    A_1("k,i") +=  - 0.5  *   conj(this->antiSymMoints["abij"]("c,d,l,k")) * this->T2_("c,d,l,i"); 
    A_1("k,i") += -fockMatrix_ta["oo"]("i,k");
    TA::TArray<MatsT> A_11(TA::get_default_world(),TA::TiledRange{orange,vrange});
    A_11.fill(0.0);
    A_11("k,c") = - conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T1_("d,l");
    A_11("k,c") += -fockMatrix_ta["ov"]("k,c");
    A_1("k,i") += A_11("k,c") * this->T1_("c,i"); 
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_2() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    A_2("a,c") = - conj(this->antiSymMoints["abci"]("d,c,a,k")) * this->T1_("d,k");
    A_2("a,c") +=   fockMatrix_ta["vv"]("a,c");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_3() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    A_3("k,c") =    conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T1_("d,l");
    A_3("k,c") +=   fockMatrix_ta["ov"]("k,c");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_4() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    A_4("k,l,c,i") = - 0.5  *   conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T1_("d,i");
    A_4("k,l,c,i") += - 0.5  *   conj(this->antiSymMoints["aijk"]("c,i,k,l"));
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::updateT2(const TArray T1_old, TArray T2_old) {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV; 
    this->T2_("a,b,i,j") =    this->antiSymMoints["abij"]("a,b,i,j");
    this->T2_("a,b,i,j") +=   B_1("k,i") * T2_old("a,b,k,j");
    this->T2_("a,b,i,j") += -B_1("k,j") * T2_old("a,b,k,i");
    this->T2_("a,b,i,j") +=   B_2("b,c") * T2_old("c,a,i,j");
    this->T2_("a,b,i,j") += -B_2("a,c") * T2_old("c,b,i,j");
    this->T2_("a,b,i,j") +=   B_3("k,a,i,j") * T1_old("b,k");
    this->T2_("a,b,i,j") += -B_3("k,b,i,j") * T1_old("a,k");
    this->T2_("a,b,i,j") += B_5("k,l,i,j") * T2_old("a,b,k,l");
    this->T2_("a,b,i,j") +=   B_6("k,a,c,j") * T2_old("c,b,k,i");
    this->T2_("a,b,i,j") += -B_6("k,b,c,j") * T2_old("c,a,k,i");
    this->T2_("a,b,i,j") +=   B_6("k,b,c,i") * T2_old("c,a,k,j");
    this->T2_("a,b,i,j") += -B_6("k,a,c,i") * T2_old("c,b,k,j");
    this->T2_("a,b,i,j") +=  this->antiSymMoints["abci"]("a,b,c,j") * T1_old("c,i");
    this->T2_("a,b,i,j") += -this->antiSymMoints["abci"]("a,b,c,i") * T1_old("c,j");
    //reuse the T2_old container
    T2_old("c,d,i,j") += T1_old("c,i") * T1_old("d,j");
    T2_old("c,d,i,j") += - T1_old("c,j") * T1_old("d,i");
    this->T2_("a,b,i,j") += 0.5 * this->antiSymMoints["abcd"]("a,b,c,d") * T2_old("c,d,i,j");
  
    TArray oneoverD_abij(TA::get_default_world(),TA::TiledRange{vrange,vrange,orange,orange});
    for(auto it = std::begin(oneoverD_abij); it != std::end(oneoverD_abij); ++it) {
      // Construct a tile
    typename TA::Array<MatsT,4>::value_type tile(
    oneoverD_abij.trange().make_tile_range(it.ordinal()));
    const auto& lobound = tile.range().lobound();
    const auto& upbound = tile.range().upbound();
    std::size_t i[] = {0,0,0,0};  
    for(i[0] = lobound[0]; i[0] != upbound[0]; ++i[0])
      for(i[1] = lobound[1]; i[1] != upbound[1]; ++i[1])
        for(i[2] = lobound[2]; i[2] != upbound[2]; ++i[2])
          for(i[3] = lobound[3]; i[3] != upbound[3]; ++i[3])
            tile[i] = 1.0/( this->ref_.eps1[i[2]] + this->ref_.eps1[i[3]] - this->ref_.eps1[i[0]+NO] - this->ref_.eps1[i[1]+NO]);
    *it = tile;
    }  
  
    this->T2_("a,b,i,j") = this->T2_("a,b,i,j") * oneoverD_abij("a,b,i,j");
  
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_1() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    B_1("k,i") =    conj(this->antiSymMoints["aijk"]("c,i,k,l")) * this->T1_("c,l");
    B_1("k,i") +=  0.5  *   conj(this->antiSymMoints["abij"]("d,c,k,l")) * this->T2_("d,c,l,i");
    B_1("k,i") += -fockMatrix_ta["oo"]("k,i");
    TA::TArray<MatsT> B_11(TA::get_default_world(),TA::TiledRange{orange,vrange});
    B_11.fill(0.0);
    B_11("k,c") = conj(this->antiSymMoints["abij"]("d,c,k,l")) * this->T1_("d,l"); 
    B_11("k,c") += -fockMatrix_ta["ov"]("k,c");
    B_1("k,i") += B_11("k,c") * this->T1_("c,i");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_2() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    B_2("b,c") =    conj(this->antiSymMoints["abci"]("d,c,b,k")) * this->T1_("d,k");
    B_2("b,c") +=  - 0.5  *   conj(this->antiSymMoints["abij"]("c,d,l,k")) * this->T2_("d,b,l,k"); 
    B_2("b,c") += -fockMatrix_ta["vv"]("b,c");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_3() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    B_3("k,a,i,j") = 0.5  *   conj(this->antiSymMoints["abci"]("c,d,a,k")) * this->T2_("d,c,i,j"); 
    B_3("k,a,i,j") +=    this->antiSymMoints["aijk"]("a,k,j,i");
    TA::TArray<MatsT> B_31(TA::get_default_world(),TA::TiledRange{orange,vrange});
    B_31.fill(0.0);
    TA::TArray<MatsT> B_32(TA::get_default_world(),TA::TiledRange{orange,orange,orange,orange});
    B_32.fill(0.0);
    TA::TArray<MatsT> B_321(TA::get_default_world(),TA::TiledRange{orange,orange,vrange,orange});
    B_321.fill(0.0);
    TA::TArray<MatsT> B_33(TA::get_default_world(),TA::TiledRange{orange,vrange,vrange,orange});
    B_33.fill(0.0);
    TA::TArray<MatsT> B_34(TA::get_default_world(),TA::TiledRange{orange,orange,vrange,orange});
    B_34.fill(0.0);
    B_321("k,l,c,j") += - 0.25   *  conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T1_("d,j");
    B_321("k,l,c,j") += - 0.5  *   conj(this->antiSymMoints["aijk"]("c,j,k,l"));
    B_31("k,c") +=    conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T1_("d,l");
    B_31("k,c") +=   fockMatrix_ta["ov"]("k,c");
    B_32("k,l,i,j") += -0.25   *  conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T2_("c,d,i,j");
    B_32("k,l,i,j") += - 0.5  *   this->antiSymMoints["ijkl"]("k,l,i,j");
    B_32("k,l,i,j") +=   B_321("k,l,c,j") * this->T1_("c,i");
    B_32("k,l,i,j") += -B_321("k,l,c,i") * this->T1_("c,j");
    B_33("k,a,c,j") +=  0.5  *   conj(this->antiSymMoints["abci"]("d,c,a,k")) * this->T1_("d,j");
    B_33("k,a,c,j") +=    this->antiSymMoints["iabj"]("k,a,c,j");
    B_34("l,k,c,j") +=   - conj(this->antiSymMoints["abij"]("c,d,l,k")) * this->T1_("d,j"); 
    B_34("l,k,c,j") +=   - conj(this->antiSymMoints["aijk"]("c,j,l,k"));  
    B_3("k,a,i,j") += B_31("k,c") * this->T2_("c,a,i,j");
    B_3("k,a,i,j") += B_32("k,l,i,j") * this->T1_("a,l");
    B_3("k,a,i,j") +=   B_33("k,a,c,j") * this->T1_("c,i");
    B_3("k,a,i,j") += -B_33("k,a,c,i") * this->T1_("c,j");
    B_3("k,a,i,j") +=   B_34("l,k,c,j") * this->T2_("c,a,l,i");
    B_3("k,a,i,j") += -B_34("l,k,c,i") * this->T2_("c,a,l,j");
  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_5() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    B_5("k,l,i,j") = 0.25   *  conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T2_("c,d,i,j");
    B_5("k,l,i,j") +=  0.5  *   this->antiSymMoints["ijkl"]("k,l,i,j");
    TA::TArray<MatsT> B_51(TA::get_default_world(),TA::TiledRange{orange,orange,vrange,orange});
    B_51.fill(0.0);
    B_51("k,l,c,i") += 0.25   *  conj(this->antiSymMoints["abij"]("d,c,k,l")) * this->T1_("d,i"); 
    B_51("k,l,c,i") += - 0.5  *   conj(this->antiSymMoints["aijk"]("c,i,k,l"));
    B_5("k,l,i,j") +=   B_51("k,l,c,i") * this->T1_("c,j");
    B_5("k,l,i,j") += -B_51("k,l,c,j") * this->T1_("c,i");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_6() {
    size_t NO = this->ref_.nO;
    size_t NV = this->ref_.nV;
    B_6("k,a,c,j") = - conj(this->antiSymMoints["abci"]("d,c,a,k")) * this->T1_("d,j");
    B_6("k,a,c,j") += - 0.5  *   conj(this->antiSymMoints["abij"]("c,d,k,l")) * this->T2_("d,a,l,j");
    B_6("k,a,c,j") += - this->antiSymMoints["iabj"]("k,a,c,j");
  }
  
};


