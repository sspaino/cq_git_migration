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
#ifdef CQ_HAS_TA
#pragma once

#include <singleslater.hpp>
#include <particleintegrals.hpp>
#include <wavefunction.hpp>
#include <chronusq_sys.hpp>
#include <singleslater/hartreefock.hpp>
#include <util/files.hpp>
#include <util/math.hpp>
#include <integrals.hpp>
#include <cerr.hpp>
#include <tiledarray.h>
#include "./coupledcluster/DiskDIIS.hpp"


namespace ChronusQ {

  enum class CC_TYPE { CCSD};

  struct CoupledClusterSettings {
    CC_TYPE cctype = CC_TYPE::CCSD;
    double eConv = 1e-8;
    double tConv = 1e-6;
    int maxiter = 1000;
  };

  struct CCBase
  {
    virtual void getCorrEnergy() = 0;
    virtual void run() = 0;
    bool useDIIS = true;
    size_t nDIIS = 8;

    CoupledClusterSettings ccSettings;
  };

  template <typename MatsT, typename IntsT>
  class CCSD : public CCBase
  {
    using TArray = TA::TArray<MatsT>;
  protected:
    TArray T1_;
    TArray T2_;
    SingleSlater<MatsT,IntsT>& ref_;

  public:    
    CCSD(SingleSlater<MatsT,IntsT>& ref): ref_(ref){}
    
    ~CCSD(){};
    
    std::map<std::string,TArray> antiSymMoints;
    std::map<std::string,TArray> fockMatrix_ta;
    
    double CorrE = 0.0;

    void getCorrEnergy() override;
    void updateT1(const TArray T1_old, const TArray T2_old);
    void updateT2(const TArray T1_old, TArray T2_old);
    void transformInts();
    void run() override ;
    void runConventional();

//Auxiliary variables and functions to help initialization of TA tensors
    TA::TiledRange1 orange;
    TA::TiledRange1 vrange;

    std::size_t blksize{4};

    void initRanges();

//Singles intermediates and tensors
    void formA_1();
    void formA_2();
    void formA_3();
    void formA_4();
    TArray A_1;
    TArray A_2; 
    TArray A_3;
    TArray A_4;
//Doubles intermediates and tensors
    void formB_1();
    void formB_2();
    void formB_3();
    void formB_5();
    void formB_6();   
    TArray B_1;
    TArray B_2;
    TArray B_3;
    TArray B_5;
    TArray B_6;      

    void initIntermediates();
    void initAmplitudes();
    void initFock();
    void doDIIS(size_t t1offset, size_t NT1blks, size_t NT2blks, TArray T1_old, TArray T2_old, std::shared_ptr<DiskDIIS<MatsT>> diis );

  };// namespace CCSD

}; // namespace ChronusQ
#endif
