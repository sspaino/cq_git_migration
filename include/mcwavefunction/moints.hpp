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

#include <mcwavefunction.hpp>
#include <mointstransformer/impl.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>

#include <util/matout.hpp>

// #define DEBUG_MCWFN_MOINTS

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::setMORanges() {
    mointsTF->setMORanges(dynamic_cast<const MCWaveFunctionBase &>(*this));
  }; // MCWaveFunction::setMORanges

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::transformInts(EMPerturbation & pert) {
    
    // TODO: expand to RI
    // build object and allocate memory
    size_t nTOrb  = this->MOPartition.nMO;
    size_t nCorrO = this->MOPartition.nCorrO;
    size_t nInact = this->MOPartition.nInact;
    size_t nInact2 = nInact * nInact;
    
    auto & mem   = this->memManager;

    /*
     * compute inactive core energy
     */
    MatsT * h1e_ii  = mem.template malloc<MatsT>(nInact);
    MatsT * GDjj_ii = mem.template malloc<MatsT>(nInact);

    double fc1C = (this->reference().nC == 1) ? 2.0 : 1.0;
    
    mointsTF->transformHCore(pert, h1e_ii, "ii", true);
    mointsTF->transformGD(pert, 'i', GDjj_ii, "jj", true, true, "WithInactive-i");  
    
    MatsT ECore = 0.;
    // compute core enenrgy
    for (auto i = 0ul; i < nInact; i++) {
      ECore += h1e_ii[i] + 0.5 * GDjj_ii[i];
    }
    
    mem.free(h1e_ii, GDjj_ii);

    this->InactEnergy = std::real(ECore) * fc1C;
    
    /*
     * compute hCore and ERI in correlated space
     */ 
    OnePInts<MatsT> hCore_tu(mem, nCorrO);
    OnePInts<MatsT> hCoreP_tu(mem, nCorrO);
    InCore4indexTPI<MatsT> ERI_tuvw(mem, nCorrO);
    
    mointsTF->transformHCore(pert, hCore_tu.pointer(), "tu", false, 'i');
    mointsTF->transformTPI(pert, ERI_tuvw.pointer(), "tuvw", this->cacheHalfTransTPI_);
    
    /*
     * compute hCoreP
     */
    hCoreP_tu.clear();
    if (this->MOPartition.scheme == RAS) { // form hCoreP for RAS case
      std::cout<<"RAS hcore prime forming starts:"<<std::endl;
      std::vector<size_t> nActO = this->MOPartition.nActOs;
      std::vector<std::pair<size_t, size_t>> rasloop;
      rasloop.push_back(std::make_pair(0, nActO[0]));
      rasloop.push_back(std::make_pair(nActO[0], nActO[0]+nActO[1]));
      rasloop.push_back(std::make_pair(nActO[0]+nActO[1], nCorrO));

      for (auto i = 0; i < nCorrO; i++)
      for (auto j = 0; j < nCorrO; j++)
        hCoreP_tu(i, j) = hCore_tu(i, j);
      double symmFc;
      for (auto pblk = 0; pblk < 3; pblk++)
      for (auto sblk = 0; sblk < 3; sblk++)
      for (auto qrblk = 0; qrblk < 3; qrblk++) {
        if ((pblk + 3*qrblk) < (qrblk + 3*sblk)) symmFc = 1.0;
        else if ((pblk + 3*qrblk) == (qrblk + 3*sblk)) symmFc = 0.5;
        else continue;

        for (auto p = rasloop[pblk].first; p < rasloop[pblk].second; p++)
        for (auto s = rasloop[sblk].first; s < rasloop[sblk].second; s++)
        for (auto qr = rasloop[qrblk].first; qr < rasloop[qrblk].second; qr++)
          if (pblk == qrblk or sblk == qrblk)
            hCoreP_tu(p, s) -= symmFc * ERI_tuvw(p, qr, qr, s);
          else if (pblk < qrblk)
            hCoreP_tu(p, s) = hCoreP_tu(p, s) - ERI_tuvw(p, qr, qr, s) + ERI_tuvw(qr, qr, p, s);
      }
    } else {

#pragma omp parallel for schedule(static) collapse(2) default(shared)       
      for (auto u = 0ul; u < nCorrO; u++) 
      for (auto t = 0ul; t < nCorrO; t++) {

        MatsT tmp = 0.;
        for (auto v = 0ul; v < nCorrO; v++)
          tmp += 0.5 * ERI_tuvw(t, v, v, u);

        hCoreP_tu(t, u) = hCore_tu(t, u) - tmp;
      }
    }

    this->moints.addIntegral("hCore_Correlated_Space", 
      std::make_shared<OnePInts<MatsT>>(hCore_tu));
    this->moints.addIntegral("hCoreP_Correlated_Space", 
      std::make_shared<OnePInts<MatsT>>(hCoreP_tu));
    this->moints.addIntegral("ERI_Correlated_Space",  
      std::make_shared<InCore4indexTPI<MatsT>>(ERI_tuvw));

  }; // MCWaveFunction::transformInts


}; // namespace ChronusQ
