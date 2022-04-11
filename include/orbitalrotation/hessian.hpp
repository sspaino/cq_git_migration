/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#include <orbitalrotation.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/ortho.hpp>

// #define DEBUG_ORBITALROTATION_HESSIAN

namespace ChronusQ {
  
  
  /*
   * Build orbital hessian diagonals block-wisely:
   *   IN -> inactive,      indexed by i, j, k, l
   *   CO -> correlated,    indexed by t, u, v, w
   *   FV -> frozen virtual indexed by a, b, c, d
   */
  template <typename MatsT, typename IntsT>  
  void OrbitalRotation<MatsT, IntsT>::computeOrbOrbHessianDiag(EMPerturbation & pert, 
    SquareMatrix<MatsT> & oneRDM, InCore4indexTPI<MatsT> & twoRDM, MatsT * H) {
    
    bool oneETermApproximation = settings.alg == ORB_ROT_APPROX_QUASI_2ND_ORDER;
    MatsT notRotatedHessian = MatsT (1 / settings.hessianDiagScale);

    auto & mopart = mcwfn_.MOPartition;
    auto & mem    = mcwfn_.memManager;

    size_t nTOrb   = mopart.nMO;
    size_t nCorrO  = mopart.nCorrO;
    size_t nInact  = mopart.nInact;
    size_t nFVirt  = mopart.nFVirt;
    
    size_t nINCO   = nInact + nCorrO;
    size_t SCRDim  = std::max(nCorrO, std::max(nInact, nFVirt));
    size_t SCRDim2 = SCRDim * SCRDim; 
    double fc   = (mcwfn_.reference().nC > 1) ? 1.0: 2.0;
    
    // Allocate SCR
    MatsT * SCR   = mem.template malloc<MatsT>(SCRDim2);
    MatsT * F1_ii = nullptr, * F1_aa = nullptr; 
    MatsT * RDM1_tt = mem.template malloc<MatsT>(nCorrO); 
    MatsT * F1_tt = mem.template malloc<MatsT>(nCorrO);
    MatsT * F2_tt = mem.template malloc<MatsT>(nCorrO);
    
    // populate One electron terms
    for(auto t = 0ul; t < nCorrO; t++) RDM1_tt[t] = oneRDM(t,t) / fc;
    //  F1_tt -> SCR
    this->formGeneralizedFock1(pert, oneRDM, F1_tt, "tt", true); 
    
    //  F2_tt -> SCR
    this->formGeneralizedFock2(pert, oneRDM, twoRDM, SCR, "t"); 
    for(auto t = 0ul; t < nCorrO; t++) F2_tt[t] = SCR[t + t * nCorrO]; 
    
    if (nInact > 0) {
      //  F1_ii -> SCR
      F1_ii = mem.template malloc<MatsT>(nInact);
      this->formGeneralizedFock1(pert, oneRDM, F1_ii, "ii", true);
    }
    
    if (nFVirt > 0) { 
      //  F1_aa -> SCR
      F1_aa = mem.template malloc<MatsT>(nFVirt);
      this->formGeneralizedFock1(pert, oneRDM, F1_aa, "aa", true); 
    }
  
    // put one to the IN-IN and FV-FV blocks
    std::fill_n(SCR, SCRDim2, notRotatedHessian);
    SetMat('N', nInact, nInact, MatsT(1.), SCR, nInact, H, nTOrb);
    SetMat('N', nFVirt, nFVirt, MatsT(1.), SCR, nFVirt, H + nINCO * (nTOrb+1), nTOrb);

    // CO-CO block
    if (settings.rotate_within_correlated) {
      CErr("rotate MOs within correlated orbitals not implemented");
    } else {
      SetMat('N', nCorrO, nCorrO, MatsT(1.), SCR, nCorrO, H + nInact * (nTOrb + 1), nTOrb);
    }
    
    // IN-CO block
    if (settings.rotate_inact_correlated) {
      // H_{it,it} = F1_tt + 1RDM_tt * F1_ii - F1_ii - F2_tt
#pragma omp parallel for schedule(static) default(shared)       
      for (auto t = 0ul; t < nCorrO; t++)
      for (auto i = 0ul; i < nInact; i++)
        SCR[i + t * nInact] = F1_tt[t] + RDM1_tt[t] * F1_ii[i] - F2_tt[t]- F1_ii[i];
    } else {
      std::fill_n(SCR, nCorrO * nInact, notRotatedHessian);
    } 
    
    SetMat('N', nInact, nCorrO, MatsT(1.), SCR, nInact, H + nTOrb * nInact, nTOrb);
    SetMat('C', nInact, nCorrO, MatsT(1.), SCR, nInact, H + nInact, nTOrb);

    // IN-FV block
    if (settings.rotate_inact_virtual) {
      // H_{ia,ia} = F1_aa - F1_ii  
      for (auto a = 0ul; a < nFVirt; a++)
      for (auto i = 0ul; i < nInact; i++)
        SCR[i + a * nInact] = F1_aa[a] - F1_ii[i];
    } else {
      std::fill_n(SCR, nInact * nFVirt, notRotatedHessian);
    }
    
    SetMat('N', nInact, nFVirt, MatsT(1.), SCR, nInact, H + nTOrb * nINCO, nTOrb);
    SetMat('C', nInact, nFVirt, MatsT(1.), SCR, nInact, H + nINCO, nTOrb);

    // CO-FV block
    if (settings.rotate_correlated_virtual) { 
      // H_{ta,ta} = 1RDM_tt * F1_aa - F2_tt
#pragma omp parallel for schedule(static) default(shared)       
      for (auto a = 0ul; a < nFVirt; a++)
      for (auto t = 0ul; t < nCorrO; t++)
        SCR[t + a * nCorrO] = RDM1_tt[t] * F1_aa[a] - F2_tt[t];
    } else {
      std::fill_n(SCR, nCorrO * nFVirt, notRotatedHessian);
    }

    SetMat('N', nCorrO, nFVirt, MatsT(1.), SCR, nCorrO, H + nTOrb * nINCO + nInact, nTOrb);
    SetMat('C', nCorrO, nFVirt, MatsT(1.), SCR, nCorrO, H + nTOrb * nInact + nINCO, nTOrb);
    
    // Free SCR
    if(SCR)     mem.free(SCR); 
    if(F1_ii)   mem.free(F1_ii);
    if(F2_tt)   mem.free(F2_tt);
    if(F1_tt)   mem.free(F1_tt);
    if(F1_aa)   mem.free(F1_aa);
    if(RDM1_tt) mem.free(RDM1_tt);
    
    // Scale Hessian
    blas::scal(nTOrb * nTOrb, MatsT(settings.hessianDiagScale), H, 1);

/*
 *  TODO: Test Damping 
 */
    // Dampping and taking care of negative or small eigenvalues 
#pragma omp parallel for schedule(static) default(shared)       
    for (auto pq = 0ul; pq < nTOrb * nTOrb; pq++) {
//      if (std::abs(OOHessianD[i]) > this->hessianDiagDampingThreshold)
//        OOHessianD[i] /= this->hessianDiagDamping;
//      else 
      if (std::abs(H[pq]) < settings.hessianDiagMinTol) H[pq]  = MatsT(1.0e6);
    }
//*/

#ifdef DEBUG_ORBITALROTATION_HESSIAN 
  prettyPrintSmart(std::cout, "OR Orbital-Orbital Hessian Approximated Diagonal",
    H, nTOrb, nTOrb, nTOrb);
  double HDNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, H, nTOrb); 
  std::cout << "Hessian Approximated Diagonal Norm = " << HDNorm << std::endl;
#endif  

  }; // OrbitalRotation<MatsT>::computeOrbOrbHessianDiag

}; // namespace ChronusQ
