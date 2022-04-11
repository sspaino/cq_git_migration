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

// #define DEBUG_ORBITALROTATION_GRADIENT

namespace ChronusQ {
  
  /*
   * Build gradient block-wisely:
   *   IN -> inactive
   *   CO -> correlated
   *   FV -> frozen virtual
   */
  template <typename MatsT, typename IntsT>  
  double OrbitalRotation<MatsT, IntsT>::computeOrbGradient(EMPerturbation & pert,
    SquareMatrix<MatsT> & oneRDM, InCore4indexTPI<MatsT> & twoRDM) {
    
    ProgramTimer::tick("Form Gradient");
    
    auto & mopart = mcwfn_.MOPartition;
    auto & mem    = mcwfn_.memManager;
    if (not orbitalGradient_) 
      orbitalGradient_ = std::make_shared<SquareMatrix<MatsT>>(mem, mopart.nMO);

    MatsT * G = orbitalGradient_->pointer();

    size_t nTOrb   = mopart.nMO;
    size_t nCorrO  = mopart.nCorrO;
    size_t nInact  = mopart.nInact;
    size_t nFVirt  = mopart.nFVirt;
    size_t nINCO   = nInact + nCorrO;
    size_t SCRDim  = std::max(nCorrO, std::max(nInact, nFVirt));
    
    // Allocate SCR
    MatsT * SCR  = mem.template malloc<MatsT>(SCRDim*SCRDim);
    MatsT * SCR2 = mem.template malloc<MatsT>(nCorrO*nInact);

    // zero out IN-IN and FV-FV block s
    SetMat('N', nInact, nInact, MatsT(0.), SCR, nInact, G, nTOrb);
    SetMat('N', nFVirt, nFVirt, MatsT(0.), SCR, nFVirt, G + nINCO * (nTOrb + 1), nTOrb);
    
    // CO-CO block
    if (settings.rotate_within_correlated) {
      CErr("rotate MOs within correlated orbitals not implemented");
    } else {
      SetMat('N', nCorrO, nCorrO, MatsT(0.), SCR, nCorrO, G + nInact * (nTOrb + 1), nTOrb);
    }
    
    MatsT fc;

    // IN-CO block
    if (settings.rotate_inact_correlated) {
      // g_it = -2 * (F1_it - F2_ti)
      // populate SCR as F1_it and SCR2 as F2_ti
      this->formGeneralizedFock1(pert, oneRDM, SCR, "it");
      this->formGeneralizedFock2(pert, oneRDM, twoRDM, SCR2, "i");
      
      IMatCopy('T', nCorrO, nInact, MatsT(1.), SCR2, nCorrO, nInact);
      MatAdd('N', 'N', nInact, nCorrO, MatsT(1.), SCR, nInact,
        MatsT(-1.), SCR2, nInact, SCR, nInact);
      //  blas::scal(nInact * nCorrO, MatsT(2.), SCR, 1);
      fc = MatsT(1.);
    } else {
      fc = MatsT(0.); 
    }
    
    SetMat('N', nInact, nCorrO, -fc, SCR, nInact, G + nTOrb * nInact, nTOrb);
    SetMat('C', nInact, nCorrO,  fc, SCR, nInact, G + nInact, nTOrb);

    // IN-FV block
    if (settings.rotate_inact_virtual) {
      // g_ia = -2 * F1_ia
      // pupulate SCR as F1_ia
      this->formGeneralizedFock1(pert, oneRDM, SCR, "ia");
      //  blas::scal(nInact * nFVirt, MatsT(2.), SCR, 1);
      fc = MatsT(1.);
    } else {
      fc = MatsT(0.); 
    }

    SetMat('N', nInact, nFVirt, -fc, SCR, nInact, G + nTOrb * nINCO, nTOrb);
    SetMat('C', nInact, nFVirt,  fc, SCR, nInact, G + nINCO, nTOrb);

    // CO-FV block
    if (settings.rotate_correlated_virtual) { 
       // g_ta = -2 * F2_ta^*
       // pupulate SCR as F2_ta
      this->formGeneralizedFock2(pert, oneRDM, twoRDM, SCR, "a");
      // blas::scal(nCorrO * nFVirt, MatsT(2.), SCR, 1); 
      fc = MatsT(1.);
    } else {
      fc = MatsT(0.); 
    }

    SetMat('R', nCorrO, nFVirt, -fc, SCR, nCorrO, G + nTOrb * nINCO + nInact, nTOrb);
    SetMat('T', nCorrO, nFVirt,  fc, SCR, nCorrO, G + nTOrb * nInact + nINCO, nTOrb);
  
    mem.free(SCR, SCR2);
    
    double orbitalGradientNorm =  lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, G, nTOrb); 

#ifdef DEBUG_ORBITALROTATION_GRADIENT
    prettyPrintSmart(std::cout, " 1RDM ", oneRDM.pointer(), nCorrO, nCorrO, nCorrO);
    prettyPrintSmart(std::cout, " OR Orbital Gradient ", G, nTOrb, nTOrb, nTOrb);
    std::cout << "OR Orbital Gradient Norm = " << orbitalGradientNorm << std::endl;
#endif  
    
    ProgramTimer::tock("Form Gradient");
    
    return orbitalGradientNorm;
  }; // OrbitalRotation<MatsT>::computeOrbGradient

}; // namespace ChronusQ
