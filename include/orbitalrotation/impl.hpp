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
#include <util/print.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/ortho.hpp>

//#define DEBUG_ORBITALROTATION_IMPL

namespace ChronusQ {
  
  void OrbitalRotationSettings::print() {
    if (alg == OrbitalRotationAlgorithm::ORB_ROT_APPROX_QUASI_2ND_ORDER) {
      FormattedLine(std::cout,"  Oribital Rotation Algorithm:",  "Approximated Quasi-2nd Order");
    } else if (alg == OrbitalRotationAlgorithm::ORB_ROT_QUASI_2ND_ORDER) {
      FormattedLine(std::cout,"  Oribital Rotation Algorithm:",  "Quasi-2nd Order");
    } else if (alg == ORB_ROT_2ND_ORDER) {
      FormattedLine(std::cout,"  Oribital Rotation Algorithm:",  "Seconnd Order");
    } else CErr("NYI Orbital Rotation Algorithm");
   
    FormattedLine(std::cout, "  Rotation Blocks:");
    FormattedLine(std::cout, "  - Within Correlated:",               rotate_within_correlated);
    FormattedLine(std::cout, "  - Between Inactive and Correlated:", rotate_inact_correlated);
    FormattedLine(std::cout, "  - Between Inactive and Virtual:",    rotate_inact_virtual);
    FormattedLine(std::cout, "  - Between Correlated and Virtual:",  rotate_correlated_virtual);

  } // OrbitalRotationSettings::print
  
  /*
   * \brief rotate MO for one step size using Newton-Raphson:
   * 
   * steps:
   *   1. compute Gradient if not computed
   *   2. compute (approximated) Hessian and inverse it
   *   3. compute rotation paramemter X = - g / H
   *   4. compute rotation matrix by matrix exponential U = exp(X)
   *   5. rotate MO for one step size as C' = C U
   */ 

  template <typename MatsT, typename IntsT>  
  void OrbitalRotation<MatsT, IntsT>::rotateMO(EMPerturbation & pert, 
    SquareMatrix<MatsT> & oneRDM, InCore4indexTPI<MatsT> & twoRDM) {
    
    auto & mopart = mcwfn_.MOPartition;
    auto & mem    = mcwfn_.memManager;
    auto & mo     = mcwfn_.reference().mo[0];
    size_t nTOrb  = mopart.nMO;
    size_t nTOrb2 = nTOrb * nTOrb;
    size_t nAO    = mo.dimension(); 
    
    // offset it by inactive core index
    auto mo_pointer = mo.pointer() +
      mcwfn_.mointsTF->parseMOType('i').first * nAO;
    
    // get gradient
    if (not orbitalGradient_) computeOrbGradient(pert, oneRDM, twoRDM); 
    MatsT * G = orbitalGradient_->pointer();

    // allocate memory    
    MatsT * X = mem.template malloc<MatsT>(std::max(nAO,nTOrb)*nTOrb);
    MatsT * U = mem.template malloc<MatsT>(nTOrb2);

    MatsT * H = nullptr; 
    if (settings.alg == ORB_ROT_2ND_ORDER) {
      H = mem.template malloc<MatsT>(nTOrb2*nTOrb2); 
    } else {
      H = mem.template malloc<MatsT>(nTOrb2); 
    }
    
    ProgramTimer::tick("Form Hessian");
    // X = - H ^{-1} * G
    if (settings.alg == ORB_ROT_2ND_ORDER) {
      // this->computeOrbOrbHessian(H);
      CErr("2nd Order orbital rotation not implemented");
    } else {
      this->computeOrbOrbHessianDiag(pert, oneRDM, twoRDM, H);
      for (auto i = 0ul; i < nTOrb2; i++) X[i] = - G[i] / H[i];
    }
    ProgramTimer::tock("Form Hessian");
    
    ProgramTimer::tick("Rotate MO");
    
    // Damp X to avoid Large Rotation
    auto X_max = std::max_element(X, X+nTOrb2, 
      [&] (MatsT a, MatsT b) { return std::norm(a) < std::norm(b);} );
    
    if(std::abs(*X_max) > settings.XDampTol) {
      double XDampFc = settings.XDampTol / std::abs(*X_max);
      auto pos = std::distance(X, X_max);
      auto j   = int(pos / nTOrb);
      auto i   = pos - j*nTOrb;
      std::cout << "    WARNING!: Large Rotation! XMax at X[" 
                << i << ", " << j << "] = " << *X_max
                << ", step scaled by " << XDampFc
                << "\n" << std::endl;
      
      blas::scal(nTOrb2, MatsT(XDampFc), X, 1);
    }
    
    // U = exp(X)
    this->MatExpT('T', nTOrb, 1.0, X, nTOrb, U, nTOrb, mem);
    
#ifdef DEBUG_ORBITALROTATION_IMPL
    double HNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, H, nTOrb); 
    prettyPrintSmart(std::cout, " OR H ", H, nTOrb, nTOrb, nTOrb);
    double XNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, X, nTOrb); 
    prettyPrintSmart(std::cout, " OR X ", X, nTOrb, nTOrb, nTOrb);
    double UNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, U, nTOrb); 
    prettyPrintSmart(std::cout, " OR U ", U, nTOrb, nTOrb, nTOrb);
    std::cout << "OR H Norm = " << HNorm << std::endl;
    std::cout << "OR X Norm = " << XNorm << std::endl;
    std::cout << "OR U Norm = " << UNorm << std::endl;
#endif

    // Orthonormalized U and disable GramSchmidt printining
    std::cout.setstate(std::ios_base::failbit);
    size_t NUOrtho = GramSchmidt(nTOrb, 0, nTOrb, U, nTOrb, mem);   
    std::cout.clear();
    
    if(NUOrtho != nTOrb) CErr("Failed at Orthonormalizing Rotation U Matix.");
    
    // Rotate MO Coefficient C' = C * U
    // use X as SCR
    
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
      nAO, nTOrb, nTOrb, MatsT(1.), mo_pointer, 
      nAO, U, nTOrb, MatsT(0.), X, nAO);
    SetMat('N', nAO, nTOrb, MatsT(1.), X, nAO, mo_pointer, nAO);  
    
    ProgramTimer::tock("Rotate MO");
    
    // free memory
    if(orbitalGradient_) orbitalGradient_ = nullptr;
    if(H) mem.free(H);
    if(X) mem.free(X);
    if(U) mem.free(U);

  }; // OrbitalRotation::rotateOrbitals()

  /* 
   * Matrix Exponential Using Taylor Expansion
   * TODO: move this to MatExp
   */
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  void OrbitalRotation<MatsT, IntsT>::MatExpT(char ALG, size_t N, double ALPHA, MatsU *A, size_t LDA,
    MatsU *ExpA, size_t LDEXPA, CQMemManager &mem) {
    
    // allocate memory
    size_t N2 = N*N;
    MatsU * OddTerms  = mem.template malloc<MatsU>(N2);
    MatsU * EvenTerms = mem.template malloc<MatsU>(N2);
   
    std::fill_n(ExpA, N2, MatsU(0.));
    // form zero order term
    for (auto i = 0ul; i < N; i++) ExpA[i + i*N] = MatsU(1.);
    
    // form 1st order term
    // TODO: add Alpha scaling
    std::copy_n(A, N2, OddTerms);
    
    double residue;
    double small_number = std::numeric_limits<double>::epsilon();
    bool converged = false;
    size_t maxIter = 200;
    for (auto iter = 0ul; iter < maxIter; iter+=2) {
      // add the odd term
      MatAdd('N', 'N', N, N, MatsU(1.), OddTerms, N, MatsU(1.), ExpA, N, ExpA, N);
      
      residue = lapack::lange(lapack::Norm::Fro, N, N, OddTerms, N); 
      if (residue <= small_number) {
        converged = true;
        break;
      } 
      
      // form and add next even term
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
        N, N, N, MatsU(1./(iter+2)), A, N, OddTerms, N, 
        MatsU(0.), EvenTerms, N);
      MatAdd('N', 'N', N, N, MatsU(1.), EvenTerms, N, MatsU(1.), ExpA, N, ExpA, N);
      
      residue = lapack::lange(lapack::Norm::Fro, N, N, EvenTerms, N); 
      if (residue <= small_number) {
        converged = true;
        break;
      }

      // form next odd term
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
        N, N, N, MatsU(1./(iter+3)), A, N, EvenTerms, N, 
        MatsU(0.), OddTerms, N);
    
    }
    
    mem.free(OddTerms, EvenTerms);

    if(not converged) CErr("MatExpT failed to converge.");

    return;
  }; // MCSCF::MatExpT

}; // namespace ChronusQ

// Other headers
#include <orbitalrotation/fock.hpp>      // fock matrix implementation 
#include <orbitalrotation/gradient.hpp>  // gradient implementation
#include <orbitalrotation/hessian.hpp>   // hessian implementation
#include <orbitalrotation/ivo.hpp>   // hessian implementation


