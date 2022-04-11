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

#include <matrix.hpp>
#include <chronusq_sys.hpp>

namespace ChronusQ {

enum ORTHO_TYPE { LOWDIN, CHOLESKY };   ///< Orthonormalization Scheme

template<typename MatsT>
class Orthogonalization {
private:
  std::shared_ptr<SquareMatrix<MatsT>> overlap;         ///< Shared_ptr to overlap Matrix
  std::shared_ptr<SquareMatrix<MatsT>> forwardTrans;    ///< Transformation from the nonorthogonal basis to orthogonal basis (S^{-1/2})
  std::shared_ptr<SquareMatrix<MatsT>> backwardTrans;   ///< Transformation from the orthogonal basis to the nonorthogonal basis (S^{1/2})
  ORTHO_TYPE orthoType = LOWDIN;                        ///< Using Lowdin Type Orthogonalization

public:
  // Constructors
  Orthogonalization() {};
  Orthogonalization(SquareMatrix<MatsT>& s) {setOverlap(s);};

  // Copy/Move Constructors
  Orthogonalization(Orthogonalization<MatsT>&)  = default;
  Orthogonalization(Orthogonalization<MatsT>&&) = default;

  template<typename MatsU>
  Orthogonalization(Orthogonalization<MatsU>& other):
    overlap(std::make_shared<SquareMatrix<MatsT>>(*(other.overlapPointer()) )),
    forwardTrans(std::make_shared<SquareMatrix<MatsT>>(*(other.forwardPointer()) )),
    backwardTrans(std::make_shared<SquareMatrix<MatsT>>(*(other.backwardPointer()) )),
    orthoType(other.getOrthoType()) 
  {};

  template<typename MatsU>
  Orthogonalization(Orthogonalization<MatsU>&& other):
    overlap(std::make_shared<SquareMatrix<MatsT>>(std::move(*(other.overlapPointer()) ))),
    forwardTrans(std::make_shared<SquareMatrix<MatsT>>(std::move(*(other.forwardPointer()) ))),
    backwardTrans(std::make_shared<SquareMatrix<MatsT>>(std::move(*(other.backwardPointer()) ))),
    orthoType(other.getOrthoType()) 
  {};

  // Destructor
  ~Orthogonalization() {};

  // Getter/Setters for overlap
  inline std::shared_ptr<SquareMatrix<MatsT>> overlapPointer() const {
    if( not overlap ) CErr("Overlap has not been initialized");
    return overlap;
  };

  inline bool hasOverlap() const { return (overlap != nullptr); };

  // Change out the Overlap Matrix
  void setOverlap(SquareMatrix<MatsT>& s) {
    if (not overlap or s.dimension() != overlap->dimension()) {
      size_t nDim   = s.dimension();
      if (overlap and nDim != overlap->dimension()) {
        std::cout << std::endl;
        std::cout << "WARNING: The overlap matrix given to setOverlap is a different dimension" << std::endl;
        std::cout << "WARNING: Changing the dimension of overlap and setting the new matrix" << std::endl << std::endl;
      }
      overlap       = std::make_shared<SquareMatrix<MatsT>>(s.memManager(), nDim);
      forwardTrans  = std::make_shared<SquareMatrix<MatsT>>(s.memManager(), nDim);
      backwardTrans = std::make_shared<SquareMatrix<MatsT>>(s.memManager(), nDim);
    }
    *overlap    = s;
    computeOrtho();
  }

  // Getter functions for Transformations
  inline void setOrthoType(ORTHO_TYPE o) { orthoType = o; };
  inline ORTHO_TYPE getOrthoType() const {return orthoType; };
  inline std::shared_ptr<SquareMatrix<MatsT>> forwardPointer() const {
    if( not overlap ) CErr("Overlap has not been initialized");
    return forwardTrans;
  };
  inline std::shared_ptr<SquareMatrix<MatsT>> backwardPointer() const {
    if( not overlap ) CErr("Overlap has not been initialized");
    return backwardTrans;
  };

  // Member Functions
  void computeOrtho();

  // Transform Operators
  SquareMatrix<MatsT> nonortho2ortho(SquareMatrix<MatsT>&) const;
  SquareMatrix<MatsT> ortho2nonortho(SquareMatrix<MatsT>&) const;
  PauliSpinorSquareMatrices<MatsT> nonortho2ortho(PauliSpinorSquareMatrices<MatsT>&) const;
  PauliSpinorSquareMatrices<MatsT> ortho2nonortho(PauliSpinorSquareMatrices<MatsT>&) const;

  // Transform the basis for Coefficients
  void nonortho2orthoCoeffs(SquareMatrix<MatsT>&) const;
  void nonortho2orthoCoeffs(std::vector<SquareMatrix<MatsT>>&) const;
  void nonortho2orthoCoeffs(std::vector<std::reference_wrapper<SquareMatrix<MatsT>>>&) const;
  void ortho2nonorthoCoeffs(SquareMatrix<MatsT>&) const;
  void ortho2nonorthoCoeffs(std::vector<SquareMatrix<MatsT>>&) const;
  void ortho2nonorthoCoeffs(std::vector<std::reference_wrapper<SquareMatrix<MatsT>>>&) const;

  // Orthogonalize States e.g. to orthogonalize the occupied MO's but virtual orbitals are not transformed
  void orthogonalizeStates(SquareMatrix<MatsT>& mo, size_t nStates, size_t disp=0) const;
    // nStates = number of States to orthogonalize
    // disp    = number of vectors to shift by e.g. NBC/2 to get to the positive energy states for 4C

  void getOrthogonalizationGradients(std::vector<SquareMatrix<MatsT>>&,
    std::vector<SquareMatrix<MatsT>>&);
};

};  // namespace ChronusQ
