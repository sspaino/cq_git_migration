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

namespace ChronusQ {

  template <typename ScalarT, typename MatsT>
  class ScaledSquareMatrix {

    template <typename ScalarT2, typename MatsT2>
    friend class ScaledSquareMatrix;

  protected:
    bool isPauli_; // Keep track of declaration type for static binding effect
    ScalarT scalar_;
    const SquareMatrix<MatsT> &mat_;

  public:

    ScaledSquareMatrix() = delete;
    ScaledSquareMatrix( const ScaledSquareMatrix& ) = default;
    ScaledSquareMatrix( ScaledSquareMatrix&& ) = default;
    ScaledSquareMatrix(ScalarT scalar, const SquareMatrix<MatsT> &ints, bool isPauli = false):
        isPauli_(isPauli), scalar_(scalar), mat_(ints) {}
    ScaledSquareMatrix(ScalarT scalar, const PauliSpinorSquareMatrices<MatsT> &ints):
        isPauli_(true), scalar_(scalar), mat_(ints) {}
    template <typename ScalarT1, typename ScalarT2>
    ScaledSquareMatrix(ScalarT1 scalar,
        const ScaledSquareMatrix<ScalarT2, MatsT> &scaled):
        isPauli_(scaled.isPauli_), scalar_(scalar * scaled.scalar_), mat_(scaled.mat_) {}

    bool isPauli() const { return isPauli_; }
    ScalarT scalar() const { return scalar_; }
    const SquareMatrix<MatsT>& matrix() const { return mat_; }

    ScaledSquareMatrix operator-() const {
      return ScaledSquareMatrix(-1.0, *this);
    }

  }; // class ScaledSquareMatrix

  template <typename ScalarT, typename MatsT>
  ScaledSquareMatrix<ScalarT, MatsT>
  operator*(ScalarT a, const SquareMatrix<MatsT> &ints) {
    return ScaledSquareMatrix<ScalarT, MatsT>(a, ints);
  }
  template <typename ScalarT, typename MatsT>
  ScaledSquareMatrix<ScalarT, MatsT>
  operator*(const SquareMatrix<MatsT> &ints, ScalarT a) {
    return ScaledSquareMatrix<ScalarT, MatsT>(a, ints);
  }
  template <typename ScalarT, typename MatsT>
  ScaledSquareMatrix<ScalarT, MatsT>
  operator*(ScalarT a, const PauliSpinorSquareMatrices<MatsT> &ints) {
    return ScaledSquareMatrix<ScalarT, MatsT>(a, ints);
  }
  template <typename ScalarT, typename MatsT>
  ScaledSquareMatrix<ScalarT, MatsT>
  operator*(const PauliSpinorSquareMatrices<MatsT> &ints, ScalarT a) {
    return ScaledSquareMatrix<ScalarT, MatsT>(a, ints);
  }
  template <typename ScalarT1, typename ScalarT2, typename MatsT>
  ScaledSquareMatrix<typename std::conditional<
  (std::is_same<ScalarT1, dcomplex>::value or
   std::is_same<ScalarT2, dcomplex>::value),
  dcomplex, double>::type, MatsT>
  operator*(ScalarT1 a, const ScaledSquareMatrix<ScalarT2, MatsT> &ints) {
    return ScaledSquareMatrix<typename std::conditional<
        (std::is_same<ScalarT1, dcomplex>::value or
         std::is_same<ScalarT2, dcomplex>::value),
        dcomplex, double>::type, MatsT>(a, ints);
  }
  template <typename ScalarT1, typename ScalarT2, typename MatsT>
  ScaledSquareMatrix<typename std::conditional<
  (std::is_same<ScalarT1, dcomplex>::value or
   std::is_same<ScalarT2, dcomplex>::value),
  dcomplex, double>::type, MatsT>
  operator*(const ScaledSquareMatrix<ScalarT1, MatsT> &ints, ScalarT2 a) {
    return ScaledSquareMatrix<typename std::conditional<
        (std::is_same<ScalarT1, dcomplex>::value or
         std::is_same<ScalarT2, dcomplex>::value),
        dcomplex, double>::type, MatsT>(a, ints);
  }

  template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<ScalarT1, dcomplex>::value or
   std::is_same<ScalarT2, dcomplex>::value or
   std::is_same<MatsT1, dcomplex>::value or
   std::is_same<MatsT2, dcomplex>::value),
  dcomplex, double>::type>
  operator+(const ScaledSquareMatrix<ScalarT1, MatsT1>&, const ScaledSquareMatrix<ScalarT2, MatsT2>&);

  template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<ScalarT1, dcomplex>::value or
   std::is_same<ScalarT2, dcomplex>::value or
   std::is_same<MatsT1, dcomplex>::value or
   std::is_same<MatsT2, dcomplex>::value),
  dcomplex, double>::type>
  operator-(const ScaledSquareMatrix<ScalarT1, MatsT1>&, const ScaledSquareMatrix<ScalarT2, MatsT2>&);

}; // namespace ChronusQ
