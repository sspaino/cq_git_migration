/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute itnd/or modify
 *  it under the terms of the GNU General Public Licenses published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option)ny later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUTNY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received copy of the GNU General Public Licenselong
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *
 */

#include <matrix.hpp>
#include <cqlinalg.hpp>

namespace ChronusQ {

  template ScaledSquareMatrix<double, double>::
      ScaledSquareMatrix(double, const ScaledSquareMatrix<double, double>&);
  template ScaledSquareMatrix<dcomplex, double>::
      ScaledSquareMatrix(dcomplex, const ScaledSquareMatrix<dcomplex, double>&);
  template ScaledSquareMatrix<dcomplex, double>::
      ScaledSquareMatrix(dcomplex, const ScaledSquareMatrix<double, double>&);
  template ScaledSquareMatrix<dcomplex, double>::
      ScaledSquareMatrix(double, const ScaledSquareMatrix<dcomplex, double>&);
  template ScaledSquareMatrix<double, dcomplex>::
      ScaledSquareMatrix(double, const ScaledSquareMatrix<double, dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>::
      ScaledSquareMatrix(dcomplex, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>::
      ScaledSquareMatrix(dcomplex, const ScaledSquareMatrix<double, dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>::
      ScaledSquareMatrix(double, const ScaledSquareMatrix<dcomplex, dcomplex>&);

  template ScaledSquareMatrix<double, double>
  operator*(double, const SquareMatrix<double>&);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(dcomplex, const SquareMatrix<double>&);
  template ScaledSquareMatrix<double, dcomplex>
  operator*(double, const SquareMatrix<dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(dcomplex, const SquareMatrix<dcomplex>&);
  template ScaledSquareMatrix<double, double>
  operator*(const SquareMatrix<double>&, double);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(const SquareMatrix<double>&, dcomplex);
  template ScaledSquareMatrix<double, dcomplex>
  operator*(const SquareMatrix<dcomplex>&, double);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(const SquareMatrix<dcomplex>&, dcomplex);

  template ScaledSquareMatrix<double, double>
  operator*(double, const PauliSpinorSquareMatrices<double>&);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(dcomplex, const PauliSpinorSquareMatrices<double>&);
  template ScaledSquareMatrix<double, dcomplex>
  operator*(double, const PauliSpinorSquareMatrices<dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(dcomplex, const PauliSpinorSquareMatrices<dcomplex>&);
  template ScaledSquareMatrix<double, double>
  operator*(const PauliSpinorSquareMatrices<double>&, double);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(const PauliSpinorSquareMatrices<double>&, dcomplex);
  template ScaledSquareMatrix<double, dcomplex>
  operator*(const PauliSpinorSquareMatrices<dcomplex>&, double);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(const PauliSpinorSquareMatrices<dcomplex>&, dcomplex);

  template ScaledSquareMatrix<double, double>
  operator*(double, const ScaledSquareMatrix<double, double>&);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(dcomplex, const ScaledSquareMatrix<double, double>&);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(double, const ScaledSquareMatrix<dcomplex, double>&);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(dcomplex, const ScaledSquareMatrix<dcomplex, double>&);
  template ScaledSquareMatrix<double, dcomplex>
  operator*(double, const ScaledSquareMatrix<double, dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(dcomplex, const ScaledSquareMatrix<double, dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(double, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(dcomplex, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template ScaledSquareMatrix<double, double>
  operator*(const ScaledSquareMatrix<double, double>&, double);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(const ScaledSquareMatrix<double, double>&, dcomplex);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(const ScaledSquareMatrix<dcomplex, double>&, double);
  template ScaledSquareMatrix<dcomplex, double>
  operator*(const ScaledSquareMatrix<dcomplex, double>&, dcomplex);
  template ScaledSquareMatrix<double, dcomplex>
  operator*(const ScaledSquareMatrix<double, dcomplex>&, double);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(const ScaledSquareMatrix<double, dcomplex>&, dcomplex);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(const ScaledSquareMatrix<dcomplex, dcomplex>&, double);
  template ScaledSquareMatrix<dcomplex, dcomplex>
  operator*(const ScaledSquareMatrix<dcomplex, dcomplex>&, dcomplex);

  template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<ScalarT1, dcomplex>::value or
   std::is_same<ScalarT2, dcomplex>::value or
   std::is_same<MatsT1, dcomplex>::value or
   std::is_same<MatsT2, dcomplex>::value),
  dcomplex, double>::type>
  operator+(const ScaledSquareMatrix<ScalarT1, MatsT1> &lhs, const ScaledSquareMatrix<ScalarT2, MatsT2> &rhs) {
    if (lhs.matrix().dimension() != rhs.matrix().dimension())
      CErr("Cannot add two SquareMatrix of different size.");
    typedef typename std::conditional<
        (std::is_same<ScalarT1, dcomplex>::value or
         std::is_same<ScalarT2, dcomplex>::value or
         std::is_same<MatsT1, dcomplex>::value or
         std::is_same<MatsT2, dcomplex>::value),
        dcomplex, double>::type RetT;
    if (rhs.isPauli()) {
      PauliSpinorSquareMatrices<RetT> result(rhs);
      result += lhs;
      return result;
    }
    PauliSpinorSquareMatrices<RetT> result(lhs);
    result += rhs;
    return result;
  }
  template PauliSpinorSquareMatrices<double>
  operator+(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator+(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);

  template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<ScalarT1, dcomplex>::value or
   std::is_same<ScalarT2, dcomplex>::value or
   std::is_same<MatsT1, dcomplex>::value or
   std::is_same<MatsT2, dcomplex>::value),
  dcomplex, double>::type>
  operator-(const ScaledSquareMatrix<ScalarT1, MatsT1>& lhs, const ScaledSquareMatrix<ScalarT2, MatsT2>& rhs) {
    return lhs + (-rhs);
  }
  template PauliSpinorSquareMatrices<double>
  operator-(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, double>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<double, dcomplex>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<double, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<double, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<dcomplex, double>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, double>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);
  template PauliSpinorSquareMatrices<dcomplex>
  operator-(const ScaledSquareMatrix<dcomplex, dcomplex>&, const ScaledSquareMatrix<dcomplex, dcomplex>&);

}; // namespace ChronusQ
