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

#include <matrix.hpp>

namespace ChronusQ {

  template class SquareMatrix<double>;
  template class SquareMatrix<dcomplex>;

  template class PauliSpinorSquareMatrices<double>;
  template class PauliSpinorSquareMatrices<dcomplex>;

  /**
   *  \brief The pointer convertor. This static function converts
   *  the underlying polymorphism correctly to hold a different
   *  type of matrices. It is called when the corresponding
   *  SingleSlater object is being converted.
   */
  template <typename IntsT>
  template <typename IntsU>
  std::shared_ptr<SquareMatrix<IntsU>>
  SquareMatrix<IntsT>::convert(const std::shared_ptr<SquareMatrix<IntsT>> &mat) {

    if (not mat) return nullptr;

    const std::type_info &tID(typeid(*mat));

    if (tID == typeid(SquareMatrix<IntsT>)) {
      return std::make_shared<SquareMatrix<IntsU>>(*mat);

    } else if (tID == typeid(PauliSpinorSquareMatrices<IntsT>)) {
      return std::make_shared<PauliSpinorSquareMatrices<IntsU>>(
          *std::dynamic_pointer_cast<PauliSpinorSquareMatrices<IntsT>>(mat));

    } else {
      std::stringstream errMsg;
      errMsg << "SquareMatrix implementation \"" << tID.name()
             << "\" not registered in convert." << std::endl;
      CErr(errMsg.str(),std::cout);
    }

    return nullptr;

  }

  template <typename MatsT>
  std::ostream& operator<<(std::ostream &out, const SquareMatrix<MatsT> &mat) {
    mat.output(out);
    return out;
  }

  template std::ostream& operator<<(std::ostream&, const SquareMatrix<double>&);
  template std::ostream& operator<<(std::ostream&, const SquareMatrix<dcomplex>&);

  template class ScaledSquareMatrix<double, double>;
  template class ScaledSquareMatrix<dcomplex, double>;
  template class ScaledSquareMatrix<double, dcomplex>;
  template class ScaledSquareMatrix<dcomplex, dcomplex>;

  template SquareMatrix<double>::SquareMatrix(const PauliSpinorSquareMatrices<double>&);
  template SquareMatrix<dcomplex>::SquareMatrix(const PauliSpinorSquareMatrices<double>&);
  template SquareMatrix<dcomplex>::SquareMatrix(const PauliSpinorSquareMatrices<dcomplex>&);

  template PauliSpinorSquareMatrices<double>::PauliSpinorSquareMatrices(const SquareMatrix<double>&, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(const SquareMatrix<double>&, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(const SquareMatrix<dcomplex>&, bool, bool);

  template PauliSpinorSquareMatrices<double>::PauliSpinorSquareMatrices(const PauliSpinorSquareMatrices<double>&, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(const PauliSpinorSquareMatrices<double>&, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(const PauliSpinorSquareMatrices<dcomplex>&, bool, bool);

}; // namespace ChronusQ
