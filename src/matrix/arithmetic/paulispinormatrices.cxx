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
#include <matrix.hpp>

namespace ChronusQ {

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<MatsT>::PauliSpinorSquareMatrices(
      const ScaledSquareMatrix<ScalarT, MatsU> &scaled, bool addXY, bool addZ ):
      SquareMatrix<MatsT>(scaled.matrix().memManager(), scaled.matrix().dimension()) {
    size_t nZYX = addXY ? 3 : addZ ? 1 : 0;
    if (scaled.isPauli()) {
      auto scaledMat = dynamic_cast<const PauliSpinorSquareMatrices<MatsU>&>(scaled.matrix());
      SquareMatrix<MatsT>::operator=(scaled.scalar() * scaledMat.S());
      nZYX = std::max(nZYX, scaledMat.components_.size());
      components_.reserve(nZYX);
      for (auto &mat : scaledMat.components_)
        components_.emplace_back(scaled.scalar() * mat);
      while (nZYX > scaledMat.components_.size()) {
        components_.emplace_back(this->memManager(), this->dimension());
        components_.back().clear();
        nZYX--;
      }
    } else {
      SquareMatrix<MatsT>::operator=(scaled);
      components_.reserve(nZYX);
      while (nZYX > 0) {
        components_.emplace_back(this->memManager(), this->dimension());
        components_.back().clear();
        nZYX--;
      }
    }
  }
  template PauliSpinorSquareMatrices<double>::PauliSpinorSquareMatrices(
      const ScaledSquareMatrix<double, double>&, bool, bool );
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(
      const ScaledSquareMatrix<double, double>&, bool, bool );
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(
      const ScaledSquareMatrix<dcomplex, double>&, bool, bool );
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(
      const ScaledSquareMatrix<double, dcomplex>&, bool, bool );
  template PauliSpinorSquareMatrices<dcomplex>::PauliSpinorSquareMatrices(
      const ScaledSquareMatrix<dcomplex, dcomplex>&, bool, bool );

  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator=(const SquareMatrix<double>&);
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=(const SquareMatrix<double>&);
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=(const SquareMatrix<dcomplex>&);
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator=(SquareMatrix<double>&&);
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=(SquareMatrix<double>&&);
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=(SquareMatrix<dcomplex>&&);

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator=( const ScaledSquareMatrix<ScalarT, MatsU> &rhs ) {
    if (rhs.isPauli()) {
      const PauliSpinorSquareMatrices<MatsU>& r =
          dynamic_cast<const PauliSpinorSquareMatrices<MatsU>&>(rhs.matrix());
      S() = rhs.scalar() * r.S();
      if (r.hasZ()) {
        if (not hasZ())
          components_.emplace_back(this->memManager(), this->dimension());
        Z() = rhs.scalar() * r.Z();
      } else if (hasZ()) {
        Z().clear();
      }
      if (r.hasXY()) {
        if (not hasXY()) {
          components_.emplace_back(this->memManager(), this->dimension());
          components_.emplace_back(this->memManager(), this->dimension());
        }
        Y() = rhs.scalar() * r.Y();
        X() = rhs.scalar() * r.X();
      } else if (hasXY()) {
        Y().clear();
        X().clear();
      }
    } else {
      S() = rhs;
      for (auto &mat : components_) mat.clear();
    }
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator=( const ScaledSquareMatrix<double, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=( const ScaledSquareMatrix<double, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=( const ScaledSquareMatrix<dcomplex, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=( const ScaledSquareMatrix<double, dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator=( const PauliSpinorSquareMatrices<MatsT> &other ) {
    if (this != &other) { // self-assignment check expected
      return operator=(1.0 * other);
    }
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator=( const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=( const PauliSpinorSquareMatrices<dcomplex>& );

  template <typename MatsT>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator=( PauliSpinorSquareMatrices<MatsT> &&other ) {
    if (this != &other) { // self-assignment check expected
      S() = std::move(other.S());
      if (other.hasZ()) {
        if (hasZ())
          Z() = std::move(other.Z());
        else
          components_.emplace_back(std::move(other.Z()));
      } else if (hasZ()) {
        Z().clear();
      }
      if (other.hasXY()) {
        if (hasXY()) {
          Y() = std::move(other.Y());
          X() = std::move(other.X());
        } else {
          components_.emplace_back(std::move(other.Y()));
          components_.emplace_back(std::move(other.X()));
        }
      } else if (hasXY()) {
        Y().clear();
        X().clear();
      }
    }
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator=( PauliSpinorSquareMatrices<double>&& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator=( PauliSpinorSquareMatrices<dcomplex>&& );

  template <typename MatsT>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator*=( MatsT scalar ) {
    S() *= scalar;
    for (SquareMatrix<MatsT> &mat : components_)
      mat *= scalar;
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator*=( double scalar );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator*=( dcomplex scalar );

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator+=( const PauliSpinorSquareMatrices<MatsU> &other ) {
    S() += other.S();
    if (hasZ() and other.hasZ())
      Z() += other.Z();
    else if (other.hasZ())
      components_.emplace_back(other.Z());
    if (hasXY() and other.hasXY()) {
      Y() += other.Y();
      X() += other.X();
    } else if (other.hasXY()) {
      components_.emplace_back(other.Y());
      components_.emplace_back(other.X());
    }
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator+=( const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const PauliSpinorSquareMatrices<dcomplex>& );

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator+=( const SquareMatrix<MatsU> &other ) {
    S() += other;
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator+=( const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const SquareMatrix<dcomplex>& );

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator-=( const SquareMatrix<MatsU> &other ) {
    S() = other;
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator-=( const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const SquareMatrix<dcomplex>& );

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type>
  PauliSpinorSquareMatrices<MatsT>::operator+( const SquareMatrix<MatsU> &other ) const {
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> sum(*this);
    sum += other;
    return sum;
  }
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::operator+( const SquareMatrix<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator+( const SquareMatrix<dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const SquareMatrix<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const SquareMatrix<dcomplex>& ) const;

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type>
  PauliSpinorSquareMatrices<MatsT>::operator-( const SquareMatrix<MatsU> &other ) const {
    return *this + (-other);
  }
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::operator-( const SquareMatrix<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator-( const SquareMatrix<dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const SquareMatrix<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const SquareMatrix<dcomplex>& ) const;

  template PauliSpinorSquareMatrices<double>
  operator+( const SquareMatrix<double>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const SquareMatrix<dcomplex>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const SquareMatrix<double>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const SquareMatrix<dcomplex>&, const PauliSpinorSquareMatrices<dcomplex>& );

  template PauliSpinorSquareMatrices<double>
  operator-( const SquareMatrix<double>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const SquareMatrix<dcomplex>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const SquareMatrix<double>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const SquareMatrix<dcomplex>&, const PauliSpinorSquareMatrices<dcomplex>& );

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value), dcomplex, double>::type>
  PauliSpinorSquareMatrices<MatsT>::operator+( const PauliSpinorSquareMatrices<MatsU> &other ) const {
    typedef typename std::conditional<
        (std::is_same<MatsT, dcomplex>::value or
         std::is_same<MatsU, dcomplex>::value),
        dcomplex, double>::type ResultsT;
    if (other.components_.size() > components_.size()) {
      PauliSpinorSquareMatrices<ResultsT> sum(other);
      sum += *this;
      return sum;
    }
    PauliSpinorSquareMatrices<ResultsT> sum(*this);
    sum += other;
    return sum;
  }
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::operator+( const PauliSpinorSquareMatrices<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const PauliSpinorSquareMatrices<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator+( const PauliSpinorSquareMatrices<dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const PauliSpinorSquareMatrices<dcomplex>& ) const;

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value), dcomplex, double>::type>
  PauliSpinorSquareMatrices<MatsT>::operator-( const PauliSpinorSquareMatrices<MatsU> &other ) const {
    return *this + (-other);
  }
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::operator-( const PauliSpinorSquareMatrices<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const PauliSpinorSquareMatrices<double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator-( const PauliSpinorSquareMatrices<dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const PauliSpinorSquareMatrices<dcomplex>& ) const;

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator-=( const PauliSpinorSquareMatrices<MatsU> &other ) {
    return operator+=(-other);
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator-=( const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const PauliSpinorSquareMatrices<dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator+=( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) {
    if (scaled.isPauli()) {
      const PauliSpinorSquareMatrices<MatsU>& scaledMat =
          dynamic_cast<const PauliSpinorSquareMatrices<MatsU>&>(scaled.matrix());
      S() += scaled.scalar() * scaledMat.S();
      if (scaledMat.hasZ()) {
        if (hasZ())
          Z() += scaled.scalar() * scaledMat.Z();
        else {
          components_.emplace_back(this->memManager(), this->dimension());
          Z() = scaled.scalar() * scaledMat.Z();
        }
      }
      if (scaledMat.hasXY()) {
        if (hasXY()) {
          Y() += scaled.scalar() * scaledMat.Y();
          X() += scaled.scalar() * scaledMat.X();
        } else {
          components_.emplace_back(this->memManager(), this->dimension());
          Y() = scaled.scalar() * scaledMat.Y();
          components_.emplace_back(this->memManager(), this->dimension());
          X() = scaled.scalar() * scaledMat.X();
        }
      }
    } else {
      S() += scaled;
    }
    return *this;
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator+=( const ScaledSquareMatrix<double, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const ScaledSquareMatrix<double, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const ScaledSquareMatrix<dcomplex, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const ScaledSquareMatrix<double, dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator+=( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<MatsT>&
  PauliSpinorSquareMatrices<MatsT>::operator-=( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) {
    return operator+=(-scaled);
  }
  template PauliSpinorSquareMatrices<double>&
  PauliSpinorSquareMatrices<double>::operator-=( const ScaledSquareMatrix<double, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const ScaledSquareMatrix<double, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const ScaledSquareMatrix<dcomplex, double>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const ScaledSquareMatrix<double, dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>&
  PauliSpinorSquareMatrices<dcomplex>::operator-=( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type>
  PauliSpinorSquareMatrices<MatsT>::operator+( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) const {
    PauliSpinorSquareMatrices<typename std::conditional<
        (std::is_same<MatsT, dcomplex>::value or
         std::is_same<MatsU, dcomplex>::value or
         std::is_same<ScalarT, dcomplex>::value),
        dcomplex, double>::type> sum(*this);
    sum += scaled;
    return sum;
  }
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::operator+( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator+( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator+( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator+( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator+( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type>
  PauliSpinorSquareMatrices<MatsT>::operator-( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) const {
    return *this + (-scaled);
  }
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::operator-( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator-( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator-( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<double>::operator-( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::operator-( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;

  template PauliSpinorSquareMatrices<double>
  operator+( const ScaledSquareMatrix<double, double>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<double, dcomplex>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<double, double>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<double, dcomplex>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, double>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, dcomplex>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, double>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, dcomplex>&, const PauliSpinorSquareMatrices<dcomplex>& );

  template PauliSpinorSquareMatrices<double>
  operator-( const ScaledSquareMatrix<double, double>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<double, dcomplex>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<double, double>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<double, dcomplex>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, double>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, dcomplex>&, const PauliSpinorSquareMatrices<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, double>&, const PauliSpinorSquareMatrices<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, dcomplex>&, const PauliSpinorSquareMatrices<dcomplex>& );

}; // namespace ChronusQ
