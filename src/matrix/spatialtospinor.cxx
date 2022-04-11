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
#include <cqlinalg.hpp>

namespace ChronusQ {

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsU>
  SquareMatrix<MatsT>::spinScatter(bool hasXY, bool hasZ) const {
    size_t N = dimension() / 2;
    PauliSpinorSquareMatrices<MatsU> pauli(memManager(), N, hasXY, hasZ);

    MatsU *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
    if (hasZ) Z = pauli.Z().pointer();
    if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }

    SpinScatter(N, pointer(), dimension(), S, N, Z, N, Y, N, X, N);
    return pauli;
  }

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>
  PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(
      const SquareMatrix<MatsU> &AA, bool hasXY, bool hasZ) {
    size_t N = AA.dimension();
    PauliSpinorSquareMatrices<MatsT> pauli(AA.memManager(), N, hasXY, hasZ);
    MatsT *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
    if (hasZ) Z = pauli.Z().pointer();
    if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }
    SpinScatter(N, N, AA.pointer(), N, reinterpret_cast<MatsU*>(NULL), N,
        reinterpret_cast<MatsU*>(NULL), N, reinterpret_cast<MatsU*>(NULL), N,
        S, N, Z, N, Y, N, X, N, true, true);
    return pauli;
  }

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>
  PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(
      const SquareMatrix<MatsU> &AA, const SquareMatrix<MatsU> &BB,
      bool hasXY, bool hasZ) {
    size_t N = AA.dimension();
    PauliSpinorSquareMatrices<MatsT> pauli(AA.memManager(), N, hasXY, hasZ);
    MatsT *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
    if (hasZ) Z = pauli.Z().pointer();
    if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }
    SpinScatter(N, N, AA.pointer(), N, reinterpret_cast<MatsU*>(NULL), N,
        reinterpret_cast<MatsU*>(NULL), N, BB.pointer(), N,
        S, N, Z, N, Y, N, X, N, true, false);
    return pauli;
  }

  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT>
  PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(
      const SquareMatrix<MatsU> &AA, const SquareMatrix<MatsU> &AB,
      const SquareMatrix<MatsU> &BA, const SquareMatrix<MatsU> &BB,
      bool hasXY, bool hasZ) {
    size_t N = AA.dimension();
    PauliSpinorSquareMatrices<MatsT> pauli(AA.memManager(), N, hasXY, hasZ);
    MatsT *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
    if (hasZ) Z = pauli.Z().pointer();
    if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }
    SpinScatter(N, N, AA.pointer(), N, AB.pointer(), N,
                      BA.pointer(), N, BB.pointer(), N,
                      S, N, Z, N, Y, N, X, N, false, false);
    return pauli;
  }

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<MatsU> PauliSpinorSquareMatrices<MatsT>::spinGather() const {
    size_t N = this->dimension();
    SquareMatrix<MatsU> mat(this->memManager(), 2*N);

    const MatsT *AS = S().pointer(), *AZ = nullptr, *AY = nullptr, *AX = nullptr;
    if (hasZ()) AZ = Z().pointer();
    if (hasXY()) { AY = Y().pointer(); AX = X().pointer(); }

    SpinGather(N, mat.pointer(), 2*N, AS, N, AZ, N, AY, N, AX, N, not hasXY(), not hasZ());
    return mat;
  }

  template <typename MatsT>
  template <typename MatsU>
  std::vector<SquareMatrix<MatsU>>
  PauliSpinorSquareMatrices<MatsT>::spinGatherToBlocks(
      bool genABBA, bool genBB) const {
    size_t N = this->dimension();
    std::vector<SquareMatrix<MatsU>> blocks;
    blocks.reserve(1 + (genABBA ? 2 : 0) + (genBB ? 1 : 0));
    blocks.emplace_back(this->memManager(), N);
    MatsU *AA = nullptr, *AB = nullptr, *BA = nullptr, *BB = nullptr;
    if (genABBA) {
      blocks.emplace_back(this->memManager(), N);
      blocks.emplace_back(this->memManager(), N);
    }
    if (genBB) { blocks.emplace_back(this->memManager(), N); BB = blocks.back().pointer(); }
    AA = blocks[0].pointer();
    if (genABBA) { AB = blocks[1].pointer(); BA = blocks[2].pointer(); }

    const MatsT *AS = S().pointer(), *AZ = nullptr, *AY = nullptr, *AX = nullptr;
    if (hasZ()) AZ = Z().pointer();
    if (hasXY()) { AY = Y().pointer(); AX = X().pointer(); }

    SpinGather(N, N, AA, N, AB, N, BA, N, BB, N,
        AS, N, AZ, N, AY, N, AX, N, not hasXY(), not hasZ());
    return blocks;
  }

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<MatsU> SquareMatrix<MatsT>::spatialToSpinBlock() const {
    SquareMatrix<MatsU> spinor(memManager(), 2*dimension());
/*
    for ( auto sp = 0ul; sp < 2; sp++)
    for ( auto nu = 0ul; nu < N_; nu++)
    for ( auto mu = 0ul; mu < N_; mu++) {
      spinor(sp*N_ + mu, sp*N_ + nu) = (*this)(mu, nu);
    }
*/
    SetMatDiag(N_, N_, pointer(), N_, spinor.pointer(), 2*N_);
    return spinor;
  }

  template PauliSpinorSquareMatrices<double> SquareMatrix<double>::spinScatter(bool,bool) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::spinScatter(bool,bool) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::spinScatter(bool,bool) const;

  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::spinBlockScatterBuild(const SquareMatrix<double> &AA, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::spinBlockScatterBuild(const SquareMatrix<double> &AA, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::spinBlockScatterBuild(const SquareMatrix<dcomplex> &AA, bool, bool);
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::spinBlockScatterBuild(
      const SquareMatrix<double> &AA, const SquareMatrix<double> &BB, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::spinBlockScatterBuild(
      const SquareMatrix<double> &AA, const SquareMatrix<double> &BB, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::spinBlockScatterBuild(
      const SquareMatrix<dcomplex> &AA, const SquareMatrix<dcomplex> &BB, bool, bool);
  template PauliSpinorSquareMatrices<double>
  PauliSpinorSquareMatrices<double>::spinBlockScatterBuild(
      const SquareMatrix<double> &AA, const SquareMatrix<double> &AB,
      const SquareMatrix<double> &BA, const SquareMatrix<double> &BB, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::spinBlockScatterBuild(
      const SquareMatrix<double> &AA, const SquareMatrix<double> &AB,
      const SquareMatrix<double> &BA, const SquareMatrix<double> &BB, bool, bool);
  template PauliSpinorSquareMatrices<dcomplex>
  PauliSpinorSquareMatrices<dcomplex>::spinBlockScatterBuild(
      const SquareMatrix<dcomplex> &AA, const SquareMatrix<dcomplex> &AB,
      const SquareMatrix<dcomplex> &BA, const SquareMatrix<dcomplex> &BB, bool, bool);

  template SquareMatrix<double> PauliSpinorSquareMatrices<double>::spinGather() const;
  template SquareMatrix<dcomplex> PauliSpinorSquareMatrices<double>::spinGather() const;
  template SquareMatrix<dcomplex> PauliSpinorSquareMatrices<dcomplex>::spinGather() const;

  template std::vector<SquareMatrix<double>>
  PauliSpinorSquareMatrices<double>::spinGatherToBlocks(bool genABBA, bool genBB) const;
  template std::vector<SquareMatrix<dcomplex>>
  PauliSpinorSquareMatrices<double>::spinGatherToBlocks(bool genABBA, bool genBB) const;
  template std::vector<SquareMatrix<dcomplex>>
  PauliSpinorSquareMatrices<dcomplex>::spinGatherToBlocks(bool genABBA, bool genBB) const;

  template SquareMatrix<double> SquareMatrix<double>::spatialToSpinBlock() const;
  template SquareMatrix<dcomplex> SquareMatrix<double>::spatialToSpinBlock() const;
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::spatialToSpinBlock() const;

}; // namespace ChronusQ
