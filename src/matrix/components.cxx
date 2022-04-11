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
  void SquareMatrix<MatsT>::componentScatter(SquareMatrix<MatsU> & LL,
    SquareMatrix<MatsU> & LS, SquareMatrix<MatsU> & SL, SquareMatrix<MatsU> & SS, 
    bool increment) const {
      
      // dimension check
      size_t N = this->dimension();
      if (N % 2 == 1) CErr("componentGather is only supported for SqaureMatrix with even size"); 
      N /= 2; 

      if( (LL.pointer() and N != LL.dimension()) or
          (LS.pointer() and N != LS.dimension()) or 
          (SL.pointer() and N != SL.dimension()) or
          (SS.pointer() and N != SS.dimension()))
        CErr("Mismatched dimensioins in componentGather!");
      
      ComponentScatter(N, N, this->pointer(), this->dimension(),
        LL.pointer(), N, LS.pointer(), N, SL.pointer(), N, SS.pointer(), N, increment); 
  
  }
  
  template <typename MatsT>
  template <typename MatsU>
  void SquareMatrix<MatsT>::componentGather(const SquareMatrix<MatsU> & LL,
    const SquareMatrix<MatsU> & LS, const SquareMatrix<MatsU> & SL, const SquareMatrix<MatsU> & SS, 
    bool increment) {
      
      // dimension check
      size_t N = this->dimension() / 2;

      if( (LL.pointer() and N != LL.dimension()) or
          (LS.pointer() and N != LS.dimension()) or 
          (SL.pointer() and N != SL.dimension()) or
          (SS.pointer() and N != SS.dimension()))
        CErr("Mismatched dimensioins in componentGather!");

      ComponentGather(N, N, this->pointer(), this->dimension(),
        LL.pointer(), N, LS.pointer(), N, SL.pointer(), N, SS.pointer(), N, increment); 

  }

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<MatsT> SquareMatrix<MatsT>::componentGatherBuild(const SquareMatrix<MatsU> & LL,
    const SquareMatrix<MatsU> & LS, const SquareMatrix<MatsU> & SL, const SquareMatrix<MatsU> & SS) { 
    
    // get dimension
    size_t N;
    if (LL.pointer()) N = LL.dimension();
    else if (LS.pointer()) N = LS.dimension();
    else if (SL.pointer()) N = SL.dimension();
    else if (SS.pointer()) N = SS.dimension();
    else CErr("Nothing to build in componentGatherBuild");
    
    SquareMatrix<MatsT> mat(LL.memManager(), N*2);
    mat.componentGather(LL, LS, SL, SS, false);
    
    return mat;
  
  }    
  
  template void SquareMatrix<double>::componentScatter(SquareMatrix<double> & LL,
    SquareMatrix<double> & LS, SquareMatrix<double> & SL, SquareMatrix<double> & SS, bool increment) const;
  template void SquareMatrix<double>::componentScatter(SquareMatrix<dcomplex> & LL,
    SquareMatrix<dcomplex> & LS, SquareMatrix<dcomplex> & SL, SquareMatrix<dcomplex> & SS, bool increment) const;
  template void SquareMatrix<dcomplex>::componentScatter(SquareMatrix<dcomplex> & LL,
    SquareMatrix<dcomplex> & LS, SquareMatrix<dcomplex> & SL, SquareMatrix<dcomplex> & SS, bool increment) const;

  template void SquareMatrix<double>::componentGather(const SquareMatrix<double> & LL,
    const SquareMatrix<double> & LS, const SquareMatrix<double> & SL, const SquareMatrix<double> & SS, bool increment);
  template void SquareMatrix<dcomplex>::componentGather(const SquareMatrix<double> & LL,
    const SquareMatrix<double> & LS, const SquareMatrix<double> & SL, const SquareMatrix<double> & SS, bool increment);
  template void SquareMatrix<dcomplex>::componentGather(const SquareMatrix<dcomplex> & LL,
    const SquareMatrix<dcomplex> & LS, const SquareMatrix<dcomplex> & SL, const SquareMatrix<dcomplex> & SS, bool increment);

  template SquareMatrix<double> SquareMatrix<double>::componentGatherBuild(const SquareMatrix<double> & LL,
    const SquareMatrix<double> & LS, const SquareMatrix<double> & SL, const SquareMatrix<double> & SS);
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::componentGatherBuild(const SquareMatrix<double> & LL,
    const SquareMatrix<double> & LS, const SquareMatrix<double> & SL, const SquareMatrix<double> & SS);
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::componentGatherBuild(const SquareMatrix<dcomplex> & LL,
    const SquareMatrix<dcomplex> & LS, const SquareMatrix<dcomplex> & SL, const SquareMatrix<dcomplex> & SS);
  
  template <typename MatsT>
  template <typename MatsU>
  void PauliSpinorSquareMatrices<MatsT>::componentScatter(
    PauliSpinorSquareMatrices<MatsU> & LL, 
    PauliSpinorSquareMatrices<MatsU> & LS, 
    PauliSpinorSquareMatrices<MatsU> & SL, 
    PauliSpinorSquareMatrices<MatsU> & SS, 
    bool increment) const {
      
      this->S().componentScatter(LL.S(), LS.S(), SL.S(), SS.S(), increment); 
      
      SquareMatrix<MatsU> dummy(this->memManager(), 0);
      if (this->hasZ()) {
        SquareMatrix<MatsU> & LLZ = LL.hasZ() ? LL.Z(): dummy;
        SquareMatrix<MatsU> & LSZ = LS.hasZ() ? LS.Z(): dummy;
        SquareMatrix<MatsU> & SLZ = SL.hasZ() ? SL.Z(): dummy;
        SquareMatrix<MatsU> & SSZ = SS.hasZ() ? SS.Z(): dummy;
        this->Z().componentScatter(LLZ, LSZ, SLZ, SSZ, increment);
      }

      if (this->hasXY()) {
        SquareMatrix<MatsU> & LLX = LL.hasXY() ? LL.X(): dummy;
        SquareMatrix<MatsU> & LLY = LL.hasXY() ? LL.Y(): dummy;
        SquareMatrix<MatsU> & LSX = LS.hasXY() ? LS.X(): dummy;
        SquareMatrix<MatsU> & LSY = LS.hasXY() ? LS.Y(): dummy;
        SquareMatrix<MatsU> & SLX = SL.hasXY() ? SL.X(): dummy;
        SquareMatrix<MatsU> & SLY = SL.hasXY() ? SL.Y(): dummy;
        SquareMatrix<MatsU> & SSX = SS.hasXY() ? SS.X(): dummy;
        SquareMatrix<MatsU> & SSY = SS.hasXY() ? SS.Y(): dummy;
        this->X().componentScatter(LLX, LSX, SLX, SSX, increment);
        this->Y().componentScatter(LLY, LSY, SLY, SSY, increment);
      }

  }

  template <typename MatsT>
  template <typename MatsU>
  void PauliSpinorSquareMatrices<MatsT>::componentGather(
    const PauliSpinorSquareMatrices<MatsU> & LL, 
    const PauliSpinorSquareMatrices<MatsU> & LS, 
    const PauliSpinorSquareMatrices<MatsU> & SL, 
    const PauliSpinorSquareMatrices<MatsU> & SS, 
    bool increment) {
    
      this->S().componentGather(LL.S(), LS.S(), SL.S(), SS.S(), increment); 
      
      SquareMatrix<MatsU> dummy(this->memManager(), 0);
      if (this->hasZ()) {
        const SquareMatrix<MatsU> & LLZ = LL.hasZ() ? LL.Z(): dummy;
        const SquareMatrix<MatsU> & LSZ = LS.hasZ() ? LS.Z(): dummy;
        const SquareMatrix<MatsU> & SLZ = SL.hasZ() ? SL.Z(): dummy;
        const SquareMatrix<MatsU> & SSZ = SS.hasZ() ? SS.Z(): dummy;
        this->Z().componentGather(LLZ, LSZ, SLZ, SSZ, increment);
      }
      
      if (this->hasXY()) {
        const SquareMatrix<MatsU> & LLX = LL.hasXY() ? LL.X(): dummy;
        const SquareMatrix<MatsU> & LLY = LL.hasXY() ? LL.Y(): dummy;
        const SquareMatrix<MatsU> & LSX = LS.hasXY() ? LS.X(): dummy;
        const SquareMatrix<MatsU> & LSY = LS.hasXY() ? LS.Y(): dummy;
        const SquareMatrix<MatsU> & SLX = SL.hasXY() ? SL.X(): dummy;
        const SquareMatrix<MatsU> & SLY = SL.hasXY() ? SL.Y(): dummy;
        const SquareMatrix<MatsU> & SSX = SS.hasXY() ? SS.X(): dummy;
        const SquareMatrix<MatsU> & SSY = SS.hasXY() ? SS.Y(): dummy;
        this->X().componentGather(LLX, LSX, SLX, SSX, increment);
        this->Y().componentGather(LLY, LSY, SLY, SSY, increment);
      }
  }
  
  template <typename MatsT>
  template <typename MatsU>
  PauliSpinorSquareMatrices<MatsT> 
  PauliSpinorSquareMatrices<MatsT>::componentGatherBuild(
    const PauliSpinorSquareMatrices<MatsU> & LL, 
    const PauliSpinorSquareMatrices<MatsU> & LS, 
    const PauliSpinorSquareMatrices<MatsU> & SL, 
    const PauliSpinorSquareMatrices<MatsU> & SS) {
      
      size_t N;
      if (LL.S().pointer()) N = LL.dimension();
      else if (LS.S().pointer()) N = LS.dimension();
      else if (SL.S().pointer()) N = SL.dimension();
      else if (SS.S().pointer()) N = SS.dimension();
      else CErr("Nothing to build in componentGatherBuild");
      
      bool any_hasZ  = LL.hasZ()  or LS.hasZ()  or SL.hasZ()  or SS.hasZ();
      bool any_hasXY = LL.hasXY() or LS.hasXY() or SL.hasXY() or SS.hasXY();
      
      PauliSpinorSquareMatrices<MatsT> pauli(LL.memManager(), 2*N, any_hasXY, any_hasZ);

      pauli.S() = SquareMatrix<MatsT>::componentGatherBuild(LL.S(), LS.S(), SL.S(), SS.S());
      
      if(any_hasZ)
        pauli.Z() = SquareMatrix<MatsT>::componentGatherBuild(LL.Z(), LS.Z(), SL.Z(), SS.Z());
      
      if(any_hasXY) {
        pauli.X() = SquareMatrix<MatsT>::componentGatherBuild(LL.X(), LS.X(), SL.X(), SS.X());
        pauli.Y() = SquareMatrix<MatsT>::componentGatherBuild(LL.Y(), LS.Y(), SL.Y(), SS.Y());
      }

      return pauli;
  }
  

  template void PauliSpinorSquareMatrices<double>::componentScatter(
    PauliSpinorSquareMatrices<double> & LL, PauliSpinorSquareMatrices<double> & LS, 
    PauliSpinorSquareMatrices<double> & SL, PauliSpinorSquareMatrices<double> & SS, bool increment) const;
  template void PauliSpinorSquareMatrices<double>::componentScatter(
    PauliSpinorSquareMatrices<dcomplex> & LL, PauliSpinorSquareMatrices<dcomplex> & LS, 
    PauliSpinorSquareMatrices<dcomplex> & SL, PauliSpinorSquareMatrices<dcomplex> & SS, bool increment) const;
  template void PauliSpinorSquareMatrices<dcomplex>::componentScatter(
    PauliSpinorSquareMatrices<dcomplex> & LL, PauliSpinorSquareMatrices<dcomplex> & LS, 
    PauliSpinorSquareMatrices<dcomplex> & SL, PauliSpinorSquareMatrices<dcomplex> & SS, bool increment) const;
  
  template void PauliSpinorSquareMatrices<double>::componentGather(
    const PauliSpinorSquareMatrices<double> & LL, const PauliSpinorSquareMatrices<double> & LS, 
    const PauliSpinorSquareMatrices<double> & SL, const PauliSpinorSquareMatrices<double> & SS, bool increment);
  template void PauliSpinorSquareMatrices<dcomplex>::componentGather(
    const PauliSpinorSquareMatrices<double> & LL, const PauliSpinorSquareMatrices<double> & LS, 
    const PauliSpinorSquareMatrices<double> & SL, const PauliSpinorSquareMatrices<double> & SS, bool increment);
  template void PauliSpinorSquareMatrices<dcomplex>::componentGather(
    const PauliSpinorSquareMatrices<dcomplex> & LL, const PauliSpinorSquareMatrices<dcomplex> & LS, 
    const PauliSpinorSquareMatrices<dcomplex> & SL, const PauliSpinorSquareMatrices<dcomplex> & SS, bool increment);

  template PauliSpinorSquareMatrices<double> PauliSpinorSquareMatrices<double>::componentGatherBuild(
    const PauliSpinorSquareMatrices<double> & LL, const PauliSpinorSquareMatrices<double> & LS, 
    const PauliSpinorSquareMatrices<double> & SL, const PauliSpinorSquareMatrices<double> & SS);
  template PauliSpinorSquareMatrices<dcomplex> PauliSpinorSquareMatrices<dcomplex>::componentGatherBuild(
    const PauliSpinorSquareMatrices<double> & LL, const PauliSpinorSquareMatrices<double> & LS, 
    const PauliSpinorSquareMatrices<double> & SL, const PauliSpinorSquareMatrices<double> & SS);
  template PauliSpinorSquareMatrices<dcomplex> PauliSpinorSquareMatrices<dcomplex>::componentGatherBuild(
    const PauliSpinorSquareMatrices<dcomplex> & LL, const PauliSpinorSquareMatrices<dcomplex> & LS, 
    const PauliSpinorSquareMatrices<dcomplex> & SL, const PauliSpinorSquareMatrices<dcomplex> & SS);
    
  template <typename MatsT>
  template <typename MatsU>
  void PauliSpinorSquareMatrices<MatsT>::componentAdd(
    const char TRANS, MatsU scale, const std::string & comp,
    const PauliSpinorSquareMatrices<MatsU> & pauli) {
    
     if (this->dimension() != pauli.dimension()*2) 
       CErr("Dimension mismatch in componentAdd");
     
     size_t N = pauli.dimension();
     size_t twoN = 2*N;
     size_t offset;
     if (not comp.compare("LL")) offset = 0;
     else if (not comp.compare("LS")) offset = N*twoN;
     else if (not comp.compare("SL")) offset = N;
     else if (not comp.compare("SS")) offset = N + N*twoN;
     else CErr("Unsupported component types in componentAdd"); 
     
     MatAdd(TRANS,'N',N,N,scale,pauli.S().pointer(),N,
       MatsT(1.),this->S().pointer()+offset,twoN,this->S().pointer()+offset,twoN); 
     
     if (this->hasZ() and pauli.hasZ()) {
       MatAdd(TRANS,'N',N,N,scale,pauli.Z().pointer(),N,
         MatsT(1.),this->Z().pointer()+offset,twoN,this->Z().pointer()+offset,twoN); 
     }
     
     if (this->hasXY() and pauli.hasXY()) {
       MatAdd(TRANS,'N',N,N,scale,pauli.X().pointer(),N,
         MatsT(1.),this->X().pointer()+offset,twoN,this->X().pointer()+offset,twoN); 
       MatAdd(TRANS,'N',N,N,scale,pauli.Y().pointer(),N,
         MatsT(1.),this->Y().pointer()+offset,twoN,this->Y().pointer()+offset,twoN); 
     }

  }
  
  template void PauliSpinorSquareMatrices<double>::componentAdd(const char TRANS, 
    double scale, const std::string & comp, const PauliSpinorSquareMatrices<double> & pauli);
  template void PauliSpinorSquareMatrices<dcomplex>::componentAdd(const char TRANS,
    double scale, const std::string & comp, const PauliSpinorSquareMatrices<double> & pauli);
  template void PauliSpinorSquareMatrices<dcomplex>::componentAdd(const char TRANS,
    dcomplex scale, const std::string & comp, const PauliSpinorSquareMatrices<dcomplex> & pauli);
 
  template <typename MatsT>
  void PauliSpinorSquareMatrices<MatsT>::symmetrizeLSSL(char TRANS, bool get_SL_from_LS) {
      
     size_t twoN = this->dimension();
     size_t N = twoN / 2;
     size_t P1, P2;

     if (get_SL_from_LS) {
       P1 = N*twoN; // LS
       P2 = N;      // SL
     } else { 
       P1 = N;      // SL
       P2 = N*twoN; // LS
     }
     
     SetMat(TRANS,N,N,MatsT(1.),this->S().pointer()+P1,twoN,this->S().pointer()+P2,twoN); 
     if(hasZ()) SetMat(TRANS,N,N,MatsT(1.),this->Z().pointer()+P1,twoN,this->Z().pointer()+P2,twoN); 
     if(hasXY()) {
       SetMat(TRANS,N,N,MatsT(1.),this->Y().pointer()+P1,twoN,this->Y().pointer()+P2,twoN);
       SetMat(TRANS,N,N,MatsT(1.),this->X().pointer()+P1,twoN,this->X().pointer()+P2,twoN);
     }
  }

  template void PauliSpinorSquareMatrices<double>::symmetrizeLSSL(char TRANS, bool get_SL_from_LS);
  template void PauliSpinorSquareMatrices<dcomplex>::symmetrizeLSSL(char TRANS, bool get_SL_from_LS);

}; // namespace ChronusQ
