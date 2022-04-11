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

// Dummy macros
#define DUMMY(X)
#define DUMMY2(X,Y)
#define DUMMY3(X,Y,Z)


// Macros for the deallocation of raw memory using the 
// CQMemManager
#define DEALLOC_OP(mem,PTR) if(PTR != nullptr) mem.free(PTR);
#define DEALLOC_VEC_OP(mem, VEC_PTR) \
  for(auto i = 0; i < VEC_PTR.size(); i++) \
    if(VEC_PTR[i] != nullptr) DEALLOC_OP(mem,VEC_PTR[i]); \
  VEC_PTR.clear();

// Some common wrappers around the deallocation macros for
// consistant calling signatures

// For consistant calling signiture with MOVE and COPY macros
#define DEALLOC_OP_5(X,Y,Z,mem,PTR) DEALLOC_OP(mem,PTR);
#define DEALLOC_VEC_OP_5(X,Y,Z,mem,PTR) DEALLOC_VEC_OP(mem,PTR);


// COPY preprocessor macros
#define COPY_OTHER_MEMBER(this,other,X) this->X = other.X;
#define COPY_OTHER_MEMBER_OP(T,this,other,mem,PTR) \
  if(other.PTR != nullptr) { \
    size_t OPSZ = mem.getSize(other.PTR); \
    this->PTR    = mem.template malloc<T>(OPSZ); \
    std::copy_n(other.PTR, OPSZ, this->PTR); \
  } else this->PTR = nullptr;

#define COPY_OTHER_MEMBER_VEC_OP(T,this,other,mem,VEC_PTR) \
  this->VEC_PTR.clear(); \
  for(auto i = 0; i < other.VEC_PTR.size(); i++) { \
    this->VEC_PTR.push_back(nullptr); \
    COPY_OTHER_MEMBER_OP(T,this,other,mem,VEC_PTR[i]); \
  }


// MOVE preprocessor macros
#define MOVE_OTHER_MEMBER(this,other,X) this->X = std::move(other.X);
#define MOVE_OTHER_MEMBER_OP(T,this,other,mem,PTR) \
  if(other.PTR != nullptr) { \
    size_t OPSZ = mem.getSize(other.PTR); \
    this->PTR    = mem.template malloc<T>(OPSZ); \
    std::copy_n(other.PTR, OPSZ, this->PTR); \
    DEALLOC_OP(mem,other.PTR); \
  } else this->PTR = nullptr;

#define MOVE_OTHER_MEMBER_VEC_OP(T,this,other,mem,VEC_PTR) \
  this->VEC_PTR.clear(); \
  for(auto i = 0; i < other.VEC_PTR.size(); i++) { \
    this->VEC_PTR.push_back(nullptr); \
    MOVE_OTHER_MEMBER_OP(T,this,other,mem,VEC_PTR[i]); \
  } \
  other.VEC_PTR.clear();

