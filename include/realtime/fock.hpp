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

#include <realtime.hpp>


namespace ChronusQ {

/*
  void operator<<(std::array<double,3> &lhs, std::valarray<double> &rhs){
    assert(rhs.size() == 3);
    std::copy_n(&rhs[0],3,lhs.begin());
  }
*/

  template <size_t N, typename T>
  std::array<T,N> valarray2array(const std::valarray<T> &x) {
    assert( x.size() == N );

    std::array<T,N> arr;
    std::copy_n(&x[0],N,arr.begin());
    return arr;
   
  };

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::formFock(bool increment, double t, size_t idx) {

    ProgramTimer::timeOp("Form Fock", [&]() {

      // Get perturbation for the current time and build a Fock matrix
      EMPerturbation pert_t = pert.getPert(t);

      // Add on the SCF Perturbation
      if ( intScheme.includeSCFField )
        for( auto& field : scfPert.fields )
          pert_t.addField( field );

      systems_[idx]->formFock(pert_t,increment);

    });

  };

}; // namespace ChronusQ


