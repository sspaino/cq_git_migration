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

#include "scf.hpp"


// O2 6-31G(d) test
TEST( UHF, O2_631Gd ) {

  CQSCFTEST( "scf/serial/uhf/oxygen_6-31Gd", "oxygen_6-31Gd.bin.ref" );

};

// MnHe sto-3g test
TEST( UHF, MnHe_sto3G ) {

  CQSCFTEST( "scf/serial/uhf/MnHe_sto-3G", "MnHe_sto-3G.bin.ref" );

};

// RbHe sto-3g test
TEST( UHF, RbHe_sto3G ) {

  CQSCFTEST( "scf/serial/uhf/RbHe_sto-3G", "RbHe_sto-3G.bin.ref" );

};

// KCaKrK sto-3g test for fchk parsing
TEST( UHF, KCaKrK_sto3G ) {

  CQSCFTEST( "scf/serial/uhf/KCaKrK_sto-3G", "KCaKrK_sto-3G.bin.ref", 1e-8,
              true, true, true, true, true, true,
              false, "KCaKrK_sto-3G.fchk");

};

// H atom 3-21g post-scf test
TEST( UHF, H_321G ) {

  CQSCFTEST( "scf/serial/uhf/H_3-21G", "H_3-21G.bin.ref", 1e-7,
              true, true, true, true, true, true);

};

// Ne(Z=10.5)He(Z=1.5)H(Z=0.5) sto-3g test
TEST( UHF, NeHeH_fracZ_sto3G ) {

  CQSCFTEST( "scf/serial/uhf/NeHeH_fracZ_uhf_sto-3G", "NeHeH_fracZ_uhf_sto-3G.bin.ref" );

};

// B sto-3g test for MO swapping
TEST( UHF, B_swap_UHF_sto3G ) {

  CQSCFTEST( "scf/serial/uhf/B_swap_UHF_sto-3g", "B_swap_UHF_sto-3g.bin.ref",
    1e-8, true, true, true, true, true, true,
    true );

};

#ifdef _CQ_DO_PARTESTS

// SMP O2 6-31G(d) test
TEST( UHF, PAR_O2_631Gd ) {

  CQSCFTEST( "scf/parallel/uhf/oxygen_6-31Gd", "oxygen_6-31Gd.bin.ref" );

};

#endif




