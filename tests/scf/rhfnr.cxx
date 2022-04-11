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

#ifndef _CQ_GENERATE_TESTS
#include "scf.hpp"

// Water 6-31G(d) test
TEST( RHF_NR, Water_631Gd_NR ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_nr", "water_6-31Gd.bin.ref" );
 
};

#ifdef _CQ_DO_PARTESTS

// SMP Water 6-31G(d) test
TEST( RHF_NR, PAR_Water_631Gd_NR ) {

  CQSCFTEST( "scf/parallel/rhf/water_6-31Gd_nr", "water_6-31Gd.bin.ref" );
 
};

#endif

#endif
