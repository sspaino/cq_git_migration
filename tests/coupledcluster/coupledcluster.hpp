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

#include <ut.hpp>

#include <cxxapi/procedural.hpp>
#include <util/files.hpp>
#include <util/mpi.hpp>

#define CC_TEST_REF TEST_ROOT "/coupledcluster/reference/"

using namespace ChronusQ;



static void CQCCTEST( std::string in, std::string ref ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    CC_TEST_REF + ref, TEST_OUT + in + ".scr");

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin",TEST_OUT + in + ".scr");


  SafeFile refFile(CC_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);

  double xDummy, yDummy;

  resFile.readData("/CC/CORRELATION_ENERGY",&xDummy);
  refFile.readData("/CC/CORRELATION_ENERGY",&yDummy);

  EXPECT_NEAR( xDummy, yDummy, 1e-8 ) <<
    "CORRELATION ENERGY TEST FAILED "; 

  double xDummy1, yDummy1;

  resFile.readData("/CC/T_NORM",&xDummy1);
  refFile.readData("/CC/T_NORM",&yDummy1);

  EXPECT_NEAR( xDummy1, yDummy1, 1e-6 ) <<
    "CORRELATION ENERGY NORM TEST FAILED "; 

#endif


};
