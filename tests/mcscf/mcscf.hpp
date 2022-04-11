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
//#include <util/mpi.hpp>

// Directory containing reference files
#define MCSCF_TEST_REF TEST_ROOT "/mcscf/reference/"

using namespace ChronusQ;


inline void CQMCSCFTEST( std::string in, std::string ref,
  double tol        = 1e-8,
  bool checkOsc     = false,
  bool checkEne     = true,
  std::string scr   = "") {


  std::string inputSCR;
  if (scr.compare("")) inputSCR = MCSCF_TEST_REF + scr + ".scr";

#ifdef _CQ_GENERATE_TESTS

//  MPI_Barrier(MPI_COMM_WORLD);

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    MCSCF_TEST_REF + ref, inputSCR);

//  MPI_Barrier(MPI_COMM_WORLD);
//  if(MPIRank(MPI_COMM_WORLD) != 0) return;

#else

//  MPI_Barrier(MPI_COMM_WORLD);
  
  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin", inputSCR);

//  MPI_Barrier(MPI_COMM_WORLD);
//  if(MPIRank(MPI_COMM_WORLD) != 0) return;

  SafeFile refFile(MCSCF_TEST_REF + ref,  true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
  
  double xNS, yNS;
  std::cout << " * PERFORMING MSWFN ENERGY CHECK " << std::endl;
  std::cout << "MCSCF_TEST_REF=" << MCSCF_TEST_REF <<std::endl;
  
  refFile.readData("MCWFN/NSTATES", &xNS);
  resFile.readData("MCWFN/NSTATES", &yNS);

  EXPECT_NEAR( xNS, yNS, tol ) << "NUMBER OF STATES TEST FAILED ";

  std::vector<double> xStateEnergy(xNS);
  std::vector<double> yStateEnergy(yNS);
//  std::array<double,3> xDummy3, yDummy3; 
//  std::array<std::array<double,3>,3> xDummy33, yDummy33;
//  std::array<std::array<std::array<double,3>,3>,3> xDummy333, yDummy333;
  
  /* Check Energy */ 
  if( checkEne ) {

    refFile.readData("MCWFN/STATE_ENERGY", &xStateEnergy[0]);
    resFile.readData("MCWFN/STATE_ENERGY", &yStateEnergy[0]);

    for(auto i = 0; i < xNS; i++)
      EXPECT_NEAR(xStateEnergy[i], yStateEnergy[i], tol ) << 
        "ENERGY TEST FAILED ISTATE = " << i; 
  }
  
  /* Check Oscillator Strength */
  if( checkOsc ) {

    std::cout << " * PERFORMING MCSCF OSCILLATOR STRENGTH CHECK " << std::endl;

    auto oscDim     = resFile.getDims("MCWFN/OSC_STR");
    ASSERT_EQ(oscDim.size(),2);
    std::vector<double> xDummy, yDummy;
    xDummy.resize(oscDim[0] * oscDim[1]);
    yDummy.resize(oscDim[0] * oscDim[1]);
    refFile.readData("MCWFN/OSC_STR", &xDummy[0]);
    resFile.readData("MCWFN/OSC_STR", &yDummy[0]);

    for(auto i = 0; i < oscDim[0]; i++)
    for(auto j = 0; j < oscDim[1]; j++)
      EXPECT_NEAR(yDummy[i*oscDim[1]+j], xDummy[i*oscDim[1]+j], tol ) <<
        "OSCILLATOR STRENGTH TEST FAILED STATE 1 = " << i
                                     << "STATE 2 = " << j;

  }

  /* Check Multipoles */ 

//  if( checkDipLen ) {
//
//    std::cout << " * PERFORMING MCSCF DIPOLE (LEN) CHECK " << std::endl;
//
//    refFile.readData("MCSCF/LEN_ELECTRIC_DIPOLE",&xDummy3[0]);
//    resFile.readData("SCF/LEN_ELECTRIC_DIPOLE",&yDummy3[0]);
//    for(auto i = 0; i < 3; i++)
//      EXPECT_NEAR(yDummy3[i], xDummy3[i], tol ) << 
//        "DIPOLE TEST FAILED IXYZ = " << i; 
//
//
//  }
//  
//
//  if( checkQuadLen ) {
//
//    std::cout << " * PERFORMING SCF QUADRUPOLE (LEN) CHECK " << std::endl;
//
//    refFile.readData("SCF/LEN_ELECTRIC_QUADRUPOLE",&xDummy33[0][0]);
//    resFile.readData("SCF/LEN_ELECTRIC_QUADRUPOLE",&yDummy33[0][0]);
//    for(auto i = 0; i < 3; i++)
//    for(auto j = 0; j < 3; j++)
//      EXPECT_NEAR(yDummy33[i][j], xDummy33[i][j],  tol) << 
//        "QUADRUPOLE TEST FAILED IXYZ = " << i 
//                           << " JXYZ = " << j; 
//
//  }
//  
//  if( checkOctLen ) {
//
//    std::cout << " * PERFORMING SCF OCTUPOLE (LEN) CHECK " << std::endl;
//
//    refFile.readData("SCF/LEN_ELECTRIC_OCTUPOLE",&xDummy333[0][0][0]);
//    resFile.readData("SCF/LEN_ELECTRIC_OCTUPOLE",&yDummy333[0][0][0]);
//    for(auto i = 0; i < 3; i++)
//    for(auto j = 0; j < 3; j++)
//    for(auto k = 0; k < 3; k++)
//      EXPECT_NEAR(yDummy333[i][j][k],  xDummy333[i][j][k],  tol) <<
//        "OCTUPOLE TEST FAILED IXYZ = " << i 
//                                       << " JXYZ = " << j 
//                                       << " KXYZ = " << k; 
//
//  }
//  
//  /* Check Spin */
//
//  if( checkSEXP ) {
//
//    std::cout << " * PERFORMING SCF <S> CHECK " << std::endl;
//
//    refFile.readData("SCF/S_EXPECT",&xDummy3[0]);
//    resFile.readData("SCF/S_EXPECT",&yDummy3[0]);
//    for(auto i = 0; i < 3; i++)
//      EXPECT_NEAR(yDummy3[i], xDummy3[i], tol ) << 
//        "<S> TEST FAILED IXYZ = " << i; 
//
//  }
//  
//  if( checkSSq ) {
//
//    std::cout << " * PERFORMING SCF <S^2> CHECK " << std::endl;
//
//    refFile.readData("SCF/S_SQUARED",&xDummy);
//    resFile.readData("SCF/S_SQUARED",&yDummy);
//    EXPECT_NEAR(yDummy,  xDummy, tol) << "<S^2> TEST FAILED " ;
//
//  }

#endif

}

inline void CQMCSCFREFTEST( std::string in, std::string ref,
  double tol        = 1e-8,
  bool checkOsc     = false,
  bool checkEne     = true,
  std::string scr   = "" ) {

 
  std::string inputSCR;
  if (scr.compare("")) inputSCR = MCSCF_TEST_REF + scr + ".scr";

#ifdef _CQ_GENERATE_TESTS

//  MPI_Barrier(MPI_COMM_WORLD);

  // Assumes user added bin to RESP_TEST_REF
  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MCSCF_TEST_REF + ref, inputSCR);

//  MPI_Barrier(MPI_COMM_WORLD);
//  if(MPIRank(MPI_COMM_WORLD) != 0) return;

#else

//  MPI_Barrier(MPI_COMM_WORLD);

  std::ifstream  src(MCSCF_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(TEST_OUT + in + ".bin", std::ios::binary);
  dst << src.rdbuf();
  dst.flush();
  
  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin", inputSCR);

//  MPI_Barrier(MPI_COMM_WORLD);
//  if(MPIRank(MPI_COMM_WORLD) != 0) return;

  SafeFile refFile(MCSCF_TEST_REF + ref,  true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);

  double xNS, yNS;
  std::cout << " * PERFORMING MSWFN ENERGY CHECK " << std::endl;
  std::cout << "MCSCF_TEST_REF=" << MCSCF_TEST_REF <<std::endl;

  refFile.readData("MCWFN/NSTATES", &xNS);
  resFile.readData("MCWFN/NSTATES", &yNS);

  EXPECT_NEAR( xNS, yNS, tol ) << "NUMBER OF STATES TEST FAILED ";

  std::vector<double> xStateEnergy(xNS);
  std::vector<double> yStateEnergy(yNS);

  /* Check Energy */
  if( checkEne ) {

    refFile.readData("MCWFN/STATE_ENERGY", &xStateEnergy[0]);
    resFile.readData("MCWFN/STATE_ENERGY", &yStateEnergy[0]);

    for(auto i = 0; i < xNS; i++)
      EXPECT_NEAR(xStateEnergy[i], yStateEnergy[i], tol ) <<
        "ENERGY TEST FAILED ISTATE = " << i;
  }

  /* Check Oscillator Strength */
  if( checkOsc ) {

    std::cout << " * PERFORMING MCSCF OSCILLATOR STRENGTH CHECK " << std::endl;

    auto oscDim     = resFile.getDims("MCWFN/OSC_STR");
    ASSERT_EQ(oscDim.size(),2);
    std::vector<double> xDummy, yDummy;
    xDummy.resize(oscDim[0] * oscDim[1]);
    yDummy.resize(oscDim[0] * oscDim[1]);
    refFile.readData("MCWFN/OSC_STR", &xDummy[0]);
    resFile.readData("MCWFN/OSC_STR", &yDummy[0]);

    for(auto i = 0; i < oscDim[0]; i++)
    for(auto j = 0; j < oscDim[1]; j++)
      EXPECT_NEAR(yDummy[i*oscDim[1]+j], xDummy[i*oscDim[1]+j], tol ) <<
        "OSCILLATOR STRENGTH TEST FAILED STATE 1 = " << i
                                     << "STATE 2 = " << j;

  }




#endif

}





