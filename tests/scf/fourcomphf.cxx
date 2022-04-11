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

/* 
 * Note: Currently only SCF energies are tested. Properties need
 *       to be implemented.
 */

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, UPu_184_plus_P_NR_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_NR_pointnuc",
    "UPu_184+_P_NR_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, UPu_184_plus_P_DC_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_pointnuc",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

/*
 * FIXME: Different compilers/LA libraries give different solutions to this
 *        test
// Two electron U-Pu 184+ test Dirac-Couloumb with X2CHF guess (readden)
TEST( FOURCHF, UPu_184_plus_P_DC_X2CGuess_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_X2CGuess_pointnuc",
    "UPu_184+_P_DC_X2CGuess_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true, false, "UPu_184+_X2C_pointnuc.scr.bin" );

};
*/

#ifndef _CQ_GENERATE_TESTS
// Two electron U-Pu 184+ test Dirac-Couloumb with readmo
TEST( FOURCHF, UPu_184_plus_P_DC_pointnuc_readmo ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_pointnuc_readmo",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true, true);

};
#endif

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, UPu_184_plus_P_DCG_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCG_pointnuc",
    "UPu_184+_P_DCG_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, UPu_184_plus_P_NR_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_NR_finitenuc",
    "UPu_184+_P_NR_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, UPu_184_plus_P_DC_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_finitenuc",
    "UPu_184+_P_DC_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb AO Direct
TEST( FOURCHF, UPu_184_plus_P_DC_finitenuc_direct ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_finitenuc_direct",
    "UPu_184+_P_DC_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, UPu_184_plus_P_DCG_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCG_finitenuc",
    "UPu_184+_P_DCG_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt-Gauge-SSSS AO Direct
TEST( FOURCHF, UPu_184_plus_P_DCGGS_finitenuc_direct ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCGGS_finitenuc_direct",
    "UPu_184+_P_DCGGS_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt-Gauge-SSSS
TEST( FOURCHF, UPu_184_plus_P_DCGGS_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCGGS_finitenuc",
    "UPu_184+_P_DCGGS_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Distorted CuH3 to test C1 symmetry Dirac-Coulomb-Gaunt
TEST( FOURCHF, CuH3_321g_DCG ) {

  CQSCFTEST( "scf/serial/fourcomp/CuH3_321g_DCG",
    "CuH3_321g_DCG.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Ag Neutral doublet Dirac-Coulomb-Gaunt
TEST( FOURCHF, Ag_sapporoDZ_DCG ) {

  CQSCFTEST( "scf/serial/fourcomp/Ag_sapporoDZ_DCG",
    "Ag_sapporoDZ_DCG.bin.ref",1e-8,
    false, false, false, false, false, true);

};



#ifdef _CQ_DO_PARTESTS

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, PAR_UPu_184_plus_P_NR_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_NR_pointnuc",
    "UPu_184+_P_NR_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_pointnuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_pointnuc",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, PAR_UPu_184_plus_P_DCG_pointnuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCG_pointnuc",
    "UPu_184+_P_DCG_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#ifndef _CQ_GENERATE_TESTS
// Two electron U-Pu 184+ test Dirac-Couloumb with readmo
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_pointnuc_readmo ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_pointnuc_readmo",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true, true);

};
#endif

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, PAR_UPu_184_plus_P_NR_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_NR_finitenuc",
    "UPu_184+_P_NR_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_finitenuc",
    "UPu_184+_P_DC_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb AO Direct
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_finitenuc_direct ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_finitenuc_direct",
    "UPu_184+_P_DC_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, PAR_UPu_184_plus_P_DCG_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCG_finitenuc",
    "UPu_184+_P_DCG_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt-Gauge-SSSS AO Direct
TEST( FOURCHF, PAR_UPu_184_plus_P_DCGGS_finitenuc_direct ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCGGS_finitenuc_direct",
    "UPu_184+_P_DCGGS_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt-Gauge-SSSS
TEST( FOURCHF, PAR_UPu_184_plus_P_DCGGS_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCGGS_finitenuc",
    "UPu_184+_P_DCGGS_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#endif


