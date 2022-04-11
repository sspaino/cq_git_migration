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

#include "rt.hpp"



// Water 6-31G(d) Delta Spike (along Y)
TEST( RHF_RT, Water_631Gd_Delta_Y ) {

  CQRTTEST( rt/serial/rrt/water_6-31Gd_rhf_delta_y,
    water_6-31Gd_rhf_delta_y.bin.ref );

}

// Magnus 2 delta electric field
TEST( RHF_RT, Water_631Gd_Magnus2 ) {

  CQRTTEST( rt/serial/rrt/water_6-31Gd_rhf_magnus2,
    water_6-31Gd_rhf_magnus2.bin.ref );

}

// MMUT w/ Magnus 2 restart non-delta electric field
TEST( RHF_RT, Water_631Gd_MMUT_Magnus2 ) {

  CQRTTEST( rt/serial/rrt/water_6-31Gd_rhf_mmut_magnus2,
    water_6-31Gd_rhf_mmut_magnus2.bin.ref );

}

#ifdef _CQ_DO_PARTESTS

// SMP Water 6-31G(d) Delta Spike (along Y)
TEST( RHF_RT, PAR_Water_631Gd_Delta_Y ) {

  CQRTTEST( rt/parallel/rrt/water_6-31Gd_rhf_delta_y,
    water_6-31Gd_rhf_delta_y.bin.ref );

}

#endif














// Water 6-31G(d) B3LYP Delta Spike (along Y)
TEST( RKS_RT, Water_631Gd_B3LYP_Delta_Y ) {

  CQRTTEST( rt/serial/rrt/water_6-31Gd_rb3lyp_delta_y,
    water_6-31Gd_rb3lyp_delta_y.bin.ref );

}


// MMUT w/ Magnus 2 restart non-delta electric field
TEST( RKS_RT, Water_631Gd_MMUT_Magnus2 ) {

  CQRTTEST( rt/serial/rrt/water_6-31Gd_rb3lyp_mmut_magnus2,
    water_6-31Gd_rb3lyp_mmut_magnus2.bin.ref );

}


#ifdef _CQ_DO_PARTESTS

// SMP Water 6-31G(d) B3LYP Delta Spike (along Y)
TEST( RKS_RT, PAR_Water_631Gd_B3LYP_Delta_Y ) {

  CQRTTEST( rt/parallel/rrt/water_6-31Gd_rb3lyp_delta_y,
    water_6-31Gd_rb3lyp_delta_y.bin.ref );

}

#endif


//
// Restart testing
//

TEST( RESTART_RT, Restart_Water_631Gd_B3LYP_Delta_Y ) {

  CQRTRESTARTTEST( rt/serial/rrt/water_6-31Gd_rb3lyp_delta_y_restart_mid,
    water_6-31Gd_rb3lyp_delta_y_restart_mid.bin,
    rt/serial/rrt/water_6-31Gd_rb3lyp_delta_y_restart,
    water_6-31Gd_rb3lyp_delta_y_restart.bin.ref );

}

#ifdef _CQ_DO_PARTESTS

TEST( RESTART_RT, PAR_Restart_Water_631Gd_B3LYP_Delta_Y ) {

  CQRTRESTARTTEST( rt/parallel/rrt/water_6-31Gd_rb3lyp_delta_y_restart_mid,
    water_6-31Gd_rb3lyp_delta_y_restart_mid.bin,
    rt/parallel/rrt/water_6-31Gd_rb3lyp_delta_y_restart,
    water_6-31Gd_rb3lyp_delta_y_restart.bin.ref );

}

#endif
