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



// Oxygen 6-31G(d) Delta Spike (along Y)
TEST( UHF_RT, O2_631Gd_Delta_Y ) {

  CQRTTEST( rt/serial/urt/oxygen_6-31Gd_uhf_delta_y,
    oxygen_6-31Gd_uhf_delta_y.bin.ref );

}

// Magnus 2 delta electric field
TEST( UHF_RT, O2_631Gd_Magnus2 ) {

  CQRTTEST( rt/serial/urt/oxygen_6-31Gd_uhf_magnus2,
    oxygen_6-31Gd_uhf_magnus2.bin.ref );

}

// MMUT w/ Magnus 2 restart non-delta electric field
TEST( UHF_RT, O2_631Gd_MMUT_Magnus2 ) {

  CQRTTEST( rt/serial/urt/oxygen_6-31Gd_uhf_mmut_magnus2,
    oxygen_6-31Gd_uhf_mmut_magnus2.bin.ref );

}


#ifdef _CQ_DO_PARTESTS

// SMP Oxygen 6-31G(d) Delta Spike (along Y)
TEST( UHF_RT, PAR_O2_631Gd_Delta_Y ) {

  CQRTTEST( rt/parallel/urt/oxygen_6-31Gd_uhf_delta_y,
    oxygen_6-31Gd_uhf_delta_y.bin.ref );

}

#endif









// Oxygen 6-31G(d) B3LYP Delta Spike (along Y)
TEST( UKS_RT, O2_631Gd_B3LYP_Delta_Y ) {

  CQRTTEST( rt/serial/urt/oxygen_6-31Gd_ub3lyp_delta_y,
    oxygen_6-31Gd_ub3lyp_delta_y.bin.ref );

}

// MMUT w/ Magnus 2 restart non-delta electric field
TEST( UKS_RT, O2_631Gd_MMUT_Magnus2 ) {

  CQRTTEST( rt/serial/urt/oxygen_6-31Gd_ub3lyp_mmut_magnus2,
    oxygen_6-31Gd_ub3lyp_mmut_magnus2.bin.ref );

}

#ifdef _CQ_DO_PARTESTS

// SMP Oxygen 6-31G(d) B3LYP Delta Spike (along Y)
TEST( UKS_RT, PAR_O2_631Gd_B3LYP_Delta_Y ) {

  CQRTTEST( rt/parallel/urt/oxygen_6-31Gd_ub3lyp_delta_y,
    oxygen_6-31Gd_ub3lyp_delta_y.bin.ref );

}

#endif


TEST( RESTART_RT, Restart_O2_631Gd_B3LYP_Delta_Y ) {

  CQRTRESTARTTEST( rt/serial/urt/oxygen_6-31Gd_ub3lyp_delta_y_restart_mid,
    oxygen_6-31Gd_ub3lyp_delta_y_restart_mid.bin,
    rt/serial/urt/oxygen_6-31Gd_ub3lyp_delta_y_restart,
    oxygen_6-31Gd_ub3lyp_delta_y_restart.bin.ref );

}

