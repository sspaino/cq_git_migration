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



// Water 6-311+G(d,p) X2CDelta Spike (along Y)
TEST( X2CHF_RT, Water_6311pGdp_X2C_Delta_Y ) {

  CQRTTEST( rt/serial/grt/water_6-311pGdp_x2c_delta_y,
    water_6-311pGdp_x2c_delta_y.bin.ref );

}

// Water 6-311+G(d,p) X2CDelta Spike (along Y) w/magnus2 restart
TEST( X2CHF_RT, Water_6311pGdp_X2C_MMUT_Magnus2 ) {

  CQRTTEST( rt/serial/grt/water_6-311pGdp_x2c_mmut_magnus2,
    water_6-311pGdp_x2c_mmut_magnus2.bin.ref );

}

#ifdef _CQ_DO_PARTESTS

// SMP Water 6-311+G(d,p) X2CDelta Spike (along Y)
TEST( X2CHF_RT, PAR_Water_6311pGdp_X2C_Delta_Y ) {

  CQRTTEST( rt/parallel/grt/water_6-311pGdp_x2c_delta_y,
    water_6-311pGdp_x2c_delta_y.bin.ref );

}

#endif



TEST( RESTART_RT, Restart_Water_6311pGdp_X2C_Delta_Y ) {

  CQRTRESTARTTEST( rt/serial/grt/water_6-311pGdp_x2c_delta_y_restart_mid,
    water_6-311pGdp_x2c_delta_y_restart_mid.bin,
    rt/serial/grt/water_6-311pGdp_x2c_delta_y_restart,
    water_6-311pGdp_x2c_delta_y_restart.bin.ref );

}

