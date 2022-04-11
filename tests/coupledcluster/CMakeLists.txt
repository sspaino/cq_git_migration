#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2022 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#
if(CQ_HAS_TA)
# Directories for Coupled-Cluster tests
set(CC_TEST_SOURCE_ROOT "${TEST_ROOT}/coupledcluster" )
set(CC_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/coupledcluster" )


# Set up compilation of Coupled-Cluster test exe
add_executable(cctest ../ut.cxx x2c_cc.cxx)


target_include_directories(cctest PUBLIC ${CC_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(cctest PUBLIC ChronusQ::Dependencies ${CQEX_LINK})
target_link_libraries(cctest PUBLIC ChronusQ::gtest)


if(CQEX_DEP)
  add_dependencies(cctest ${CQEX_DEP})
endif()


file(MAKE_DIRECTORY ${CC_TEST_BINARY_ROOT}/ccsd)

include( CQTestGeneration )


add_cq_test( X2C_CCSD cctest "X2C_CCSD*" )
add_cq_test( GHF_CCSD cctest "GHF_CCSD*" )
endif()

