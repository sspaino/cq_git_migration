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

# Directories for RT tests
set(MCSCF_TEST_SOURCE_ROOT "${TEST_ROOT}/mcscf" )
set(MCSCF_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/mcscf" )
  

# Set up compilation of MCSCF test exe
add_executable(mcscftest ../ut.cxx ras.cxx cas.cxx)

target_include_directories(mcscftest PUBLIC ${MCSCF_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(mcscftest PUBLIC ChronusQ::Dependencies ${CQEX_LINK})
target_link_libraries(mcscftest PUBLIC ChronusQ::gtest)

if(CQEX_DEP)
  add_dependencies(mcscftest ${CQEX_DEP})
endif()

# Generate directories
file(MAKE_DIRECTORY ${MCSCF_TEST_BINARY_ROOT}/serial/cas)
file(MAKE_DIRECTORY ${MCSCF_TEST_BINARY_ROOT}/serial/ras)

file(MAKE_DIRECTORY ${MCSCF_TEST_BINARY_ROOT}/parallel/cas)
file(MAKE_DIRECTORY ${MCSCF_TEST_BINARY_ROOT}/parallel/ras)


include( CQTestGeneration )

# CAS test
add_cq_test(OneC_CASSCF_FULLMATRIX  mcscftest "OneC_CASSCF_FULLMATRIX.*")
add_cq_test(X2C_CASSCF_FULLMATRIX   mcscftest "X2C_CASSCF_FULLMATRIX.*")
add_cq_test(FourC_CASSCF_FULLMATRIX mcscftest "FourC_CASSCF_FULLMATRIX.*")
add_cq_test(CASCI_READMO_SKIPSCF mcscftest "CASCI_READMO_SKIPSCF.*")
add_cq_test(CASCI_DAVIDSON       mcscftest "CASCI_DAVISON.*")
add_cq_test(OneC_CAS_SWAP        mcscftest "OneC_CAS_SWAP.*")
add_cq_test(TwoC_CAS_SWAP        mcscftest "TwoC_CAS_SWAP.*")
add_cq_test(GHF_CAS_OSC          mcscftest "GHF_CAS_OSC.*")

# RAS test
add_cq_test(RAS_DAVIDSON         mcscftest "RAS_DAVIDSON.*")

if( CQ_ENABLE_MPI )

endif()

