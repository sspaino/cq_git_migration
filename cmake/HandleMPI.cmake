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

include(ExternalProject)

if(CQ_ENABLE_MPI)

  message( "" )
  
  # FindMPI
  find_package(MPI REQUIRED)
  target_link_libraries( ChronusQ::Dependencies INTERFACE MPI::MPI_CXX )
  copy_header_properties( MPI::MPI_CXX ChronusQ::DepHeaders )
  
  message( "" )
  
  # Print out extraneous information
  message( STATUS "MPIEXEC found to be: ${MPIEXEC}" )
  message( STATUS "MPIEXEC_NUMPROC_FLAG found to be: ${MPIEXEC_NUMPROC_FLAG}" )
  message( STATUS "MPI_INCLUDE_PATH found to be: ${MPI_INCLUDE_PATH}" )
  
  
  message( "" )
  
  
  # MXX
  
  message( STATUS "Adding CMake Target for MXX" )
  set( MXX_PREFIX     ${PROJECT_SOURCE_DIR}/external/mxx )
  set( MXX_INCLUDEDIR ${MXX_PREFIX}/src/mxx/include )
  
  ExternalProject_Add(mxx
    PREFIX ${MXX_PREFIX}
    GIT_REPOSITORY https://github.com/wavefunction91/mxx.git
    CONFIGURE_COMMAND echo 'No MXX Configure'
    UPDATE_COMMAND echo 'No ScaLAPACK MXX Command'
    BUILD_COMMAND echo 'No MXX Build'
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND echo 'No MXX Install'
  )
  
  list(APPEND CQEX_DEP mxx)
  
  # MXX Includes
  file( MAKE_DIRECTORY ${MXX_INCLUDEDIR} )
  set_property( TARGET ChronusQ::Dependencies APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${MXX_INCLUDEDIR}
  )
  set_property( TARGET ChronusQ::DepHeaders APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${MXX_INCLUDEDIR}
  )

endif()
