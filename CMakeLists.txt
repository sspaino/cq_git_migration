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

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)


project(ChronusQ
    VERSION 0.6.0
    DESCRIPTION
    "Chronus Quantum is an ab initio computational chemistry \
software package that focuses on explicitly time-dependent and unconventional \
post-SCF quantum mechanical methods."
    LANGUAGES CXX C Fortran)


# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CQHelperFunctions)

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

option(BLAS_EXTERNAL       "Search for external BLAS"        ON)
option(CQ_ENABLE_MPI       "Enable MPI parallelism"          OFF)
option(ENABLE_COVERAGE     "Enable coverage and profiling"   OFF)
option(CQ_LINALG_USESYSTEM "Use System defaults for LA libs" OFF)
option(CQ_EXTERNAL_OPENMP  "Force linking to an external OpenMP library" OFF)
option(CQ_ENABLE_TA  "Enable TA"  OFF)
option(CQ_EXTERNAL_IN_BUILD_DIR "Build external dependencies in the build directory" OFF)


# Required package building options
#   These set the variables CQ_BUILD_${PKGNAME}_TYPE as ALLOW, FORCE, or NONE
#   e.g.the option for libint is CQ_BUILD_LIBINT_TYPE
build_pkg_option( LIBINT ALLOW )
build_pkg_option( LIBCINT ALLOW )

# -----------------------------------------------------------------------------


# Move FetchContent base to either external in project root or build root
#   Can be overridden by specifying FETCHCONTENT_BASE_DIR
#   CMake's default for this is _deps in build root
if ( NOT FETCHCONTENT_BASE_DIR AND NOT CQ_EXTERNAL_IN_BUILD_DIR )
  set ( FETCHCONTENT_BASE_DIR "${PROJECT_SOURCE_DIR}/external"
        CACHE STRING "" FORCE )
elseif( NOT FETCHCONTENT_BASE_DIR )
  set ( FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/external"
        CACHE STRING "" FORCE )
endif()

# Add specified paths to CMAKE_PREFIX_PATH
if( MPI_PREFIX )
  set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${MPI_PREFIX} )
endif()

# C++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)


# Dummy target to collect general libraries to link
add_library( ChronusQ::Dependencies INTERFACE IMPORTED )
# Dummy target to collect header/compile definitions to link
add_library( ChronusQ::DepHeaders   INTERFACE IMPORTED )


# Parallelism

message( "\n\n" )
message( "ChronusQ Parallelism Settings" )

# OpenMP
message( "" )
find_package(OpenMP REQUIRED)
target_link_libraries( ChronusQ::Dependencies INTERFACE OpenMP::OpenMP_CXX )
copy_header_properties( OpenMP::OpenMP_CXX ChronusQ::DepHeaders )

# MPI
if( MPI_PREFIX )
  # MPI_PREFIX paths to CMAKE_PREFIX_PATH
  message( "" )
  message( STATUS "Adding ${MPI_PREFIX} to CMAKE_PREFIX_PATH" )
  set( CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${MPI_PREFIX} )
endif()
include(HandleMPI)


message( "\n\n" )


# Define include directories
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_BINARY_DIR}/include")



# Non-trivial dependencies
message( "ChronusQ Dependencies" )

# Linear Algebra
include(HandleLA)


include(HandleLibint)
include(HandleLibcint)
include(FindLibXC)

set(CQEX_LINK ${CQ_EXT_LINK})

# HDF5
# HDF5 targets aren't until CMake 3.19, so mock them here
find_package(HDF5 REQUIRED CXX)
add_library( ChronusQ::HDF5 INTERFACE IMPORTED )
set_target_properties( ChronusQ::HDF5 PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${HDF5_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES      "${HDF5_CXX_LIBRARIES}"
)
target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::HDF5 )
list(APPEND CQEX_LINK ChronusQ::HDF5 )


message( "\n\n" )


include(HandleTiledArray)


# Compiling with coverage report
if(ENABLE_COVERAGE)

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "ENABLE_COVERAGE Only Compatible with GNU Compilers!")
  else()
    message(STATUS "Enabling Coverage and Profiling -- Warning! May impede performance")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  endif()
endif()



set(BASIS_PATH     "${PROJECT_SOURCE_DIR}/basis"     )
set(AUTHOR_PATH    "${PROJECT_SOURCE_DIR}/AUTHORS"   )

configure_file(
  "${PROJECT_SOURCE_DIR}/include/chronusq_config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/chronusq_config.hpp"
)

# Recurse into subdirectories
add_subdirectory(src)

enable_testing()
add_subdirectory(tests)
