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
include(FetchContent)

message( "\n\n" )
message( "ChronusQ Linear Algebra Settings:\n" )

# Eigen3
find_package(Eigen3 REQUIRED)
set_property( TARGET ChronusQ::Dependencies APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR}
)
set_property( TARGET ChronusQ::DepHeaders APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR}
)


# Better BLAS discovery
FetchContent_Declare( la_cmake
  GIT_REPOSITORY https://github.com/ajaypanyala/linalg-cmake-modules.git
  GIT_TAG 4f7bc30697f0351012356ddc2505924940482f98
)
FetchContent_GetProperties( la_cmake )
if( NOT la_cmake_POPULATED )
  FetchContent_Populate( la_cmake )
endif()
list( APPEND CMAKE_MODULE_PATH ${la_cmake_SOURCE_DIR} )

################################# BLAS ########################################

# Find external BLAS installation
if( NOT TARGET ChronusQ::BLAS AND BLAS_EXTERNAL )
  message("")
  message( STATUS "Searching for external BLAS\n" )
  find_package( BLAS QUIET )

  if( TARGET BLAS::BLAS )
    add_library( ChronusQ::BLAS ALIAS BLAS::BLAS )
  endif()
endif()

# Try to find BLAS already compiled by CQ
if (NOT TARGET ChronusQ::BLAS )
  message("")
  message( STATUS "Searching for BLAS built for CQ\n" )
  set( OPENBLAS_PREFIX ${PROJECT_SOURCE_DIR}/external/openblas )

  # Workaround to not having NO_DEFAULT_PATH in module mode
  set( temp_storage ${CMAKE_PREFIX_PATH} )
  set( CMAKE_PREFIX_PATH ${OPENBLAS_PREFIX} )

  find_package( BLAS QUIET )
  if( TARGET BLAS::BLAS )
    add_library( ChronusQ::BLAS ALIAS BLAS::BLAS )
  endif()

  # If this was built for CQ, we will need to install it to somewhere else
  set( BLAS_INCLUDE_DIRS
    $<BUILD_INTERFACE:${BLAS_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
  )

  set_property( TARGET BLAS::BLAS PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES
      ${BLAS_INCLUDE_DIRS}
  )

  set( CMAKE_PREFIX_PATH ${temp_storage} ${OPENBLAS_PREFIX} )

endif() 
    
# A BLAS needs to be built externally
if (NOT TARGET ChronusQ::BLAS ) 
  message(FATAL_ERROR "
No BLAS/LAPACK Libraries Have Been Found!
You can install a BLAS library by running ./bin/buildblas from ${PROJECT_SOURCE_DIR}
")
else()
  message("")
  message( STATUS "Found BLAS library: ${BLAS_LIBRARIES}" )
endif()

############################### LAPACK ########################################

# Find external BLAS installation
if( NOT TARGET ChronusQ::LAPACK AND BLAS_EXTERNAL )
  message("")
  message( STATUS "Searching for external LAPACK\n" )
  find_package( LAPACK QUIET )

  if( TARGET LAPACK::LAPACK )
      add_library( ChronusQ::LAPACK ALIAS LAPACK::LAPACK )
  endif()
endif()

# Try to find BLAS already compiled by CQ
if (NOT TARGET ChronusQ::LAPACK )
  message("")
  message( STATUS "Searching for LAPACK built for CQ\n" )
  set( OPENBLAS_PREFIX ${PROJECT_SOURCE_DIR}/external/openblas )

  # Workaround to not having NO_DEFAULT_PATH in module mode
  set( temp_storage ${CMAKE_PREFIX_PATH} )
  set( CMAKE_PREFIX_PATH ${OPENBLAS_PREFIX} )

  find_package( LAPACK QUIET )
  if( TARGET LAPACK::LAPACK )
    add_library( ChronusQ::LAPACK ALIAS LAPACK::LAPACK )
  endif()

  # If this was installed for CQ, we will need to install it somewhere else
  set( LAPACK_INCLUDE_DIRS
    $<BUILD_INTERFACE:${LAPACK_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
  )

  set_property( TARGET LAPACK::LAPACK PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES
      ${LAPACK_INCLUDE_DIRS}
  )


  set( CMAKE_PREFIX_PATH ${temp_storage} ${OPENBLAS_PREFIX} )

endif() 
    
# A BLAS needs to be built externally
if (NOT TARGET ChronusQ::LAPACK ) 
  message(FATAL_ERROR "
No LAPACK Libraries Have Been Found!
You can install a LAPACK library by running ./bin/buildblas from ${PROJECT_SOURCE_DIR}
")
else()
  message("")
  message( STATUS "Found LAPACK library: ${BLAS_LIBRARIES}" )
endif()

set( CQ_LINALG_LIBRARIES ChronusQ::BLAS ChronusQ::LAPACK )

target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::BLAS )
copy_header_properties( ChronusQ::BLAS ChronusQ::DepHeaders )

target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::LAPACK )
copy_header_properties( ChronusQ::LAPACK ChronusQ::DepHeaders )

# TODO: Switch this to BLAS_VENDOR
if( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
  set( _CQ_MKL 1 )
endif()

####################### BLACS + ScaLAPACK LIBRARIES ###########################
if( CQ_ENABLE_MPI )


  message( STATUS "CQ_ENABLE_MPI Triggers Search for ScaLAPACK/BLACS" )

  # CXXBLACS
  set( CQ_ENABLE_CXXBLACS TRUE )
  message( STATUS "---> Creating CMake Target for CXXBLACS" )

  set( CXXBLACS_PREFIX     ${PROJECT_SOURCE_DIR}/external/cxxblacs )
  set( CXXBLACS_INCLUDEDIR ${CXXBLACS_PREFIX}/src/cxxblacs/include )
  
  ExternalProject_Add(cxxblacs
    PREFIX ${CXXBLACS_PREFIX}
    GIT_REPOSITORY https://github.com/wavefunction91/CXXBLACS.git
    CONFIGURE_COMMAND echo 'No CXXBLACS Configure'
    UPDATE_COMMAND echo 'No CXXBLACS Update Command'
    BUILD_COMMAND echo 'No CXXBLACS Build'
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND echo 'No CXXBLACS Install'
  )

  list(APPEND CQEX_DEP cxxblacs)
  
  # CXXBLACS Includes
  file( MAKE_DIRECTORY ${CXXBLACS_INCLUDEDIR} )
  set_property( TARGET ChronusQ::Dependencies APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CXXBLACS_INCLUDEDIR}
  )
  set_property( TARGET ChronusQ::DepHeaders APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CXXBLACS_INCLUDEDIR}
  )



  if( NOT TARGET ScaLAPACK )
    # Externally installed scalapack
    find_package( ScaLAPACK QUIET )

    if( TARGET ScaLAPACK::ScaLAPACK )
      message( STATUS "Found external ScaLAPACK installation" )
      add_library( ChronusQ::ScaLAPACK ALIAS ScaLAPACK::ScaLAPACK )
    else()
      list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/scalapack-build )
      find_package( ScaLAPACK QUIET )

      if( TARGET ScaLAPACK::ScaLAPACK )
        message( STATUS "Found previously installed ScaLAPACK" )
        add_library( ChronusQ::ScaLAPACK ALIAS ScaLAPACK::ScaLAPACK )
      else()
        message( STATUS "Building ScaLAPACK" )
        FetchContent_Declare( ScaLAPACK
          GIT_REPOSITORY https://github.com/Reference-ScaLAPACK/scalapack.git
          GIT_TAG bc6cad585362aa58e05186bb85d4b619080c45a9
        )
        FetchContent_MakeAvailable( ScaLAPACK )
	set_property( TARGET scalapack APPEND PROPERTY COMPILE_FLAGS "-Wno-implicit-function-declaration" )
	set_property( TARGET scalapack APPEND PROPERTY INTERFACE_COMPILE_FLAGS "-Wno-implicit-function-declaration" )
        add_library( ChronusQ::ScaLAPACK ALIAS scalapack )
      endif()
    endif()
  endif()





  # Append ScaLAPACK / BLACS to linker
  list( APPEND CQ_LINALG_LIBRARIES ChronusQ::ScaLAPACK )
  target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::ScaLAPACK )
  copy_header_properties( ChronusQ::ScaLAPACK ChronusQ::DepHeaders )

  if( CQ_BLACS_LIBRARIES )
    set( CQ_LINALG_LIBRARIES ${CQ_LINALG_LIBRARIES} "${CQ_BLACS_LIBRARIES}")
  endif()

endif( CQ_ENABLE_MPI )










# Add CQ_LINALG_LIBRARIES to linker
if( CQ_LINALG_LIBRARIES )

  list(APPEND CQ_EXT_LINK ${CQ_LINALG_LIBRARIES})

  message(STATUS "CQ_LINALG_LIBRARIES = ${CQ_LINALG_LIBRARIES}")

endif( CQ_LINALG_LIBRARIES )


############################# BLAS++ ##########################################

# If a dependency has already pulled in BLAS++, use theirs
if( NOT TARGET blaspp )
  # Find external BLAS++ installation
  find_package( blaspp QUIET )
  
  if( TARGET blaspp )
    message( STATUS "Found external BLAS++ installation" )
    add_library( ChronusQ::blaspp ALIAS blaspp )
  else() 
  
    set( blaspp_PREFIX ${PROJECT_SOURCE_DIR}/external/blaspp )
    list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/blaspp-build )
    find_package( blaspp QUIET )
  
    if( TARGET blaspp )
      message( STATUS "Found previously installed BLAS++" )
      add_library( ChronusQ::blaspp ALIAS blaspp )
    else()
      message( STATUS "Could not find BLAS++. Compiling BLAS++." )
      FetchContent_Declare( blaspp 
        GIT_REPOSITORY https://bitbucket.org/icl/blaspp.git
        GIT_TAG ed392fe
      )
      FetchContent_MakeAvailable( blaspp )
    endif()
  endif()
endif()
add_library( ChronusQ::blaspp ALIAS blaspp )

############################# LAPACK++ ########################################

# If a dependency has already pulled in LAPACK++, use theirs
if( NOT TARGET lapackpp )
  # Find external LAPACK++ installation
  find_package( lapackpp QUIET )
  
  if( TARGET lapackpp )
    message( STATUS "Found external LAPACK++ installation" )
  else() 
  
    set( lapackpp_PREFIX ${PROJECT_SOURCE_DIR}/external/lapackpp )
    list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/lapackpp-build )
    find_package( lapackpp QUIET )
  
    if( TARGET lapackpp )
      message( STATUS "Found previously installed LAPACK++" )
    else()
      message( STATUS "Could not find LAPACK++. Compiling LAPACK++." )
      FetchContent_Declare( lapackpp
        GIT_REPOSITORY https://bitbucket.org/icl/lapackpp.git
        GIT_TAG dbcf60f
        PATCH_COMMAND cd ${FETCHCONTENT_BASE_DIR}/lapackpp-src
        && patch -p0 < ${lapackpp_PREFIX}/patch/fortran.patch || exit 0
      )
      # && patch -p0 < ${lapackpp_PREFIX}/patch/config.patch || exit 0
      FetchContent_MakeAvailable( lapackpp )
    endif()
  endif()
endif()
add_library( ChronusQ::lapackpp ALIAS lapackpp )


target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::blaspp )
target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::lapackpp )

copy_header_properties( ChronusQ::blaspp ChronusQ::DepHeaders )
copy_header_properties( ChronusQ::lapackpp ChronusQ::DepHeaders )

list(APPEND CQEX_DEP ChronusQ::blaspp ChronusQ::lapackpp )

message( "\n\n\n" )
