# 
#  This file is part of the Chronus Quantum (ChronusQ) software package
#  
#  Copyright (C) 2014-2022 Li Research Group (University of Washington)
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#  
#  Contact the Developers:
#    E-Mail: xsli@uw.edu
#  
#

include( FetchContent )

# Do nothing if TA is disabled
if( CQ_ENABLE_TA )

  if( NOT CQ_ENABLE_MPI )
    message(FATAL_ERROR "TA requires MPI. Please rerun with CQ_ENABLE_MPI=On" )
  endif()

  # Try to load all builds that tiledarray may have recieved by fetchcontent
  set( ta_deps "MADNESS" "BTAS" )
  foreach( ta_dep IN LISTS ta_deps )
    string( TOLOWER ${ta_dep} ta_dep_lower )
    list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/${ta_dep_lower}-build )
    find_package( ${ta_dep} QUIET )
    if( ${ta_dep}_FOUND )
      message( "Found TA dependency ${ta_dep}" )
    endif()
  endforeach()

  
  find_package( tiledarray QUIET )
  
  if( TARGET tiledarray )
    message( STATUS "Found external TiledArrray installation" )
    add_library( ChronusQ::TiledArray ALIAS tiledarray )
  else()

    list( APPEND CMAKE_PREFIX_PATH ${FETCHCONTENT_BASE_DIR}/tiledarray-build )
    find_package( tiledarray QUIET )
    
    if( TARGET tiledarray )
      message( STATUS "Found previously built TiledArray" )
      add_library( ChronusQ::TiledArray INTERFACE IMPORTED )

      # Mock the target instead of alias because the headers aren't correctly set by tiledarray
      set_property( TARGET ChronusQ::TiledArray PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES
	  $<BUILD_INTERFACE:${TILEDARRAY_BUILD_INCLUDE_DIRS}>
          $<INSTALL_INTERFACE:${TILEDARRAY_INSTALL_INCLUDE_DIRS}>
      )
      set_property( TARGET ChronusQ::TiledArray PROPERTY
        INTERFACE_LINK_LIBRARIES
	  tiledarray
      )

    else()

      message( STATUS "Building TiledArray" )
      # Workaround for tiledarray using an old version of the LA discovery modules
      if( DEFINED BLA_VENDOR )
        set( temp_BLA_VENDOR ${BLA_VENDOR} )
        unset( BLA_VENDOR )
      endif()
      if( DEFINED CACHE{BLA_VENDOR} )
        set( temp_BLA_VENDOR_CACHE $CACHE{BLA_VENDOR} )
        unset( BLA_VENDOR CACHE )
      endif()
      
      # Variables for tiledarray to figure out if this is a built or installed version afterwards
      set( TILEDARRAY_BINARY_DIR "${FETCHCONTENT_BASE_DIR}/tiledarray-build" )
      set( TILEDARRAY_SOURCE_DIR "${FETCHCONTENT_BASE_DIR}/tiledarray-src" )
      
      FetchContent_Declare( tiledarray
        GIT_REPOSITORY https://github.com/ValeevGroup/tiledarray.git
        GIT_TAG b3b25af5f541f579c149dabbf5cd5646a8847cb2
        PATCH_COMMAND git am ${PROJECT_SOURCE_DIR}/external/tiledarray/patch/external-umpire.patch
      )
      
      FetchContent_MakeAvailable( tiledarray )
      
      add_library( ChronusQ::TiledArray ALIAS tiledarray )
      list( APPEND CQEX_DEP ChronusQ::TiledArray )
      list( APPEND CQEX_LINK ChronusQ::TiledArray )

      if( TARGET Umpire )
        file( MAKE_DIRECTORY ${FETCHCONTENT_BASE_DIR}/umpire-src/src/umpire/tpl/camp/include )
        file( MAKE_DIRECTORY ${FETCHCONTENT_BASE_DIR}/umpire-build/include )
      endif()
      
      if( DEFINED temp_BLA_VENDOR )
        set( BLA_VENDOR ${temp_BLA_VENDOR} )
        unset( temp_BLA_VENDOR )
      endif()
      if( DEFINED temp_BLA_VENDOR_CACHE )
        set( BLA_VENDOR ${temp_BLA_VENDOR_CACHE} CACHE STRING "" FORCE )
      endif()

    endif()
  endif()
  
  if( TARGET ChronusQ::TiledArray )
    set( CQ_HAS_TA ON )
  endif()

  target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::TiledArray )
  copy_header_properties( ChronusQ::TiledArray ChronusQ::DepHeaders )

else( CQ_ENABLE_TA )
  message( STATUS "TiledArray not enabled; Coupled-cluster functonality in CQ is disabled." )
endif( CQ_ENABLE_TA )
