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

include(CheckCXXSourceCompiles)
include(FetchContent)

function ( check_libint_am )

   set ( LIBINT_MAX_AM_SOURCE 
   "
   #include <libint2/config.h>

   #if LIBINT_MAX_AM < 5
   #error Libint Angular Momentum must be >= 5
   #endif

   int main() { return 0; }
   " )

   set ( CMAKE_REQUIRED_INCLUDES ${LIBINT2_INCLUDE_DIRS} )
   check_cxx_source_compiles("${LIBINT_MAX_AM_SOURCE}" LIBINT_MAX_AM_GOOD )
   set ( LIBINT_MAX_AM_GOOD PARENT_SCOPE )

endfunction()


message ( "\n == Libint ==" )

#
#  Find preinstalled Libint unless turned off
#

if ( NOT CQ_BUILD_LIBINT_TYPE STREQUAL "FORCE" )

  find_package ( Libint2 QUIET )


  # Prefer CMake installed
  if ( TARGET Libint2::int2 )

    check_libint_am()

    if ( LIBINT_MAX_AM_GOOD )

      get_target_property ( libint_alias Libint2::int2 ALIASED_TARGET )

      if ( libint_alias )
        add_library ( ChronusQ::Libint2 ALIAS ${libint_alias} )
        get_target_property ( libint_location ${libint_alias} LOCATION )
      else()
        add_library ( ChronusQ::Libint2 ALIAS Libint2::int2 )
        get_target_property ( libint_location Libint2::int2 LOCATION )
      endif()

    endif() 

  # Otherwise create dummy target
  elseif ( DEFINED Libint2_ROOT )

    find_library ( LIBINT2_LIBRARIES
                   NAMES int2 libint2
                   PATHS "${Libint2_ROOT}/lib"
                 )

    if ( LIBINT2_LIBRARIES AND
         EXISTS "${Libint2_ROOT}/include/libint2.hpp" )
      
      set ( LIBINT2_INCLUDE_DIRS
          ${Libint2_ROOT}/include
          ${Libint2_ROOT}/include/libint2
      )
      
      check_libint_am()

      if ( LIBINT_MAX_AM_GOOD )

        add_library ( ChronusQ::Libint2 INTERFACE IMPORTED )
        set_target_properties ( ChronusQ::Libint2 PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${LIBINT2_INCLUDE_DIRS}"
          INTERFACE_LINK_LIBRARIES      "${LIBINT2_LIBRARIES}"
        )

        set ( libint_location ${LIBINT2_LIBRARIES} )
        set ( LIBINT2_FOUND TRUE )

      else()

        message( STATUS
"Found Libint installation ${LIBINT2_LIBRARIES}, but the maximum angular \
momentum is below that required by ChronusQuantum (5)" )

      endif()

    endif()

  endif()

  if ( LIBINT2_FOUND )
    message ( STATUS "Found External Libint Installation: ${libint_location}" )
  endif()

endif()


#
#  Build Libint if a suitable libint hasn't been found
#

if ( NOT TARGET ChronusQ::Libint2 )

  if ( NOT CQ_BUILD_LIBINT_TYPE STREQUAL "NONE" )

    if ( NOT CQ_BUILD_LIBINT_TYPE STREQUAL "FORCE" )
 
      message ( STATUS "Checking for a previously built Libint for CQ" )
 
      find_library ( LIBINT2_LIBRARIES
                     NAMES int2 libint2
                     PATHS "${FETCHCONTENT_BASE_DIR}/libint2-build"
                     NO_DEFAULT_PATH )
 
      if ( LIBINT2_LIBRARIES AND
           EXISTS "${FETCHCONTENT_BASE_DIR}/libint2-build/include/libint2/config.h" AND
           EXISTS "${FETCHCONTENT_BASE_DIR}/libint2-src/include/libint2.hpp" )
 
 
        set ( LIBINT2_INCLUDE_DIRS
              "${FETCHCONTENT_BASE_DIR}/libint2-build/include"
              "${FETCHCONTENT_BASE_DIR}/libint2-build/include/libint2"
              "${FETCHCONTENT_BASE_DIR}/libint2-src/include" 
              "${FETCHCONTENT_BASE_DIR}/libint2-src/include/libint2"
        )
 
        add_library ( ChronusQ::Libint2 INTERFACE IMPORTED )
        set_target_properties ( ChronusQ::Libint2 PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${LIBINT2_INCLUDE_DIRS}"
          INTERFACE_LINK_LIBRARIES      "${LIBINT2_LIBRARIES}"
          INTERFACE_COMPILE_DEFINITIONS "__COMPILING_LIBINT2"
        )
        
        message ( STATUS "Found CQ Libint Installation: ${LIBINT2_LIBRARIES}" )

      endif()

    endif()

    # If we're forcing building or failed finding a previously built version of
    #   libint for CQ
    if ( NOT TARGET ChronusQ::Libint2 )

      message ( STATUS "Opting to build a copy of Libint" )

      # Update policy for Libint
      set ( CMAKE_POLICY_DEFAULT_CMP0074 NEW )

  
    FetchContent_Declare (
      Libint2
      PREFIX ${CUSTOM_LIBINT_PREFIX}
      GIT_REPOSITORY "https://urania.chem.washington.edu/chronusq/libint-cq.git"
      GIT_TAG "2.7.0-beta.6"
    )
  
      FetchContent_GetProperties ( Libint2 )

      if ( NOT Libint2_POPULATED )

        message ( STATUS "Downloading Libint... (This may take several minutes)" )
        FetchContent_Populate ( Libint2 )
        message ( STATUS "Downloading Libint - Done" )

        add_subdirectory ( ${libint2_SOURCE_DIR} ${libint2_BINARY_DIR} )

        # ALWAYS build libint using unity build by setting UNITY_BUILD for libint2_obj
        # this saves time and works around issues with Libint library being too big for Ninja on MacOS
        if (NOT CMAKE_UNITY_BUILD)
            set_target_properties(libint2_obj PROPERTIES UNITY_BUILD ON)
            message(STATUS "Will unity-build Libint2")
        endif()

      endif()
  
      install(TARGETS libint2  
              ARCHIVE DESTINATION  "lib/libint2"
              LIBRARY DESTINATION  "lib/libint2"
              INCLUDES DESTINATION "include/libint2")
    
      add_library ( ChronusQ::Libint2 ALIAS libint2 )

    endif()
  
  else()
  
    message ( FATAL_ERROR "Suitable Libint installation could not be found! \
Set Libint2_ROOT to the prefix of the Libint installation or set \
CQ_BUILD_LIBINT to ALLOW or FORCE."
    )
  
  endif()

endif()

target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::Libint2 )
copy_header_properties( ChronusQ::Libint2 ChronusQ::DepHeaders )
list(APPEND CQ_EXT_LINK ChronusQ::Libint2)

# Workaround to intel compiler bug breaking libint2 ERI evaluation
if( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
  set_property( TARGET ChronusQ::Dependencies APPEND PROPERTY
    INTERFACE_COMPILE_DEFINITIONS "LIBINT2_CONSTEXPR_STATICS=0"
  )
endif()

message ( " == End Libint ==\n" )
