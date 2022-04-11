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


# Sets a string option with a limited number of values and checks those values
#   Note: This is case insensitive for VARNAME and requires OPTLIST in uppercase
#
#   VARNAME: Name of the option you want to set
#   DEFAULT: Default value for VARNAME
#   OPTLIST: List of possibilities for VARNAME in uppercase
#   DESCRIPTION: Description of what the option controls
#
macro( string_limited_option VARNAME DEFAULT OPTLIST DESCRIPTION )

  set ( ${VARNAME} ${DEFAULT} CACHE STRING ${DESCRIPTION} )
  set_property ( CACHE ${VARNAME} PROPERTY STRINGS ${${OPTLIST}} )

  string(TOUPPER ${${VARNAME}} ${VARNAME})
  if ( NOT ${VARNAME} IN_LIST ${OPTLIST} )
    message ( FATAL_ERROR "${VARNAME} must be one of ${${OPTLIST}}" )
  endif()

endmacro()


# Specific implementation of string_limited_option that specifies ALLOW, FORCE,
#   or NONE as the build type of a package. Sets the option
#   CQ_BUILD_${PKGNAME}_TYPE.
#
#   PKGNAME: Name of the package for the option
#   DEFAULT: What the default build type should be (ALLOW|FORCE|NONE)
#
macro ( build_pkg_option PKGNAME DEFAULT )

  set ( BUILD_${PKGNAME}_TYPES ALLOW FORCE NONE )
  string ( TOLOWER ${PKGNAME} ${PKGNAME}_lower )

  string_limited_option( 
    CQ_BUILD_${PKGNAME}_TYPE
    ${DEFAULT}
    BUILD_${PKGNAME}_TYPES
    "Disable, enable, or force CQ to build ${${PKGNAME}_lower}"
  )

endmacro()

# Adds the given target properties from the first target to the second 
macro( copy_property TARGET1 TARGET2 PROP )
  set_property( TARGET ${TARGET2} APPEND PROPERTY
    ${PROP}
    $<TARGET_PROPERTY:${TARGET1},${PROP}>
  )
endmacro()

function( copy_header_properties TARGET1 TARGET2 )
  copy_property( ${TARGET1} ${TARGET2} INTERFACE_INCLUDE_DIRECTORIES )
  copy_property( ${TARGET1} ${TARGET2} INTERFACE_COMPILE_DEFINITIONS )
  copy_property( ${TARGET1} ${TARGET2} INTERFACE_COMPILE_OPTIONS )
endfunction()
