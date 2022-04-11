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

include(ExternalProject)
set( LIBXC_PREFIX ${PROJECT_SOURCE_DIR}/external/libxc )
set( LIBXC_INCLUDEDIR ${LIBXC_PREFIX}/include )
set( LIBXC_LIBDIR ${LIBXC_PREFIX}/lib )

if( NOT EXISTS "${LIBXC_PREFIX}/include/xc.h" )

  ExternalProject_Add(libxc
    PREFIX ${LIBXC_PREFIX}
    URL "${LIBXC_PREFIX}/libxc-4.0.4.tar.gz"
    CONFIGURE_COMMAND ./configure 
      --prefix=${LIBXC_PREFIX} 
      #CC=gcc
      #CFLAGS=${CMAKE_C_FLAGS} 
      #FC=gfortran
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
  )

  list(APPEND CQEX_DEP libxc)
  file( MAKE_DIRECTORY ${LIBXC_INCLUDEDIR} )
  message( STATUS "Opting to build a copy of LibXC" )
else()
  message( STATUS "Found LibXC installation!" )
endif()

set_property( TARGET ChronusQ::Dependencies APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${LIBXC_INCLUDEDIR}
)
set_property( TARGET ChronusQ::DepHeaders APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${LIBXC_INCLUDEDIR}
)

set_property( TARGET ChronusQ::Dependencies APPEND PROPERTY
  INTERFACE_LINK_LIBRARIES ${LIBXC_LIBDIR}/libxc.a
)
