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
add_library(particleintegrals STATIC
transform.cxx
spatialtospinor.cxx
aointegrals/aointegrals_onee_drivers.cxx
aointegrals/aointegrals_onee_inhouse.cxx
aointegrals/aointegrals_twoe_libint.cxx
aointegrals/aointegrals_twoe_libcint.cxx
aointegrals/aointegrals_twoe_inhouse.cxx
aointegrals/aointegrals_twoe_ri.cxx
aointegrals/aointegrals_utilities.cxx
aointegrals/giaointegrals_twoe_drivers.cxx
aointegrals/giaointegrals_onee_drivers.cxx
aointegrals/giaointegrals_onee_inhouse.cxx
aointegrals/giaointegrals_twoe_inhouse.cxx
impl.cxx)

target_link_libraries( particleintegrals PUBLIC ChronusQ::DepHeaders )

if(TARGET libint)
  add_dependencies(particleintegrals libint)
endif()

# Append particleintegrals to executable link
list(APPEND CQEX_LINK particleintegrals)
set(CQEX_LINK ${CQEX_LINK} PARENT_SCOPE)


if(CQEX_DEP)
  add_dependencies(particleintegrals ${CQEX_DEP})
endif()
