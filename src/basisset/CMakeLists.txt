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
add_library(basisset STATIC reference.cxx basisset.cxx basisset_eval.cxx)
target_link_libraries ( basisset PUBLIC ChronusQ::DepHeaders )

if(TARGET libint)
  add_dependencies(basisset libint)
endif()

# Append basisset to executable link
list(APPEND CQEX_LINK basisset)
set(CQEX_LINK ${CQEX_LINK} PARENT_SCOPE)

if(CQEX_DEP)
  add_dependencies(basisset ${CQEX_DEP})
endif()
