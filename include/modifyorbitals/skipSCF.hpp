/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *
 */
#pragma once

#include <modifyorbitals/conventionalSCF.hpp>

namespace ChronusQ {

template <typename MatsT>
class SkipSCF : public ConventionalSCF<MatsT> {
public:
  // Constructor
  SkipSCF() = delete;
  SkipSCF(SCFControls scfControls, MPI_Comm comm, ModifyOrbitalsOptions<MatsT> mod, CQMemManager& mem):
      ConventionalSCF<MatsT>(scfControls, comm, mod, mem) {}

  // Functions

  // RunModifyOrbitals is a stripped down version from OptimizeOrbitals
  void runModifyOrbitals(EMPerturbation& pert, VecMORef<MatsT>& mo, VecEPtr& eps) {
    ProgramTimer::tick("SCF Total");

    // Form the Fock matrix D(k) -> F(k)
    ProgramTimer::timeOp("Form Fock", [&]() { this->modOrbOpt.formFock(pert); });
    this->computeEigenvalues(mo,eps);

    // Compute initial properties
    this->modOrbOpt.computeProperties(pert);

    if( this->scfControls.printLevel > 0 and MPIRank(this->comm) == 0 ) {
      this->printRunHeader(std::cout, pert);
      this->printIteration(std::cout, false);
    }

    // Save current state of the wave function (method specific)
    this->modOrbOpt.saveCurrentState();

    // printSCFFooter(isConverged);
    if( this->scfControls.printLevel > 0 ) {
      std::cout << std::endl
                << "SCF Completed: E(" << this->scfControls.refShortName_ << ") = " << std::fixed << std::setprecision(10)
                << this->modOrbOpt.getTotalEnergy() << " Eh after " << this->scfConv.nSCFIter << " SCF Iterations" << std::endl;
    }

    if( this->scfControls.printLevel > 0 ) std::cout << BannerEnd << std::endl;

    if( this->scfControls.printLevel > 1 ) { this->modOrbOpt.printProperties(); }

    ProgramTimer::tock("SCF Total");

  };   // SkipSCF<MatsT>::RunModifyOrbitals()
};
};   // namespace ChronusQ
