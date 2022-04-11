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

#include <realtime.hpp>
#include <cxxapi/output.hpp>
#include <physcon.hpp>


namespace ChronusQ {

  void RTFormattedLine(std::ostream &out, std::string s) {
    out << std::setw(38) << "  " + s << std::endl;
  }

  template <typename T>
  void RTFormattedLine(std::ostream &out, std::string s, T v) {
    out << std::setw(38) << "  " + s << v << std::endl;
  }

  template <typename T, typename U>
  void RTFormattedLine(std::ostream &out, std::string s, T v, U u) {
    out << std::setw(38) << "  " + s << v << u << std::endl;
  }


  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::printRTHeader() {
    
    // No printing if silent
    if( this->printLevel == 0 ) return;

    std::cout << BannerTop << std::endl;
    std::cout << "Real-Time Propagation Settings:" << std::endl << std::endl;

    std::cout << std::left << std::setprecision(7);
    std::string AUTime = " \u0127 / Eh";


    RTFormattedLine(std::cout,"* Simulation Parameters:");


    int nSteps = intScheme.tMax / intScheme.deltaT;
    RTFormattedLine(std::cout,"Simulation Time:",intScheme.tMax,AUTime);
    RTFormattedLine(std::cout,"Number of Steps:",nSteps);
    RTFormattedLine(std::cout,"Step Size:",intScheme.deltaT,AUTime);
    RTFormattedLine(std::cout," ",intScheme.deltaT * FSPerAUTime ," fs");




    std::cout << std::endl;
    RTFormattedLine(std::cout,"* Integration Parameters:");

    std::string methString;
    if( intScheme.intAlg == MMUT ) 
      methString = "Modified Midpoint Unitary Transformation (MMUT)"; 
    else if(intScheme.intAlg == ExpMagnus2) 
      methString = "Explicit 2nd Order Magnus"; 

    RTFormattedLine(std::cout,"Electronic Integration:",methString); 

    if( intScheme.intAlg == MMUT ) {

      std::string rstString;
      if( intScheme.rstStep == ForwardEuler )
        rstString = "Forward Euler";
      else if( intScheme.rstStep == ExplicitMagnus2 )
        rstString = "Explicit 2nd Order Magnus";

       RTFormattedLine(std::cout,"Restarting MMUT every ",
         intScheme.iRstrt," steps with a(n) " + rstString + " step");

    }


    if( pert.fields.size() > 0 ) {
      std::cout << std::endl;
      RTFormattedLine(std::cout,"* Perturbation:\n");
  
      for(auto &field : pert.fields) {
        std::cout << std::setw(4) << " ";
        std::cout << "Field " << std::distance(&field,&pert.fields[0]) +1
                  << ":  ";
        
        auto amp = field->getAmp(0);
        if( dynamic_cast<TDDipoleField&>(*field).emFieldTyp == Electric )
          std::cout << "Electric";
        else
          std::cout << "Magnetic";

        std::cout << " ";

        if( amp.size() == 3 ) std::cout << "Dipole";

        std::cout << " Field\n";


        std::cout << std::setw(4) << " ";
        std::cout << std::setw(20) << " * Amplitude (AU)" << "{ ";
        for(auto i = 0; i < amp.size(); i++) {
          std::cout << amp[i]; if(i != amp.size() - 1) std::cout << ", ";
        }
        std::cout << " }\n";
          

        std::cout << std::setw(4) << " ";
        try {
          StepField &env = dynamic_cast<StepField&>(*field->envelope);
          std::cout << std::setw(20) << " * Step Field";
          std::cout << std::setw(9) << "TON = "  << std::setw(10) << env.tOn;
          std::cout << "   ";
          std::cout << std::setw(9) << "TOFF = " << std::setw(10) << env.tOff;
          std::cout << std::endl;
        } catch(...) { }


        


      }

    }


    std::cout << std::endl;
    RTFormattedLine(std::cout,"* Misc Parameters:");
 
    std::string expString;
    if( intScheme.prpAlg == Diagonalization )
      expString = "Eigen Decomposition";
    else if( intScheme.prpAlg == TaylorExpansion )
      expString = "Taylor Expansion";

    RTFormattedLine(std::cout,"Matrix Exponential Method:",expString);
    
    std::cout << std::endl << BannerTop << std::endl;

    std::cout << std::endl << std::fixed << std::right;

    if( this->printLevel == 1 ) {
      std::cout << std::setprecision(4);
      std::cout << std::setw(11) << "Time (AU)" << " ";

      std::cout << std::setprecision(10);
      std::cout << std::setw(16) << "Energy (Eh)" << " ";

      std::cout << std::setprecision(8);
      std::cout << std::setw(16) << "Dipole (X)" << " ";
      std::cout << std::setw(16) << "Dipole (Y)" << " ";
      std::cout << std::setw(16) << "Dipole (Z)" << " ";

      std::cout << std::endl << bannerTop << std::endl << std::endl;
    }

    if( this->printLevel == -1 )
      this->printLevel = 0;

  };

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::printRTStep() { 
    if(this->printLevel == 1) {
      printStepSummary();
    }
    else if(this->printLevel > 1) {
      printStepDetail();
    }
  };

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::printStepSummary() { 
    std::cout << std::fixed << std::right;

    std::cout << std::setprecision(4);
    std::cout << std::setw(11) << curState.xTime << " ";

    std::cout << std::setprecision(10);
    std::cout << std::setw(16) << propagator_.totalEnergy << " ";

    std::cout << std::setprecision(8);
    std::cout << std::setw(16) << propagator_.elecDipole[0] << " ";
    std::cout << std::setw(16) << propagator_.elecDipole[1] << " ";
    std::cout << std::setw(16) << propagator_.elecDipole[2] << " ";

    std::cout << std::endl;

  }; 

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::printStepDetail() { 
    std::cout << bannerTop << "\n\n";
    std::cout << std::fixed << std::right;
    std::cout << "Step: " << std::setw(7) << curState.iStep << '\n';

    std::cout << std::setprecision(5) << "Time: ";
    std::cout << std::setw(11) << curState.xTime << " (au) | ";
    std::cout << std::setw(11) << curState.xTime * FSPerAUTime << " (fs)\n";

    std::cout << std::setprecision(12) << "Energy: ";
    std::cout << std::setw(24) << propagator_.totalEnergy << " (Hartree)\n";

    std::cout << std::setprecision(8) << "Dipole: ";
    std::cout << std::setw(16) << propagator_.elecDipole[0]*EBohrPerDebye << " ";
    std::cout << std::setw(16) << propagator_.elecDipole[1]*EBohrPerDebye << " ";
    std::cout << std::setw(16) << propagator_.elecDipole[2]*EBohrPerDebye << " (Debye)";
    std::cout << std::endl;
  };

}; // namespace ChronusQ

