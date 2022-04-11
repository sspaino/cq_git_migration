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

#include <cxxapi/options.hpp>
#include <cerr.hpp>
#include <typeinfo>
#include <stdlib.h>

namespace ChronusQ {
  /**
   *
   *  Check valid keywords in the section.
   *
  */
 void CQCC_VALID( std::ostream& out, CQInputFile& input){
   //Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "TYPE",
      "USEDIIS",
      "NDIIS",
      "ETOL",
      "TTOL",
      "MAXITER"
    };
      // Specified keywords
    std::vector<std::string> ccKeywords = input.getDataInSection("CC"); 
     // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : ccKeywords ) {
    auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
    if( ipos == allowedKeywords.end() ) 
      CErr("Keyword CC." + keyword + " is not recognized",std::cout);// Error
    }
  // Check for disallowed combinations (if any)
  }

#ifdef CQ_HAS_TA
  //Construct a CCBase object from input file
  std::shared_ptr<CCBase> CQCCOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<SingleSlaterBase>& ss) {
     
    if( not input.containsSection("CC") )
      CErr("CC Section must be specified for CC job",out);
    //check if ALG is set to INCORE. If not, abort.
    std::string ALG;
    OPTOPT( ALG = input.getData<std::string>("INTS.ALG"); )
    trim(ALG);
    if(  ALG.compare("INCORE") ) {
      CErr("To use CCSD code, you must set ALG=INCORE in the [INTS] section.");
    }
        
    out << "\n  *** Parsing CC options ***\n";

    std::shared_ptr<CCBase> cc;


    // Determine  reference and construct CC object

    bool isCCSD = false;
    bool found = false;
    #define CONSTRUCT_CC(_MBTYPE,_TYPEMATST,_TYPEINTST,_CLASS) \
      if(not found) try { \
      cc = std::dynamic_pointer_cast<CCBase>(\
          std::make_shared<_MBTYPE<_TYPEMATST,_TYPEINTST>>(\
            dynamic_cast<_CLASS<_TYPEMATST,_TYPEINTST>&>(*ss)\
          )\
        );\
        found = true; \
      } catch(...) { }
      
    OPTOPT(
      std::string ccopts = input.getData<std::string>("CC.TYPE");

      if( not ccopts.compare("CCSD") ) {isCCSD = true; }
      else CErr(ccopts + " NOT RECOGNIZED CC.TYPE");

    );

    if (isCCSD){
      CONSTRUCT_CC(CCSD, dcomplex, double, HartreeFock);
      if(input.containsData("CC.USEDIIS")){
        OPTOPT(cc->useDIIS = input.getData<bool>("CC.USEDIIS");)
      }

      if(input.containsData("CC.NDIIS")){
        if(not cc->useDIIS){
          std::cout << "You must enable DIIS to specify nDIIS. " << std::endl;
        }
        else{
          OPTOPT(cc->nDIIS = input.getData<size_t>("CC.NDIIS");)
        }  
      }

      if(input.containsData("CC.ETOL")){
        OPTOPT(cc->ccSettings.eConv = input.getData<double>("CC.ETOL");)
      }

      if(input.containsData("CC.TTOL")){
        OPTOPT(cc->ccSettings.tConv = input.getData<double>("CC.TTOL");)
      }

      if(input.containsData("CC.MAXITER")){
        OPTOPT(cc->ccSettings.maxiter = input.getData<int>("CC.MAXITER");)
      }
    }
    else {
      CErr("NYI");
    }

    return cc;
  }
#endif
};

