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

#include <unordered_set>

#include <cxxapi/input.hpp>
#include <cerr.hpp>

namespace ChronusQ {


  /**
   *  \brief Parses a ChronusQ input file
   *
   *  Parses the file and populates the dict_ map which holds the
   *  input data fields to control the ChronusQ calculation
   */
  void CQInputFile::parse() {
  
    // Check if file actually exists
    if(not inFile_->good()) CErr("Input File Couldn't Be Found!",std::cout);


    bool parseSection(false);
    bool prevLineData(false);
    size_t prevIndent(0);
  
    std::string sectionHeader;
    std::string dataHeader;

    // Keywords that are case sensitive (do *not* transform data to UPPER)
    std::unordered_map<std::string, std::unordered_set<std::string> > caseSens;

    // Add case sensitive data keywords here
    caseSens["BASIS"] = {"BASIS"};
        
  
    // Loop over all lines of the file
    while( not inFile_->eof() ) {
  
      std::string line;
      std::getline(*inFile_,line);
  
      // Skip blank lines
      if(line.length() < 1) {
        prevLineData = false;
        continue;
      }
        
  
      // Determine position of first and last non-space character
      size_t firstNonSpace = line.find_first_not_of(" ");
      size_t lastNonSpace  = line.find_last_not_of(" ");
  
      size_t comPos = line.find("#");
  
      // Skip lines in which the first non-space character is #
      // (Comment line)
      if(comPos == firstNonSpace) continue;
  
      // Remove comment portion of the line if it exists
      //  - This is general to when # does not appear in the line
      line = line.substr(0,comPos); 
  
      // Strip trailing spaces
      trim_right(line);
  
  
      auto strToUpper = [](std::string& s){
        std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){ return std::toupper(c);});
      };
  
  
  
      size_t lBrckPos = line.find('[');
      size_t rBrckPos = line.find(']');
  
      size_t eqPos  = line.find('=');
      size_t colPos = line.find(':');
  
      // Determine if this is a line with a section header
      bool sectionLine = 
        lBrckPos == firstNonSpace and rBrckPos == lastNonSpace;
  
      // Determine if this is a line that contains a data field
      bool dataLine    = 
        eqPos != std::string::npos or 
        colPos != std::string::npos;
  
      // Determine if this is a line continuation of a previous data field
      bool multiLine   = prevLineData and firstNonSpace > prevIndent;
  
      // Section line
      if(sectionLine) {
  
        // Strip first spaces
        line = line.substr(firstNonSpace,line.length());
  
        // Obtain the section header name
        sectionHeader = line.substr(1,line.length()-2);
        
        // Convert to UPPER
        strToUpper(sectionHeader);
  
        // Create a dictionary entry for the section header
        dict_[sectionHeader] = 
          std::unordered_map<std::string,std::string>();
  
        // XXX: Possibly check if the section is already defined?
  
        parseSection = true;
        prevLineData = false;
        continue;
  
      }
  
  
      // Data line
      if(parseSection and dataLine) {
  
        line = 
          line.substr(firstNonSpace,line.length()-firstNonSpace);
  
        // Split the line into tokens, trim spaces
        std::vector<std::string> tokens;
        split(tokens,line,"=:");
        for(auto &X : tokens) { trim(X); }
  
        dataHeader = tokens[0];
        strToUpper(dataHeader);

        // Capitalize data if not case sensitive
        if ( not (
          caseSens.find(sectionHeader) != caseSens.end() &&
          caseSens[sectionHeader].find(dataHeader) != caseSens[sectionHeader].end()) &&
          tokens.size() > 1)
          strToUpper(tokens[1]);

        // Create a dictionary entry for the data field in the current
        // section header
        if(tokens.size() > 1) 
          dict_[sectionHeader][dataHeader] = tokens[1];
        else 
          dict_[sectionHeader][dataHeader] = " ";
  
        prevLineData = true;
        prevIndent = firstNonSpace;
      }
  
      // Multiline data
      else if(parseSection and multiLine) {
        // Capitalize data if not case sensitive
        if ( not (
          caseSens.find(sectionHeader) != caseSens.end() &&
          caseSens[sectionHeader].find(dataHeader) != caseSens[sectionHeader].end()))
          strToUpper(line);
 
        line = 
          line.substr(firstNonSpace,line.length()-firstNonSpace);
        dict_[sectionHeader][dataHeader] += "\n" + line;
      }
      
    };
  
  /* Debug code which prints out the contents of the dict_ map
    for(auto &sec : dict_) {
      std::cout << "Section: " << sec.first << std::endl;
      for(auto &data : sec.second) {
        std::cout << "  DATA: " << data.first << " ; " << data.second << std::endl;
      }
    }
  */
  
  }; // CQInputFile::parse
  
  
  
  /** 
   *  \brief Splits a query string on a period "."
   * 
   *  This is a helpder function for the getData function which takes a 
   *  formatted string and splits it into a section and data field.
   *
   *  i.e.  "QM.REFERENCE" -> { "QM", "REFERENCE" }
   *
   *  \param [in] query Query string to be split
   *  \return     std::pair containing the two fields separated by a "."
   */
  std::pair<std::string,std::string> CQInputFile::splitQuery(
    const std::string &query) {
  
    std::vector<std::string> tokens;
  
    // Make sure that the query contains a period
  //assert( query.find(".") != query.end() );
  
    split(tokens,query,".");
    for(auto &X : tokens) {
      trim(X);
      std::transform(X.begin(),X.end(),X.begin(),
        [](unsigned char c){ return std::toupper(c);} );
    }
  
    return 
      std::pair<std::string,std::string>(tokens[0],tokens[1]);

  }; // CQInputFile::splitQuery
  
  
  /**
   *  \brief Custom exception type for handeling the case when
   *  a data field is not found for a query
   */
  class data_not_found : public std::exception {
  
    std::string msg; ///< Error message
  
  public:
  
    // Disable default constructor
    data_not_found() = delete;
  
    /**
     *  Exception constructor. Creates a useful error message
     *  which specifies the failed query
     */ 
    data_not_found(std::string x) { 
      msg = "Data ";
      msg += x; 
      msg += " Not Found\n";
    };
  
    /**
     *  Specialization of std::exception::what. Outputs the error message
     */ 
    virtual const char* what() const throw() {
      return msg.c_str();
    }
  
  }; // data_not_found class
  
  /**
   *  \brief Custom exception type for handeling the case when
   *  a section header is not found for a query
   */
  class section_not_found : public std::exception {
  
    std::string msg; ///< Error message
  
  public:
  
    // Disable default constructor
    section_not_found() = delete;
  
    /**
     *  Exception constructor. Creates a useful error message
     *  which specifies the failed query
     */ 
    section_not_found(std::string x) { 
      msg = "Section ";
      msg += x; 
      msg += " Not Found\n";
    };
  
    /**
     *  Specialization of std::exception::what. Outputs the error message
     */ 
    virtual const char* what() const throw() {
      return msg.c_str();
    }
  
  }; // section_not_found class
  
  
  /**
   *  \brief Specialization of getData to return std::string of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a std::string
   */
  template<>
  std::string CQInputFile::getData(std::string query) {
  
    auto tokenPair = splitQuery(query);
    auto hasSection = dict_.find(tokenPair.first);
  
    if(hasSection != dict_.end()) {
      auto hasData = 
        dict_[tokenPair.first].find(tokenPair.second);
  
      if(hasData != dict_[tokenPair.first].end())
        return 
          dict_[tokenPair.first][tokenPair.second];
      else throw data_not_found(query);
    } else throw section_not_found(tokenPair.first);
  
  }; // CQInputFile::getData<std::string>
  
  /**
   *  \brief Specialization of getData to return int of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a int
   */
  template<>
  int CQInputFile::getData(std::string query) {
  
    return std::stoi(getData<std::string>(query));
  
  }; // CQInputFile::getData<int>
  
  /**
   *  \brief Specialization of getData to return bool of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a bool
   */
  template<>
  bool CQInputFile::getData(std::string query) {
  
    query = getData<std::string>(query);
    if (not query.compare("TRUE") or not query.compare("ON")){
      return true;
    }
      
    if (not query.compare("FALSE") or not query.compare("OFF")){
      return false;
    }
  
    CErr("Invalid Input For Boolean-Type Keyword!");

    return false;
  }; // CQInputFile::getData<bool>
  
  /**
   *  \brief Specialization of getData to return size_t of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a size_t
   */
  template<>
  size_t CQInputFile::getData(std::string query) {
  
    return std::stoul(getData<std::string>(query));
  
  }; // CQInputFile::getData<size_t>
  
  /**
   *  \brief Specialization of getData to return double of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a double
   */
  template<>
  double CQInputFile::getData(std::string query) {
  
    return std::stod(getData<std::string>(query));
  
  }; // CQInputFile::getData<double>

}; // namespace ChronusQ

