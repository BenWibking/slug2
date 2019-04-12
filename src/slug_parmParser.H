/*********************************************************************
Copyright (C) 2014 Robert da Silva, Michele Fumagalli, Mark Krumholz
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#ifndef _slug_parmParser_H_
#define _slug_parmParser_H_

#include <boost/filesystem.hpp>
#include <boost/variant.hpp>
#include <cassert>
#include "slug.H"
#include "slug_IO.H"
#include "slug_MPI.H"
#include "slug_sim.H"
#include "filters/slug_filter_set.H"
#include "specsyn/slug_specsyn.H"
#include "yields/slug_yields.H"

// Typedef to hold keyword data, listing all the possible types it
// could have; note that we only list here types that are not
// convertible to one another; in particular, though some keywords
// hold things we treat as enum's, an enum is the same as an int, so
// we do not declare enum types separately here
typedef 
boost::variant<int, double, std::string, std::vector<double>,
	       std::vector<std::string> >
slug_parmParser_keyword;

////////////////////////////////////////////////////////////////////////
// class slug_parmParser
//
// This class handles reading slug2 parameter files, and storing their
// contents in a form that can then be read elsewhere in the code.
////////////////////////////////////////////////////////////////////////

class slug_parmParser {

public:
  slug_parmParser(int argc, char **argv, slug_ostreams &ostreams_
#ifdef ENABLE_MPI
		  , MPI_Comm comm
#endif
		  ); // The constructor
  ~slug_parmParser();                     // The destructor

  // Function to return the value of a keyword
  template <class T> const T query(const std::string& name) const {
    return boost::get<T>(keywords.at(name)); }

  // Function to get the path to a directory or file
  const char *fpath(const std::string& name) const {
    return filepaths.at(name).c_str(); }

  // Functions related to checkpointing and restarting
  unsigned int get_checkpoint_ctr() const { return checkpoint_ctr; }
  unsigned int get_checkpoint_trials() const { return checkpoint_trials; }

  // Function to write out the parameters being used to a file
  void writeParams() const;

private:
  
  // Internal functions
  void printUsage();                      // Print a usage message
  void checkParams();                     // Ensure valid parameter values
  void parseFile(std::ifstream &paramFile); // File parsing function
  void setDefaults();                     // Set parameters to default values
  void restartSetup();                    // Set up restart runs

  // Error routines
  [[noreturn]] void parseError(std::string line); 
  [[noreturn]] void valueError(std::string line);

  // Buffer classes to handle I/O
  slug_ostreams &ostreams;

  // Structure to hold the list of key-value pairs
  std::map<const std::string, slug_parmParser_keyword> keywords;

  // List of file paths and directories; most of these are taken
  // directly from keywords in the keywords list, but they may be
  // modified by pre-pending default paths. Thus we keep a separate
  // list of variables for them.
  std::map<const std::string, boost::filesystem::path> filepaths;
  
  // Storage related to checkpointing and restarting; these are stored
  // separately from the keywords because they're not part of the
  // parameter file
  unsigned int checkpoint_ctr;
  unsigned int checkpoint_trials;

#ifdef ENABLE_MPI
  // MPI stuff
  MPI_Comm comm;
  int rank;
#endif
};

#endif
// _slug_parmParser_H_
