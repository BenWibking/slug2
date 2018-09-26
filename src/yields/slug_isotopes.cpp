/*********************************************************************
Copyright (C) 2014-6 Robert da Silva, Michele Fumagalli, Mark
Krumholz, Evan Demers
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

#include "slug_isotopes.H"
#include "../slug_MPI.H"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// Isotope sorter function
////////////////////////////////////////////////////////////////////////
bool slug_isotopes::isotope_sort(const isotope_data *iso1,
				 const isotope_data *iso2) {
  if (iso1->num() == iso2->num()) {
    return iso1->wgt() < iso2->wgt();
  } else {
    return iso1->num() < iso2->num();
  }
}


////////////////////////////////////////////////////////////////////////
// Isotope table class constructor
////////////////////////////////////////////////////////////////////////
isotope_table::isotope_table(const char *data_dir,
			     slug_ostreams &ostreams) {

  // Open the isotopic data file
  path dirname(data_dir);
  path datafile = dirname / path("isotope_data.txt");
  std::ifstream isofile;
  isofile.open(datafile.c_str());
  if (!(isofile.is_open())) {
    // Bail out if file open failed
    ostreams.slug_err_one << "unable to open isotopic data file "
			  << datafile.string() << endl;
    bailout(1);
  }

  // Process file line by line
  string line;
  while (getline(isofile, line)) {

    // Skip blank and comment lines
    trim(line);
    if (line.length() == 0) continue;
    if (line.compare(0, 1, "#") == 0) continue;

    // Break up line
    vector<string> tokens;
    split(tokens, line, is_any_of("\t "), token_compress_on);
    if (tokens.size() < 3) {
      ostreams.slug_err_one
	<< "badly formatted line '" << line << "' in isotope data file "
	<< datafile.string() << endl;
      bailout(1);
    }
    
    // First three entries are Z, A, lifetime
    unsigned int Z, A;
    double lifetime;
    try {
      Z = lexical_cast<unsigned int>(tokens[0]);
      A = lexical_cast<unsigned int>(tokens[1]);
      lifetime = lexical_cast<double>(tokens[2]);
    } catch (const bad_lexical_cast& ia) {
      (void) ia; // No-op to suppress compiler warning
      ostreams.slug_err_one
	<< "badly formatted line '" << line << "' in isotope data file "
	<< datafile.string() << endl;
      bailout(1);
    }

    // If lifetime < 0, isotope is stable, so push it onto list and
    // continue to next isotope
    if (lifetime < 0) {
      isotope_data *iso = new isotope_data(periodic_table::symbols[Z-1], Z, A);
      all_isotopes_.push_back(iso);
      stable_isotopes_.push_back(iso);
      continue;
    }
    
    // If lifetime is > 0, isotope is unstable, so see if we have
    // daughter nuclides and branching ratios. The format of the rest
    // of the line is
    // branch1  branch2  ...  Z1  A1  Z2  A2 ...
    // so number of daughter products can be deduced from number of
    // entries
    unsigned int ndaughter = (tokens.size()-3) / 3;
    vector<double> branch(ndaughter);
    vector<unsigned int> daughters_Z(ndaughter), daughters_A(ndaughter);
    try {
      for (unsigned int i=0; i<ndaughter; i++) {
	branch[i] = lexical_cast<double>(tokens[i+3]);
	daughters_Z[i] = lexical_cast<unsigned int>(tokens[2*i+3+ndaughter]);
	daughters_A[i] = lexical_cast<unsigned int>(tokens[2*i+4+ndaughter]);
      }
    } catch (const bad_lexical_cast& ia) {
      (void) ia; // No-op to suppress compiler warning
      ostreams.slug_err_one
	<< "badly formatted line '" << line << "' in isotope data file "
	<< datafile.string() << endl;
      bailout(1);
    }

    // Now create the isotope object and push it onto our lists
    isotope_data *iso = new isotope_data(periodic_table::symbols[Z-1], Z, A,
					 lifetime, daughters_Z, daughters_A,
					 branch);
    all_isotopes_.push_back(iso);
    unstable_isotopes_.push_back(iso);
  }

  // Close file
  isofile.close();

  // Build maps that go from both symbols and (Z,A) pairs to isotope
  // data
  for (auto it = all_isotopes_.begin(); it != all_isotopes_.end(); ++it) {
    string symbol = (*it)->symbol() + 
      to_string((unsigned long long) (*it)->wgt());
    isotope_map_symbol[symbol] = *it;
    isotope_map_za[make_pair((*it)->num(), (*it)->wgt())] = *it;
  }

  // Last step: for each unstable isotope, provide the pointers to its
  // daughters
  for (auto it = unstable_isotopes_.begin(); it != unstable_isotopes_.end();
       ++ it) {
    // Loop over daughters
    for (vector<unsigned int>::size_type i=0; i<(*it)->daughters.size(); i++) {
      (*it)->daughters[i] =
	isotope_map_za.at(make_pair((*it)->daughters_Z[i],
				    (*it)->daughters_A[i]));
      assert((*it)->daughters[i] != nullptr);
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Isotope table destructor
////////////////////////////////////////////////////////////////////////
isotope_table::~isotope_table() {
  // Delete all isotopic data
  for (auto it = all_isotopes_.begin(); it != all_isotopes_.end(); ++it) {
    delete *it;
  }
}

