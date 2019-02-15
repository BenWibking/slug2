/*********************************************************************
Copyright (C) 2014-7 Robert da Silva, Michele Fumagalli, Mark Krumholz
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

#include "slug_yields_karakas16_doherty14.H"
#include "../constants.H"
#include "../slug_MPI.H"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/algorithm/sort.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// A little utility typedef to facilitate sorting yield files by mass
////////////////////////////////////////////////////////////////////////
namespace yields {

  // A data type for holding names of yield files and their
  // corresponding masses, and for sorting them by mass, partial
  // mixing zone mass, and overshoot
  struct karakas16_file_data {
    string fname;
    double mass;
    double pmz_mass;
    AGB_Karakas16_overshoot overshoot;
    friend bool operator<(const karakas16_file_data& f1, 
			  const karakas16_file_data& f2) {
      if (f1.mass != f2.mass)
	return f1.mass < f2.mass;
      else if (f1.pmz_mass != f2.pmz_mass)
	return f1.pmz_mass < f2.pmz_mass;
      else
	return f1.overshoot < f2.overshoot;
    }
  };
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////

slug_yields_karakas16_doherty14::~slug_yields_karakas16_doherty14() {
  // De-allocate all splines and accelerateors
  for (vector<double>::size_type i=0; i<isotopes.size(); i++) {
    gsl_spline_free(yield_interp[i]);
    gsl_interp_accel_free(yield_accel[i]);
  }
}


////////////////////////////////////////////////////////////////////////
// Initializatio routine
////////////////////////////////////////////////////////////////////////
void slug_yields_karakas16_doherty14::
init(const char *yield_dir,
     const double pmz_mass,
     const AGB_Karakas16_overshoot overshoot,
     const bool include_all) {

  // Read the Karakas & Lugaro tables
  read_karakas16_tables(yield_dir, pmz_mass, overshoot);
  
  // Read the Doherty tables
  read_doherty14_tables(yield_dir);

  // Construct the table of isotopes combining the two tables. We do
  // this in one of two ways: either including all isotopes even if
  // only one table has them, or only including isotopes that are
  // present in both tables; if including all isotopes, issue warning
  // about isotopes only found in one list
  if (include_all) {

    // All isotopes included
    set<const isotope_data *> isotope_set;
    for (vector<int>::size_type i = 0; i < isotopes_k16.size(); i++)
      isotope_set.insert(isotopes_k16[i]);
    for (vector<int>::size_type i = 0; i < isotopes_d14.size(); i++)
      isotope_set.insert(isotopes_d14[i]);
    isotopes.assign(isotope_set.begin(), isotope_set.end());

    // Issue warning
    stringstream ss;
    for (vector<int>::size_type i = 0; i < isotopes.size(); i++) {
      if (find(isotopes_d14.begin(), isotopes_d14.end(), isotopes[i]) ==
	  isotopes_d14.end())
	ss << " " << isotopes[i]->symbol() << isotopes[i]->wgt();
    }
    ostreams.slug_warn_one
      << "the following isotopes are present in the AGB yields from "
      << "Karakas & Lugaro (2016), but not the super-AGB yields from "
      << "Doherty+ (2014); their yield will be set to zero for "
      << "super-AGB stars:"
      << ss.str()
      << endl;

  } else {

    // Only isotopes in common included
    set_intersection(isotopes_k16.begin(), isotopes_k16.end(),
		     isotopes_d14.begin(), isotopes_d14.end(),
		     back_inserter(isotopes),
		     slug_isotopes::isotope_sort);
    
  }    

  // Now that we have all the isotopes set up, call the parent
  // isotope_init function
  isotope_init();

  // Record mass range
  agb_mass.push_back(min(mass_k16.front(), mass_d14.front()));
  agb_mass.push_back(max(mass_k16.back(), mass_d14.back()));
  mmin = agb_mass.front();
  mmax = agb_mass.back();

  // Get common list of masses
  set<double> mass_set;
  for (vector<int>::size_type i = 0; i < mass_k16.size(); i++)
      mass_set.insert(mass_k16[i]);
  for (vector<int>::size_type i = 0; i < mass_d14.size(); i++)
      mass_set.insert(mass_d14[i]);
  mass.assign(mass_set.begin(), mass_set.end());
  sort(mass.begin(), mass.end());

  // Construct combined yield table; our strategy is as follows: we
  // construct a union of the mass elements present in either
  // table. In the mass range covered by only one or the other table,
  // we just use the values from that table. In the mass range that
  // overlaps, if the isotope in question is present in both tables,
  // we use a weighted average of the table values, with the weighting
  // chosen so that the weights from from 0 to 1 over the overlapping
  // mass range. For isotopes present in only one table, we just use
  // that tables values.
  array2d::extent_gen extent;
  yield_tab.resize(extent[isotopes.size()][mass.size()]);
  for (vector<int>::size_type j=0; j<mass.size(); j++) {
    for (vector<int>::size_type i=0; i<isotopes.size(); i++) {
      if (mass[j] < mass_d14.front()) {
	
	// In mass range just covered by K16 table
	vector<const isotope_data *>::iterator isotope_it
	  = find(isotopes_k16.begin(), isotopes_k16.end(),
		 isotopes[i]);
	if (isotope_it == isotopes_k16.end()) continue; // Not in table
	vector<double>::iterator mass_it
	  = find(mass_k16.begin(), mass_k16.end(), mass[j]);
	vector<int>::size_type iso_idx
	  = std::distance(isotopes_k16.begin(), isotope_it);
	vector<int>::size_type mass_idx
	  = std::distance(mass_k16.begin(), mass_it);
	yield_tab[i][j] = yield_tab_k16[iso_idx][mass_idx];
	
      } else if (mass[j] > mass_k16.back()) {

	// In mass range just covered by D14 table
	vector<const isotope_data *>::iterator isotope_it
	  = find(isotopes_d14.begin(), isotopes_d14.end(),
		 isotopes[i]);
	if (isotope_it == isotopes_d14.end()) continue; // Not in table
	vector<double>::iterator mass_it
	  = find(mass_d14.begin(), mass_d14.end(), mass[j]);
	vector<int>::size_type iso_idx
	  = std::distance(isotopes_d14.begin(), isotope_it);
	vector<int>::size_type mass_idx
	  = std::distance(mass_d14.begin(), mass_it);
	yield_tab[i][j] = yield_tab_d14[iso_idx][mass_idx];

      } else {

	// In the overlapping mass range; get yields from each table;
	// this may require iterpolation, since the mass grids from
	// the two tables don't line up perfectly
	vector<const isotope_data *>::iterator iso_it_k16, iso_it_d14;
	iso_it_k16 = find(isotopes_k16.begin(), isotopes_k16.end(),
			  isotopes[i]);
	iso_it_d14 = find(isotopes_d14.begin(), isotopes_d14.end(),
			  isotopes[i]);
	double yld_k16, yld_d14;
	if (iso_it_k16 == isotopes_k16.end()) {
	  yld_k16 = -1.0;
	} else {
	  vector<double>::size_type iidx =
	    std::distance(isotopes_k16.begin(), iso_it_k16);
	  vector<double>::size_type midx = 0;
	  while (mass_k16[midx+1] <= mass[j]) {
	    midx++;
	    if (midx == mass_k16.size()) {
	      midx--;
	      break;
	    }
	  }
	  if (mass_k16[midx] == mass[j])
	    yld_k16 = yield_tab_k16[iidx][midx];
	  else {
	    double wgt = (mass_k16[midx+1] - mass[j]) /
	      (mass_k16[midx+1]  - mass_k16[midx]);
	    yld_k16 = wgt * yield_tab_k16[iidx][midx] +
	      (1.0 - wgt) * yield_tab_k16[iidx][midx+1];
	  }
	}
	if (iso_it_d14 == isotopes_d14.end()) {
	  yld_d14 = -1.0;
	} else {
	  vector<double>::size_type iidx =
	    std::distance(isotopes_d14.begin(), iso_it_d14);
	  vector<double>::size_type midx = 0;
	  while (mass_d14[midx+1] <= mass[j]) {
	    midx++;
	    if (midx == mass_d14.size()) {
	      midx--;
	      break;
	    }
	  }
	  if (mass_d14[midx] == mass[j])
	    yld_d14 = yield_tab_d14[iidx][midx];
	  else {
	    double wgt = (mass_d14[midx+1] - mass[j]) /
	      (mass_d14[midx+1]  - mass_d14[midx]);
	    yld_d14 = wgt * yield_tab_d14[iidx][midx] +
	      (1.0 - wgt) * yield_tab_d14[iidx][midx+1];
	  }
	}
      
	// Now construct averaged yield from two tables if both are
	// available; if not, just use whichever is available
	if (yld_k16 < 0.0) {
	  yield_tab[i][j] = yld_d14;
	} else if (yld_d14 < 0.0) {
	  yield_tab[i][j] = yld_k16;
	} else {
	  double wgt = (mass_k16.back() - mass[j]) /
	    (mass_k16.back() - mass_d14.front());
	  yield_tab[i][j] = wgt * yld_k16 +
	    (1.0 - wgt) * yld_d14;
	}
      }
    }
  }

  // Now build interpolation functions for each isotope
  yield_interp.resize(niso);
  yield_accel.resize(niso);
  for (vector<int>::size_type i=0; i<niso; i++) {
    yield_interp[i] = gsl_spline_alloc(
#if GSLVERSION == 2
				       gsl_interp_steffen,
#else
				       gsl_interp_akima,
#endif
				       mass.size());
    yield_accel[i] = gsl_interp_accel_alloc();
    vector<double> tmp(mass.size());
    for (vector<double>::size_type j=0; j<mass.size(); j++) {
      tmp[j] = yield_tab[i][j];
    }
    gsl_spline_init(yield_interp[i], mass.data(), tmp.data(), mass.size());
  }

}

////////////////////////////////////////////////////////////////////////
// Routine to read the Karakas & Lugaro 2016 yield tables
////////////////////////////////////////////////////////////////////////
void slug_yields_karakas16_doherty14::
read_karakas16_tables(const char *yield_dir,
		      const double pmz_mass,
		      const AGB_Karakas16_overshoot overshoot) {

  // Metallicities and file names available, hardcoded; note that
  // these are relative to Solar, while the file names are absolute
  const vector<double> z_file = { 0.5, 1.0, 0.03/0.014 };
  const vector<string> z_dir = { "z007models", "z014models", "z03models" };

  // Find nearest metallicity to input value, and get directory name
  // for it
  vector<double>::size_type file_idx = 0;
  for (vector<double>::size_type i=1; i<z_file.size(); i++)
    if (abs(log(metallicity/z_file[i])) <
	abs(log(metallicity/z_file[file_idx]))) file_idx = i;
  path yield_path(yield_dir);
  yield_path /= path("AGB_Karakas16");
  yield_path /= path(z_dir[file_idx]);
  
  // Issue warning if metallicity does not match one of inputs exactly
  if (metallicity != z_file[file_idx]) {
    ostreams.slug_warn_one
      << "Solar-normalized metallicity is "
      << metallicity
      << ", which does not exactly match any of the metallicities "
      << "computed by Karakas & Lugaro (2016); using closest match, "
      << "Z/Zsun = " << z_file[file_idx]
      << ". Computation will continue."
      << endl;
  }
  
  // Now read list of files available in the relevant metallicity
  // directory
  vector<yields::karakas16_file_data> yield_files_avail;
  for (directory_iterator it(yield_path); it != directory_iterator(); ++it) {
    try {

      // Record information about this file
      yields::karakas16_file_data yld_file;
      yld_file.fname = it->path().filename().string();
      
      // Get the mass of this file
      size_t ptr1, ptr2;
      ptr1 = yld_file.fname.find('z');
      if (ptr1 == string::npos) {
	ostreams.slug_err_one
	  << "unexpected file name "
	  << yld_file.fname
	  << " in directory " << yield_dir
	  << endl;
	bailout(1);
      }
      yld_file.mass = lexical_cast<double>(yld_file.fname.substr(1, ptr1-1));

      // Get the partial mixing zone mass
      ptr1 = yld_file.fname.find("pmz=");
      if (ptr1 == string::npos)
	yld_file.pmz_mass = 0.0;
      else {
	ptr2 = yld_file.fname.find(".", ptr1);
	if (ptr2 == string::npos) {
	  ostreams.slug_err_one
	    << "unexpected file name "
	    << yld_file.fname
	    << " in directory " << yield_dir
	    << endl;
	  bailout(1);
	}
	yld_file.pmz_mass
	  = lexical_cast<double>(yld_file.fname.substr(ptr1+4, ptr2-ptr1-4));
      }

      // Get overshooting type
      ptr1 = yld_file.fname.find("overshoot");
      ptr2 = yld_file.fname.find("noovershoot");
      if (ptr1 == string::npos)
	yld_file.overshoot = NO;
      else if (ptr2 == string::npos)
	yld_file.overshoot = YES;
      else
	yld_file.overshoot = NO;

      // Save data
      yield_files_avail.push_back(yld_file);
      
    } catch (const bad_lexical_cast& ia) {
      (void) ia; // No-op to suppress compiler warning
      ostreams.slug_err_one
	<< "unexpected file name "
	<< it->path().filename()
	<< " in directory " << yield_dir
	<< endl;
      bailout(1);
    }
  }

  // Sort the list of files by mass
  sort(yield_files_avail.begin(), yield_files_avail.end());

  // Now construct the list of which yield files to use at each mass
  vector<yields::karakas16_file_data> yield_files;
  vector<double>::size_type ptr1, ptr2;
  ptr1 = ptr2 = 0;
  while (ptr1 < yield_files_avail.size()) {

    // Find all the files that have the same mass
    ptr2 = ptr1 + 1;
    while (yield_files_avail[ptr1].mass == yield_files_avail[ptr2].mass) {
      ptr2++;
      if (ptr2 == yield_files_avail.size()) break;
    }

    // Record the mass
    mass_k16.push_back(yield_files_avail[ptr1].mass);

    // Is there more than one yield file available for this mass
    if (ptr2 == ptr1+1) {

      // No, so use the one available
      yield_files.push_back(yield_files_avail[ptr1]);
      
    } else {
      
      // Yes, so we need to decide which file to use; A. Karakas
      // recommends the following defaults:
      // pmz = 2e-3 Msun for M <= 3 Msun
      // pmz = 1e-3 Msun for 3 Msun < M <= 4 Msun
      // pmz = 1e-4 Msun for 4 Msun < M < 5 Msun
      // pmz = 0 for M >= 5 Msun
      // Overshooting turned on if available
      // However, we allow the user to override these defaults
      double pmz_mass_target;
      if (pmz_mass >= 0.0) pmz_mass_target = pmz_mass;
      else {
	if (mass_k16.back() <= 3.0) pmz_mass_target = 2.0e-3;
	else if (mass_k16.back() <= 4.0) pmz_mass_target = 1.0e-3;
	else if (mass_k16.back() < 5.0) pmz_mass_target = 1.0e-4;
	else pmz_mass_target = 0.0;
      }

      // Find all the files where the pmz_mass equals the target;
      // there may be more than one
      vector<yields::karakas16_file_data> candidates;
      double pmz_mass_diff = 1e10;
      for (vector<int>::size_type i=ptr1; i<ptr2; i++)
	pmz_mass_diff = min(pmz_mass_diff,
			    abs(yield_files_avail[i].pmz_mass -
				pmz_mass_target));
      for (vector<int>::size_type i=ptr1; i<ptr2; i++)
	if (pmz_mass_diff == abs(yield_files_avail[i].pmz_mass -
				 pmz_mass_target))
	  candidates.push_back(yield_files_avail[i]);	  

      // Now see if we have more than one candidate left
      if (candidates.size() == 1) {
	// One candidate, so that's it
	yield_files.push_back(candidates[0]);
      } else {
	// More than one candidate left, so we must have some
	// overshooting and some non-overshooting; choose overshooting
	// unless instructed otherwise
	vector<double>::size_type i;
	for (i=0; i<candidates.size(); i++) {
	  if ((candidates[i].overshoot == YES && overshoot == DEFAULT) ||
	      (candidates[i].overshoot == overshoot)) {
	    yield_files.push_back(candidates[i]);
	    break;
	  }
	}
	// We can get here only if we have a user-specified PMZ mass
	// that happens to be halfway between two available values,
	// neither of which has overshooting; in this case, just
	// choose one
	if (i == candidates.size())
	  yield_files.push_back(candidates[0]);
      }
    }

    // If the use requested a specific pmz mass or overshoot value
    // and we can't match it exactly, issue a warning
    if (pmz_mass >= 0.0 && yield_files.back().pmz_mass != pmz_mass) {
      ostreams.slug_warn_one
	<< "slug_yields_karakas16_doherty14: requested pmz mass = "
	<< pmz_mass
	<< ", but for AGB star of initial mass "
	<< yield_files.back().mass
	<< " closest pmz mass found was "
	<< yield_files.back().pmz_mass
	<< "; computation will continue with closest match"
	<< endl;
    }
    if (overshoot != DEFAULT &&
	yield_files.back().overshoot != overshoot) {
      ostreams.slug_warn_one
	<< "slug_yields_karakas16_doherty14: requested overshoot = "
	<< (overshoot == YES)
	<< ", but for AGB star of initial mass "
	<< yield_files.back().mass
	<< " only found a model with overshoot = "
	<< (yield_files.back().overshoot == YES)
	<< "; computation will continue with that model"
	<< endl;
    }

    // Advance pointers
    ptr1 = ptr2;
  }

  // We now have a list of masses and a list of files containing their
  // yields; now we can actually begin to read the data

  // Go through first file just to grab the list of isotopes
  std::ifstream yf_stream;
  path file_path = yield_path / path(yield_files[0].fname.c_str());
  yf_stream.open(file_path.c_str());
  if (!yf_stream.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open yield file " 
			  << yield_path.string() << endl;
    bailout(1);
  }

  // Toss out the first two lines; the first is a header, and the
  // second is the yield of pure neutrons
  vector<string> tokens;
  string line;
  getline(yf_stream, line);
  getline(yf_stream, line);

  // Loop over entries
  while (getline(yf_stream, line)) {

    // Extract isotope name
    trim(line);
    split(tokens, line, is_any_of("\t "), token_compress_on);

    // Grab isotope data for this isotope; handle special cases
    // where Karakas's naming convention differs from ours
    const isotope_data *iso;
    if (tokens[0] == "p") {
      iso = iso_table->data("h1");
    } else if (tokens[0] == "d") {
      iso = iso_table->data("h2");
    } else if (tokens[0] == "al-6") {
      iso = iso_table->data("al26");
    } else if (tokens[0] == "al*6") {
      continue;
    } else if (tokens[0] == "kr-5") {
      iso = iso_table->data("kr85");
    } else if (tokens[0] == "kr*5") {
      continue;
    } else {
      try {
	iso = iso_table->data(tokens[0]);
      } catch (const out_of_range& oor) {
	ostreams.slug_err_one
	  << "unrecognized isotope " << tokens[0]
	  << " in " << yield_files[0].fname
	  << endl;
	bailout(1);
      }
    }

    // Is this isotope name already in our list? If so, skip it and
    // do nothing. Otherwise push it onto our isotope list.
    vector<double>::size_type j;
    for (j=0; j<isotopes_k16.size(); j++) {
      if (iso == isotopes_k16[j])
	break;
    }
    if (j == isotopes_k16.size()) isotopes_k16.push_back(iso);
  }

  // Close file
  yf_stream.close();

  // Now that we know how many isotopes we have, build arrays to hold
  // yield data
  array2d::extent_gen extent;
  yield_tab_k16.resize(extent[isotopes_k16.size()][mass_k16.size()]);
  for (vector<double>::size_type i=0; i<isotopes_k16.size(); i++) {
    for (vector<double>::size_type j=0; j<mass_k16.size(); j++) {
      yield_tab_k16[i][j] = 0.0;
    }
  }

  // Build a map from isotrope string to index; we'll use this to fill
  // the table below
  sort(isotopes_k16.begin(), isotopes_k16.end(),
       slug_isotopes::isotope_sort);
  map<string, vector<double>::size_type> iso_map_str;
  for (vector<double>::size_type i=0; i<isotopes_k16.size(); i++) {
    stringstream ss;
    ss << isotopes_k16[i]->symbol() << isotopes_k16[i]->wgt();
    iso_map_str[ss.str()] = i;
  }
  
  // Now go through all files, grabbing yield data
  for (vector<double>::size_type i=0; i<mass_k16.size(); i++) {

    // Open file
    file_path = yield_path / path(yield_files[i].fname.c_str());
    yf_stream.open(file_path.c_str());
    if (!yf_stream.is_open()) {
      // Couldn't open file, so bail out
      ostreams.slug_err_one << "unable to open yield file " 
			    << yield_path.string() << endl;
      bailout(1);
    }

    // Burn 2 lines
    getline(yf_stream, line);
    getline(yf_stream, line);

    // Loop over entries
    while (getline(yf_stream, line)) {

      // Extract isotope name
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);

      // Handle special cases
      string iso_str;
      if (tokens[0] == "p") {
	iso_str = "h1";
      } else if (tokens[0] == "d") {
	iso_str = "h2";
      } else if (tokens[0] == "al-6") {
	iso_str = "al26";
      } else if (tokens[0] == "al*6") {
	continue;
      } else if (tokens[0] == "kr-5") {
	iso_str = "kr85";
      } else if (tokens[0] == "kr*5") {
	continue;
      } else {
	iso_str = tokens[0];
      }

      // Record yield; note that we put a safety in here to ensure
      // values are non-negative, because some entries in Karakas's
      // tables are due to numerical issues
      yield_tab_k16[iso_map_str.at(iso_str)][i]
	= max(0.0, lexical_cast<double>(tokens[3]));
    }

    // Close file
    yf_stream.close();
  }
}

////////////////////////////////////////////////////////////////////////
// Routine to read the Doherty+ 2014 yield tables
////////////////////////////////////////////////////////////////////////
void slug_yields_karakas16_doherty14::
read_doherty14_tables(const char *yield_dir) {

  // Construct name of file to read
  path yield_path(yield_dir);
  yield_path /= path("superAGB_Doherty14");
  yield_path /= path("doherty14a_table1.txt");

  // Open file
  std::ifstream yf_stream;
  yf_stream.open(yield_path.c_str());
  if (!yf_stream.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open yield file " 
			  << yield_path.string() << endl;
    bailout(1);
  }

  // First pass: read list of available masses, metallicities, and
  // isotopes
  string line;
  vector<string> tokens;
  vector<double> Z;
  vector< vector<double> > m;
  while (getline(yf_stream, line)) {

    // Break line into tokens
    trim(line);
    split(tokens, line, is_any_of("\t "), token_compress_on);

    // Does this line start with the word "species"? If so, it's a
    // header line, and can be discarded
    if (tokens[0] == "Species") continue;

    // Is the first token of the form numberM? If so, this line
    // declares start of a new mass and metallicity block, so extract
    // that information; note that we treat Z = 0.02 as "solar"
    if (tokens[0][tokens[0].length()-1] == 'M') {
      
      // Get m and Z of this block
      double m_block
	= lexical_cast<double>(tokens[0].substr(0, tokens[0].size()-1));
      double Z_block = lexical_cast<double>(tokens[1].substr(2)) / 0.02;

      // Find index in metallicity list, creating a new entry if one
      // corresponding to this metallicity does not yet exist
      vector<int>::size_type idx;
      if (Z.size() == 0) {
	Z.push_back(Z_block);
	vector<double> m_z;
	m_z.push_back(m_block);
	m.push_back(m_z);
	idx = 0;
      } else {
	idx = -1;
	for (vector<double>::size_type i=0; i<Z.size(); i++) {
	  if (Z[i] == Z_block) {
	    idx = i;
	    break;
	  } else if (Z[i] > Z_block) {
	    Z.insert(Z.begin()+i, Z_block);
	    vector<double> m_z;
	    m_z.push_back(m_block);
	    m.insert(m.begin()+i, m_z);
	    idx = i;
	    break;
	  }
	}
	if (idx == -1) {
	  Z.push_back(Z_block);
	  vector<double> m_z;
	  m_z.push_back(m_block);
	  m.push_back(m_z);
	  idx = Z.size()-1;
	} else {
	  m[idx].push_back(m_block);
	}
      }

    } else {

      // Any other line is a line declaring an isotope
      const isotope_data *iso;
      if (tokens[0] == "p") {
	iso = iso_table->data("h1");
      } else if (tokens[0] == "g") {
	continue;
      } else {
	try {
	  iso = iso_table->data(tokens[0]);
	} catch (const out_of_range& oor) {
	  ostreams.slug_err_one
	    << "unrecognized isotope " << tokens[0]
	    << " in " << yield_path.string()
	    << endl;
	  bailout(1);
	}
      }

      // Is this isotope name already in our list? If so, skip it and
      // do nothing. Otherwise push it onto our isotope list.
      vector<double>::size_type j;
      for (j=0; j<isotopes_d14.size(); j++) {
	if (iso == isotopes_d14[j])
	  break;
      }
      if (j == isotopes_d14.size()) isotopes_d14.push_back(iso);
    }

  }

  // Re-set file to beginning; clear EOF flag
  yf_stream.clear();
  yf_stream.seekg(0, ios::beg);

  // See if we have the requested metallicity; if not, choose closest
  // match
  vector<double>::size_type Z_idx;
  for (Z_idx = 0; Z_idx < Z.size(); Z_idx++) {
    if (metallicity == Z[Z_idx]) {
      break;
    } else if (metallicity < Z[Z_idx]) {
      if (Z_idx != 0) {
	if (log(metallicity/Z[Z_idx-1]) < log(Z[Z_idx]/metallicity)) {
	  Z_idx--;
	}
      }
      break;
    }
  }
  if (Z_idx == Z.size()) Z_idx--;
  if (Z[Z_idx] != metallicity) {
    ostreams.slug_warn_one
      << "Solar-normalized metallicity is "
      << metallicity
      << ", which does not exactly match any of the metallicities "
      << "computed by Doherty et al. (2014); using closest match, "
      << "Z/Zsun = " << Z[Z_idx]
      << ". Computation will continue."
      << endl;
  }

  // Load masses for this metallicity into mass list
  mass_d14 = m[Z_idx];

  // Initialize storage for yield table data
  array2d::extent_gen extent;
  yield_tab_d14.resize(extent[isotopes_d14.size()][mass_d14.size()]);
  for (vector<double>::size_type i=0; i<isotopes_d14.size(); i++) {
    for (vector<double>::size_type j=0; j<mass_d14.size(); j++) {
      yield_tab_d14[i][j] = 0.0;
    }
  }

  // Now go through file a second time, reading yields
  vector<double>::size_type m_idx = 0;
  while (getline(yf_stream, line)) {

    // Read metallicity of this block
    trim(line);
    split(tokens, line, is_any_of("\t "), token_compress_on);
    double Z_block = lexical_cast<double>(tokens[1].substr(2)) / 0.02;

    // If this is not the metallicity we want, skip rest of block
    if (Z_block != Z[Z_idx]) {
      for (vector<int>::size_type i=0; i<isotopes_d14.size()+1; i++)
	getline(yf_stream, line);
      continue;
    }

    // If we're here, this is a metallicity we want, so ingest the
    // data
    getline(yf_stream, line); // Burn one header line
    for (vector<int>::size_type i=0; i<isotopes_d14.size(); i++) {
      getline(yf_stream, line);
      trim(line);
      split(tokens, line, is_any_of("\t "), token_compress_on);
      yield_tab_d14[i][m_idx] = lexical_cast<double>(tokens[2]);
    }

    // Burn a line
    getline(yf_stream, line);

    // Increment mass index
    m_idx++;
  }

  // Close file
  yf_stream.close();
}

////////////////////////////////////////////////////////////////////////
// Functions to return yields of all isotopes or of a single isotope
////////////////////////////////////////////////////////////////////////

vector<double>
slug_yields_karakas16_doherty14::get_yield(const double m) const {
  vector<double> yld(niso);
  for (vector<double>::size_type i=0; i<niso; i++) {
    yld[i] = gsl_spline_eval(yield_interp[i], m, yield_accel[i]);
  }
  return yld;
}

double
slug_yields_karakas16_doherty14::get_yield(const double m,
				 const vector<double>::size_type i) const {
  return gsl_spline_eval(yield_interp[i], m, yield_accel[i]);
}
