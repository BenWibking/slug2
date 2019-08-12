/*********************************************************************
Copyright (C) 2017 Michele Fumagalli, Mark Krumholz
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

#ifdef ENABLE_FITS

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif

#include <cmath>
#include <cstring>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include "../constants.H"
#include "../slug_MPI.H"
#include "slug_tracks_mist.H"
extern "C" {
#include "fitsio.h"
}

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace boost::multi_array_types;


////////////////////////////////////////////////////////////////////////
// Constructor from file name; note that parent class does the actual
// file reading here; we just get the metallicity data from the file
// name.
////////////////////////////////////////////////////////////////////////

slug_tracks_mist::
slug_tracks_mist(const char *fname, slug_ostreams& ostreams_,
		 const bool force_mcur_monotonic_) :
  slug_tracks_slime(fname,
		    ostreams_,
		    force_mcur_monotonic_) {

  // Read metallicity information from file name
  const regex Z_pattern("feh_[pm][0-9].[0-9][0-9]");
  match_results<std::basic_string<char>::iterator> Z_match;
  string fname_str(fname);
  if (regex_search(fname_str.begin(),
		   fname_str.end(),
		   Z_match, Z_pattern, match_posix)) {
    string Z_str(Z_match[0].first, Z_match[0].second);
    Z_str = Z_str.substr(4);
    if (Z_str[0] == 'p') Z_str[0] = '+';
    else Z_str[0] = '-';
    metallicity = pow(10.0, lexical_cast<double>(Z_str));
  } else {
    ostreams.slug_err_one << "slug_tracks_mist: could not determine "
			  << "metallicity from specified file name"
			  << endl;
    bailout(1);
  }
  filenames.push_back(fname_str);
  Z_files.push_back(metallicity);
}

////////////////////////////////////////////////////////////////////////
// Constructor from track set
////////////////////////////////////////////////////////////////////////

slug_tracks_mist::
slug_tracks_mist(const trackSet tr_set,
		 const double metallicity_,
		 const char *track_dir,
		 slug_ostreams& ostreams_,
		 const ZInterpMethod Z_int_meth_,
		 const bool force_mcur_monotonic_) :
  slug_tracks_slime(ostreams_, metallicity_, Z_int_meth_,
		    force_mcur_monotonic_) {
  
  // Select the subdirectory based on the input track set
  path subdir(track_dir);
  subdir /= path("mist");
  switch (tr_set) {
  case MIST_2016_VVCRIT_00: {
    subdir /= path("vvcrit000");
    break;
  }
  case MIST_2016_VVCRIT_40: {
    subdir /= path("vvcrit040");
    break;
  }
  default: {
    ostreams.slug_err_one
      << "slug_tracks_mist constructor invoked with "
      << "non-MIST track set!" << endl;
    bailout(1);
  }
  }

  // Iterate through files in the subdirectory to figure out the
  // metallicities available
  vector<path> files;
  vector<double> log_Z;
  try {
    for (directory_iterator dit(subdir); dit != directory_iterator();
	 dit++) {
      files.push_back(dit->path());
      const regex Z_pattern("feh_[pm][0-9].[0-9][0-9]");
      match_results<std::basic_string<char>::iterator> Z_match;
      string fname = dit->path().filename().string();
      if (regex_search(fname.begin(),
		       fname.end(),
		       Z_match, Z_pattern, match_posix)) {
	string Z_str(Z_match[0].first, Z_match[0].second);
	Z_str = Z_str.substr(4);
	if (Z_str[0] == 'p') Z_str[0] = '+';
	else Z_str[0] = '-';
	log_Z.push_back(lexical_cast<double>(Z_str));
      }
    }
  } catch (const filesystem_error& ex) {
    ostreams.slug_err_one << "could not read directory "
			  << subdir
			  << "; boost::filesystem says "
			  << ex.what()
			  << endl;
    bailout(1);
  }

  // Call the data initializer from parent to interpolate the tracks
  // in metallicity for us
  init_from_files(files, log_Z);
}


////////////////////////////////////////////////////////////////////////
// Methods to determine if stars are WR stars, and, if so, what type;
// these methods assume that the isochrone is properly set for this
// age.
////////////////////////////////////////////////////////////////////////
void slug_tracks_mist::
set_WR_type(const double m, const double t, slug_stardata& star) const {

  // Check phase
  double logm = log(m);
  double logt = log(t);
  double phase = (*interp)(logt, logm, idx_phase);
  if (phase < 8.5) {
    star.WR = NONE;
    return;
  }

  // Star is a WR star; determine type based on surface abundance
  // fractions
  double H_frac = (*interp)(logt, logm, idx_h_surf);
  if (H_frac > 0.1) {
    star.WR = WN;
    return;
  }
  double C_frac = (*interp)(logt, logm, idx_c_surf);
  double N_frac = (*interp)(logt, logm, idx_n_surf);
  if (C_frac/(N_frac+constants::small) < 10.0) {
    star.WR = WN;
  } else {
    star.WR = WC;
  }  
}

void
slug_tracks_mist::set_WR_type(const double m, 
			      spl_arr_view_1d& isochrone_,
			      acc_arr_view_1d& isochrone_acc_,
			      slug_stardata& star) const {

  // Check phase
  double logm = log(m);
  double phase = gsl_spline_eval(isochrone_[idx_phase],
				 logm,
				 isochrone_acc_[idx_phase]);
  if (phase < 8.5) {
    star.WR = NONE;
    return;
  }

  // If phase is WR, decide type based on surface abundances
  double H_frac = gsl_spline_eval(isochrone_[idx_h_surf],
				  logm,
				  isochrone_acc_[idx_h_surf]);
  if (H_frac > 0.1) {
    star.WR = WN;
    return;
  }
  double C_frac = gsl_spline_eval(isochrone_[idx_c_surf],
				  logm,
				  isochrone_acc_[idx_c_surf]);
  double N_frac = gsl_spline_eval(isochrone_[idx_n_surf],
				  logm,
				  isochrone_acc_[idx_n_surf]);
  if (C_frac/(N_frac+constants::small) < 10.0) {
    star.WR = WN;
  } else {
    star.WR = WC;
  }
}


#endif
// ENABLE_FITS
