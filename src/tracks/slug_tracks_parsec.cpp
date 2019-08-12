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
#include "slug_tracks_parsec.H"
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

slug_tracks_parsec::
slug_tracks_parsec(const char *fname, slug_ostreams& ostreams_,
		   const bool force_mcur_monotonic_) :
  slug_tracks_slime(ostreams_, tracks::null_metallicity, Z_UNSPECIFIED,
		    force_mcur_monotonic_) {

  // Get metallicity information from file name
  const regex Z_pattern("(?<=Z)(.*?)(?=Y)");
  match_results<std::basic_string<char>::iterator> Z_match;
  string fname_str(fname);
  if (regex_search(fname_str.begin(),
		   fname_str.end(),
		   Z_match, Z_pattern, match_posix)) {
    string Z_str(Z_match[0].first, Z_match[0].second);
    metallicity = lexical_cast<double>(Z_str);
  } else {
    ostreams.slug_err_one << "slug_tracks_parsec: could not determine "
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

slug_tracks_parsec::
slug_tracks_parsec(const double metallicity_,
		   const char *track_dir,
		   slug_ostreams& ostreams_,
		   const ZInterpMethod Z_int_meth_,
		   const bool force_mcur_monotonic_) :
  slug_tracks_slime(ostreams_, metallicity_, Z_int_meth_,
		    force_mcur_monotonic_) {

  // Set directory to where PARSEC tracks are located
  path subdir(track_dir);
  subdir /= path("parsec") / path("v1.2");

  // Iterate through files in the subdirectory to figure out the
  // metallicities available
  vector<path> files;
  vector<double> log_Z;
  try {
    for (directory_iterator dit(subdir); dit != directory_iterator();
	 dit++) {
      files.push_back(dit->path());
      const regex Z_pattern("(?<=Z)(.*?)(?=Y)");
      match_results<std::basic_string<char>::iterator> Z_match;
      string fname = dit->path().filename().string();
      if (regex_search(fname.begin(),
		       fname.end(),
		       Z_match, Z_pattern, match_posix)) {
	string Z_str(Z_match[0].first, Z_match[0].second);
	// N.B. We take "Solar metallicity" for this track set to be 0.014
	log_Z.push_back(log10(lexical_cast<double>(Z_str) / 0.014));
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
void slug_tracks_parsec::
set_WR_type(const double m, const double t, slug_stardata& star) const {
  // Placeholder
  star.WR = NONE;
}

void
slug_tracks_parsec::set_WR_type(const double m, 
				spl_arr_view_1d& isochrone_,
				acc_arr_view_1d& isochrone_acc_,
				slug_stardata& star) const {
  // Placeholder
  star.WR = NONE;
}  
#endif
