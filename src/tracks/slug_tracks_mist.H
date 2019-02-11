/*********************************************************************
Copyright (C) 2017 Robert da Silva, Michele Fumagalli, Mark Krumholz
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

////////////////////////////////////////////////////////////////////////
// class slug_tracks_mist
//
// This class holds a set of evolutionary tracks structured following
// the format used in the MIST database. It is a derived class
// of slug_tracks, which defines a general interface to evolutionary
// tracks.
//
// Note that this class makes use of FITS capability, and is disabled
// if FITS is not included
////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_FITS

#ifndef _slug_tracks_mist_H_
#define _slug_tracks_mist_H_

#include <cstdio>
#include "slug_tracks_2d.H"

class slug_tracks_mist : public slug_tracks_2d {

public:

  // Construct a set of tracks from a specified file name
  slug_tracks_mist(const char *fname, slug_ostreams& ostreams_,
		   const bool force_mcur_monotonic_ = false);

  // Construct a set of tracks from a specified track set
  slug_tracks_mist(const trackSet tr_set,
		   const double metallicity_,
		   const char *track_dir,
		   slug_ostreams& ostreams_,
		   const ZInterpMethod Z_int_meth_ = Z_STEFFEN,
		   const bool force_mcur_monotonic_= false);

  // Destructor
  virtual ~slug_tracks_mist() { };

  // Method to decide of a star is a Wolf-Rayet star, and which type
  // of WR star, for these tracks
  virtual void set_WR_type(const double m, const double t,
			   slug_stardata& star) const;
  virtual void set_WR_type(const double m,
			   spl_arr_view_1d& isochrone_,
			   acc_arr_view_1d& isochrone_acc_,
			   slug_stardata& star) const;

  
protected:

  // Number of properties in these tracks, and indices of variables
  static const boost::multi_array_types::size_type nprop = 10,
    idx_h_surf = 4,
    idx_he_surf = 5,
    idx_c_surf = 6,
    idx_n_surf = 7,
    idx_o_surf = 8,
    idx_phase = 9;

  // Names of columns in the FITS file that go with each property;
  // initialized in the cpp file
  static const char* const fits_colnames[];

private:

  // Method to read the header of a track file; this just extracts the
  // list of initial masses and time points in contains, then closes
  // the file
  void read_trackfile_header(const char *fname,
			     array1d& logm,
			     boost::multi_array_types::size_type& ntime);

  // Method to read a track file; templated so that it can accept
  // either arrays or array views -- S must be a 2D arraylike object,
  // and T must be a 3D arraylike object, and both must be the correct
  // shape at the time they are passed to this routine
  template<typename S, typename T>
  void read_trackfile(const char *fname, S& logt, T& trackdata);

};

#endif
// _slug_tracks_mist_H_

#endif
// _ENABLE_FITS_
