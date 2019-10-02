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
#include "slug_tracks_slime.H"

class slug_tracks_mist : public slug_tracks_slime {

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

#ifdef WINDS_ON
  // Method to set stellar wind properties; this is a pure virtual
  // that must be implemented in a derived class
  virtual void set_wind(const double m, const double t,
			slug_stardata& star) const;
  virtual void set_wind(const double m,
			spl_arr_view_1d& isochrone_,
			acc_arr_view_1d& isochrone_acc_,
			slug_stardata& star) const;
#endif
  
};

#endif
// _slug_tracks_mist_H_

#endif
// _ENABLE_FITS_
