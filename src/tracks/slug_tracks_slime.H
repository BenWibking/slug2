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
// This class holds a set of evolutionary tracks in the format
// produced by the SLIME (SLug Ingestion of MESA Evolution)
// script. This script was originally written to process MESA outputs,
// but has since been generalised to process other track sets. The
// parts of the parsing that are generic to all track sets are located
// in this class, and further classes that handle individual track
// sets are derived from it. This is a pure abstract class.
//
// Note that this class makes use of FITS capability, and is disabled
// if FITS is not included
////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_FITS

#ifndef _slug_tracks_slime_H_
#define _slug_tracks_slime_H_

#include <cstdio>
#include <boost/filesystem.hpp>
#include "slug_tracks_2d.H"
#include "../constants.H"
#include "../slug_MPI.H"
extern "C" {
#include "fitsio.h"
}

class slug_tracks_slime : public slug_tracks_2d {

public:
  
  // Construct a set of tracks from a specified file name
  slug_tracks_slime(const char *fname, slug_ostreams& ostreams_,
		    const bool force_mcur_monotonic_ = false);

  // Placeholder constructor -- this is a generic constructor that
  // just initiliazes the parent slug_tracks_2d class with
  // placeholders, so that we can fill them in later
  slug_tracks_slime(slug_ostreams& ostreams_,
		    double metallicity_ = tracks::null_metallicity,
		    ZInterpMethod Z_int_meth_ = Z_UNSPECIFIED,
		    bool force_mcur_monotonic_ = false) :
    slug_tracks_2d(ostreams_, metallicity_, Z_int_meth_,
		   force_mcur_monotonic_) { };

  // Destructor
  virtual ~slug_tracks_slime() {};

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
  const char* const fits_colnames[10] =
    { "mass", "log_L", "log_Teff", "mdot",
      "h_surf", "he_surf", "c_surf", "n_surf",
      "o_surf", "phase" };

  // Method to obtain track data by interpolating in metallicity; this
  // is implemented as a protected method rather than a constructor
  // because we usually need to do the construction of the derived
  // class first in order to get the list of file names and
  // metallicities, and only then do we have the information we need
  // to invoke this method. The inputs are the list of files and their
  // corresponding log metallicities. By the time this method is
  // invoked, the desired metallicity and interpolation types should
  // already be stored. A side effect of invoking this method is that
  // the files and log_Z values will be sorted in order of increasing
  // metallicity.
  void init_from_files(std::vector<boost::filesystem::path> &files,
		       std::vector<double>& log_Z);
  
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


////////////////////////////////////////////////////////////////////////
//
// We define the methods here, rather than in a separate cpp file,
// because read_trackfile is a template method, which needs to be
// included full in the source of every derived class that is going to
// make use of it in order not to cause linking issues
//
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Method to read the data from a track file
////////////////////////////////////////////////////////////////////////

template<typename S, typename T> void
slug_tracks_slime::read_trackfile(const char *fname, S& logt,
				  T& trackdata) {

  // Store number of masses and tracks for convenience
  boost::multi_array_types::size_type ntime = trackdata.shape()[0]-1;
  boost::multi_array_types::size_type ntrack = trackdata.shape()[1];

  // Open file
  fitsfile *fptr;
  int fits_status = 0;
  fits_open_table(&fptr, fname, READONLY, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to open file "
			  << fname << std::endl;
    bailout(1);
  }
  
  // Loop over masses
  double *buf = new double[ntime];
  for (boost::multi_array_types::size_type i=0; i<ntrack; i++) {

    // Move to next HDU
    int hdutype;
    fits_movrel_hdu(fptr, 1, &hdutype, &fits_status);
    if (fits_status) {
      ostreams.slug_err_one << "unable to read data from file "
			    << fname << std::endl;
      bailout(1);
    }

    // Read times; note that times stored in the file are linear, but
    // for internal purposes we use log
    int anynul, colnum;
    char colname[] = "age";
    fits_get_colnum(fptr, CASEINSEN, colname, &colnum,
		    &fits_status);
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ntime, nullptr, buf, &anynul,
		  &fits_status);
    if (fits_status) {
      ostreams.slug_err_one << "unable to read data from file "
			    << fname << std::endl;
      bailout(1);
    }
    for (boost::multi_array_types::size_type j=0; j<ntime; j++)
      logt[j+1][i] = log(buf[j]);

    // Read data; again, note that mass and mdot are stored as linear
    // values in the MIST tracks, but we use log values internally in
    // slug
    for (boost::multi_array_types::size_type n=0; n<nprop; n++) {
      char tmp_colname[10];
      strcpy(tmp_colname, fits_colnames[n]); // Needed to avoid const issues
      fits_get_colnum(fptr, CASEINSEN, tmp_colname,
		      &colnum, &fits_status);
      fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ntime, nullptr, buf,
		    &anynul, &fits_status);
      if (fits_status) {
	ostreams.slug_err_one << "unable to read data from file "
			      << fname << std::endl;
	bailout(1);
      }
      if (n == idx_log_cur_mass || n == idx_log_mDot) {
	for (boost::multi_array_types::size_type j=0; j<ntime; j++)
	  trackdata[j+1][i][n] = log(buf[j]);
      } else {
	for (boost::multi_array_types::size_type j=0; j<ntime; j++)
	  trackdata[j+1][i][n] = buf[j];
      }
    }
  }
  delete[] buf;

  // Close file
  fits_close_file(fptr, &fits_status);

  // Populate the dummy row at time 0. We add this row to avoid
  // running into problems if we try to interpolate to very young
  // ages.
  for (boost::multi_array_types::size_type i=0; i<ntrack; i++) {
    logt[0][i] = -log10(constants::big);
    for (boost::multi_array_types::size_type j=0;
	 j<trackdata.shape()[2]; j++) {
      trackdata[0][i][j] = trackdata[1][i][j];
    }
  }
}

#endif
// _slug_tracks_slime_H_

#endif
// _ENABLE_FITS_
