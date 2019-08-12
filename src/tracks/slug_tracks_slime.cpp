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

#include "../constants.H"
#include "slug_tracks_slime.H"
extern "C" {
#include "fitsio.h"
}

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace boost::multi_array_types;

////////////////////////////////////////////////////////////////////////
// Little utility function to sort file names by metallicity
////////////////////////////////////////////////////////////////////////
namespace tracks {
  typedef struct {
    path fname;
    double log_Z;
  } Z_file_info;
  bool slime_filesort(const Z_file_info &f1, const Z_file_info &f2) {
    return f1.log_Z < f2.log_Z;
  }
}

////////////////////////////////////////////////////////////////////////
// Method to construct the tracks from a file name
////////////////////////////////////////////////////////////////////////
slug_tracks_slime::
slug_tracks_slime(const char *fname, slug_ostreams& ostreams_,
		 const bool force_mcur_monotonic_) :
  slug_tracks_2d(ostreams_, tracks::null_metallicity, Z_UNSPECIFIED,
		 force_mcur_monotonic_) {

  // Read header information from file so that we can size our data
  // arrays
  array1d logm;
  array2d logt;
  array3d trackdata;
  size_type ntime, ntrack;
  read_trackfile_header(fname, logm, ntime);
  ntrack = logm.size();
  logt.resize(boost::extents[ntime+1][ntrack]);
  trackdata.resize(boost::extents[ntime+1][ntrack][nprop]);

  // Read data
  read_trackfile(fname, logt, trackdata);

  // Set interpolation type
  vector<const gsl_interp_type *> interp_type(nprop);
  for (vector<int>::size_type i=0; i<nprop; i++) {
    if ((i == idx_log_cur_mass && force_mcur_monotonic) || i == idx_phase)
      interp_type[i] = gsl_interp_linear;
    else
      interp_type[i] = slug_default_interpolator;
  }
  
  // Build the interpolation class that will interpolate on the tracks
  interp = new
    slug_mesh2d_interpolator_vec(logt, logm, trackdata,
				 interp_type, force_mcur_monotonic);

}


////////////////////////////////////////////////////////////////////////
// Method to initialize track data by interpolating in metallicity
////////////////////////////////////////////////////////////////////////

void slug_tracks_slime::init_from_files(vector<path>& files,
					vector<double>& log_Z) {
  
  // Sort input files by metallicity
  vector<tracks::Z_file_info> file_info(log_Z.size());
  for (vector<double>::size_type i=0; i<log_Z.size(); i++) {
    file_info[i].fname = files[i];
    file_info[i].log_Z = log_Z[i];
  }
  std::sort(file_info.begin(), file_info.end(), tracks::slime_filesort);

  // Store file names and metallicities for future reference
  filenames.resize(log_Z.size());
  Z_files.resize(log_Z.size());
  for (vector<double>::size_type i=0; i<log_Z.size(); i++) {
    filenames[i] = file_info[i].fname.string();
    Z_files[i] = pow(10., file_info[i].log_Z);
    files[i] = file_info[i].fname;
    log_Z[i] = file_info[i].log_Z;
  }

  // Make sure metallicity is in the covered range; if not, bail out
  double log_Z_in = log10(metallicity);
  if (log_Z_in < log_Z.front() ||
      log_Z_in > log_Z.back()) {
    ostreams.slug_err_one
      << "slug_tracks_slime: requested metallicity " << metallicity
      << " is outside range of Z = " << pow(10.0, log_Z.front())
      << " - " << pow(10.0, log_Z.back())
      << " covered by requested track set"
      << endl;
    bailout(1);
  }

  // Get index of best metallicity match to requested value
  double Z_dist = constants::big;
  vector<double>::size_type Z_idx;
  for (vector<double>::size_type i=0; i<log_Z.size(); i++) {
    double dlogZ = abs(log_Z_in - log_Z[i]);
    if (dlogZ < Z_dist) {
      Z_dist = dlogZ;
      Z_idx = i;
    }
  }

  // Data holders
  array1d logm;
  array2d logt;
  array3d trackdata;
  size_type ntime, ntrack;

  // Handle separately cases where we need to read only a single data
  // file (because interpolation method is nearest neighbor of input
  // exactly matches available value) and cases where we need to read
  // and interpolate.
  if (Z_int_meth == Z_NEAR_NEIGHBOR || Z_dist == 0.0) {

    // Exact match or nearest neighbor interpolation; just read the
    // appropriate file

    // Read file header and size data holders
    read_trackfile_header(files[Z_idx].c_str(), logm, ntime);
    ntrack = logm.size();
    logt.resize(boost::extents[ntime+1][ntrack]);
    trackdata.resize(boost::extents[ntime+1][ntrack][nprop]);

    // Read file data
    read_trackfile(files[Z_idx].c_str(), logt, trackdata);
    
  } else {

    // In this case we need to read all available metallicities so
    // that we can interpolate on them. In principle we could reduce
    // this to just reading two files in the case of linear
    // interpolation, but since this is all a one-time operation that
    // will be carried out during startup, it's not worth optimizing
    // to that level.

    // Read through all headers to get number of times in each file
    size_type nZ = files.size();
    vector<size_type> ntime_Z(nZ);
    ntime = 0;
    for (vector<double>::size_type i=0; i<nZ; i++) {
      read_trackfile_header(files[i].c_str(), logm, ntime_Z[i]);
      ntime = max(ntime, ntime_Z[i]);
    }

    // Make appropriately-sized arrays to hold all data from all
    // metallicities
    array3d logt_Z;
    array4d trackdata_Z;
    ntrack = logm.size();
    logt_Z.resize(boost::extents[nZ][ntime+1][ntrack]);
    trackdata_Z.resize(boost::extents[nZ][ntime+1][ntrack][nprop]);

    // Loop over metallicities and read into array
    for (vector<double>::size_type i=0; i<nZ; i++) {
      view2d logt_tmp
	= logt_Z[indices[i][range_t(0,ntime_Z[i]+1)]
		 [range_t(0,ntrack)]];
      view3d trackdata_tmp
	= trackdata_Z[indices[i][range_t(0,ntime_Z[i]+1)]
		      [range_t(0,ntrack)][range_t(0,nprop)]];
      read_trackfile(files[i].c_str(), logt_tmp, trackdata_tmp);
    }

    // For metallicities that have fewer than the maximum number of
    // data points, pad by replicating the final time to all further
    // points
    for (vector<double>::size_type i=0; i<nZ; i++) {
      for (size_type j=ntime_Z[i]; j<ntime+1; j++) {
	for (size_type k=0; k<ntrack; k++) {
	  logt_Z[i][j][k] = logt_Z[i][ntime_Z[i]][k];
	  for (size_type n=0; n<nprop; n++) {
	    trackdata_Z[i][j][k][n] = trackdata_Z[i][ntime_Z[i]][k][n];
	  }
	}
      }
    }

    // Set the interpolation type
    const gsl_interp_type* Z_interp_type;
    switch (Z_int_meth) {
    case Z_LINEAR: {
      Z_interp_type = gsl_interp_linear;
      break;
    }
    case Z_AKIMA: {
      if (nZ >= gsl_interp_type_min_size(gsl_interp_akima))
	Z_interp_type = gsl_interp_akima;
      else
	Z_interp_type = gsl_interp_linear;
      break;
    }
#if GSLVERSION >= 2
    case Z_STEFFEN: {
      if (nZ >= gsl_interp_type_min_size(gsl_interp_steffen))
	Z_interp_type = gsl_interp_steffen;
      else
	Z_interp_type = gsl_interp_linear;
      break;
    }
#endif
    default: {
      ostreams.slug_err_one << "slug_tracks_mist constructor invoked with"
			    << " invalid Z interpolation method" << endl;
      bailout(1);
    }
    }

    // Get weight for linear interpolations; we use linear
    // interpolation for times, masses, and phases, because we require
    // that the interpolated grid be strictly monotone in these
    // quantities, and only linear interpolation guarantees this
    if (log_Z_in < log_Z[Z_idx]) Z_idx--;
    double wgt =
      (log_Z[Z_idx+1] - log_Z_in) / (log_Z[Z_idx+1] - log_Z[Z_idx]);
    
    // Allocate storage to hold final interpolated tracks
    logt.resize(boost::extents[ntime+1][ntrack]);
    trackdata.resize(boost::extents[ntime+1][ntrack][nprop]);
    
    // Interpolate times and all track data; use linear interpolation
    // for times, current mass, and phase, since these must be
    // preserve strict ordering, and whatever method was specified by
    // the command line argument for interpolating all other
    // quantities
    for (size_type i=0; i<ntime+1; i++) {
      for (size_type j=0; j<ntrack; j++) {
	logt[i][j] = wgt*logt_Z[Z_idx][i][j] +
	  (1.0-wgt)*logt_Z[Z_idx+1][i][j];
	for (size_type n=0; n<nprop; n++) {
	  if (n == idx_log_cur_mass || n == idx_phase) {
	    trackdata[i][j][n] = wgt * trackdata_Z[Z_idx][i][j][n] +
	      (1.0-wgt) * trackdata_Z[Z_idx+1][i][j][n];
	  } else {	
	    vector<double> data_tmp(nZ);
	    for (size_type k=0; k<nZ; k++)
	      data_tmp[k] = trackdata_Z[k][i][j][n];
	    gsl_interp *interp_tmp = gsl_interp_alloc(Z_interp_type, nZ);
	    gsl_interp_init(interp_tmp, log_Z.data(), data_tmp.data(), nZ);
	    trackdata[i][j][n] =
	      gsl_interp_eval(interp_tmp, log_Z.data(), data_tmp.data(),
			      log_Z_in, nullptr);
	    gsl_interp_free(interp_tmp);
	  }
	}
      }
    }
  }

  // Set interpolation type
  vector<const gsl_interp_type *> interp_type(nprop);
  for (vector<int>::size_type i=0; i<nprop; i++) {
    if ((i == idx_log_cur_mass && force_mcur_monotonic) || i == idx_phase)
      interp_type[i] = gsl_interp_linear;
    else
      interp_type[i] = slug_default_interpolator;
  }
  
  // Build the interpolation class that will interpolate on the tracks
  interp = new
    slug_mesh2d_interpolator_vec(logt, logm, trackdata,
				 interp_type, force_mcur_monotonic);
}

////////////////////////////////////////////////////////////////////////
// Method to read the header of a track file
////////////////////////////////////////////////////////////////////////

void
slug_tracks_slime::read_trackfile_header(const char *fname,
					 array1d& logm,
					 size_type& ntime) {

  // Open file
  fitsfile *fptr;
  int fits_status = 0;
  fits_open_table(&fptr, fname, READONLY, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to open file " << fname << std::endl;
    bailout(1);
  }

  // Read the number of masses in the file
  long ntrack;
  fits_get_num_rows(fptr, &ntrack, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to read data from file "
			  << fname << std::endl;
    bailout(1);
  }

  // Resize the logm array to hold the data
  logm.resize(boost::extents[ntrack]);

  // Read the masses; we read this into a flat c array first, then
  // transfer them to the logm array
  double *m = new double[ntrack];
  int anynul;
  fits_read_col(fptr, TDOUBLE, 1, 1, 1, ntrack, nullptr, m, &anynul,
		&fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "read data from file "
			  << fname << std::endl;
    bailout(1);
  }
  for (long i=0; i<ntrack; i++) logm[i] = log(m[i]);
  delete[] m;

  // Move to the next HDU so we can see how many times there are
  int hdutype;
  fits_movrel_hdu(fptr, 1, &hdutype, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "read data from file " << fname << std::endl;
    bailout(1);
  }

  // Read number of times
  long ntime_fits;  // Read into a long, since that is what FITS will
		    // give us, and this is not the same type as the
		    // ntimes we want to return
  fits_get_num_rows(fptr, &ntime_fits, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to read data from file " << fname
			  << std::endl;
    bailout(1);
  }
  ntime = ntime_fits;

  // Close file
  fits_close_file(fptr, &fits_status);
}


////////////////////////////////////////////////////////////////////////
//
// We need to force the compiler to instantiate the template methods
// here, so that they will be available at link time. We do this just
// by having a dummy function that instantiates the class with the
// required template methods.
//
////////////////////////////////////////////////////////////////////////



#endif
// ENABLE_FITS
