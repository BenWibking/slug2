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
// Little utility function to sort file names by metallicity
////////////////////////////////////////////////////////////////////////
namespace tracks {
  typedef struct {
    path fname;
    double log_Z;
  } Z_file_info;
  bool mist_filesort(const Z_file_info &f1, const Z_file_info &f2) {
    return f1.log_Z < f2.log_Z;
  }
}

////////////////////////////////////////////////////////////////////////
// Column names in FITS files
////////////////////////////////////////////////////////////////////////
const char* const slug_tracks_mist::fits_colnames[] =
  { "mass", "log_L", "log_Teff", "mdot",
    "h_surf", "he_surf", "c_surf", "n_surf",
    "o_surf", "phase" };

////////////////////////////////////////////////////////////////////////
// Constructor from file name
////////////////////////////////////////////////////////////////////////

slug_tracks_mist::
slug_tracks_mist(const char *fname, slug_ostreams& ostreams_,
		 const bool force_mcur_monotonic_) :
  slug_tracks_2d(ostreams_, tracks::null_metallicity, Z_UNSPECIFIED,
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
// Constructor from track set
////////////////////////////////////////////////////////////////////////

slug_tracks_mist::
slug_tracks_mist(const trackSet tr_set,
		 const double metallicity_,
		 const char *track_dir,
		 slug_ostreams& ostreams_,
		 const ZInterpMethod Z_int_meth_,
		 const bool force_mcur_monotonic_) :
  slug_tracks_2d(ostreams_, metallicity_, Z_int_meth_,
		 force_mcur_monotonic_) {

  // Set interpolation type
  vector<const gsl_interp_type *> interp_type(nprop);
  for (vector<int>::size_type i=0; i<nprop; i++) {
    if ((i == idx_log_cur_mass && force_mcur_monotonic) || i == idx_phase)
      interp_type[i] = gsl_interp_linear;
    else
      interp_type[i] = slug_default_interpolator;
  }
  
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
      filenames.push_back(fname);
      if (regex_search(fname.begin(),
		       fname.end(),
		       Z_match, Z_pattern, match_posix)) {
	string Z_str(Z_match[0].first, Z_match[0].second);
	Z_str = Z_str.substr(4);
	if (Z_str[0] == 'p') Z_str[0] = '+';
	else Z_str[0] = '-';
	log_Z.push_back(lexical_cast<double>(Z_str));
	Z_files.push_back(pow(10.0, log_Z.back()));
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

  // Sort files by metallicity
  vector<tracks::Z_file_info> file_info(log_Z.size());
  for (vector<double>::size_type i=0; i<log_Z.size(); i++) {
    file_info[i].fname = files[i];
    file_info[i].log_Z = log_Z[i];
  }
  std::sort(file_info.begin(), file_info.end(), tracks::mist_filesort);
  for (vector<double>::size_type i=0; i<log_Z.size(); i++) {
    files[i] = file_info[i].fname;
    log_Z[i] = file_info[i].log_Z;
    filenames[i] = file_info[i].fname.string();
    Z_files[i] = pow(10.0, log_Z[i]);
  }

  // Make sure metallicity is in the covered range; if not, bail out
  double log_Z_in = log10(metallicity);
  if (log_Z_in < log_Z.front() || log_Z_in > log_Z.back()) {
    ostreams.slug_err_one
      << "slug_tracks_mist: requested metallicity " << metallicity
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
    size_type nZ = log_Z.size();
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

  // Build the interpolation class that will interpolate on the tracks
  interp = new
    slug_mesh2d_interpolator_vec(logt, logm, trackdata,
				 interp_type, force_mcur_monotonic);
}

////////////////////////////////////////////////////////////////////////
// Method to read the header of a track file
////////////////////////////////////////////////////////////////////////

void
slug_tracks_mist::read_trackfile_header(const char *fname,
					array1d& logm,
					size_type& ntime) {

  // Open file
  fitsfile *fptr;
  int fits_status = 0;
  fits_open_table(&fptr, fname, READONLY, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to open file " << fname << endl;
    bailout(1);
  }

  // Read the number of masses in the file
  long ntrack;
  fits_get_num_rows(fptr, &ntrack, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to read data from file "
			  << fname << endl;
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
			  << fname << endl;
    bailout(1);
  }
  for (long i=0; i<ntrack; i++) logm[i] = log(m[i]);
  delete[] m;

  // Move to the next HDU so we can see how many times there are
  int hdutype;
  fits_movrel_hdu(fptr, 1, &hdutype, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "read data from file " << fname << endl;
    bailout(1);
  }

  // Read number of times
  long ntime_fits;  // Read into a long, since that is what FITS will
		    // give us, and this is not the same type as the
		    // ntimes we want to return
  fits_get_num_rows(fptr, &ntime_fits, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to read data from file " << fname << endl;
    bailout(1);
  }
  ntime = ntime_fits;

  // Close file
  fits_close_file(fptr, &fits_status);
}

////////////////////////////////////////////////////////////////////////
// Method to read the data from a track file
////////////////////////////////////////////////////////////////////////

template<typename S, typename T> void
slug_tracks_mist::read_trackfile(const char *fname, S& logt,
				 T& trackdata) {

  // Store number of masses and tracks for convenience
  size_type ntime = trackdata.shape()[0]-1;
  size_type ntrack = trackdata.shape()[1];

  // Open file
  fitsfile *fptr;
  int fits_status = 0;
  fits_open_table(&fptr, fname, READONLY, &fits_status);
  if (fits_status) {
    ostreams.slug_err_one << "unable to open file "
			  << fname << endl;
    bailout(1);
  }
  
  // Loop over masses
  double *buf = new double[ntime];
  for (size_type i=0; i<ntrack; i++) {

    // Move to next HDU
    int hdutype;
    fits_movrel_hdu(fptr, 1, &hdutype, &fits_status);
    if (fits_status) {
      ostreams.slug_err_one << "unable to read data from file "
			    << fname << endl;
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
			    << fname << endl;
      bailout(1);
    }
    for (size_type j=0; j<ntime; j++)
      logt[j+1][i] = log(buf[j]);

    // Read data; again, note that mass and mdot are stored as linear
    // values in the MIST tracks, but we use log values internally in
    // slug
    for (size_type n=0; n<nprop; n++) {
      char tmp_colname[10];
      strcpy(tmp_colname, fits_colnames[n]); // Needed to avoid const issues
      fits_get_colnum(fptr, CASEINSEN, tmp_colname,
		      &colnum, &fits_status);
      fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ntime, nullptr, buf,
		    &anynul, &fits_status);
      if (fits_status) {
	ostreams.slug_err_one << "unable to read data from file "
			      << fname << endl;
	bailout(1);
      }
      if (n == idx_log_cur_mass || n == idx_log_mDot) {
	for (size_type j=0; j<ntime; j++)
	  trackdata[j+1][i][n] = log(buf[j]);
      } else {
	for (size_type j=0; j<ntime; j++)
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
  for (size_type i=0; i<ntrack; i++) {
    logt[0][i] = -log10(constants::big);
    for (size_type j=0; j<trackdata.shape()[2]; j++) {
      trackdata[0][i][j] = trackdata[1][i][j];
    }
  }
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
