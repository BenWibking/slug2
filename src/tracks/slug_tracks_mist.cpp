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

#ifdef WINDS_ON
////////////////////////////////////////////////////////////////////////
// Methods to set the stellar wind velocity. The models here use the
// same models used to derive the mass loss rates in MIST for
// consistency as far as possible, but for some stellar types MIST's
// mass loss rates are purely empirical fits that make no predictions
// about velocity. See Roy+ (2019) for a discussion of the various
// assumptions we make to handle these cases.
////////////////////////////////////////////////////////////////////////
void slug_tracks_mist::
set_wind(const double m, const double t,
	 slug_stardata& star) const {

  // Set wind mass loss rate
  double logm = log(m);
  double logt = log(t);
  star.mDot = pow(10., (*interp)(logt, logm, idx_log_mDot));

  // Get phase and effecitve temperature
  double phase = (*interp)(logt, logm, idx_phase);
  double T_eff = pow(10., star.logTeff);

  // Set wind speed based on type of star
  if (m < 10.0 && phase <= 4.0) {

    // Low mass main sequence stars: we treat these as having a wind
    // speed equal to the surface escape speed; this is very, very
    // crude, but, given the very small mass loss rates and power, this
    // basically doesn't matter
    star.vWind = 1.0e-5 * sqrt(2.) *
      pow(10., 0.5 * (constants::logG + constants::logMsun -
		      constants::logRsun + star.logM - star.logR));
  
  } else if ((m < 10.0 && phase > 4.0) ||
	     (m >= 10.0 && T_eff <= 1.0e4 && star.WR == NONE)) {

    // AGB and RGB stars: MIST uses Reimers (1975) + Blocker (1995),
    // but these are empirical calibrations that make no
    // predictions for the wind velocity. We use the Elitzur & Ivezic
    // (2001, MNRAS, 327, 403) scaling v ~ R_gd^-1/2 L^1/4, where R_gd
    // = gas to dust ratio. We set the proportionality constant using
    // the empirical calibration given in Goldman et al. (2017, MNRAS,
    // 465, 403), which gives v = 9.4 km/s at L = 10^4 Lsun, R_gd =
    // Milky Way / Solar value. Note that the scaling is not what 
    // Goldman et al. find empirically, but I'm using the Elitzur 
    // scaling with the Goldman normalisation by recommendation of
    // Jacco van Loon. Jacco also notes that the metallicity scaling
    // almost certainly fails for metal poor stars because their winds
    // are not dust-driven, but was not able to give me a better 
    // estimate to use in that case.
    star.vWind = 9.4 * (pow(10., 0.5*star.logL) - 4.0) /
      sqrt(metallicity);

  } else if (m >= 10.0 && T_eff > 1.0e4 && star.WR == NONE) {

    // Massive stars

    // Decide what to do based on WR status and T_eff
    if (T_eff >= 1.1e4 && star.WR == NONE) {

      // O star regime: use the Vink et al. (2000, 2001) scheme

      // Escape speed in km/s
      double v_esc = 1.0e-5 * sqrt(2.) *
	pow(10., 0.5 * (constants::logG + constants::logMsun -
			constants::logRsun + star.logM - star.logR));
      if (T_eff && T_eff >= 2.75e4) {
        // Hot regime
        star.vWind = 2.6 * pow(metallicity, 0.13) * v_esc;
      } else if (T_eff > 1.1e4 && T_eff < 2.25e4) {
        // Cold regime
        star.vWind = 1.3 * pow(metallicity, 0.13) * v_esc;
      } else if (T_eff < 1.1e4) {
	// Very cold regime; here we interpolate linearly between the
	// O star and AGB star cases
	double vWind_O = 1.3 * pow(metallicity, 0.13) * v_esc;
	double vWind_AGB = 9.4 * (pow(10., 0.5*star.logL) - 4.0) /
	  sqrt(metallicity);
	star.vWind = vWind_O * (T_eff - 1.0e4) / 1.0e3 +
	  vWind_AGB * (1.1e4 - T_eff) / 1.0e3;
      } else {
        // Bistable regime; use equation 15 of Vink+ 2001
        double logrho = -13.636 + 0.889 * log(metallicity);
        double T_eff_jump = 6.1e4 + 2.59e3 * logrho;
        if (T_eff < T_eff_jump)
          star.vWind = 1.3 * pow(metallicity, 0.13) * v_esc;
        else
          star.vWind = 2.6 * pow(metallicity, 0.13) * v_esc;
      }

    }

  } else if (star.WR != NONE) {

    // Wolf-Rayet star; we set the velocities in these stars by
    // setting the wind momentum to L/c. While more detailed fits are
    // given in Nugis & Lamers (2000), which is where the mass loss
    // rates come from, extrapolating their fits outside the limited
    // range in metallicities over which Nugis & Lamers had data leads
    // to nonsense results. Rather than trying to fix up their
    // formulae, we'll just get something right to within a factor of
    // ~2, but which will be non-crazy over all parameter regimes.
    double L = pow(10., star.logL + constants::logLsun);
    double mDot = star.mDot * constants::Msun / constants::yr;
    star.vWind = L / (mDot * constants::c) / 1.0e5; // convert to km/s

  }

}

void slug_tracks_mist::
set_wind(const double m, 
	 spl_arr_view_1d& isochrone_,
	 acc_arr_view_1d& isochrone_acc_,
	 slug_stardata& star) const {

  // Set mass loss rate
  double logm = log(m);
  star.mDot =
    pow(10.0, gsl_spline_eval(isochrone_[idx_log_mDot],
			      logm,
			      isochrone_acc_[idx_log_mDot]));

  // Get phase and effecitve temperature
  double phase = gsl_spline_eval(isochrone_[idx_phase],
				 logm,
				 isochrone_acc_[idx_phase]);
  double T_eff = pow(10., star.logTeff);
  
  // Set wind speed based on type of star
  if (m < 10.0 && phase <= 4.0) {

    // Low mass main sequence stars: we treat these as having a wind
    // speed equal to the surface escape speed; this is very, very
    // crude, but, given the very small mass loss rates and power, this
    // basically doesn't matter
    star.vWind = 1.0e-5 * sqrt(2.) *
      pow(10., 0.5 * (constants::logG + constants::logMsun -
		      constants::logRsun + star.logM - star.logR));
  
  } else if ((m < 10.0 && phase > 4.0) ||
	     (m >= 10.0 && T_eff <= 1.0e4 && star.WR == NONE)) {

    // AGB and RGB stars: MIST uses Reimers (1975) + Blocker (1995),
    // but these are empirical calibrations that make no
    // predictions for the wind velocity. We use the Elitzur & Ivezic
    // (2001, MNRAS, 327, 403) scaling v ~ R_gd^-1/2 L^1/4, where R_gd
    // = gas to dust ratio. We set the proportionality constant using
    // the empirical calibration given in Goldman et al. (2017, MNRAS,
    // 465, 403), which gives v = 9.4 km/s at L = 10^4 Lsun, R_gd =
    // Milky Way / Solar value. Note that the scaling is not what 
    // Goldman et al. find empirically, but I'm using the Elitzur 
    // scaling with the Goldman normalisation by recommendation of
    // Jacco van Loon. Jacco also notes that the metallicity scaling
    // almost certainly fails for metal poor stars because their winds
    // are not dust-driven, but was not able to give me a better 
    // estimate to use in that case.
    star.vWind = 9.4 * (pow(10., 0.5*star.logL) - 4.0) /
      sqrt(metallicity);

  } else if (m >= 10.0 && T_eff > 1.0e4 && star.WR == NONE) {

    // Massive stars

    // Decide what to do based on WR status and T_eff
    if (T_eff >= 1.1e4 && star.WR == NONE) {

      // O star regime: use the Vink et al. (2000, 2001) scheme

      // Escape speed in km/s
      double v_esc = 1.0e-5 * sqrt(2.) *
	pow(10., 0.5 * (constants::logG + constants::logMsun -
			constants::logRsun + star.logM - star.logR));
      if (T_eff && T_eff >= 2.75e4) {
        // Hot regime
        star.vWind = 2.6 * pow(metallicity, 0.13) * v_esc;
      } else if (T_eff > 1.1e4 && T_eff < 2.25e4) {
        // Cold regime
        star.vWind = 1.3 * pow(metallicity, 0.13) * v_esc;
      } else if (T_eff < 1.1e4) {
	// Very cold regime; here we interpolate linearly between the
	// O star and AGB star cases
	double vWind_O = 1.3 * pow(metallicity, 0.13) * v_esc;
	double vWind_AGB = 9.4 * (pow(10., 0.5*star.logL) - 4.0) /
	  sqrt(metallicity);
	star.vWind = vWind_O * (T_eff - 1.0e4) / 1.0e3 +
	  vWind_AGB * (1.1e4 - T_eff) / 1.0e3;
      } else {
        // Bistable regime; use equation 15 of Vink+ 2001
        double logrho = -13.636 + 0.889 * log(metallicity);
        double T_eff_jump = 6.1e4 + 2.59e3 * logrho;
        if (T_eff < T_eff_jump)
          star.vWind = 1.3 * pow(metallicity, 0.13) * v_esc;
        else
          star.vWind = 2.6 * pow(metallicity, 0.13) * v_esc;
      }

    }

  } else if (star.WR != NONE) {

    // Wolf-Rayet star; we set the velocities in these stars by
    // setting the wind momentum to L/c. While more detailed fits are
    // given in Nugis & Lamers (2000), which is where the mass loss
    // rates come from, extrapolating their fits outside the limited
    // range in metallicities over which Nugis & Lamers had data leads
    // to nonsense results. Rather than trying to fix up their
    // formulae, we'll just get something right to within a factor of
    // ~2, but which will be non-crazy over all parameter regimes.
    double L = pow(10., star.logL + constants::logLsun);
    double mDot = star.mDot * constants::Msun / constants::yr;
    star.vWind = L / (mDot * constants::c) / 1.0e5; // convert to km/s

  }
  
}

#endif


#endif
// ENABLE_FITS
