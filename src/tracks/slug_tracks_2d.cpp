/*********************************************************************
Copyright (C) 2014 Robert da Silva, Michele Fumagalli, Mark Krumholz
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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_errno.h>
#include "slug_tracks_2d.H"
#include "../slug_MPI.H"
#include "../constants.H"

using namespace std;
using namespace boost;
using namespace boost::multi_array_types;

////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_tracks_2d::~slug_tracks_2d() {

  // De-allocate the cached isochrone
  free_isochrone();

  // Destroy the interpolation object
  delete interp;
}


////////////////////////////////////////////////////////////////////////
// Mass of star dying at a particular time; this form of the function
// only works with monotonic tracks
////////////////////////////////////////////////////////////////////////
double slug_tracks_2d::death_mass(const double time,
				  const double Z) const {

  // Safety assertion: input metallicity should always be the null
  // value for the 2d case
  assert(Z == tracks::null_metallicity);

  // Safety assertion: this routine should only be called on tracks
  // that are monotonic
  assert(monotonic);

  // Get log of time
  double logt;
  if (time > exp(interp->x_min())) logt = log(time);
  else logt = interp->x_min();
  
  // Is this time less than the smallest death time we have? If so,
  // return a big number.
  if (logt < interp->x_max(interp->y_max()))
    return(constants::big);

  // Is this time bigger than the biggest death time we have? If so,
  // return 0.
  if (logt > interp->x_max(interp->y_min()))
    return(0.0);

  // If the time is between the largest and smallest death times, get
  // it from the tracks
  vector<double> logm_lim = interp->y_lim(logt);
  return exp(logm_lim[0]);
}


////////////////////////////////////////////////////////////////////////
// Masses of stars dying at a particular time; this form of the
// function works with both monotonic and non-monotonic tracks
////////////////////////////////////////////////////////////////////////
void slug_tracks_2d::death_mass(const double time,
				vector<double>& m_death,
				const double Z) const {

  // Safety assertion: input metallicity should always be the null
  // value for the 2d case
  assert(Z == tracks::null_metallicity);

  // Get log of time
  double logt;
  if (time > exp(interp->x_min())) logt = log(time);
  else logt = interp->x_min();
  
  // Is this time less than the smallest death time we have? If so,
  // set the death mass vector to be empty, then return.
  if (logt < interp->x_max(interp->y_max())) {
    m_death.resize(0);
    return;
  }
 
  // Is this time bigger than the biggest death time we have? If so,
  // set the death mass vector to be empty, then return.
  if (logt > interp->x_max(interp->y_min())) {
    m_death.resize(0);
    return;
  }

  // If the time is between the largest and smallest death times, get
  // it from the tracks
  vector<double> logm_lim = interp->y_lim(logt);

  // The limits returned by the tracks may include the minimum and
  // maximum mass in the tracks; if so, eliminate these
  if (logm_lim.size() > 0) {
    if (logm_lim[0] == interp->y_min()) logm_lim.erase(logm_lim.begin());
  }
  if (logm_lim.size() > 0) {
    if (logm_lim.back() == interp->y_max()) logm_lim.pop_back();
  }

  // Convert from log mass to mass
  m_death.resize(logm_lim.size());
  for (vector<double>::size_type i=0; i<logm_lim.size(); i++)
    m_death[i] = exp(logm_lim[i]);
}


////////////////////////////////////////////////////////////////////////
// Calculate what range of stellar masses are alive at the specified
// time
////////////////////////////////////////////////////////////////////////
vector<double>
slug_tracks_2d::live_mass_range(const double time,
				const double Z) const {

  // Safety assertion: input metallicity should always be the null
  // value for the 2d case
  assert(Z == tracks::null_metallicity);

  // Get log of time
  double logt = log(time);
  logt = max(logt, interp->x_min());

  // If time is larger than largest time, return empty vector;
  // otherwise get answer from interpolation object; note that the
  // interpolation object works in log mass, so we need to exponentiate
  // the answer we get back from it; we also adjust the first entry
  // downward to zero, under the assumption that stars with masses
  // below those included in the tracks always have an infinite
  // lifetime
  vector<double> m;
  if (logt <= interp->x_max()) {
    m = interp->y_lim(logt);
    m[0] = 0.0;
    for (vector<double>::size_type i=1; i<m.size(); i++)
      m[i] = exp(m[i]);
  }
  return m;
}

////////////////////////////////////////////////////////////////////////
// Mass of remnant left by a star of a given mass; for now this is
// hardcoded to the compilation of Kruijssen (2009). Versions
// including an age argument return 0 if the star has not yet
// died. Note that for very massive stars the Kruijssen (2009) formula
// returns a value larger than the final mass of the star in the
// tracks; this is just an inconsistency between the tracks and
// Kruijssen's adopted initial-final mass relation. We test for this
// case and reduce the remnant mass if necessary, to avoid having a
// situation where the predicted remnant mass exceeds the mass of the
// star at death.
////////////////////////////////////////////////////////////////////////
double
slug_tracks_2d::remnant_mass(const double mass,
			     const double age,
			     const double Z) const {
  // Safety assertion: input metallicity should always be the null value
  assert(Z == tracks::null_metallicity);

  // If age > 0, check lifetime, and return 0 if star is still alive
  if (age > 0)
    if (age < star_lifetime(mass)) return 0.0;

  // Return remnant mass
  double rem_mass;
  if (mass < 8.0) rem_mass = 0.109*mass + 0.394;
  else if (mass < 30) rem_mass = 0.03636*(mass-8.0) + 1.02;
  else rem_mass = 0.06*(mass-30.0) + 8.3;
  double final_mass = star_mass_at_death(mass);
  if (rem_mass < final_mass) return rem_mass;
  else return final_mass;
}


vector<double>
slug_tracks_2d::remnant_mass(const vector<double>& mass,
			     const double age,
			     const double Z) const {
  // Safety assertion: input metallicity should always be the null value
  assert(Z == tracks::null_metallicity);
  vector<double> rmass(mass.size(), 0.0);
  for (std::vector<double>::size_type i=0; i<mass.size(); i++) {
    rmass[i] = remnant_mass(mass[i], age);
  }
  return rmass;
}


////////////////////////////////////////////////////////////////////////
// Method to free an isochrone
////////////////////////////////////////////////////////////////////////
void slug_tracks_2d::free_isochrone() const {

  // Do nothing if isochrone is not set
  if (isochrone_age < 0.0) return;

  // Free
  for (size_type i=0; i<isochrone.shape()[0]; i++) {
    for (size_type j=0; j<isochrone.shape()[1]; j++) {
      gsl_spline_free(isochrone[i][j]);
      gsl_interp_accel_free(isochrone_acc[i][j]);
    }
  }

  // Flag that isochrone is not set by setting age to a negative value
  isochrone_age = -1.0;
}

////////////////////////////////////////////////////////////////////////
// Method to set the properties of a single star
////////////////////////////////////////////////////////////////////////
slug_stardata
slug_tracks_2d::get_star(const double m, const double t,
			 const double Z) const {
  
  // Get log m and log t
  double logm = log(m);
  double logt = log(t);

  // Use interpolator to interpolate to requested values
  slug_stardata star;
  star.logM = constants::loge * (*interp)(logt, logm,
					  idx_log_cur_mass);
  star.logL = (*interp)(logt, logm, idx_log_L);
  star.logTeff = (*interp)(logt, logm, idx_log_Teff);
  star.logR = 0.5*(star.logL+constants::logLsun) 
    - 0.5*log10(4.0*M_PI) 
    - 0.5*constants::logsigmaSB - 2.0*star.logTeff
    - constants::logRsun;
  star.logg = constants::logG + star.logM +
    constants::logMsun - 2.0*star.logR - 
    2.0*constants::logRsun;
  set_WR_type(m, t, star);
  return star;
}

////////////////////////////////////////////////////////////////////////
// Method to get data for a set of stars on an isochrone.
////////////////////////////////////////////////////////////////////////

// Notes: 
//
// (1) array of input masses m must be sorted in increasing order
//
// (2) m may contain stars whose mass is below the minimum track mass
// or above the maximum track mass. These will simply be omitted from
// the output vectors, so that these vectors may have fewer elements
// than the input mass vector.

vector<slug_stardata>
slug_tracks_2d::
get_isochrone(const double t, const vector<double> &m,
	      const double Z) const {

  // Get log time to feed to interpolator; if time is past end of
  // tracks, return empty vector
  vector<slug_stardata> stars;
  double logt;
  if (t < exp(interp->x_min())) logt = interp->x_min();
  else if (t > exp(interp->x_max())) return stars;
  else logt = log(t);

  // Build an isochrone for this time if necessary
  if (t != isochrone_age) {
    free_isochrone();
    isochrone_age = t;
    interp->build_interp_const_x(logt, isochrone, isochrone_acc,
				 isochrone_logm_lim);
  }

  // Initialize interval pointer and alive bit
  bool alive = false;
  size_type ptr = -1;

  // Turn off GSL error handling so that we can handle this on our own
  // and print out more informative error messages
  gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

  // Loop over input stellar masses
  for (vector<double>::size_type i=0; i<m.size(); i++) {

    // Log of initial mass
    double logm = log(m[i]);

    // Check if the alive bit should flip; check for running past
    // alive mass interval
    while (logm >= isochrone_logm_lim[ptr+1]) {
      alive = !alive;
      ptr++;
      if (ptr == isochrone_logm_lim.size()-1) return stars;
    }

    // If this star is not alive, do nothing
    if (!alive) continue;

    // This star is alive; interpolate to get all its physical
    // properties
    slug_stardata star;
    int gsl_errstat =
      gsl_spline_eval_e(isochrone[ptr/2][idx_log_L], logm,
			isochrone_acc[ptr/2][idx_log_L], &star.logL);
    if (gsl_errstat == GSL_EDOM) {
      ostreams.slug_err
	<< "slug_tracks_2d::get_isochrone: interpolation error at time " 
	<< setprecision(20) << exp(logt) << ", mass = "
	<< setprecision(20) << exp(logm)
	<< "; at this time, living mass intervals are: ";
      for (size_type j=0; j<isochrone_logm_lim.size()/2; j++) {
	if (j != 0) ostreams.slug_err << ", ";
	ostreams.slug_err
	  << setprecision(20) << exp(isochrone_logm_lim[2*j]) << " - " 
	  << setprecision(20) << exp(isochrone_logm_lim[2*j+1]);
      }
      ostreams.slug_err << endl;
      bailout(1);
    }

    // Teff
    star.logTeff
      = gsl_spline_eval(isochrone[ptr/2][idx_log_Teff], logm,
			isochrone_acc[ptr/2][idx_log_Teff]);

    // R from L_bol and T
    star.logR = 0.5*(star.logL+constants::logLsun) 
      - 0.5*log10(4.0*M_PI) 
      - 0.5*constants::logsigmaSB - 2.0*star.logTeff
      - constants::logRsun;

    // Current mass; safety check to ensure that mass is strictly
    // non-increasing, avoiding mass gain due to roundoff error
    star.logM = constants::loge *
      gsl_spline_eval(isochrone[ptr/2][idx_log_cur_mass], logm, 
		      isochrone_acc[ptr/2][idx_log_cur_mass]);
    if (pow(10.0, star.logM) > m[i]) {
      star.logM = log10(m[i]);
    }
    //if (m[i] > 8.0)
    //  cout << "i = " << i
    //	   << ", m_init = " << m[i]
    //	   << ", m_cur = " << exp(star.logM) << endl;

    // log g
    star.logg = constants::logG + star.logM +
      constants::logMsun - 2.0*star.logR - 
      2.0*constants::logRsun;

    // WR status
    spl_arr_view_1d isochrone_view
      = isochrone[indices[ptr/2][range_t(0,interp->shape()[2])]];
    acc_arr_view_1d isochrone_acc_view
      = isochrone_acc[indices[ptr/2][range_t(0,interp->shape()[2])]];
    set_WR_type(m[i], isochrone_view, isochrone_acc_view, star);

    // Push onto output list
    stars.push_back(star);
  }

  // Restore the GSL error handler
  gsl_set_error_handler(gsl_err);
  
  // Return
  return stars;
}
