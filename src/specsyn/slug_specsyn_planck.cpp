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
#include <cmath>
#include <cassert>
#include "slug_specsyn_planck.H"
#include "../constants.H"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Constructor with default wavelengths. Default is 1001 wavelengths,
// logarithmically spaced from 9.1 x 10^1 to 9.1 x 10^5 Angstrom.
////////////////////////////////////////////////////////////////////////

slug_specsyn_planck::
slug_specsyn_planck(const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf,
		    const slug_PDF *my_sfh,
		    slug_ostreams& ostreams_,
		    const double z_in) : 
  slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in)
{
  for (unsigned int i=0; i<1001; i++) {
    lambda_rest.push_back(91.0 * 
			  pow(10., i/1000.0*4));
    lambda_obs.push_back(lambda_rest[i]*(1.0+z));
  }
  v_integ.set_nvec(lambda_rest.size()+1);
}


////////////////////////////////////////////////////////////////////////
// Constructor from given min, max, number of table entries
////////////////////////////////////////////////////////////////////////
slug_specsyn_planck::
slug_specsyn_planck(const double lambda_min, const double lambda_max, 
		    const unsigned int nlambda, 
		    const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf, const slug_PDF *my_sfh,
		    slug_ostreams& ostreams_,
		    const double z_in, const bool rest) : 
  slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in) {
  if (rest) {
    for (unsigned int i=0; i<nlambda; i++) {
      lambda_rest.push_back(lambda_min * 
			    pow(lambda_max/lambda_min, 
				((double) i)/(nlambda - 1)));
      lambda_obs.push_back(lambda_rest[i]*(1.0+z));
    }
  } else {
    for (unsigned int i=0; i<nlambda; i++) {
      lambda_obs.push_back(lambda_min * 
			    pow(lambda_max/lambda_min, 
				((double) i)/(nlambda - 1)));
      lambda_rest.push_back(lambda_obs[i]/(1.0+z));
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Constructor from specified wavelength vector
////////////////////////////////////////////////////////////////////////
slug_specsyn_planck::
slug_specsyn_planck(const vector<double>& lambda_in,
		    const slug_tracks *my_tracks, 
		    const slug_PDF *my_imf, 
		    const slug_PDF *my_sfh,
		    slug_ostreams& ostreams_,
		    const double z_in, 
		    const bool rest) :
  slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in) {
  if (rest) {
    lambda_rest = lambda_in;
    lambda_obs.resize(lambda_rest.size());
    for (unsigned int i=0; i<lambda_rest.size(); i++)
      lambda_obs[i] = lambda_rest[i] * (1.0+z);
  } else {
    lambda_obs = lambda_in;
    lambda_rest.resize(lambda_obs.size());
    for (unsigned int i=0; i<lambda_rest.size(); i++)
      lambda_rest[i] = lambda_obs[i] / (1.0+z);
  }
}

////////////////////////////////////////////////////////////////////////
// Method to change the wavelength table
////////////////////////////////////////////////////////////////////////
void
slug_specsyn_planck::
set_lambda(const vector<double>& lambda_in, bool rest) {
  if (rest) {
    lambda_rest = lambda_in;
    lambda_obs.resize(lambda_rest.size());
    for (unsigned int i=0; i<lambda_rest.size(); i++)
      lambda_obs[i] = lambda_rest[i] * (1.0+z);
  } else {
    lambda_obs = lambda_in;
    lambda_rest.resize(lambda_obs.size());
    for (unsigned int i=0; i<lambda_rest.size(); i++)
      lambda_rest[i] = lambda_obs[i] / (1.0+z);
  }
}

////////////////////////////////////////////////////////////////////////
// The spectral synthesizer. This just evaluates the Planck function
// to get the spectrum for each star, then sums.
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_planck::
get_spectrum_const(const vector<slug_stardata>& stardata) const {

  // Build output vector
  vector<double> L_lambda(lambda_rest.size(), 0.0);

  // Loop over stars
  for (unsigned int i=0; i<stardata.size(); i++) {

    // Add 4 pi R^2 B_lambda(T), in units of erg/s/Angstrom
    double surf_area = 4.0 * M_PI * 
      pow(10.0, 2.0*(stardata[i].logR+constants::logRsun));
    double Teff = pow(10.0, stardata[i].logTeff);
    for (unsigned int j=0; j<lambda_rest.size(); j++) {
      double lcgs = lambda_rest[j] * constants::Angstrom;
      double x = constants::hcOverkB / (lcgs * Teff);
      double b_lambda = 2.0 * constants::hc2 /
	(pow(lcgs, 5.0) * (exp(x) - 1));
      L_lambda[j] += surf_area * b_lambda * constants::Angstrom;
    }
  }

  // Return
  return L_lambda;
}



////////////////////////////////////////////////////////////////////////
// Spectral synthesizer for a single star
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_planck::
get_spectrum(const slug_stardata& stardata) const {

  // Build output bector
  vector<double> L_lambda(lambda_rest.size());

  // Surface area, Teff
  double surf_area = 4.0 * M_PI * 
    pow(10.0, 2.0*(stardata.logR+constants::logRsun));
  double Teff = pow(10.0, stardata.logTeff);

  // Add 4 pi R^2 (pi B_lambda(T)), in units of erg/s/Angstrom
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    double lcgs = lambda_rest[j] * constants::Angstrom;
    double x = constants::hcOverkB / (lcgs * Teff);
    double b_lambda = 2.0 * constants::hc2 /
      (pow(lcgs, 5.0) * (exp(x) - 1));
    L_lambda[j] = surf_area * M_PI * b_lambda * constants::Angstrom;
  }

  // Return
  return L_lambda;
}
