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
#include "slug_specsyn_hillier.H"
#include "../constants.H"
#include "../slug_MPI.H"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;


////////////////////////////////////////////////////////////////////////
// Convenient comparator to sort stars. The sorting rule is that WN's
// come before WC's, and within each class sort by Teff. Note that we
// stick this in its own namesapce rather than attaching it to the
// class because std::sort requires that sorting functions be static.
////////////////////////////////////////////////////////////////////////
namespace hillier {
  bool starsort(const slug_stardata& star1, 
		const slug_stardata& star2) {
    if (star1.WR == star2.WR) {
      return (star1.logTeff < star2.logTeff);
    } else {
      return (star1.WR == WNE || star1.WR == WNL);
    }
  }
}

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_hillier::
slug_specsyn_hillier(const char *dirname, const slug_tracks *my_tracks, 
		     const slug_PDF *my_imf, const slug_PDF *my_sfh,
		     slug_ostreams& ostreams_, const double z_in,
		     const bool check_data_in) :
  slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in),
  Teff_wn(12), Teff_wc(12), logT_break_wn(11), logT_break_wc(11),
  check_data(check_data_in)
{

  // Set arrays to correct size
  array2d::extent_gen extent;
  F_lam_wn.resize(extent[12][1221]);
  F_lam_wc.resize(extent[12][1221]);
  lambda_rest.resize(1221);

  // Choose the atmosphere file closest in metallicity to the
  // metallicity we are using for the tracks
  const string extensions[] = {"001", "004", "008", "020", "040"};
  const double zrel[] = { 0.05, 0.2, 0.4, 1.0, 2.0 }; // Z relative to Solar
  const int nfile = 5;
  double zdiff = constants::big;
  int idx = -1;
  for (int i=0; i<nfile; i++) {
    double zdiff1 = abs(log10(zrel[i]) - log10(my_tracks->get_metallicity()));
    if (zdiff > zdiff1) {
      idx = i;
      zdiff = zdiff1;
    }
  }

  // Print warning if this is a big extrapolation in metallicity
  if (zdiff > 0.1) {
    ostreams.slug_warn_one
      << "stellar track metallicity is [Z/H] = "
      << log10(my_tracks->get_metallicity())
      << ", closest Hillier atmosphere metallicity available is [Z/H] = "
      << log10(zrel[idx])
      << "; calculation will proceed" << endl;
  }

  // Construct the file name for the WC models and try to open the file
  string fname = "CMFGEN_WC_Z"+extensions[idx]+".dat";
  std::ifstream atmos_file;
  path dirname_path(dirname);
  path atmos_path = dirname_path / path(fname.c_str());
  atmos_file.open(atmos_path.c_str());
  if (!atmos_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open atmosphere file " 
      << atmos_path.string() << endl;
    bailout(1);
  }

  // Save file name
  wc_file_name = atmos_path.string();

  // Read the WC models
  string modelhdr;
  double modelnum;
  for (unsigned int i=0; i<11; i++) 
    getline(atmos_file, modelhdr); // Burn 11 lines
  for (unsigned int i=0; i<12; i++) {
    getline(atmos_file, modelhdr);   // Burn a line
    getline(atmos_file, modelhdr);   // Burn a line 
    atmos_file >> modelnum;
    atmos_file >> Teff_wc[i];
    for (unsigned int j=0; j<1221; j++) {
      atmos_file >> lambda_rest[j] >> F_lam_wc[i][j];
    }
  }

  // Close the file
  atmos_file.close();

  // Repeat for WN models
  fname = "CMFGEN_WN_Z"+extensions[idx]+".dat";
  atmos_path = dirname_path / path(fname.c_str());
  atmos_file.open(atmos_path.c_str());
  if (!atmos_file.is_open()) {
    // Couldn't open file, so bail out
    ostreams.slug_err_one
      << "unable to open atmosphere file " 
      << atmos_path.string() << endl;
    bailout(1);
  }
  wn_file_name = atmos_path.string();
  for (unsigned int i=0; i<11; i++) getline(atmos_file, modelhdr);
  for (unsigned int i=0; i<12; i++) {
    getline(atmos_file, modelhdr);
    getline(atmos_file, modelhdr);
    atmos_file >> modelnum;
    atmos_file >> Teff_wn[i];
    for (unsigned int j=0; j<1221; j++) {
      atmos_file >> lambda_rest[j] >> F_lam_wn[i][j];
    }
  }
  atmos_file.close();

  // Compute observed frame wavelengths
  lambda_obs.resize(lambda_rest.size());
  for (unsigned int i=0; i<1221; i++) 
    lambda_obs[i] = (1.0+z)*lambda_rest[i];

  // Renormalize the Hillier models so that they give fluxes at
  // stellar surfaces, rather than at 1 kpc. Radii of stars are
  // hardcoded into SB99, and are hardcoded here too -- uggh this is
  // ugly -- what we do for backwards compatibility...
  const double wcradii[] = {9.31, 8.04, 7.04, 5.95, 4.94, 4.14, 3.05, 
			    2.33, 1.84, 1.50, 1.03, 0.7};
  const double wnradii[] = {20.30, 17.24, 14.90, 11.40, 9.00, 7.30, 
			    5.07, 3.72, 2.84, 2.26, 1.82, 1.27};
  for (unsigned int i=0; i<12; i++) {
    for (unsigned int j=0; j<1221; j++) {
      F_lam_wc[i][j] *= 
	pow(constants::kpc/(constants::Rsun*wcradii[i]), 2);
      F_lam_wn[i][j] *= 
	pow(constants::kpc/(constants::Rsun*wnradii[i]), 2);
    }
  }

  // Generate breakpoints between models in log T, since that's what
  // will be passed in
  for (unsigned int i=0; i<11; i++) {
    logT_break_wn[i] = log10(0.5*(Teff_wn[i]+Teff_wn[i+1]));
    logT_break_wc[i] = log10(0.5*(Teff_wc[i]+Teff_wc[i+1]));
  }

  // If using data range checking, set up the Kurucz and Planck
  // synthesizers to fall back on
  if (check_data) {
    kurucz = new slug_specsyn_kurucz(dirname, tracks, imf, sfh,
				     ostreams, z, false);
    planck = new slug_specsyn_planck(lambda_obs, tracks, imf, sfh,
				     ostreams, z);
  } else {
    kurucz = NULL;
    planck = NULL;
  }

  // Set the size of the integrator for continuous IMFs to #
  // wavelengths + 1 (for the bolometric luminosity)
  v_integ.set_nvec(lambda_rest.size()+1);
}

////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_specsyn_hillier::~slug_specsyn_hillier() {
  if (planck != NULL) delete planck;
  if (kurucz != NULL) delete kurucz;
}


////////////////////////////////////////////////////////////////////////
// Function to turn data checking on or off
////////////////////////////////////////////////////////////////////////
void
slug_specsyn_hillier::set_check_data(bool val) {
  if ((planck == NULL) && (val == true)) 
    planck = new slug_specsyn_planck(lambda_obs, tracks, imf, sfh,
				     ostreams, z);
  if ((kurucz == NULL) && (val == true)) {
    path filename(wc_file_name);
    path dirname = filename.parent_path();
    kurucz = new slug_specsyn_kurucz(dirname.string().c_str(), 
				     tracks, imf, sfh, ostreams,
				     z, false);
  }
  check_data = val;
}

////////////////////////////////////////////////////////////////////////
// Wrapper function to get stellar spectra that first sorts stars by
// whether they should use the Hillier WN/WC model atmospheres, Kurucz
// atmospheres, or Planck (blackbody) atmospheres.
////////////////////////////////////////////////////////////////////////

vector<double>
slug_specsyn_hillier::
get_spectrum(vector<slug_stardata>& stars) const {

  // If not doing any safety checking, just call the function that
  // operates on the full list of stars, then return. This will
  // probably not crash the code if there are non-WR stars in the
  // list, but the results will be junk.
  if (!check_data) return get_spectrum_clean(stars);

  // Initialize
  vector<double> L_lambda(lambda_rest.size());

  // Separate list into WR, non-WR stars within the temperature range
  // allowed by the Kurucz atmospheres, and non-WR stars outside the
  // Kurucz range (which will use Planck)
  vector<slug_stardata> stars_WR, stars_ku, stars_pl;
  for (unsigned int i=0; i<stars.size(); i++) {
    if (stars[i].WR != NONE) {
      stars_WR.push_back(stars[i]);
    } else if ((stars[i].logTeff >= kurucz->get_logTeff_min()) &&
	       (stars[i].logTeff <= kurucz->get_logTeff_max())) {
      stars_ku.push_back(stars[i]);
    } else {
      stars_pl.push_back(stars[i]);
    }
  }

  // Handle the WR stars here
  if (stars_WR.size() > 0)
    L_lambda = get_spectrum_clean(stars_WR);

  // Call Kurucz synthesizer for stars it can handle
  if (stars_ku.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      kurucz->get_spectrum(stars_ku);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Call Planck synthesizer for remaining stars
  if (stars_pl.size() > 0) {
    const vector<double> & L_lambda_tmp = 
      planck->get_spectrum(stars_pl);
    for (unsigned int i=0; i<L_lambda.size(); i++)
      L_lambda[i] += L_lambda_tmp[i];
  }

  // Return
  return L_lambda;
}


////////////////////////////////////////////////////////////////////////
// Function to compute spectrum for WR stars using Hillier model
// atmospheres. Closely modeled on starburst99's hillier subroutine in
// terms of functionality, but it should be significantly faster for
// large group of stars because the process of searching the model
// grid has been optimized significantly.
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_hillier::
get_spectrum_clean(vector<slug_stardata>& stars) const {

  // Initialize
  vector<double> L_lambda(lambda_rest.size(), 0.0);

  // Sort stars so that WN's come first, then WC's, and within each
  // class things are sorted by Teff
  sort(stars.begin(), stars.end(), hillier::starsort);

  // Compute surface areas of stars
  vector<double> surf_area(stars.size());
  for (unsigned int i=0; i<stars.size(); i++) {
    assert(stars[i].WR != NONE); // Safety check in debug mode
    surf_area[i] = 4.0 * M_PI * 
      pow(10.0, 2.0 * (stars[i].logR + constants::logRsun));
  }

  // Figure out where the break between WN and WC occurs in the star
  // list; note that, for this purpose, we treat WNL's and WNE's as
  // the same
  unsigned int wcptr = 0;
  while (stars[wcptr].WR == WNL || stars[wcptr].WR == WNE) {
    wcptr++;
    if (wcptr == stars.size()) break;
  }

  // Go through the WN stars and assign spectra based on model
  // temperatures
  unsigned int Tptr = 0;
  vector<double>::size_type ptr1 = 0, ptr2 = 0;
  while (ptr2 < wcptr) {

    // Move pointer 2 until it hits the next Teff value or the end of
    // the star list; exception: if Tptr is already has high as it
    // goes, then use the last model for the rest of the list
    if (Tptr == Teff_wn.size()-1) {
      ptr2 = wcptr;
    } else {
      while (1) {
	if (ptr2 == wcptr) break;
	if (stars[ptr2].logTeff > logT_break_wn[Tptr]) break;
	ptr2++;
      }
    }

    // Add contribution from stars in this temperature block, with a
    // normalization factor to account for the difference in Teff
    // between the star and the grid model.
    for (vector<double>::size_type i=ptr1; i<ptr2; i++) {
      double Tfac = pow(10.0, 4.0*stars[i].logTeff) 
	/ pow(Teff_wn[Tptr], 4);
      for (unsigned int j=0; j<L_lambda.size(); j++) {
	L_lambda[j] += surf_area[i] * F_lam_wn[Tptr][j] * Tfac;
      }
    }

    // Move ptr1 and go to next temperature
    ptr1 = ptr2;
    Tptr++;
  }

  // Repeat the procedure for the WC stars; for this purpose we also
  // treat WO's like WC's
  Tptr = 0;
  ptr2 = wcptr;
  while (ptr2 < stars.size()) {

    // Move pointer 2 until it hits the next Teff value or the end of
    // the star list; exception: if Tptr is already has high as it
    // goes, then use the last model for the rest of the list
    if (Tptr == Teff_wc.size()-1) {
      ptr2 = stars.size();
    } else {
      while (1) {
	if (ptr2 == stars.size()) break;
	if (stars[ptr2].logTeff > logT_break_wc[Tptr]) break;
	ptr2++;
      }
    }

    // Add contribution from stars in this temperature block, with a
    // normalization factor to account for the difference in Teff
    // between the star and the grid model
    for (vector<double>::size_type i=ptr1; i<ptr2; i++) {
      double Tfac = pow(10.0, 4.0*stars[i].logTeff) 
	/ pow(Teff_wc[Tptr], 4);
      for (unsigned int j=0; j<L_lambda.size(); j++) {
	L_lambda[j] += surf_area[i] * F_lam_wc[Tptr][j] * Tfac;
      }
    }

    // Go to next temperature
    Tptr++;
  }

  // Return
  return L_lambda;
}


////////////////////////////////////////////////////////////////////////
// Function to return spectrum from a single star
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_hillier::
get_spectrum(const slug_stardata& stardata) const {

  // Safety check if requested
  if (check_data && (stardata.WR == NONE)) {
    if ((stardata.logTeff >= kurucz->get_logTeff_min()) &&
	(stardata.logTeff <= kurucz->get_logTeff_max()))
      return kurucz->get_spectrum(stardata);
    else
      return planck->get_spectrum(stardata);
  }

  // Initialize
  vector<double> L_lambda(lambda_rest.size());

  // Get surface area
  double surf_area = 4.0 * M_PI * 
    pow(10.0, 2.0 * (stardata.logR + constants::logRsun));

  // Is this a WC or WN star? For this purpose, WN includes both WNE
  // and WNL, and WC includes both true WC's and WO's
  if (stardata.WR == WNE || stardata.WR == WNL) {

    // WN star. Now find the right effective temperature.
    unsigned int Tptr = 0;
    while (stardata.logTeff > logT_break_wn[Tptr]) {
      Tptr++;
      if (Tptr == Teff_wn.size()-1) break;
    }

    // Compute spectrum
    double Tfac = pow(10.0, 4.0*stardata.logTeff) 
      / pow(Teff_wn[Tptr], 4);
    for (unsigned int j=0; j<L_lambda.size(); j++)
      L_lambda[j] = surf_area * F_lam_wn[Tptr][j] * Tfac;

  } else {

    // WC star. Same procedure as for WN.
    unsigned int Tptr = 0;
    while (stardata.logTeff > logT_break_wc[Tptr]) {
      Tptr++;
      if (Tptr == Teff_wc.size()-1) break;
    }
    double Tfac = pow(10.0, 4.0*stardata.logTeff) 
      / pow(Teff_wc[Tptr], 4);
    for (unsigned int j=0; j<L_lambda.size(); j++)
      L_lambda[j] = surf_area * F_lam_wc[Tptr][j] * Tfac;

  }

  // Return
  return L_lambda;
}


