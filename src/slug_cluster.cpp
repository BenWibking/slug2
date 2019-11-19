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
#include "constants.H"
#include "slug_cluster.H"
#include "utils/int_tabulated.H"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <boost/bind.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Trivial little helper functions
////////////////////////////////////////////////////////////////////////
namespace cluster {
  double curMass(const slug_stardata &data) {
    return exp(data.logM/constants::loge);
  }
  double yield(const double m, const slug_yields *yields,
	       const vector<double>::size_type i) {
    return yields->yield(m, i);
  }
#ifdef WINDS_ON
  double star_mDot(const slug_stardata &data) { return data.mDot; }
  double star_pDot(const slug_stardata &data) {
    return data.mDot * data.vWind;
  }
  double star_wind_lum(const slug_stardata &data) {
    return 0.5 * (data.mDot * constants::Msun/constants::yr) *
      data.vWind*data.vWind*1e10 / constants::Lsun;
}
#endif
}
////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_cluster::slug_cluster(const unsigned long id_, 
			   const double mass_, 
			   const double time, const slug_PDF *imf_, 
			   const slug_tracks *tracks_, 
			   const slug_specsyn *specsyn_, 
			   const slug_filter_set *filters_,
			   const slug_extinction *extinct_,
			   const slug_nebular *nebular_,
			   const slug_yields *yields_,
			   const slug_line_list *lines_,
			   slug_ostreams& ostreams_,
			   const slug_PDF *clf_,
			   const bool stoch_contrib_only_) :
  ostreams(ostreams_),
  targetMass(mass_),
  imf(imf_),
  clf(clf_),
  tracks(tracks_), 
  specsyn(specsyn_),
  filters(filters_),
  extinct(extinct_),
  nebular(nebular_),
  yields(yields_),
  lines(lines_),
  integ(tracks_, imf_, nullptr, ostreams_),
  id(id_),
  formationTime(time),
  curTime(time),
  last_yield_time(time),
  stoch_contrib_only(stoch_contrib_only_)
{

  // Initialize to non-disrupted
  is_disrupted = false;

  // Populate with stars
  stochBirthMass = stochAliveMass = stochStellarMass =
    imf->drawPopulation(targetMass, stars);

  // If the population only represents part of the mass range due to
  // restrictions on what range is being treated stochastically, be
  // sure to account for that
  nonStochBirthMass = nonStochAliveMass = nonStochStellarMass =
     targetMass * (1.0 - imf->mass_frac_restrict());

  // Set the various mass tallies
  birthMass = stochBirthMass + nonStochBirthMass;
  aliveMass = stochAliveMass + nonStochAliveMass;
  stellarMass = stochStellarMass + nonStochStellarMass;
  stochRemnantMass = nonStochRemnantMass = 0.0;

  // Sort the stars
  sort(stars.begin(), stars.end());

  // If we were given a lifetime function, use it to draw a lifetime
  if (clf != NULL) {
    lifetime = clf->draw();
  } else {
    lifetime = constants::big;
  }

  // If we are using extinction, draw one
  if (extinct != NULL) {
    A_V = extinct->draw_AV();
    A_Vneb = A_V * extinct->draw_neb_extinct_fac();
  }

  // Initialize supernova counts
  tot_sn = 0.0;
  stoch_sn = 0;

#ifdef WINDS_ON
  // Initialize wind information
  wind_mDot = wind_pDot = wind_lum = 0.0;
#endif

  // Initialize yields
  if (yields) {
    all_yields.resize(yields->get_niso());
    stoch_yields.resize(yields->get_niso());
    if (nonStochBirthMass > 0) non_stoch_yields.resize(yields->get_niso());
  }

  // Initialize status flags for what data has been stored
  spec_set = Lbol_set = data_set = phot_set = yield_set = ew_set = false;
#ifdef WINDS_ON
  winds_set = false;
#endif
}

////////////////////////////////////////////////////////////////////////
// Constructor to build an object from a serialized buffer
////////////////////////////////////////////////////////////////////////
slug_cluster::slug_cluster(const slug_cluster_buffer *buf,
			   const slug_PDF *imf_, 
			   const slug_tracks *tracks_, 
			   const slug_specsyn *specsyn_,
			   const slug_filter_set *filters_,
			   const slug_extinction *extinct_,
			   const slug_nebular *nebular_,
			   const slug_yields *yields_,
			   const slug_line_list *lines_,
			   slug_ostreams& ostreams_,
			   const slug_PDF *clf_,
			   const bool stoch_contrib_only_) :
  ostreams(ostreams_),
  targetMass(((double *) buf)[0]),
  imf(imf_),
  clf(clf_),
  tracks(tracks_), 
  specsyn(specsyn_),
  filters(filters_),
  extinct(extinct_),
  nebular(nebular_),
  yields(yields_),
  lines(lines_),
  integ(tracks_, imf_, nullptr, ostreams_),
  stoch_contrib_only(stoch_contrib_only_)
{
  // Pull out the basic data first
  // Doubles
  const double *buf_dbl = (const double *) buf;
#ifndef WINDS_ON
  const int n_dbl = 23;
#else
  const int n_dbl = 26;
#endif
  birthMass = buf_dbl[1];
  aliveMass = buf_dbl[2];
  stochBirthMass = buf_dbl[3];
  stochAliveMass = buf_dbl[4];
  nonStochBirthMass = buf_dbl[5];
  nonStochAliveMass = buf_dbl[6];
  stochRemnantMass = buf_dbl[7];
  nonStochRemnantMass = buf_dbl[8];
  stellarMass = buf_dbl[9];
  stochStellarMass = buf_dbl[10];
  nonStochStellarMass = buf_dbl[11];
  formationTime = buf_dbl[12];
  curTime = buf_dbl[13];
  clusterAge = buf_dbl[14];
  lifetime = buf_dbl[15];
  stellarDeathMass = buf_dbl[16];
  A_V = buf_dbl[17];
  A_Vneb = buf_dbl[18];
  Lbol = buf_dbl[19];
  Lbol_ext = buf_dbl[20];
  tot_sn = buf_dbl[21];
  last_yield_time = buf_dbl[22];
#ifdef WINDS_ON
  wind_mDot = buf_dbl[23];
  wind_pDot = buf_dbl[24];
  wind_lum = buf_dbl[25];
#endif

  // Unsigned longs
  const unsigned long *buf_ul = (const unsigned long *) (buf_dbl+n_dbl);
  id = buf_ul[0];
  stoch_sn = buf_ul[1];

  // Bools; could pack these into an integer, but this is easier to
  // read, the the storage space and communication cost isn't going to
  // be affected significantly one way or the other
  const bool *buf_bool = (const bool *) (buf_ul+2);
#ifdef WINDS_ON
  const int n_bool = 8;
#else
  const int n_bool = 7;
#endif
  is_disrupted = buf_bool[0];
  data_set = buf_bool[1];
  Lbol_set = buf_bool[2];
  spec_set = buf_bool[3];
  phot_set = buf_bool[4];
  yield_set = buf_bool[5];
  ew_set = buf_bool[6];
#ifdef WINDS_ON
  winds_set = buf_bool[7];
#endif
  
  // Sizes of various vectors
  const vector<double>::size_type *buf_sz =
    (const vector<double>::size_type *) (buf_bool+n_bool);
  stars.resize(buf_sz[0]);
  dead_stars.resize(buf_sz[1]);
  L_lambda.resize(buf_sz[2]);
  phot.resize(buf_sz[3]);
  L_lambda_ext.resize(buf_sz[4]);
  phot_ext.resize(buf_sz[5]);
  L_lambda_neb.resize(buf_sz[6]);
  phot_neb.resize(buf_sz[7]);
  L_lambda_neb_ext.resize(buf_sz[8]);
  phot_neb_ext.resize(buf_sz[9]);
  all_yields.resize(buf_sz[10]);
  stoch_yields.resize(buf_sz[11]);
  non_stoch_yields.resize(buf_sz[12]);
  stardata.resize(buf_sz[13]);
  ew.resize(buf_sz[14]);

  // Contents of the vectors
  const double *buf_vec = (const double *) (buf_sz+15);
  vector<double>::size_type offset, offset0;
  for (offset = 0; offset < stars.size(); offset++)
    stars[offset] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < dead_stars.size(); offset++)
    dead_stars[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < L_lambda.size(); offset++)
    L_lambda[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < phot.size(); offset++)
    phot[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < L_lambda_ext.size(); offset++)
    L_lambda_ext[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < phot_ext.size(); offset++)
    phot_ext[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < L_lambda_neb.size(); offset++)
    L_lambda_neb[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < phot_neb.size(); offset++)
    phot_neb[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < L_lambda_neb_ext.size(); offset++)
    L_lambda_neb_ext[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < phot_neb_ext.size(); offset++)
    phot_neb_ext[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < all_yields.size(); offset++)
    all_yields[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < stoch_yields.size(); offset++)
    stoch_yields[offset-offset0] = buf_vec[offset];
  for (offset0 = offset; offset-offset0 < non_stoch_yields.size(); offset++)
    non_stoch_yields[offset-offset0] = buf_vec[offset];
  const slug_stardata *buf_stardata =
    (const slug_stardata *) (buf_vec+offset);
  for (offset = 0; offset < stardata.size(); offset++)
    stardata[offset] = buf_stardata[offset];
  for (offset =0; offset < ew.size(); offset++)
  {
    ew[offset] = buf_vec[offset];
  }
}

////////////////////////////////////////////////////////////////////////
// Copy constructor
////////////////////////////////////////////////////////////////////////
slug_cluster::slug_cluster(const slug_cluster &obj,
			   const unsigned long id_) :
  ostreams(obj.ostreams),
  targetMass(obj.targetMass),
  imf(obj.imf),
  clf(obj.clf),
  tracks(obj.tracks),
  specsyn(obj.specsyn),
  filters(obj.filters),
  extinct(obj.extinct),
  nebular(obj.nebular),
  yields(obj.yields),
  lines(obj.lines),
  integ(obj.integ)
{

  // Copy ID or set new one
  if (id_) id = id_;
  else id = obj.id;

  // Copy data
  birthMass = obj.birthMass;
  aliveMass = obj.aliveMass;
  stochBirthMass = obj.stochBirthMass;
  stochAliveMass = obj.stochAliveMass;
  stochRemnantMass = obj.stochRemnantMass;
  nonStochBirthMass = obj.nonStochBirthMass;
  nonStochAliveMass = obj.nonStochAliveMass;
  nonStochRemnantMass = obj.nonStochRemnantMass;
  formationTime = obj.formationTime;
  curTime = obj.curTime;
  clusterAge = obj.clusterAge;
  lifetime = obj.lifetime;
  stellarDeathMass = obj.stellarDeathMass;
  A_V = obj.A_V;
  A_V = obj.A_Vneb;
  Lbol = obj.Lbol;
  Lbol_ext = obj.Lbol_ext;
  tot_sn = obj.tot_sn;
  last_yield_time = obj.last_yield_time;
#ifdef WINDS_ON
  wind_mDot = obj.wind_mDot;
  wind_pDot = obj.wind_pDot;
  wind_lum = obj.wind_lum;
#endif
  stoch_sn = obj.stoch_sn;
  is_disrupted = obj.is_disrupted;
  data_set = obj.data_set;
  Lbol_set = obj.Lbol_set;
  spec_set = obj.spec_set;
  phot_set = obj.phot_set;
  yield_set = obj.yield_set;
  ew_set = obj.ew_set;
#ifdef WINDS_ON
  winds_set = obj.winds_set;
#endif
  
  // Copy vectors
  stars = obj.stars;
  dead_stars = obj.dead_stars;
  L_lambda = obj.L_lambda;
  phot = obj.phot;
  L_lambda_ext = obj.L_lambda_ext;
  phot_ext = obj.phot_ext;
  L_lambda_neb = obj.L_lambda_neb;
  phot_neb = obj.phot_neb;
  L_lambda_neb_ext = obj.L_lambda_neb_ext;
  phot_neb_ext = obj.phot_neb_ext;
  all_yields = obj.all_yields;
  stoch_yields = obj.stoch_yields;
  non_stoch_yields = obj.non_stoch_yields;
  stardata = obj.stardata;
  ew = obj.ew;
}

////////////////////////////////////////////////////////////////////////
// Routines to build and manipulate serialized buffers
////////////////////////////////////////////////////////////////////////
size_t
slug_cluster::buffer_size() const {
  // Add up the storage needed for the buffer
#ifdef WINDS_ON
  const int n_dbl = 26;
#else
  const int n_dbl = 23;
#endif
  size_t bufsize = n_dbl*sizeof(double) + 2*sizeof(unsigned long) +
    7*sizeof(bool) + 15*sizeof(vector<double>::size_type) +
    sizeof(double) * (stars.size() + dead_stars.size() +
		      L_lambda.size() + phot.size() +
		      L_lambda_ext.size() + phot_ext.size() +
		      L_lambda_neb.size() + phot_neb.size() +
		      L_lambda_neb_ext.size() + phot_neb_ext.size() +
		      all_yields.size() + stoch_yields.size() +
		      non_stoch_yields.size()) +
    sizeof(slug_stardata) * stardata.size();
  return bufsize;
}

slug_cluster_buffer *
slug_cluster::make_buffer() const {

  // Allocate the buffer
  slug_cluster_buffer *buf = malloc(buffer_size());

  // Pack the buffer
  pack_buffer(buf);
  
  // Return the buffer
  return(buf);
}

void
slug_cluster::pack_buffer(slug_cluster_buffer *buf) const {

  // Doubles
  double *buf_dbl = (double *) buf;
#ifdef WINDS_ON
  const int n_dbl = 26;
#else
  const int n_dbl = 23;
#endif
  buf_dbl[0] = targetMass;
  buf_dbl[1] = birthMass;
  buf_dbl[2] = aliveMass;
  buf_dbl[3] = stochBirthMass;
  buf_dbl[4] = stochAliveMass;
  buf_dbl[5] = nonStochBirthMass;
  buf_dbl[6] = nonStochAliveMass;
  buf_dbl[7] = stochRemnantMass;
  buf_dbl[8] = nonStochRemnantMass;
  buf_dbl[9] = stellarMass;
  buf_dbl[10] = stochStellarMass;
  buf_dbl[11] = nonStochStellarMass;
  buf_dbl[12] = formationTime;
  buf_dbl[13] = curTime;
  buf_dbl[14] = clusterAge;
  buf_dbl[15] = lifetime;
  buf_dbl[16] = stellarDeathMass;
  buf_dbl[17] = A_V;
  buf_dbl[18] = A_Vneb;
  buf_dbl[19] = Lbol;
  buf_dbl[20] = Lbol_ext;
  buf_dbl[21] = tot_sn;
  buf_dbl[22] = last_yield_time;
#ifdef WINDS_ON
  buf_dbl[23] = wind_mDot;
  buf_dbl[24] = wind_pDot;
  buf_dbl[25] = wind_lum;
#endif

  // Unsigned longs
  unsigned long *buf_ul = (unsigned long *) (buf_dbl+n_dbl);
  buf_ul[0] = id;
  buf_ul[1] = stoch_sn;

  // Bools; could pack these into an integer, but this is much easier
  // to read, and the overall storage space and cost isn't going to be
  // affected in any significant way
  bool *buf_bool = (bool *) (buf_ul+2);
#ifdef WINDS_ON
  const int n_bool = 8;
#else
  const int n_bool = 7;
#endif
  buf_bool[0] = is_disrupted;
  buf_bool[1] = data_set;
  buf_bool[2] = Lbol_set;
  buf_bool[3] = spec_set;
  buf_bool[4] = phot_set;
  buf_bool[5] = yield_set;
  buf_bool[6] = ew_set;
#ifdef WINDS_ON
  buf_bool[7] = winds_set;
#endif
  
  // Sizes of the various vectors
  vector<double>::size_type *buf_sz =
    (vector<double>::size_type *) (buf_bool+n_bool);
  buf_sz[0] = stars.size();
  buf_sz[1] = dead_stars.size();
  buf_sz[2] = L_lambda.size();
  buf_sz[3] = phot.size();
  buf_sz[4] = L_lambda_ext.size();
  buf_sz[5] = phot_ext.size();
  buf_sz[6] = L_lambda_neb.size();
  buf_sz[7] = phot_neb.size();
  buf_sz[8] = L_lambda_neb_ext.size();
  buf_sz[9] = phot_neb_ext.size();
  buf_sz[10] = all_yields.size();
  buf_sz[11] = stoch_yields.size();
  buf_sz[12] = non_stoch_yields.size();
  buf_sz[13] = stardata.size();
  buf_sz[14] = ew.size();

  // Data in the vectors
  buf_dbl = (double *) (buf_sz+15);
  vector<double>::size_type offset, offset0;
  for (offset = 0; offset < stars.size(); offset++)
    buf_dbl[offset] = stars[offset];
  for (offset0 = offset; offset-offset0 < dead_stars.size(); offset++)
    buf_dbl[offset] = dead_stars[offset-offset0];
  for (offset0 = offset; offset-offset0 < L_lambda.size(); offset++)
    buf_dbl[offset] = L_lambda[offset-offset0];
  for (offset0 = offset; offset-offset0 < phot.size(); offset++)
    buf_dbl[offset] = phot[offset-offset0];
  for (offset0 = offset; offset-offset0 < L_lambda_ext.size(); offset++)
    buf_dbl[offset] = L_lambda_ext[offset-offset0];
  for (offset0 = offset; offset-offset0 < phot_ext.size(); offset++)
    buf_dbl[offset] = phot_ext[offset-offset0];
  for (offset0 = offset; offset-offset0 < L_lambda_neb.size(); offset++)
    buf_dbl[offset] = L_lambda_neb[offset-offset0];
  for (offset0 = offset; offset-offset0 < phot_neb.size(); offset++)
    buf_dbl[offset] = phot_neb[offset-offset0];
  for (offset0 = offset; offset-offset0 < L_lambda_neb_ext.size(); offset++)
    buf_dbl[offset] = L_lambda_neb_ext[offset-offset0];
  for (offset0 = offset; offset-offset0 < phot_neb_ext.size(); offset++)
    buf_dbl[offset] = phot_neb_ext[offset-offset0];
  for (offset0 = offset; offset-offset0 < all_yields.size(); offset++)
    buf_dbl[offset] = all_yields[offset-offset0];
  for (offset0 = offset; offset-offset0 < stoch_yields.size(); offset++)
    buf_dbl[offset] = stoch_yields[offset-offset0];
  for (offset0 = offset; offset-offset0 < non_stoch_yields.size(); offset++)
    buf_dbl[offset] = non_stoch_yields[offset-offset0];
  slug_stardata *buf_stardata =
    (slug_stardata *) (buf_dbl+offset);
  for (offset = 0; offset < stardata.size(); offset++) {
    buf_stardata[offset] = stardata[offset];
  }
  
  for (offset = 0; offset < ew.size(); offset++)
  {
    buf_dbl[offset] = ew[offset-offset0]; //Check
  }
}

void
slug_cluster::free_buffer(slug_cluster_buffer *buffer) const {
  free(buffer);
}


////////////////////////////////////////////////////////////////////////
// Routine to reset the cluster
////////////////////////////////////////////////////////////////////////
void
slug_cluster::reset(bool keep_id) {

  // Get new ID
  if (!keep_id) id++;

  // Reset the time, the disruption state, and all flags
  curTime = last_yield_time = 0.0;
  is_disrupted = false;
  data_set = Lbol_set = spec_set = phot_set = yield_set = ew_set =  false;
#ifdef WINDS_ON
  winds_set = false;
#endif

  // Delete current stellar masses and data
  stars.resize(0);
  dead_stars.resize(0);
  stardata.resize(0);
#ifndef __INTEL_COMPILER
  // At this time, intel does not support this c++ function
  stars.shrink_to_fit();
  dead_stars.shrink_to_fit();
  stardata.shrink_to_fit();
#endif

  // Re-populate with stars
  stochBirthMass = stochAliveMass = stochStellarMass
    = imf->drawPopulation(targetMass, stars);

  // If the population only represents part of the mass range due to
  // restrictions on what range is being treated stochastically, be
  // sure to account for that
  nonStochBirthMass = nonStochAliveMass = nonStochStellarMass
     = targetMass * (1.0 - imf->mass_frac_restrict());

  // Reset the various mass tallies
  birthMass = stochBirthMass + nonStochBirthMass;
  aliveMass = stochAliveMass + nonStochAliveMass;
  stellarMass = stochStellarMass + nonStochStellarMass;
  stochRemnantMass = nonStochRemnantMass = 0.0;

  // Sort the stars
  sort(stars.begin(), stars.end());

  // If we were given a lifetime function, use it to draw a lifetime
  if (clf != NULL) {
    lifetime = clf->draw();
  } else {
    lifetime = constants::big;
  }

  // If we are using extinction, draw one
  if (extinct != NULL) {
    A_V = extinct->draw_AV();
    A_Vneb = A_V * extinct->draw_neb_extinct_fac();
  }

  // Reset supernova counts
  tot_sn = 0.0;
  stoch_sn = 0;

#ifdef WINDS_ON
  // Reset wind quantities
  wind_mDot = wind_pDot = wind_lum = 0.0;
#endif

  // Reset yields
  all_yields.assign(all_yields.size(), 0.0);
  stoch_yields.assign(stoch_yields.size(), 0.0);
  if (nonStochBirthMass > 0.0) 
    non_stoch_yields.assign(non_stoch_yields.size(), 0.0);
}


////////////////////////////////////////////////////////////////////////
// Advance routine. All we do is destroy all the stars whose lifetime
// is less than the current time, and set the flag to disrupted if
// we've exceeded the survival time.
////////////////////////////////////////////////////////////////////////
void
slug_cluster::advance(double time) {

  // Make sure we're not trying to go back into the past
  assert(time >= curTime);

  // Do nothing if the time hasn't changed
  if (time == curTime) return;

  // Get current age
  clusterAge = time - formationTime;

  // Reset the list of stars that died this time step
  dead_stars.resize(0);

  // Handle cases of monotonic and non-monotonic tracks differently
  if (tracks->check_monotonic()) {

    // Monotonic track case

    // Get stellar death mass corresponding to this age; save the
    // current one as a temporary in case we need it below
    stellarDeathMass = tracks->death_mass(clusterAge);

    // Traverse the list, popping off stars that have died, and adding
    // to the remnant mass tally as we go; also increment the SN count
    while (stars.size() > 0) {
      if (stars.back() > stellarDeathMass) {
	stochRemnantMass += tracks->remnant_mass(stars.back());
	if (yields)
	  if (yields->produces_sn(stars.back())) stoch_sn++;
	dead_stars.push_back(stars.back());
	stars.pop_back();
      }
      else break;
    }

  } else {

    // Non-monotonic track case
    vector<double> mass_cuts = tracks->live_mass_range(clusterAge);

    // Traverse the list to pop off stars that are more massive than
    // the upper end of the most massive "alive mass" interval; add to
    // the remnant mass total. Be careful to properly handle
    // pathological case where all stars are dead.
    while (stars.size() > 0) {
      double max_alive_mass = 0.0;
      if (mass_cuts.size() > 0) max_alive_mass = mass_cuts.back();
      if (stars.back() > max_alive_mass) {
	stochRemnantMass += tracks->remnant_mass(stars.back());
	if (yields)
	  if (yields->produces_sn(stars.back())) stoch_sn++;
	dead_stars.push_back(stars.back());
	stars.pop_back();
      }
      else break;
    }

    // Now kill off stars outside every other alive mass interval
    int interval_ptr = mass_cuts.size()-2;
    int starptr = stars.size()-1;
    while (interval_ptr >= 0) {

      // Stop if we reach the bottom of the star list
      if (starptr < 0) break;

      // Find the lowest mass star that is still alive
      while (stars[starptr] > mass_cuts[interval_ptr]) {
	starptr--;
	if (starptr < 0) break;
      }

      // If we've reached the bottom of the star list, and this star
      // is still alive, we're done
      if (starptr == -1) break;

      // If this star is below the minimum mass star we're killing
      // off, go to next iteration
      if (interval_ptr > 0) {
	if (stars[starptr] < mass_cuts[interval_ptr-1]) {
	  interval_ptr -= 2;
	  continue;
	}
      }

      // If we're here, starptr points to a star that should be dead
      // at this age. We now proceed through the star list to find
      // the minimum mass stars that is still alive.
      int starptr2 = starptr;
      for (; starptr2 > 0; starptr2--)
	if (stars[starptr2-1] < mass_cuts[interval_ptr-1]) break;

      // The interval from starptr2 to starptr represents stars whose
      // masses are such that they are dead. First add their
      // contribution to the remnant mass, then remove them from the
      // list of stars. Note that we need to add 1 to starptr because
      // the c++ vector erase method excludes the last element.
      for (int i=starptr2; i<=starptr; i++) {
	stochRemnantMass += tracks->remnant_mass(stars[i]);
	if (yields)
	  if (yields->produces_sn(stars[i])) stoch_sn++;
	dead_stars.push_back(stars[i]);
      }
      stars.erase(stars.begin()+starptr2, stars.begin()+starptr+1);

      // Move the death mass point and star pointer
      interval_ptr -= 2;
      starptr = starptr2 - 1;

    }

    // Last step: if the minimum alive mass is non-zero, then we need
    // to kill all stars smaller than the smallest alive mass
    if (mass_cuts[0] > 0.0) {

      // Find all stars smaller than the minimum alive mass
      for (starptr=0; (unsigned int) starptr<stars.size(); starptr++)
	if (stars[starptr] > mass_cuts[0]) break;

      // Kill those stars, adding to the remnant mass
      for (int i=0; i<starptr; i++) {
	stochRemnantMass += tracks->remnant_mass(stars[i]);
	if (yields)
	  if (yields->produces_sn(stars[i])) stoch_sn++;
	dead_stars.push_back(stars[i]);
      }
      stars.erase(stars.begin(), stars.begin()+starptr);

    }
  }

  // Flag if we're disrupted
  if (clusterAge > lifetime) is_disrupted = true;

  // Set current time
  curTime = time;

  // Mark that data are not current
  data_set = spec_set = Lbol_set = phot_set = yield_set = ew_set = false;
#ifdef WINDS_ON
  winds_set = false;
#endif

  // Update all the stellar data to the new isochrone
  set_isochrone();
  
  // Recompute the present-day masses of the surviving stochastic
  // stars using the new isochrone. We sum the masses of stars below
  // the minimum track mass assuming they have no mass loss, and for
  // all other stars we use the present-day mass from the tracks.
  stochAliveMass = 0.0;
  for (vector<double>::size_type i=0; i<stars.size(); i++) { 
    if (stars[i] >= tracks->min_mass()) {
      break;
    }
    stochAliveMass += stars[i];
  }
  for (vector<slug_stardata>::size_type i=0; i<stardata.size(); i++)
    stochAliveMass += pow(10.0, stardata[i].logM);

  // Now do the same calculation for the non-stochastic stars
  nonStochAliveMass = 0.0;
  if (imf->get_xStochMin() > imf->get_xMin()) {
    if (imf->get_xMin() < tracks->min_mass()) {
      nonStochAliveMass += targetMass * 
	imf->mass_frac(imf->get_xMin(), 
		       min(tracks->min_mass(), imf->get_xStochMin()));
    }
    nonStochAliveMass += integ.integrate(targetMass, curTime-formationTime,
					 boost::bind(cluster::curMass, _1));
  }

  // Recompute the remnant mass from the non-stochastic stars; note a
  // bit of fancy footwork here: the slug_tracks::remnant_mass
  // function is overloaded, so we need to static_cast to the version
  // of it we want before passing to boost::bind
  nonStochRemnantMass = 
    integ.integrate_nt(targetMass, 
		       curTime-formationTime,
		       boost::bind(static_cast<double (slug_tracks::*)
				   (const double, const double,
				    const double) const> 
				   (&slug_tracks::remnant_mass), 
				   tracks, _1, _2, tracks::null_metallicity));

  // Recompute the number of non-stochastic supernovae and thus total
  // supernovae; only do this if we have yields, because this is where
  // the information on SN production is stored
  if (yields && !stoch_contrib_only) {
    tot_sn = stoch_sn;
    if (imf->has_stoch_lim()) {

      // Some stuff we'll need
      const vector<double>& sn_mass_range = yields->sn_mass_range();
      double mbar = imf->expectationVal();

      if (tracks->check_monotonic()) {

	// Monotonic tracks, so just integrate the non-stochastic part
	// of the IMF over the range of masses that has supernovae and
	// is above the larger of the stellar death mass and below the
	// maximum non-stochastic mass
	double m_stop = min(imf->get_xStochMin(), sn_mass_range.back());
	double m = stellarDeathMass;
	bool has_sn = yields->produces_sn(m);
	double m_next = min(imf->get_xStochMin(), sn_mass_range.front());
	vector<double>::size_type ptr = 0;
	while (m < m_stop) {
	  if (has_sn) tot_sn += targetMass * imf->integral(m, m_next) / mbar;
	  has_sn = !has_sn;
	  m = m_next;
	  m_next = min(m_stop, sn_mass_range[ptr+1]);
	  ptr++;
	}

      } else {

	// Non-monotonic tracks; we need to find the intersection
	// between the mass intervals that have died, the mass
	// intervals that have supernovae, and mass intervals that are
	// being treated non-stochastically, so that we can integrate
	// over them. The variable sn_nonstoch_range will be a vector
	// with an even number of elements that defines the ranges
	// over which to integrate.
	vector<double> sn_nonstoch_range;

	// First step: get the live mass range
	vector<double> live_mass_range;
	live_mass_range = tracks->live_mass_range(curTime-formationTime);

	// Second step: get the non-stochastic mass range
	vector<double> non_stoch_range;
	if (imf->get_xMin() != imf->get_xStochMin()) {
	  non_stoch_range.push_back(imf->get_xMin());
	  non_stoch_range.push_back(imf->get_xStochMin());
	}
	if (imf->get_xMax() != imf->get_xStochMax()) {
	  non_stoch_range.push_back(imf->get_xStochMax());
	  non_stoch_range.push_back(imf->get_xMax());
	}

	// Third step: construct the intervals that result from the
	// intersection of the inverse of the live mass range, the SN
	// mass range, and the non-stochastic range
	vector<double>::size_type live_ptr, ns_ptr, sn_ptr;
	bool alive, has_sn, non_stoch;
	if (live_mass_range[0] < non_stoch_range[0]) {
	  if (live_mass_range[0] < sn_mass_range[0]) {
	    alive = true;
	    non_stoch = false;
	    has_sn = false;
	    live_ptr = 1;
	    ns_ptr = sn_ptr = 0;
	  } else {
	    alive = false;
	    non_stoch = false;
	    has_sn = true;
	    sn_ptr = 1;
	    live_ptr = ns_ptr = 0;
	  }
	} else {
	  if (non_stoch_range[0] < sn_mass_range[0]) {
	    alive = false;
	    non_stoch = true;
	    has_sn = false;
	    ns_ptr = 1;
	    live_ptr = sn_ptr = 0;
	  } else {
	    alive = false;
	    non_stoch = false;
	    has_sn = true;
	    sn_ptr = 1;
	    live_ptr = ns_ptr = 0;
	  }
	}
	while (ns_ptr != non_stoch_range.size()) {
	  // Figure out which pointer to advance
	  if (live_ptr == live_mass_range.size()) {
	    // live_ptr is at end of list
	    if (sn_ptr == sn_mass_range.size()) {
	      // sn_ptr is at end of list
	      
	      // advance ns_ptr; if we are inside an interval of dead,
	      // sn-producing, non-stochastic mass, end the
	      // interval; flip the stochasticity bit; if we are at the
	      // start of a dead, sn-producing, non-stochastic interval,
	      // start it
	      if (!alive && has_sn)
		sn_nonstoch_range.push_back(non_stoch_range[ns_ptr]);
	      ns_ptr++;
	      non_stoch = !non_stoch;
	    } else {
	      // decide if we should advance sn_ptr or ns_ptr
	      if (sn_mass_range[sn_ptr] < non_stoch_range[ns_ptr]) {
		// advance sn_ptr
		if (!alive && non_stoch)
		  sn_nonstoch_range.push_back(sn_mass_range[sn_ptr]);
		sn_ptr++;
		has_sn = !has_sn;
	      } else {
		// advance ns_ptr
		if (!alive && has_sn)
		  sn_nonstoch_range.push_back(non_stoch_range[ns_ptr]);
		ns_ptr++;
		non_stoch = !non_stoch;
	      }
	    }
	  } else {
	    // live_ptr is not at end of list
	    if (sn_ptr == sn_mass_range.size()) {
	      // sn_ptr is at end of list; decie if we should advance
	      // live_ptr or ns_ptr
	      if (live_mass_range[live_ptr] < non_stoch_range[ns_ptr]) {
		// advance live_ptr
		if (non_stoch && has_sn)
		  sn_nonstoch_range.push_back(live_mass_range[live_ptr]);
		live_ptr++;
		alive = !alive;
	      } else {
		// advance ns_ptr
		if (!alive && has_sn)
		  sn_nonstoch_range.push_back(non_stoch_range[ns_ptr]);
		ns_ptr++;
		non_stoch = !non_stoch;
	      }
	    } else {
	      // no pointers are at end of list, so compare all three
	      if (live_mass_range[live_ptr] < non_stoch_range[ns_ptr]) {
		if (live_mass_range[live_ptr] < sn_mass_range[sn_ptr]) {
		  // live_ptr is next to advance
		  if (non_stoch && has_sn)
		    sn_nonstoch_range.push_back(live_mass_range[live_ptr]);
		  live_ptr++;
		  alive = !alive;
		} else {
		  // sn_ptr is next to advance
		  if (!alive && non_stoch)
		    sn_nonstoch_range.push_back(sn_mass_range[sn_ptr]);
		  sn_ptr++;
		  has_sn = !has_sn;
		}
	      } else {
		if (non_stoch_range[ns_ptr] < sn_mass_range[sn_ptr]) {
		  // ns_ptr is next to advance
		  if (!alive && has_sn)
		    sn_nonstoch_range.push_back(non_stoch_range[ns_ptr]);
		  ns_ptr++;
		  non_stoch = !non_stoch;
		} else {
		  // sn_ptr is next to advance
		  if (!alive && non_stoch)
		    sn_nonstoch_range.push_back(sn_mass_range[sn_ptr]);
		  sn_ptr++;
		  has_sn = !has_sn;
		}
	      }
	    }
	  }
	}

	// Fourth step: integrate over the mass intervals we have found
	// to decide how many non-stochastic supernovae there are
	for (vector<double>::size_type i=0;
	     i<sn_nonstoch_range.size(); i+=2) {
	  tot_sn += targetMass *
	    imf->integral(sn_nonstoch_range[i],
			  sn_nonstoch_range[i+1]) / mbar;
	}
      }
    }
  }

  // Compute the new alive and total stellar masses
  aliveMass = nonStochAliveMass + stochAliveMass;
  stellarMass = aliveMass + stochRemnantMass + nonStochRemnantMass;
}


////////////////////////////////////////////////////////////////////////
// Routine to get stellar data at this time
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_isochrone() {

  // Do nothing if already set; if not set, refresh stellar data and
  // flag that it is now current
  if (data_set) return;
  stardata = tracks->get_isochrone(curTime-formationTime, stars);
  data_set = true;
}


////////////////////////////////////////////////////////////////////////
// Routine to get bolometric luminosity at this time
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_Lbol() {

  // Do nothing if already set
  if (Lbol_set) return;

  // Initialize
  Lbol = 0.0;

  // Stochastic stars part
  if (stars.size() > 0) {

    // Refresh the stellar data
    set_isochrone();

    // Add bolometric luminosity from stochastic stars
    for (unsigned int i=0; i<stardata.size(); i++)
      Lbol += pow(10.0, stardata[i].logL);
  }

  // Non-stochastic part
  if (imf->has_stoch_lim() && !stoch_contrib_only)
    Lbol += specsyn->get_Lbol_cts(birthMass, curTime-formationTime);

  // If using extinction, to compute Lbol we need to first compute the
  // spectrum, then integrate it
  if (extinct != NULL) {
    set_spectrum();
    Lbol_ext = int_tabulated::
      integrate(extinct->lambda(), L_lambda_ext) / constants::Lsun;
  }

  // Flag that things are set
  Lbol_set = true;
}



////////////////////////////////////////////////////////////////////////
// Spectral synthesis routine. Note that this routine also sets Lbol
// in the process because the extra cost of computing it is negligible.
////////////////////////////////////////////////////////////////////////
void
slug_cluster::set_spectrum() {

  // Do nothing if already set
  if (spec_set) return;

  // Initialize
  vector<double>::size_type nl = specsyn->n_lambda();
  L_lambda.assign(nl, 0.0);
  Lbol = 0.0;

  // Stochastic stars part
  if (stars.size() > 0) {

    // Refresh the stellar data
    set_isochrone();

    // Get spectrum for stochastic stars
    // Obtain rectified spectrum as well if it is present
    if (specsyn->get_rectify())
    {
      // Get wavelengths for rectified spectrum
      recspec_wl = specsyn->get_recspec_wl();
      
      // Resize recspec to wl size to store the spectrum
      recspec.resize(recspec_wl.size());
      // Return the spectrum
      L_lambda = specsyn->get_spectrum(stardata,recspec);         

      
      
    }
    else
    {
      L_lambda = specsyn->get_spectrum(stardata);
    }
    // Add bolometric luminosity from stochastic stars
    for (unsigned int i=0; i<stardata.size(); i++)
      Lbol += pow(10.0, stardata[i].logL);
  }

  // Non-stochastic part
  if (imf->has_stoch_lim() && !stoch_contrib_only) {
    double Lbol_tmp;
    vector<double> spec;
    specsyn->get_spectrum_cts(birthMass, curTime-formationTime, spec, 
			      Lbol_tmp);
    for (unsigned int i=0; i<nl; i++) L_lambda[i] += spec[i];
    Lbol += Lbol_tmp;
  }

  // If using nebular emission, compute the stellar+nebular spectrum
  vector<double> L_lambda_neb_only;
  if (nebular != NULL) {
    L_lambda_neb_only = nebular->get_neb_spec(L_lambda, this->get_age());
    L_lambda_neb = nebular->add_stellar_nebular_spec(L_lambda,
						     L_lambda_neb_only);
  }
  
  // If using extinction, compute the extincted spectrum and the
  // bolometric luminosity after extinction is applied
  if (extinct != NULL) {
    L_lambda_ext = extinct->spec_extinct(A_V, L_lambda);
    Lbol_ext = int_tabulated::
      integrate(extinct->lambda(), L_lambda_ext) 
      / constants::Lsun;
    if (nebular != NULL) {
      // Procedure depends on if nebular and stellar extinctions
      // differ
      if (extinct->excess_neb_extinct()) {
	// Apply nebular extinction to nebular portion
	vector<double> L_lambda_neb_only_ext =
	  extinct->spec_extinct_neb(A_Vneb, L_lambda_neb_only);
	// Add extincted nebular and stellar spectra; note that we
	// must provide the offsets, because the the extincted
	// spectrum may be truncated relative to the full stellar and
	// nebular grids
	L_lambda_neb_ext =
	  nebular->add_stellar_nebular_spec(L_lambda_ext,
					    L_lambda_neb_only_ext,
					    extinct->off(),
					    extinct->off_neb());
      } else {
	// We are using the same extinction for nebular stellar
	// emission; just apply that extinction to both
	L_lambda_neb_ext = 
	  extinct->spec_extinct_neb(A_V, L_lambda_neb);
      }
    }
  }

  // Flag that things are set
  spec_set = Lbol_set = true;
}


////////////////////////////////////////////////////////////////////////
// Photometry calculation routine
////////////////////////////////////////////////////////////////////////
void
slug_cluster::set_photometry() {

  // Do nothing if already set
  if (phot_set) return;

  // Compute the spectrum
  set_spectrum();

  // Compute photometry
  phot = filters->compute_phot(specsyn->lambda(), L_lambda);

  // If any of the photometric values are -big, that indicates that we
  // want the bolometric luminosity, so insert that
  for (vector<double>::size_type i=0; i<phot.size(); i++)
    if (phot[i] == -constants::big) phot[i] = Lbol;

  // Repeat for stellar+nebular spectrum
  if (nebular != NULL) {
    phot_neb = filters->compute_phot(nebular->lambda(), L_lambda_neb);
    // Special cases: force ionizing luminosity to be zero exactly for
    // this spectrum, and bolometric luminosity to be exactly the same
    // as for the non-nebular case
    for (vector<double>::size_type i=0; i<phot_neb.size(); i++) {
      if (phot_neb[i] == -constants::big) phot_neb[i] = Lbol;
      if (filters->get_filter(i)->photon_filter() &&
	  (filters->get_filter(i)->get_wavelength_max()
	   <= constants::lambdaHI))
	phot_neb[i] = 0.0;
    }
  }

  // If using extinction, compute photometry on the extincted
  // spectrum; in the process, be careful to mask filters whose
  // response curve doesn't overlap with the extincted spectrum
  if (extinct != NULL) {
    phot_ext = filters->compute_phot(extinct->lambda(), 
				     L_lambda_ext);
    for (vector<double>::size_type i=0; i<phot_ext.size(); i++) {
      if (phot_ext[i] == -constants::big) {
	phot_ext[i] = Lbol_ext;
      } else if (filters->get_filter(i)->photon_filter() &&
		 (filters->get_filter(i)->get_wavelength_min() >
		  extinct->lambda().back())) {
	phot_ext[i] = nan("");
      } else if ((filters->get_filter(i)->get_wavelength_min() <
		  extinct->lambda().front()) ||
		 (filters->get_filter(i)->get_wavelength_max() >
		  extinct->lambda().back())) {
	phot_ext[i] = nan("");
      }
    }

    // Repeat for stellar+nebular spectrum
    if (nebular != NULL) {
      phot_neb_ext = filters->compute_phot(extinct->lambda_neb(), 
					   L_lambda_neb_ext);
      for (vector<double>::size_type i=0; i<phot_neb_ext.size(); i++) {
	if (phot_neb_ext[i] == -constants::big) {
	  phot_neb_ext[i] = Lbol_ext;
	} else if (filters->get_filter(i)->photon_filter() &&
		   (filters->get_filter(i)->get_wavelength_min() >
		    extinct->lambda_neb().back())) {
	  phot_neb_ext[i] = nan("");
	} else if ((filters->get_filter(i)->get_wavelength_min() <
		    extinct->lambda_neb().front()) ||
		   (filters->get_filter(i)->get_wavelength_max() >
		    extinct->lambda_neb().back())) {
	  phot_neb_ext[i] = nan("");
	} else if (filters->get_filter(i)->photon_filter() &&
		   (filters->get_filter(i)->get_wavelength_max()
		    <= constants::lambdaHI)) {
	  phot_neb_ext[i] = 0.0;
	}
      }
    }
  }
  
  // Flag that the photometry is set
  phot_set = true;
}  

////////////////////////////////////////////////////////////////////////
// Routine to compute yields at this time
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_yield() {

  // Do nothing if already set
  if (yield_set) return;

  // Make unstable isotopes decay
  if (!yields->no_decay) {
    const vector<const isotope_data *>& isotopes = yields->get_isotopes();
    for (vector<double>::size_type i=0; i<stoch_yields.size(); i++) {
      if (!isotopes[i]->stable()) {
	stoch_yields[i] *= exp(-(curTime-last_yield_time)/isotopes[i]->ltime());
      }
    }
    last_yield_time = curTime;
  }
  
  // Add yield from stars that died; for unstable isotopes, be sure to
  // apply the correct amount of radioactive decay between the star's
  // death time and the current time
  if (dead_stars.size() > 0) {
    vector<double> star_yields;
    if (!yields->no_decay) {
      vector<double> decay_time(dead_stars.size());
      for (vector<double>::size_type i=0; i<dead_stars.size(); i++)
	decay_time[i] = curTime - formationTime -
	  tracks->star_lifetime(dead_stars[i]);
      star_yields = yields->yield(dead_stars, decay_time);
    } else {
      star_yields = yields->yield(dead_stars);
    }
    for (vector<double>::size_type i=0; i<star_yields.size(); i++)
      stoch_yields[i] += star_yields[i];
  }

  // Do we have non-stochastic stars?
  if (!imf->has_stoch_lim() || stoch_contrib_only) {

    // No, so just copy stochastic yields to total yields
    for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      all_yields[i] = stoch_yields[i];
    
  } else {

    // We do have non-stochastic stars
    
    // Reset the non-stochastic yield
    non_stoch_yields.assign(non_stoch_yields.size(), 0.0);

    // The next step is to find the intersection between the mass
    // intervals that have died and mass intervals that are being
    // treated non-stochastically, so that we can integrate over
    // them. The variable dead_nonstoch_range will be a vector with an
    // even number of elements that defines the ranges over which to
    // integrate.
    vector<double> dead_nonstoch_range;

    // First step: get the live mass range
    vector<double> live_mass_range;
    if (tracks->check_monotonic()) {
      live_mass_range.push_back(0.0);
      live_mass_range.push_back(stellarDeathMass);
    } else {
      live_mass_range = tracks->live_mass_range(curTime-formationTime);
    }

    // Second step: get the non-stochastic mass range
    vector<double> non_stoch_range;
    if (imf->get_xMin() != imf->get_xStochMin()) {
      non_stoch_range.push_back(imf->get_xMin());
      non_stoch_range.push_back(imf->get_xStochMin());
    }
    if (imf->get_xMax() != imf->get_xStochMax()) {
      non_stoch_range.push_back(imf->get_xStochMax());
      non_stoch_range.push_back(imf->get_xMax());
    }

    // Third step: construct the intervals that result from the
    // intersection of the inverse of the first set (the live mass
    // range) and the second set (the non-stochastic range)
    vector<double>::size_type live_ptr, ns_ptr;
    bool alive, non_stoch;
    if (live_mass_range[0] == non_stoch_range[0]) {
      alive = true;
      non_stoch = true;
      live_ptr = ns_ptr = 1;
    } else if (live_mass_range[0] > non_stoch_range[0]) {
      alive = false;
      non_stoch = true;
      dead_nonstoch_range.push_back(non_stoch_range[0]);
      live_ptr = 0;
      ns_ptr = 1;
    } else {
      alive = true;
      non_stoch = false;
      live_ptr = 1;
      ns_ptr = 0;
    }
    while (ns_ptr != non_stoch_range.size()) {
      // Figure out which pointer to advance
      if (live_ptr == live_mass_range.size()) {
	// live_ptr is at end of list, so advance ns_ptr; if we are
	// inside an interval of dead, non-stochastic mass, end the
	// interval; flip the stochasticity bit; if we are at the
	// start of a dead, non-stochastic interval, start it
	if (!alive)
	  dead_nonstoch_range.push_back(non_stoch_range[ns_ptr]);
	ns_ptr++;
	non_stoch = !non_stoch;
      } else if (ns_ptr == non_stoch_range.size()) {
	// ns_ptr is at end of list, so advance live_ptr; if we are
	// inside an interval of dead, non-stochastic mass, end the
	// interval; flip the alive bit
	if (non_stoch)
	  dead_nonstoch_range.push_back(live_mass_range[live_ptr]);
	live_ptr++;
	alive = !alive;
      } else if (live_mass_range[live_ptr] == non_stoch_range[ns_ptr]) {
	// Neither point is at end of list, and the next points are
	// identical; advance both pointers, flip both bits, and, if we are
	// inside an interval of dead, non-stochastic mass, end the
	// interval
	if (alive != non_stoch)
	  dead_nonstoch_range.push_back(live_mass_range[live_ptr]);
	live_ptr++;
	ns_ptr++;
	alive = !alive;
	non_stoch = !non_stoch;
      } else if (live_mass_range[live_ptr] < non_stoch_range[ns_ptr]) {
	// Neither pointer is at end of list, but the next mass point
	// we hit is a transition from alive to not alive; advance the
	// alive pointer, flip the alive bit, and start and end
	// intervals
	if (non_stoch)
	  dead_nonstoch_range.push_back(live_mass_range[live_ptr]);
	live_ptr++;
	alive = !alive;
      } else {
	// The next mass point we hit is a transition from
	// non-stochastic to stochastic
	if (!alive)
	  dead_nonstoch_range.push_back(non_stoch_range[ns_ptr]);
	ns_ptr++;
	non_stoch = !non_stoch;
      }
    }

    // Now integrate the IMF-weighted yield over the non-stochastic
    // dead star mass range. Note that this doesn't quite treat decay
    // of unstable isotopes correctly. To be fixed later!
    for (vector<double>::size_type i=0;
	 i<dead_nonstoch_range.size(); i+=2) {
      for (vector<double>::size_type j=0;
	   j<non_stoch_yields.size(); j++) {
	non_stoch_yields[j] +=
	  integ.integrate_nt_lim(targetMass,
				 dead_nonstoch_range[i],
				 dead_nonstoch_range[i+1],
				 boost::bind(cluster::yield, _1,
					     yields, j));
      }
    }

    // Sum stochastic and non-stochastic yields
    for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      all_yields[i] = stoch_yields[i] + non_stoch_yields[i];
  }

  // Record that yields are now current
  yield_set = true;
}

#ifdef WINDS_ON
////////////////////////////////////////////////////////////////////////
// Routine to compute winds at this time
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_winds() {
  
  // Do nothing if already set
  if (winds_set) return;

  // Initialize
  wind_mDot = wind_pDot = wind_lum = 0.0;

  // Stochastic stars part
  if (stars.size() > 0) {

    // Refresh the stellar data
    set_isochrone();

    // Add contribution from stochastic stars
    for (vector<double>::size_type i=0; i<stardata.size(); i++) {
      wind_mDot += stardata[i].mDot;
      wind_pDot += stardata[i].mDot * stardata[i].vWind;
      wind_lum += 0.5 * (stardata[i].mDot*constants::Msun / constants::yr) *
	(stardata[i].vWind*1e5)*(stardata[i].vWind*1e5) / constants::Lsun;
    }
  }

  // Non-stochastic part
  if (imf->has_stoch_lim() && !stoch_contrib_only) {
    wind_mDot += integ.integrate(targetMass, curTime-formationTime,
				 boost::bind(cluster::star_mDot, _1));
    wind_pDot += integ.integrate(targetMass, curTime-formationTime,
				 boost::bind(cluster::star_pDot, _1));
    wind_lum += integ.integrate(targetMass, curTime-formationTime,
				 boost::bind(cluster::star_wind_lum, _1));
  }

  // Flag that things are set
  winds_set = true;

}
#endif

////////////////////////////////////////////////////////////////////////
// Routines to return the spectrum, without or without nebular
// contributions and extinction, and with our without the wavelength
// table
////////////////////////////////////////////////////////////////////////
const vector<double> &slug_cluster::get_spectrum() {
  set_spectrum();
  return L_lambda;
}

const vector<double> &slug_cluster::get_spectrum_neb() {
  set_spectrum();
  return L_lambda_neb;
}

const vector<double> &slug_cluster::get_spectrum_extinct() {
  set_spectrum();
  return L_lambda_ext;
}

const vector<double> &slug_cluster::get_spectrum_neb_extinct() { 
  set_spectrum(); 
  return L_lambda_neb_ext;
}

void slug_cluster::get_spectrum(vector<double> &lambda_out, 
				vector<double> &L_lambda_out,
				bool rest) { 
  set_spectrum(); 
  lambda_out = specsyn->lambda(rest);
  L_lambda_out = L_lambda;
}

void slug_cluster::get_spectrum_neb(vector<double> &lambda_out, 
				    vector<double> &L_lambda_out,
				    bool rest) {
  set_spectrum(); 
  lambda_out = nebular->lambda(rest); 
  L_lambda_out = L_lambda_neb;
}

void slug_cluster::get_spectrum_extinct(vector<double> &lambda_out, 
					vector<double> &L_lambda_out,
					bool rest) { 
  set_spectrum(); 
  lambda_out = extinct->lambda(rest); 
  L_lambda_out = L_lambda_ext;
}

void slug_cluster::get_spectrum_neb_extinct(vector<double> &lambda_out, 
					    vector<double> &L_lambda_out,
					    bool rest) { 
  set_spectrum(); 
  lambda_out = extinct->lambda_neb(rest); 
  L_lambda_out = L_lambda_neb_ext;
}


////////////////////////////////////////////////////////////////////////
// Routines to return the photometry
////////////////////////////////////////////////////////////////////////

const vector<double> &slug_cluster::get_photometry() { 
  set_photometry();
  return phot;
}

const vector<double> &slug_cluster::get_photometry_neb() {
  set_photometry();
  return phot_neb;
}

const vector<double> &slug_cluster::get_photometry_extinct() {
  set_photometry();
  return phot_ext;
}

const vector<double> &slug_cluster::get_photometry_neb_extinct() {
  set_photometry();
  return phot_neb_ext;
}

////////////////////////////////////////////////////////////////////////
// Routine to get the yield
////////////////////////////////////////////////////////////////////////
const vector<double> &slug_cluster::get_yield() {
  set_yield();
  return all_yields;
}

////////////////////////////////////////////////////////////////////////
// Routine to set the equivalent widths
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_ew()
{
  // Do nothing if already set
  if (ew_set) 
  {
    return;
  }
  
  // Set spectrum
  set_spectrum();
  
  // Safety check - if we have no rectified spectrum to operate
  // on, kill slug
  if (recspec.size() == 0)
  {
    ostreams.slug_err_one << "Rectified spectrum missing. Bailing out."
                          << endl;
    //bailout(1);  
    exit(1);                        
  }
  
  // Compute equivalent widths
  ew = lines->compute_ew(recspec_wl, recspec);  
  
  // Set the flag to say we have calculated the equivalent widths
  ew_set = true;
}
////////////////////////////////////////////////////////////////////////
// Routine to get the equivalent widths
////////////////////////////////////////////////////////////////////////
const vector<double> &slug_cluster::get_ew()
{
  set_ew();
  return ew;
}
////////////////////////////////////////////////////////////////////////
// Routines to clear data
////////////////////////////////////////////////////////////////////////
void slug_cluster::clear_spectrum() {
  L_lambda.resize(0); 
  L_lambda_ext.resize(0); 
  L_lambda_neb.resize(0);
  L_lambda_neb_ext.resize(0);
  phot.resize(0);
  phot_neb.resize(0);
  phot_ext.resize(0);
  phot_neb_ext.resize(0);
  spec_set = false;
  phot_set = false;
  ew_set = false;
  ew.resize(0); 
  recspec.resize(0);
  recspec_wl.resize(0);
}

////////////////////////////////////////////////////////////////////////
// Output physical properties
////////////////////////////////////////////////////////////////////////
void
slug_cluster::write_prop(ofstream& outfile, const outputMode out_mode,
			 const unsigned long trial,
			 bool cluster_only, const std::vector<double>& imfvp) const {

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific
	    << setw(11) << right << id << "   "
	    << setw(11) << right << curTime << "   "
	    << setw(11) << right << formationTime << "   "
	    << setw(11) << right << lifetime << "   "
	    << setw(11) << right << targetMass << "   "
	    << setw(11) << right << birthMass << "   "
	    << setw(11) << right << aliveMass << "   "
	    << setw(11) << right << stellarMass << "   "
	    << setw(11) << right << stars.size() << "   ";
    if (stars.size() > 0)
      outfile << setw(11) << right << stars[stars.size()-1];
    else
      outfile << setw(11) << right << 0.0;
    if (extinct != NULL) {
      outfile << "   " << setw(11) << right << A_V;
      if (extinct->excess_neb_extinct()) {
	outfile << "   " << setw(11) << right << A_Vneb;
      }
    }

    // Variable parameter block
    if (imfvp.size()>0) {
      // Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imfvp.size();p++)
      {
        outfile << "   " << setw(11) << right << imfvp[p];
      }
    
    }
    
    outfile << endl;
  
  } else if (out_mode == BINARY) {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) &formationTime, sizeof formationTime);
    outfile.write((char *) &lifetime, sizeof lifetime);
    outfile.write((char *) &targetMass, sizeof targetMass);
    outfile.write((char *) &birthMass, sizeof birthMass);
    outfile.write((char *) &aliveMass, sizeof aliveMass);
    outfile.write((char *) &stellarMass, sizeof stellarMass);
    vector<double>::size_type n = stars.size();
    outfile.write((char *) &n, sizeof n);
    if (stars.size() > 0)
      outfile.write((char *) &(stars[stars.size()-1]), 
		    sizeof stars[stars.size()-1]);
    else {
      double mstar = 0.0;
      outfile.write((char *) &mstar, sizeof mstar);
    }
    if (extinct != NULL) {
      outfile.write((char *) &A_V, sizeof A_V);
      if (extinct->excess_neb_extinct()) {
	outfile.write((char *) &A_Vneb, sizeof A_Vneb);
      }
    }
    
    // Write out variable parameter values
    if (imfvp.size()>0) {
      // Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imfvp.size(); p++)
      {
        outfile.write((char *) &imfvp[p], sizeof imfvp[p]);
      }
    }   
    
  }
}


////////////////////////////////////////////////////////////////////////
// Output physical properties in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::write_prop(fitsfile *out_fits, unsigned long trial,
                          const std::vector<double>& imfvp) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write a new entry
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, 1, &formationTime,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 5, nrows+1, 1, 1, &lifetime,
		 &fits_status);
  double tmass = targetMass; // Needed to avoid compiler complaint
  fits_write_col(out_fits, TDOUBLE, 6, nrows+1, 1, 1, &tmass,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 7, nrows+1, 1, 1, &birthMass,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 8, nrows+1, 1, 1, &aliveMass,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 9, nrows+1, 1, 1, &stellarMass,
		 &fits_status);
  vector<double>::size_type n = stars.size();
  fits_write_col(out_fits, TULONG, 10, nrows+1, 1, 1, &n,
		 &fits_status);
  double mstar;
  if (n>0) mstar = stars.back();
  else mstar = 0.0;
  fits_write_col(out_fits, TDOUBLE, 11, nrows+1, 1, 1, &mstar,
		 &fits_status);
		 
  int colnum = 11;			 
  if (extinct != NULL) {
    fits_write_col(out_fits, TDOUBLE, 12, nrows+1, 1, 1, &A_V, &fits_status);
    colnum++;
    if (extinct->excess_neb_extinct()) {
      fits_write_col(out_fits, TDOUBLE, 13, nrows+1, 1, 1, &A_Vneb,
		     &fits_status);
      colnum++;
    }
  }
   
  if (imfvp.size()>0) {
    // Loop over the variable parameters
    for (vector<double>::size_type p = 0; p<imfvp.size(); p++) {
      colnum++;
      double vp_p=imfvp[p];
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1, &vp_p,
		  &fits_status);
    }
  }
}
#endif


////////////////////////////////////////////////////////////////////////
// Output spectrum
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_spectrum(ofstream& outfile, const outputMode out_mode,
	       const unsigned long trial,
	       bool cluster_only) {

  // Make sure information is current
  if (!spec_set) set_spectrum();

  if (out_mode == ASCII) {
    vector<double> lambda, L_lambda_star, L_lambda_star_ext;
    if (nebular == NULL) {
      lambda = specsyn->lambda();
      L_lambda_star = L_lambda;
      if (extinct != NULL) L_lambda_star_ext = L_lambda_ext;
    } else {
      // Need to output all data using nebular wavelength table in
      // this case, which means we need to interpolate the purely
      // stellar spectra onto it
      lambda = nebular->lambda();
      L_lambda_star = nebular->interp_stellar(L_lambda);
      if (extinct != NULL) 
	L_lambda_star_ext = nebular->interp_stellar(L_lambda_ext, 
						    extinct->off());
    }
    for (unsigned int i=0; i<lambda.size(); i++) {
      outfile << setprecision(5) << scientific 
	      << setw(11) << right << id << "   "
	      << setw(11) << right << curTime << "   "
	      << setw(11) << right << lambda[i] << "   "
	      << setw(11) << right << L_lambda_star[i];
      if (nebular != NULL)
	outfile << "   "
		<< setw(11) << right << L_lambda_neb[i];
      if (extinct != NULL) {
	int j;
	if (nebular == NULL) j = i - extinct->off();
	else j = i - extinct->off_neb();
	if ((j >= 0) && ((unsigned int) j < L_lambda_star_ext.size())) {
	  outfile << "   "
		  << setw(11) << right << L_lambda_star_ext[j];
	  if (nebular != NULL)
	    outfile << "   "
		    << setw(11) << right << L_lambda_neb_ext[j];
	}
      }
      outfile << endl;
    }
  } else {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) L_lambda.data(), 
		  L_lambda.size()*sizeof(double));
    if (nebular != NULL)
      outfile.write((char *) L_lambda_neb.data(), 
		    L_lambda_neb.size()*sizeof(double));
    if (extinct != NULL) {
      outfile.write((char *) L_lambda_ext.data(), 
		    L_lambda_ext.size()*sizeof(double));
      if (nebular != NULL)
	outfile.write((char *) L_lambda_neb_ext.data(), 
		      L_lambda_neb_ext.size()*sizeof(double));
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Output spectrum in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_spectrum(fitsfile *out_fits, unsigned long trial) {

  // Make sure information is current
  if (!spec_set) set_spectrum();

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, L_lambda.size(), 
		 L_lambda.data(), &fits_status);
  unsigned int colnum = 5;
  if (nebular != NULL) {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		   L_lambda_neb.size(), L_lambda_neb.data(), 
		   &fits_status);
    colnum++;
  }
  if (extinct != NULL) {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		   L_lambda_ext.size(), L_lambda_ext.data(), 
		   &fits_status);
    colnum++;
    if (nebular != NULL) {
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		     L_lambda_neb_ext.size(), L_lambda_neb_ext.data(), 
		     &fits_status);
      colnum++;
    }
  }
  //Output rectified spectrum if it is present
  if (specsyn->get_rectify())
  {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		     recspec.size(), recspec.data(), 
		     &fits_status);
    colnum++;
  }
}
#endif

////////////////////////////////////////////////////////////////////////
// Output equivalent in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_ew(fitsfile *out_fits, unsigned long trial) 
{

  // Make sure information is current
  if (!ew_set) set_ew();
  
  
  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);


  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		             &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		             &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		             &fits_status);
  unsigned int colnum = 4;
  
  
  for (unsigned int i=0; i<ew.size(); i++) 
  {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		               &(ew[i]), &fits_status);
    colnum++;
  }
 
}
#endif

////////////////////////////////////////////////////////////////////////
// Output photometry
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_photometry(ofstream& outfile, const outputMode out_mode,
		 const unsigned long trial,
		 bool cluster_only) {

  // Make sure information is current
  if (!phot_set) set_photometry();

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific 
	    << setw(18) << right << id << "   "
	    << setw(18) << right << curTime;
    for (vector<double>::size_type i=0; i<phot.size(); i++)
      outfile << "   " << setw(18) << right << phot[i];
    if (nebular != NULL)
      for (vector<double>::size_type i=0; i<phot_neb.size(); i++)
	outfile << "   " << setw(18) << right << phot_neb[i];
    if (extinct != NULL) {
      for (vector<double>::size_type i=0; i<phot_ext.size(); i++) {
	if (!std::isnan(phot_ext[i]))
	  outfile << "   " << setw(18) << right << phot_ext[i];
	else
	  outfile << "   " << setw(18) << right << " ";
      }
      if (nebular != NULL) {
	for (vector<double>::size_type i=0; i<phot_neb_ext.size(); i++) {
	  if (!std::isnan(phot_neb_ext[i]))
	    outfile << "   " << setw(18) << right << phot_neb_ext[i];
	  else
	    outfile << "   " << setw(18) << right << " ";
	}
      }
    }
    outfile << endl;
  } else {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) phot.data(), 
		  phot.size()*sizeof(double));
    if (nebular != NULL)
      outfile.write((char *) phot_neb.data(), 
		    phot_neb.size()*sizeof(double));
    if (extinct != NULL) {
      outfile.write((char *) phot_ext.data(), 
		    phot_ext.size()*sizeof(double));
      if (nebular != NULL)
	outfile.write((char *) phot_neb_ext.data(), 
		      phot_neb_ext.size()*sizeof(double));
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Output photometry in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_photometry(fitsfile *out_fits, unsigned long trial) {

  // Make sure information is current
  if (!phot_set) set_photometry();

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  unsigned int colnum = 4;
  for (unsigned int i=0; i<phot.size(); i++) {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		   &(phot[i]), &fits_status);
    colnum++;
  }
  if (nebular != NULL) {
    for (unsigned int i=0; i<phot_neb.size(); i++) {
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_neb[i]), &fits_status);
      colnum++;
    }
  }
  if (extinct != NULL) {
    for (unsigned int i=0; i<phot_ext.size(); i++) {
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_ext[i]), &fits_status);
      colnum++;
    }
    if (nebular != NULL) {
      for (unsigned int i=0; i<phot_neb_ext.size(); i++) {
	fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		       &(phot_neb_ext[i]), &fits_status);
	colnum++;
      }
    }
  }
}
#endif

////////////////////////////////////////////////////////////////////////
// Output supernovae
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_sn(ofstream& outfile, const outputMode out_mode,
	 const unsigned long trial,
	 bool cluster_only) {

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific 
	    << setw(11) << right << id << "   "
	    << setw(11) << right << curTime << "   "
	    << setw(11) << right << tot_sn << "   "
	    << setw(11) << right << stoch_sn << endl;
  } else {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) &tot_sn, sizeof tot_sn);
    outfile.write((char *) &stoch_sn, sizeof stoch_sn);
  }
}


////////////////////////////////////////////////////////////////////////
// Output supernova count in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_sn(fitsfile *out_fits, unsigned long trial) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, 1, &tot_sn, 
		 &fits_status);
  unsigned long tmp = stoch_sn;
  fits_write_col(out_fits, TULONG, 5, nrows+1, 1, 1, &tmp, 
		 &fits_status);
}
#endif


#ifdef WINDS_ON
////////////////////////////////////////////////////////////////////////
// Output winds
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_winds(ofstream& outfile, const outputMode out_mode,
	    const unsigned long trial,
	    bool cluster_only) {

  // Make sure information is current
  if (!winds_set) set_winds();

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific 
	    << setw(15) << right << id << "   "
	    << setw(15) << right << curTime << "   "
	    << setw(15) << right << wind_mDot << "   "
	    << setw(15) << right << wind_pDot << "   "
	    << setw(15) << right << wind_lum << endl;
  } else {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) &wind_mDot, sizeof wind_mDot);
    outfile.write((char *) &wind_pDot, sizeof wind_pDot);
    outfile.write((char *) &wind_lum, sizeof wind_lum);
  }
}


////////////////////////////////////////////////////////////////////////
// Output winds in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_winds(fitsfile *out_fits, unsigned long trial) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, 1, &wind_mDot, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 5, nrows+1, 1, 1, &wind_pDot, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 6, nrows+1, 1, 1, &wind_lum, 
		 &fits_status);
}
#endif
#endif


////////////////////////////////////////////////////////////////////////
// Output yields
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_yield(ofstream& outfile, const outputMode out_mode,
	    const unsigned long trial, const bool cluster_only) {

  // Make sure information is current
  if (!yield_set) set_yield();

  // Write
  if (out_mode == ASCII) {
    const vector<const isotope_data *>& isodata = yields->get_isotopes();
    for (vector<double>::size_type i=0; i<all_yields.size(); i++) {
      outfile << setprecision(5) << scientific
	      << setw(11) << right << id << "   "
	      << setw(11) << right << curTime << "   "
	      << setw(11) << right << isodata[i]->symbol() << "   "
	      << setw(11) << right << isodata[i]->num() << "   "
	      << setw(11) << right << isodata[i]->wgt() << "   "
	      << setw(11) << right << all_yields[i] << endl;
    }
  } else if (out_mode == BINARY) {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) all_yields.data(),
		  all_yields.size()*sizeof(double));
  }
}

////////////////////////////////////////////////////////////////////////
// Output yields in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void 
slug_cluster::write_yield(fitsfile *out_fits, unsigned long trial){

  // Make sure information is current
  if (!yield_set) set_yield();

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write a new set of entries
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial,
     &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime,
     &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, all_yields.size(), 
    all_yields.data(), &fits_status);
}
#endif
