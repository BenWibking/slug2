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
#include <boost/bind.hpp>
#include "constants.H"
#include "slug_cluster.H"
#include "slug_galaxy.H"
#include "slug_parmParser.H"
#include "specsyn/slug_specsyn.H"
#include "tracks/slug_tracks.H"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <boost/bind.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Some trivial helper functions
////////////////////////////////////////////////////////////////////////
namespace galaxy {

  // Sorts stars by death time
  bool sort_death_time_decreasing(const slug_star star1, 
				  const slug_star star2) {
    return (star1.death_time > star2.death_time);
  }

  // Returns current mass from star data
  double curMass(const slug_stardata &data) {
    return exp(data.logM/constants::loge);
  }

  // This function returns the yield of a star of mass specified mass
  // if it is old enough to have died, or a vector of 0's otherwise
  vector<double> yield(const double &m, const double &t,
		       const slug_tracks *tracks,
		       const slug_yields *yields) {
    if (tracks->star_lifetime(m) > t) {
      vector<double> yld(yields->get_niso(), 0.0);
      return yld;
    } else {
      return yields->yield(m);
    }
  }

  // This function returns int_0^{t-t_life(m)} SFR(t') dt'
  double dnSN_dm(const double m, const double t,
		 const slug_PDF* sfh, const slug_tracks* tracks) {
    double t_life = tracks->star_lifetime(m);
    if (t > t_life) {
      return sfh->integral(0, t-t_life);
    } else {
      return 0.0;
    }
  }
}

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_galaxy::slug_galaxy(const slug_parmParser& pp, 
			 const slug_PDF* imf_,
			 const slug_PDF* cmf_, 
			 const slug_PDF* clf_, 
			 const slug_PDF* sfh_, 
			 const slug_tracks* tracks_, 
			 const slug_specsyn* specsyn_,
			 const slug_filter_set* filters_,
			 const slug_extinction* extinct_,
			 const slug_nebular* nebular_,
			 const slug_yields* yields_,
			 slug_ostreams& ostreams_) :
  ostreams(ostreams_),
  imf(imf_), 
  cmf(cmf_), 
  clf(clf_),
  sfh(sfh_),
  tracks(tracks_),
  specsyn(specsyn_),
  filters(filters_),
  extinct(extinct_),
  nebular(nebular_),
  yields(yields_),
  integ(tracks_, imf_, sfh_, ostreams_),
  v_integ(tracks_, imf_, sfh_, ostreams_)
 {

  // Initialize mass and time
  curTime = 0.0;
  mass = 0.0;
  targetMass = 0.0;
  aliveMass = nonStochAliveMass = fieldAliveMass = stellarMass
    = clusterAliveMass = clusterStellarMass = fieldRemnantMass = 0.0;
  clusterMass = 0.0;
  nonStochFieldMass = 0.0;
  field_tot_sn = 0.0;
  field_stoch_sn = 0;

  // Initialize yields
  if (yields) {
    stoch_field_yields.assign(yields->get_niso(), 0.0);
    all_yields.assign(yields->get_niso(), 0.0);
    v_integ.set_nvec(yields->get_niso());
  }

  // Get fc
  fc = pp.query<double>("clust_frac");

  // Initialize the cluster ID pointer
  cluster_id = 0;

  // Initialize status flags
  Lbol_set = spec_set = field_data_set = phot_set = yield_set = false;
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_galaxy::~slug_galaxy() {

  // Destroy cluster lists
  while (disrupted_clusters.size() > 0) {
    if (disrupted_clusters.back() != nullptr)
      delete disrupted_clusters.back();
    disrupted_clusters.pop_back();
  }
  while (clusters.size() > 0) {
    if (clusters.back() != nullptr)
      delete clusters.back();
    clusters.pop_back();
  }
}


////////////////////////////////////////////////////////////////////////
// Reset function -- sets galaxy back to initial state
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::reset(bool reset_cluster_id) {
  // N.B. By default we do NOT reset the cluster_id pointer here,
  // because if we're doing multiple trials, we want the cluster IDs
  // for different trials to be distinct.
  curTime = mass = targetMass = aliveMass = nonStochAliveMass
    = fieldAliveMass = clusterMass = clusterAliveMass 
    = nonStochFieldMass = stellarMass = clusterStellarMass 
    = fieldRemnantMass = 0.0;
  Lbol_set = spec_set = field_data_set = phot_set = yield_set = false;
  field_stars.resize(0);
  field_tot_sn = 0.0;
  field_stoch_sn = 0;
  while (disrupted_clusters.size() > 0) {
    if (disrupted_clusters.back() != nullptr)
      delete disrupted_clusters.back();
    disrupted_clusters.pop_back();
  }
  while (clusters.size() > 0) {
    if (clusters.back() != nullptr)
      delete clusters.back();
    clusters.pop_back();
  }
  L_lambda.resize(0);
  if (yields) {
    stoch_field_yields.assign(yields->get_niso(), 0.0);
    all_yields.assign(yields->get_niso(), 0.0);
  }
  if (reset_cluster_id) cluster_id = 0;
}

////////////////////////////////////////////////////////////////////////
// Advance routine
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::advance(double time) {

  // Make sure we're not trying to go back into the past
  assert(time >= curTime);

  // Compute mass of new stars to be created; this is equal to the
  // expected mass of stars in the next time interval, plus or minus
  // any deficit or over-production from the previous step. This
  // second part is added to ensure that, if we're doing stop nearest
  // or something like that, we get as close as possible at each
  // time.
  double new_mass = sfh->integral(curTime, time);
  double mass_to_draw = new_mass + targetMass - mass;
  targetMass += new_mass;

  // Skip star/cluster creation if new_mass == 0, even if mass_to_draw
  // != 0, to avoid problems with trying to draw from a region of SFH
  // where no clusters should form
  if (new_mass > 0) {

    // Create new clusters
    if (fc != 0) {

      // Get masses of new clusters
      vector<double> new_cluster_masses;
      cmf->drawPopulation(fc*mass_to_draw, new_cluster_masses);

      // Get birth times of new clusters
      vector<double> birth_times = sfh->draw(curTime, time,
					     new_cluster_masses.size());

      // Create clusters of chosen masses and birth times, and push them
      // onto the master cluster list
      for (unsigned int i=0; i<new_cluster_masses.size(); i++) {
	slug_cluster *new_cluster = 
	  new slug_cluster(cluster_id++, new_cluster_masses[i],
			   birth_times[i], imf, tracks, specsyn, filters,
			   extinct, nebular, yields, lines, ostreams, clf);
	clusters.push_back(new_cluster);
	mass += new_cluster->get_birth_mass();
	clusterMass += new_cluster->get_birth_mass();
	clusterAliveMass += new_cluster->get_alive_mass();
	clusterStellarMass += new_cluster->get_stellar_mass();
      }
    }

    // Create new field stars
    if (fc != 1) {

      // Get masses of new field stars
      vector<double> new_star_masses;
      imf->drawPopulation((1.0-fc)*mass_to_draw, new_star_masses);
      int count=0;
      double mtot=0.0;
      for (vector<double>::size_type i=0; i<new_star_masses.size(); i++)
	if (new_star_masses[i] > 5.0) { count++; mtot += new_star_masses[i]; }

      // Get extinctions to field stars
      vector<double> AV, AVneb;
      if (extinct != NULL) {
	AV = extinct->draw_AV(new_star_masses.size());
	AVneb.resize(new_star_masses.size());
	for (vector<double>::size_type i = 0; i<AV.size(); i++)
	  AVneb[i] = AV[i] * extinct->draw_neb_extinct_fac();
      }

      // Push stars onto field star list; in the process, set the birth
      // time and death time for each of them
      for (unsigned int i=0; i<new_star_masses.size(); i++) {
	slug_star new_star;
	new_star.mass = new_star_masses[i];
	new_star.birth_time = sfh->draw(curTime, time);
	new_star.death_time = new_star.birth_time 
	  + tracks->star_lifetime(new_star.mass);
	field_stars.push_back(new_star);
	if (extinct != NULL) {
	  field_star_AV.push_back(AV[i]);
	  field_star_AV_neb.push_back(AVneb[i]);
	}
	mass += new_star.mass;
      }

      // Sort field star list by death time, from largest to smallest
      sort(field_stars.begin(), field_stars.end(), 
	   galaxy::sort_death_time_decreasing);

      // Increase the non-stochastic field star mass and the total
      // mass by the mass of field stars that should have formed below
      // the stochastic limit 
      double fstoch = imf->mass_frac_restrict();
      mass += (1.0-fc)*(1.0-fstoch)*mass_to_draw;
      nonStochFieldMass += (1.0-fc)*(1.0-fstoch)*mass_to_draw;
    }
  }

  // Advance all clusters to current time; recompute the clusterAliveMass
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    clusterAliveMass -= (*it)->get_alive_mass();
    clusterMass -= (*it)->get_alive_mass();
    clusterStellarMass -= (*it)->get_stellar_mass();
    (*it)->advance(time);
    clusterAliveMass += (*it)->get_alive_mass();
    clusterMass += (*it)->get_alive_mass();
    clusterStellarMass += (*it)->get_stellar_mass();
  }
  for (it = disrupted_clusters.begin(); 
       it != disrupted_clusters.end(); it++) {
    clusterAliveMass -= (*it)->get_alive_mass();
    clusterStellarMass -= (*it)->get_stellar_mass();
    (*it)->advance(time);
    clusterAliveMass += (*it)->get_alive_mass();
    clusterStellarMass += (*it)->get_stellar_mass();
  }

  // See if any clusters were disrupted over the last time step, and,
  // if so, move them to the disrupted list
  it=clusters.begin();
  while (it != clusters.end()) {
    if ((*it)->disrupted()) {
      disrupted_clusters.push_back(*it);
      clusterMass -= (*it)->get_alive_mass();
      it = clusters.erase(it);
    } else {
      ++it;
    }
  }

  // Go through the field star list and remove any field stars that
  // have died; save them so that we can compute their yields
  dead_field_stars.resize(0);
  while (field_stars.size() > 0) {
    if (field_stars.back().death_time < time) {
      fieldRemnantMass += tracks->remnant_mass(field_stars.back().mass);
      if (yields) {
	if (yields->produces_sn(field_stars.back().mass)) {
	  field_stoch_sn++;
	}
      }
      dead_field_stars.push_back(field_stars.back());
      field_stars.pop_back();
      if (extinct != NULL) {
	field_star_AV.pop_back();
	field_star_AV_neb.pop_back();
      }
    } else {
      break;
    }
  }

  // Flag that computed quantities are now out of date
  Lbol_set = spec_set = field_data_set = phot_set = yield_set = false;

  // Store new time
  curTime = time;

  // Update the field star data
  set_field_data();

  // Recompute the alive mass of the field stars; assume that field
  // stars below the lowest mass in our tracks have zero mass lass
  fieldAliveMass = 0.0;
  for (vector<double>::size_type i=0; i<field_stars.size(); i++) { 
    if (field_stars[i].mass >= tracks->min_mass()) break;
    fieldAliveMass += field_stars[i].mass;
  }
  for (vector<slug_stardata>::size_type i=0; i<field_data.size(); i++) 
    fieldAliveMass += exp(field_data[i].logM / constants::loge);

  // Compute the alive mass for non-stochastic field stars; be sure to
  // include the contribution from the part of the IMF that's below
  // the minimum track mass
  nonStochAliveMass = 0.0;
  if (imf->get_xStochMin() > imf->get_xMin()) {
    if (imf->get_xMin() < tracks->min_mass()) {
      nonStochAliveMass += targetMass * 
	imf->mass_frac(imf->get_xMin(), 
		       min(tracks->min_mass(), imf->get_xStochMin()));
    }
    nonStochAliveMass += integ.integrate_sfh(time, galaxy::curMass);
  }

  // Recompute the remnant mass for non-stochastic field stars; note a
  // bit of fancy footwork here: the slug_tracks::remnant_mass
  // function is overloaded, so we need to static_cast to the version
  // of it we want before passing to boost::bind
  nonStochRemnantMass = 0.0;
  if (imf->get_xStochMin() > imf->get_xMin())
    nonStochRemnantMass = 
      integ.integrate_sfh_nt(time, 
			     boost::bind(static_cast<double (slug_tracks::*)
					 (const double, const double,
					  const double) const> 
					 (&slug_tracks::remnant_mass), 
					 tracks, _1, _2,
					 tracks::null_metallicity));

  // Recompute the alive mass and stellar mass
  aliveMass = nonStochAliveMass + clusterAliveMass + fieldAliveMass;
  stellarMass = aliveMass + nonStochRemnantMass + fieldRemnantMass;

  // Update the total number of supernovae in the field; this is the
  // sum of the stochastic field supernovae plus the non-stochastic
  // contribution, which is obtained by integrating over the field
  // star population. The number of supernovae that have occurred
  // is given by
  //
  // N_SN = 1 / <m> int (dn/dm) Theta_SN(m)
  //     int_0^{t-t_life(m)} SFR(t') dt' dm,
  //
  // where t_life(m) is the lifetime of a star of mass m and
  // Theta_SN(m) is 0 if a particular initial mass does not end its
  // life as a supernova, and 1 if it does. Note that, because we can
  // evaluate the integral over the SFR analytically, we evaluate this
  // integral using the functionality for integration over the IMF
  // rather than that for integration over IMF and star formation
  // history.
  if (yields) {
    field_tot_sn = field_stoch_sn;
    if (imf->get_xStochMin() > imf->get_xMin()) {
      vector<double> sn_mass_range = yields->sn_mass_range();
      for (vector<double>::size_type i=0; i<=sn_mass_range.size(); i+=2) {
	double mlo = min(sn_mass_range[i], imf->get_xStochMin());
	double mhi = min(sn_mass_range[i+1], imf->get_xMax());
	if (mhi > mlo) {
	  field_tot_sn +=
	    integ.integrate_nt_lim(1.0, time, mlo, mhi,
				   boost::bind(galaxy::dnSN_dm, _1, _2,
					       sfh, tracks));
	}
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Return supernova counts
////////////////////////////////////////////////////////////////////////
double
slug_galaxy::get_sn() const {
  double sn = field_tot_sn;
  list<slug_cluster *>::const_iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    sn += (*it)->get_sn();
  }
  for (it = disrupted_clusters.begin();
       it != disrupted_clusters.end(); it++) {
    sn += (*it)->get_sn();
  }
  return sn;
}

int
slug_galaxy::get_stoch_sn() const {
  int sn = field_stoch_sn;
  list<slug_cluster *>::const_iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    sn += (*it)->get_stoch_sn();
  }
  for (it = disrupted_clusters.begin();
       it != disrupted_clusters.end(); it++) {
    sn += (*it)->get_stoch_sn();
  }
  return sn;
}

double
slug_galaxy::get_non_stoch_sn() const {
  return get_sn() - get_stoch_sn();
}

////////////////////////////////////////////////////////////////////////
// Get stellar data on all field stars
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_field_data() {

  // Do nothing if data is current
  if (field_data_set) return;

  // Initialize
  field_data.resize(0);

  // Get field star data
  for (vector<int>::size_type i=0; i<field_stars.size(); i++) {
    if (field_stars[i].mass >= tracks->min_mass() &&
	field_stars[i].mass <= tracks->max_mass()) {
      field_data.
      push_back(tracks->get_star(field_stars[i].mass,
      				   curTime - field_stars[i].birth_time));
#ifndef NDBEUG
      if (!isfinite(field_data.back().logTeff) ||
	  !isfinite(field_data.back().logL) ||
	  !isfinite(field_data.back().logg) ||
	  !isfinite(field_data.back().logM)) {
	cout << setprecision(20)
	     << "bad star: m = "
	     << field_stars[i].mass
	     << ", age = " << curTime - field_stars[i].birth_time
	     << endl;
      }
#endif
    }
  }

  // Set status flag
  field_data_set = true;
}


////////////////////////////////////////////////////////////////////////
// Compute Lbol
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_Lbol() {

  // Do nothing if already set
  if (Lbol_set) return;

  // Initialize
  Lbol = 0.0;

  // First loop over non-disrupted clusters
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++)
    Lbol += (*it)->get_Lbol();

  // Now do disrupted clusters
  for (it = disrupted_clusters.begin(); 
       it != disrupted_clusters.end(); 
       it++)
    Lbol += (*it)->get_Lbol();

  // Now do stochastic field stars
  if (!field_data_set) set_field_data();
  for (unsigned int i=0; i<field_data.size(); i++) {
    Lbol += pow(10.0, field_data[0].logL);
  }

  // Now do non-stochastic field stars
  if (imf->has_stoch_lim())
    Lbol += specsyn->get_Lbol_cts_sfh(curTime);

  // Set flag
  Lbol_set = true;
}


////////////////////////////////////////////////////////////////////////
// Compute spectrum, getting Lbol in the process
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_spectrum(const bool del_cluster) {

  // Do nothing if already set
  if (spec_set) return;

  // Initialize
  vector<double>::size_type nl = specsyn->n_lambda();
  L_lambda.assign(nl, 0.0);
  Lbol = 0.0;
  if (nebular != NULL) L_lambda_neb.assign(nebular->n_lambda(), 0.0);
  vector<double>::size_type nl_ext = 0;
  if (extinct != NULL) {
    nl_ext = extinct->n_lambda();
    L_lambda_ext.assign(nl_ext, 0.0);
    Lbol_ext = 0.0;
    if (nebular != NULL) 
      L_lambda_neb_ext.assign(extinct->n_lambda_neb(), 0.0);
  }

  // Loop over non-disrupted clusters; for each one, get spectrum and
  // bolometric luminosity and add both to global sum
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    const vector<double>& spec = (*it)->get_spectrum();
    for (vector<double>::size_type i=0; i<nl; i++) 
      L_lambda[i] += spec[i];
    Lbol += (*it)->get_Lbol();
    if (nebular != NULL) {
      const vector<double>& spec_neb = (*it)->get_spectrum_neb();
      for (vector<double>::size_type i=0; i<nebular->n_lambda(); i++) 
	L_lambda_neb[i] += spec_neb[i];
    }
    if (extinct != NULL) {
      const vector<double>& spec_ext = (*it)->get_spectrum_extinct();
      for (vector<double>::size_type i=0; i<nl_ext; i++) 
	L_lambda_ext[i] += spec_ext[i];
      Lbol_ext += (*it)->get_Lbol_extinct();
      if (nebular != NULL) {
	const vector<double>& spec_neb_ext 
	  = (*it)->get_spectrum_neb_extinct();
	for (vector<double>::size_type i=0; i<extinct->n_lambda_neb(); 
	     i++) 
	  L_lambda_neb_ext[i] += spec_neb_ext[i];
      }
    }
    if (del_cluster) {
      delete (*it);
      *it = nullptr;
    }
  }

  // Now do exactly the same thing for disrupted clusters
  for (it = disrupted_clusters.begin(); it != disrupted_clusters.end(); 
       it++) {
    const vector<double>& spec = (*it)->get_spectrum();
    for (vector<double>::size_type i=0; i<nl; i++) 
      L_lambda[i] += spec[i];
    Lbol += (*it)->get_Lbol();
    if (nebular != NULL) {
      const vector<double>& spec_neb = (*it)->get_spectrum_neb();
      for (vector<double>::size_type i=0; i<nebular->n_lambda(); i++) 
	L_lambda_neb[i] += spec_neb[i];
    }
    if (extinct != NULL) {
      const vector<double>& spec_ext = (*it)->get_spectrum_extinct();
      for (vector<double>::size_type i=0; i<nl_ext; i++) 
	L_lambda_ext[i] += spec_ext[i];
      Lbol_ext += (*it)->get_Lbol_extinct();
      if (nebular != NULL) {
	const vector<double>& spec_neb_ext 
	  = (*it)->get_spectrum_neb_extinct();
	for (vector<double>::size_type i=0; i<extinct->n_lambda_neb(); i++) 
	  L_lambda_neb_ext[i] += spec_neb_ext[i];
      }
    }
    if (del_cluster) {
      delete (*it);
      *it = nullptr;
    }
  }

  // Now do stochastic field stars
  if (!field_data_set) set_field_data();
  for (vector<slug_stardata>::size_type i=0; i<field_data.size(); i++) {
    const vector<double> &spec = specsyn->get_spectrum(field_data[i]);
    for (vector<double>::size_type j=0; j<nl; j++) 
      L_lambda[j] += spec[j];
    Lbol += pow(10.0, field_data[i].logL);
    vector<double> spec_neb, spec_neb_only;
    if (nebular != NULL) {
      spec_neb_only =
	nebular->get_neb_spec(spec, 
			      curTime-field_stars[i].birth_time);
      spec_neb = nebular->add_stellar_nebular_spec(spec, spec_neb_only);
      for (vector<double>::size_type j=0; j<nebular->n_lambda(); j++) 
	L_lambda_neb[j] += spec_neb[j];
    }
    if (extinct != NULL) {
      vector<double> spec_ext = 
	extinct->spec_extinct(field_star_AV[i], spec);
      for (vector<double>::size_type j=0; j<nl_ext; j++) 
	L_lambda_ext[j] += spec_ext[j];
      Lbol_ext += int_tabulated::
	integrate(extinct->lambda(), spec_ext) / constants::Lsun;
      if (nebular != NULL) {
	// Procedure depends on if nebular and stellar extinctions
	// differ
	vector<double> spec_neb_ext;
	if (extinct->excess_neb_extinct()) {
	  // We have differential nebular extinction; apply nebular
	  // extinction to nebular portion
	  vector<double> spec_neb_only_ext =
	    extinct->spec_extinct_neb(field_star_AV_neb[i], spec_neb_only);
	  // Add extincted stellar and nebular spectra
	  spec_neb_ext = nebular->
	    add_stellar_nebular_spec(spec_ext, spec_neb_only_ext,
				     extinct->off(),
				     extinct->off_neb());
	} else {
	  // No differential stellar and nebular spectra, so just
	  // exintct the combined spectrum we already have in hand
	  spec_neb_ext =
	    extinct->spec_extinct_neb(field_star_AV[i], spec_neb);
	}
	// Add to running total
	for (vector<double>::size_type j=0; j<extinct->n_lambda_neb(); 
	     j++) 
	  L_lambda_neb_ext[j] += spec_neb_ext[j];
      }
    }
  }

  // Finally do non-stochastic field stars; note that non-stochastic
  // field stars get extincted by the expectation value of the AV
  // distribution, and their nebular emission gets extincted by the
  // expectation value of the nebular A_V distribution.
  if (imf->has_stoch_lim()) {
    double Lbol_tmp;
    vector<double> spec;
    specsyn->get_spectrum_cts_sfh(curTime, spec, Lbol_tmp);
    for (vector<double>::size_type i=0; i<nl; i++) 
      L_lambda[i] += spec[i];
    Lbol += Lbol_tmp;
    vector<double> spec_neb, spec_neb_only;
    if (nebular != NULL) {
      spec_neb_only = nebular->get_neb_spec(spec, -1.0);
      spec_neb = nebular->add_stellar_nebular_spec(spec, spec_neb_only);
      for (vector<double>::size_type i=0; i<nebular->n_lambda(); i++) 
	L_lambda_neb[i] += spec_neb[i];
    }
    if (extinct != NULL) {
      vector<double> spec_ext =
	extinct->spec_extinct(extinct->AV_expect(), spec);
      for (vector<double>::size_type i=0; i<nl_ext; i++) 
	L_lambda_ext[i] += spec_ext[i];
      Lbol_ext += int_tabulated::
	integrate(extinct->lambda(), spec_ext) / constants::Lsun;
      if (nebular != NULL) {
	vector<double> spec_neb_ext;
	// Procedure depends on whether we have differential extinction
	if (extinct->excess_neb_extinct()) {
	  // Apply nebular extinction to nebular portion
	  vector<double> spec_neb_only_ext =
	    extinct->spec_extinct_neb(extinct->AV_neb_expect(),
				      spec_neb_only);
	  // Add extincted stellar and nebular spectra
	  spec_neb_ext = nebular->
	    add_stellar_nebular_spec(spec_ext,
				     spec_neb_only_ext,
				     extinct->off(),
				     extinct->off_neb());
	} else {
	  // Just apply a single extinction to the summed stellar plus
	  // nebular spectrum
	  spec_neb_ext = 
	    extinct->spec_extinct_neb(extinct->AV_expect(), spec_neb);
	}
	// Add to running total
	for (vector<double>::size_type i=0; i<extinct->n_lambda_neb(); 
	     i++) 
	  L_lambda_neb_ext[i] += spec_neb_ext[i];
      }
    }
  }

  // Set flags
  Lbol_set = spec_set = true;
}


////////////////////////////////////////////////////////////////////////
// Compute photometry
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_photometry(const bool del_cluster) {

  // Do nothing if already set
  if (phot_set) return;

  // Compute the spectrum
  set_spectrum(del_cluster);

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

  // Repeat for extincted values
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
// Compute yields
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_yield(const bool del_cluster) {

  // If yield is already set, do nothing
  if (yield_set) return;

  // Get yields from all clusters, disrupted and non-disrupted
  vector<double> cluster_yields(all_yields.size(), 0.0);
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    const vector<double>& ylds = (*it)->get_yield();
    for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      cluster_yields[i] += ylds[i];
  }
  for (it = disrupted_clusters.begin();
       it != disrupted_clusters.end(); it++) {
    const vector<double>& ylds = (*it)->get_yield();
    for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      cluster_yields[i] += ylds[i];
  }

  // If we don't have field stars, we're done; just copy the cluster
  // yields to the list of all yields
  if (fc == 1.0) {
    for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      all_yields[i] = cluster_yields[i];
    yield_set = true;
    return;
  }

  // If we're here, we have field stars
  
  // Make unstable isotopes produced by field stars decay; the cluster
  // update will have already taken care of this for the clusters
  if (!yields->no_decay) {
    const vector<const isotope_data *>& isotopes = yields->get_isotopes();
    for (vector<double>::size_type i=0; i<stoch_field_yields.size(); i++) {
      if (!isotopes[i]->stable()) {
	stoch_field_yields[i] *=
	  exp(-(curTime-last_yield_time)/isotopes[i]->ltime());
      }
    }
    last_yield_time = curTime;
  }
  
  // Add yields from stochastic field stars that died this time step
  if (dead_field_stars.size() > 0) {
    vector<double> yld_masses(dead_field_stars.size());
    for (vector<double>::size_type i=0; i<dead_field_stars.size(); i++)
      yld_masses[i] = dead_field_stars[i].mass;
    vector<double> star_yields;
    if (!yields->no_decay) {
      vector<double> decay_time(dead_field_stars.size());
      for (vector<double>::size_type i=0; i<dead_field_stars.size(); i++)
	decay_time[i] = curTime - dead_field_stars[i].death_time;
      star_yields = yields->yield(yld_masses, decay_time);
    } else {
      star_yields = yields->yield(yld_masses);
    }
    for (vector<double>::size_type i=0; i<star_yields.size(); i++)
      stoch_field_yields[i] += star_yields[i];
  }

  // Do we have non-stochastic field stars? If not, we're done; just
  // add the cluster yields and the stochastic field star yields, and
  // then return
  if (!imf->has_stoch_lim()) {
    for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      all_yields[i] = cluster_yields[i] + stoch_field_yields[i];
    yield_set = true;
    return;
  }

  // If we're here, we have a contribution from non-stochastic field
  // stars, which we must compute as a double-integral over the SFH
  // and IMF. This doesn't quite treat radioactive decay correctly,
  // because we should apply decay to the radioactive isotopes
  // produced during the time step, which we're not doing now. To be
  // fixed later! For now, the problem can be minimized by either not
  // using the radioactive isotopes in non-stochastic mode, or by
  // setting the time step to something much smaller than the lifetime
  // of the isotopes in question.
  vector<double> non_stoch_field_yields =
    v_integ.integrate_sfh_nt(curTime,
			     boost::bind(galaxy::yield, _1, _2,
					 tracks, yields));
  
  // Sum all yields
  for (vector<double>::size_type i=0; i<all_yields.size(); i++)
      all_yields[i] = cluster_yields[i] + stoch_field_yields[i]
	+ non_stoch_field_yields[i];
  yield_set = true;  
}


////////////////////////////////////////////////////////////////////////
// Output integrated properties
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_prop(ofstream& int_prop_file, 
				   const outputMode out_mode, 
				   const unsigned long trial,
				   const std::vector<double>& imfvp) {

  if (out_mode == ASCII) {
  
    //ASCII Output
    int_prop_file << setprecision(5) << scientific 
		  << setw(11) << right << curTime << "   "
		  << setw(11) << right << targetMass << "   "
		  << setw(11) << right << mass << "   "
		  << setw(11) << right << aliveMass << "   "
		  << setw(11) << right << stellarMass << "   "
		  << setw(11) << right << clusterMass << "   "
		  << setw(11) << right << clusters.size() << "   "
		  << setw(11) << right << disrupted_clusters.size() << "   "
		  << setw(11) << right << field_stars.size();
	  
    // Output any variable parameters  
    if (imfvp.size()>0) {
      // Loop over the variable parameters and output
      for (vector<double>::size_type p = 0; p<imfvp.size();p++) {
        int_prop_file << "   " << setw(11) << right << imfvp[p];
      }
    }

    // Close
    int_prop_file << endl;

  } else {
  
    // Binary Output
    int_prop_file.write((char *) &trial, sizeof trial);
    int_prop_file.write((char *) &curTime, sizeof curTime);
    int_prop_file.write((char *) &targetMass, sizeof targetMass);
    int_prop_file.write((char *) &mass, sizeof mass);
    int_prop_file.write((char *) &aliveMass, sizeof aliveMass);
    int_prop_file.write((char *) &stellarMass, sizeof stellarMass);
    int_prop_file.write((char *) &clusterMass, sizeof clusterMass);
    vector<slug_cluster *>::size_type n = clusters.size();
    int_prop_file.write((char *) &n, sizeof n);
    n = disrupted_clusters.size();
    int_prop_file.write((char *) &n, sizeof n);
    n = field_stars.size();
    int_prop_file.write((char *) &n, sizeof n);
    
    // Write out variable parameter values
    if (imfvp.size()>0) {
      // Loop over the variable parameters
      for (vector<double>::size_type p = 0; p<imfvp.size();p++) {
        int_prop_file.write((char *) &imfvp[p], sizeof imfvp[p]);
      }  
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Output integrated properties to fits file
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_prop(fitsfile* int_prop_fits, 
				   unsigned long trial,
				   const std::vector<double>& imfvp) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_prop_fits, &nrows, &fits_status);

  // Write a new entry
  fits_write_col(int_prop_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 3, nrows+1, 1, 1, &targetMass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 4, nrows+1, 1, 1, &mass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 5, nrows+1, 1, 1, &aliveMass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 6, nrows+1, 1, 1, &stellarMass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 7, nrows+1, 1, 1, 
		 &clusterMass, &fits_status);
  vector<slug_cluster *>::size_type n = clusters.size();
  fits_write_col(int_prop_fits, TULONG, 8, nrows+1, 1, 1, &n,	 
		 &fits_status);
  n = disrupted_clusters.size();
  fits_write_col(int_prop_fits, TULONG, 9, nrows+1, 1, 1, &n,	 
		 &fits_status);
  n = field_stars.size();
  fits_write_col(int_prop_fits, TULONG, 10, nrows+1, 1, 1, &n,	 
		 &fits_status);
		 
  if (imfvp.size()>0) {
    // Loop over the variable parameters
    int colnum=10;  //Current column number
    for (vector<double>::size_type p = 0; p<imfvp.size(); p++) {
      colnum++;
      double vp_p=imfvp[p];
      fits_write_col(int_prop_fits, TDOUBLE, colnum, nrows+1, 1, 1, &vp_p,
		  &fits_status);
    }
  }		 
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster properties
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_prop(ofstream& cluster_prop_file, 
				const outputMode out_mode,
				const unsigned long trial,
				const std::vector<double>& imfvp) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  
  if (out_mode == BINARY) {
    cluster_prop_file.write((char *) &trial, sizeof trial);
    cluster_prop_file.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    cluster_prop_file.write((char *) &n, sizeof n);
  }

  // Now write out each cluster
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_prop(cluster_prop_file, out_mode, trial,false,imfvp);
}

////////////////////////////////////////////////////////////////////////
// Output cluster properties in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_prop(fitsfile* cluster_prop_fits, 
				unsigned long trial,
				const std::vector<double>& imfvp) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_prop(cluster_prop_fits, trial,imfvp);
}
#endif


////////////////////////////////////////////////////////////////////////
// Output integrated spectra
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_spec(ofstream& int_spec_file, 
				   const outputMode out_mode,
				   const unsigned long trial,
				   const bool del_cluster) {

  // Make sure spectrum information is current. If not, compute it.
  if (!spec_set) set_spectrum(del_cluster);

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
    for (vector<double>::size_type i=0; i<lambda.size(); i++) {
      int_spec_file << setprecision(5) << scientific 
		    << setw(11) << right << curTime << "   "
		    << setw(11) << right << lambda[i] << "   "
		    << setw(11) << right << L_lambda_star[i];
      if (nebular != NULL)
	int_spec_file << "   "
		      << setw(11) << right << L_lambda_neb[i];
      if (extinct != NULL) {
	int j;
	if (nebular == NULL) j = i - extinct->off();
	else j = i - extinct->off_neb();
	if ((j >= 0) && ((unsigned int) j < L_lambda_star_ext.size())) {
	  int_spec_file << "   "
			<< setw(11) << right << L_lambda_star_ext[j];
	  if (nebular != NULL)
	    int_spec_file << "   "
			  << setw(11) << right << L_lambda_neb_ext[j];
	}
      }
      int_spec_file << endl;
    }
  } else {
    int_spec_file.write((char *) &trial, sizeof trial);
    int_spec_file.write((char *) &curTime, sizeof curTime);
    int_spec_file.write((char *) L_lambda.data(), 
			sizeof(double)*L_lambda.size());
    if (nebular != NULL)
      int_spec_file.write((char *) L_lambda_neb.data(),
			  sizeof(double)*L_lambda_neb.size());
    if (extinct != NULL) {
      int_spec_file.write((char *) L_lambda_ext.data(),
			  sizeof(double)*L_lambda_ext.size());
      if (nebular != NULL)
	int_spec_file.write((char *) L_lambda_neb_ext.data(),
			    sizeof(double)*L_lambda_neb_ext.size());
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Output integrated spectra in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_spec(fitsfile* int_spec_fits, 
				   unsigned long trial,
				   const bool del_cluster) {

  // Make sure spectrum information is current. If not, compute it.
  if (!spec_set) set_spectrum(del_cluster);

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_spec_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(int_spec_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_spec_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(int_spec_fits, TDOUBLE, 3, nrows+1, 1, 
		 L_lambda.size(), L_lambda.data(), &fits_status);
  int colnum = 4;
  if (nebular != NULL) {
    fits_write_col(int_spec_fits, TDOUBLE, colnum, nrows+1, 1, 
		   L_lambda_neb.size(), L_lambda_neb.data(),
		   &fits_status);
    colnum++;
  }
  if (extinct != NULL) {
    fits_write_col(int_spec_fits, TDOUBLE, colnum, nrows+1, 1, 
		   L_lambda_ext.size(), L_lambda_ext.data(),
		   &fits_status);
    colnum++;
    if (nebular != NULL) {
      fits_write_col(int_spec_fits, TDOUBLE, colnum, nrows+1, 1, 
		     L_lambda_neb_ext.size(), L_lambda_neb_ext.data(),
		     &fits_status);
      colnum++;
    }
  }
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster spectra
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_spec(ofstream& cluster_spec_file, 
				const outputMode out_mode,
				const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  if (out_mode == BINARY) {
    cluster_spec_file.write((char *) &trial, sizeof trial);
    cluster_spec_file.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    cluster_spec_file.write((char *) &n, sizeof n);
  }

  // Now have each cluster write
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_spectrum(cluster_spec_file, out_mode, trial);
}


////////////////////////////////////////////////////////////////////////
// Output cluster spectra in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_spec(fitsfile* cluster_spec_fits, 
				unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_spectrum(cluster_spec_fits, trial);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output integrated photometry
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_phot(ofstream& outfile, 
				   const outputMode out_mode,
				   const unsigned long trial,
				   const bool del_cluster) {

  // Make sure photometric information is current. If not, compute it.
  if (!phot_set) set_photometry(del_cluster);

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific 
	    << setw(18) << right << curTime;
    for (vector<double>::size_type i=0; i<phot.size(); i++)
      outfile << "   " << setw(18) << right << phot[i];
    if (nebular != NULL) {
      for (vector<double>::size_type i=0; i<phot_neb.size(); i++) {
	if (!std::isnan(phot_neb[i]))
	  outfile << "   " << setw(18) << right << phot_neb[i];
	else
	  outfile << "   " << setw(18) << right << " ";
      }
    }
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
    outfile.write((char *) &trial, sizeof trial);
    outfile.write((char *) &curTime, sizeof curTime);
    outfile.write((char *) phot.data(), 
			sizeof(double)*phot.size());
    if (nebular != NULL)
      outfile.write((char *) phot_neb.data(), 
		    sizeof(double)*phot_neb.size());
    if (extinct != NULL) {
      outfile.write((char *) phot_ext.data(), 
		    sizeof(double)*phot_ext.size());
      if (nebular != NULL)
	outfile.write((char *) phot_neb_ext.data(), 
		      sizeof(double)*phot_neb_ext.size());
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Output integrated photometry in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_phot(fitsfile* int_phot_fits, 
				   unsigned long trial,
				   const bool del_cluster) {

  // Make sure photometric information is current. If not, compute it.
  if (!phot_set) set_photometry(del_cluster);

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_phot_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(int_phot_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_phot_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  unsigned int colnum = 3;
  for (unsigned int i=0; i<phot.size(); i++) {
    fits_write_col(int_phot_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		   &(phot[i]), &fits_status);
    colnum++;
  }
  if (nebular != NULL) {
    for (unsigned int i=0; i<phot_neb.size(); i++) {
      fits_write_col(int_phot_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_neb[i]), &fits_status);
      colnum++;
    }
  }
  if (extinct != NULL) {
    for (unsigned int i=0; i<phot_ext.size(); i++) {
      fits_write_col(int_phot_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_ext[i]), &fits_status);
      colnum++;
    }
    if (nebular != NULL) {
      for (unsigned int i=0; i<phot_neb_ext.size(); i++) {
	fits_write_col(int_phot_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		       &(phot_neb_ext[i]), &fits_status);
	colnum++;
      }
    }
  }
}
#endif


////////////////////////////////////////////////////////////////////////
// Output integrated supernova count
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_sn(ofstream& int_sn_file, 
				 const outputMode out_mode, 
				 const unsigned long trial) {

  if (out_mode == ASCII) {
  
    //ASCII Output
    int_sn_file << setprecision(5) << scientific 
		<< setw(11) << right << curTime << "   "
		<< setw(11) << right << get_sn() << "   "
		<< setw(11) << right << get_stoch_sn() << endl;

  } else {
  
    // Binary Output
    int_sn_file.write((char *) &trial, sizeof trial);
    int_sn_file.write((char *) &curTime, sizeof curTime);
    double tot_sn = get_sn();
    int_sn_file.write((char *) &tot_sn, sizeof tot_sn);
    unsigned long stoch_sn = get_stoch_sn();
    int_sn_file.write((char *) &stoch_sn, sizeof stoch_sn);
    
  }
}

////////////////////////////////////////////////////////////////////////
// Output integrated supernova count to fits file
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_sn(fitsfile* int_sn_fits, 
				 unsigned long trial) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_sn_fits, &nrows, &fits_status);

  // Write a new entry
  fits_write_col(int_sn_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_sn_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  double totSN = get_sn();
  fits_write_col(int_sn_fits, TDOUBLE, 3, nrows+1, 1, 1, &totSN, 
		 &fits_status);
  unsigned long stochSN = get_stoch_sn();
  fits_write_col(int_sn_fits, TULONG, 4, nrows+1, 1, 1, &stochSN, 
		 &fits_status);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster supernova counts
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_sn(ofstream& cluster_sn_file, 
			      const outputMode out_mode,
			      const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data  
  if (out_mode == BINARY) {
    cluster_sn_file.write((char *) &trial, sizeof trial);
    cluster_sn_file.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    cluster_sn_file.write((char *) &n, sizeof n);
  }

  // Now write out each cluster
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_sn(cluster_sn_file, out_mode, trial, false);
}

////////////////////////////////////////////////////////////////////////
// Output cluster supernova counts in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_sn(fitsfile* cluster_sn_fits, 
			      unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_sn(cluster_sn_fits, trial);
}
#endif


////////////////////////////////////////////////////////////////////////
// Output integrated yield
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_yield(std::ofstream& outfile,
				    const outputMode out_mode,
				    const unsigned long trial,
				    const bool del_cluster) {
  
  // Make sure yield information is current. If not, compute it.
  if (!yield_set) set_yield(del_cluster);

  // Write
  if (out_mode == ASCII) {
    const vector<const isotope_data *>& isodata = yields->get_isotopes();
    for (vector<double>::size_type i=0; i<all_yields.size(); i++) {
      outfile << setprecision(5) << scientific
	      << setw(11) << right << curTime << "   "
	      << setw(11) << right << isodata[i]->symbol() << "   "
	      << setw(11) << right << isodata[i]->num() << "   "
	      << setw(11) << right << isodata[i]->wgt() << "   "
	      << setw(11) << right << all_yields[i] << endl;
    }
  } else if (out_mode == BINARY) {
    outfile.write((char *) &trial, sizeof trial);
    outfile.write((char *) &curTime, sizeof curTime);
    outfile.write((char *) all_yields.data(),
		  all_yields.size()*sizeof(double));
  }
}


////////////////////////////////////////////////////////////////////////
// Output integrated yield in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_yield(fitsfile* int_yield_fits, 
				    unsigned long trial,
				    const bool del_cluster) {

  // Make sure photometric information is current. If not, compute it.
  if (!yield_set) set_yield(del_cluster);

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_yield_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(int_yield_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_yield_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(int_yield_fits, TDOUBLE, 3, nrows+1, 1, all_yields.size(), 
		 all_yields.data(), &fits_status);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster photometry
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_phot(ofstream& outfile, 
				const outputMode out_mode,
				const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  if (out_mode == BINARY) {
    outfile.write((char *) &trial, sizeof trial);
    outfile.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    outfile.write((char *) &n, sizeof n);
  }

  // Now have each cluster write
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_photometry(outfile, out_mode, trial);
}

////////////////////////////////////////////////////////////////////////
// Output cluster photometry in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_phot(fitsfile* cluster_phot_fits, 
				unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_photometry(cluster_phot_fits, trial);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster yields
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_yield(ofstream& outfile, 
				const outputMode out_mode,
				const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  if (out_mode == BINARY) {
    outfile.write((char *) &trial, sizeof trial);
    outfile.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    outfile.write((char *) &n, sizeof n);
  }

  // Now have each cluster write
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_yield(outfile, out_mode, trial);
}

////////////////////////////////////////////////////////////////////////
// Output cluster yields in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_yield(fitsfile* cluster_yield_fits, 
				unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_yield(cluster_yield_fits, trial);
}
#endif

