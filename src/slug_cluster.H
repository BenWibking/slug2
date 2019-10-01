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

////////////////////////////////////////////////////////////////////////
// class slug_cluster
//
// This class represents a single star cluster. It knows its age, and
// it knows the masses of all stars it contains. It knows how to draw
// stars from an IMF, how to evolve the star list in time, and how to
// write its proprties to a file.
////////////////////////////////////////////////////////////////////////
#ifndef _slug_cluster_H_
#define _slug_cluster_H_

#include "slug.H"
#include "slug_IO.H"
#include "slug_extinction.H"
#include "slug_galaxy.H"
#include "slug_nebular.H"
#include "filters/slug_filter_set.H"
#include "lines/slug_line_list.H"
#include "pdfs/slug_PDF.H"
#include "specsyn/slug_specsyn.H"
#include "tracks/slug_tracks.H"
#include "utils/slug_imf_integrator.H"
#include "yields/slug_yields.H"
#include <iostream>
#include <fstream>
#ifdef ENABLE_FITS
extern "C" {
#   include "fitsio.h"
}
#endif

// slug_cluster_buffer is an alias for a general buffer, which holds
// data used to reconstruct a slug_cluster object. Note that this has
// to be declared as a general buffer instead of a struct (which would
// certainly be prettier) because we do not know in advance how many
// stars will be in the cluster.
typedef void slug_cluster_buffer;

class slug_cluster {

public:
  // Construct a cluster with specified target mass, birth time, IMF
  slug_cluster(const unsigned long id_, const double mass_, 
	       const double time, const slug_PDF *imf_, 
	       const slug_tracks *tracks_, 
	       const slug_specsyn *specsyn_,
	       const slug_filter_set *filters_,
	       const slug_extinction *extinct_,
	       const slug_nebular *nebular_,
	       const slug_yields *yields_,
	       const slug_line_list *lines_,
	       slug_ostreams &ostreams_,
	       const slug_PDF *clf_ = nullptr,
	       const bool stoch_contrib_only_ = false);

  // Constructor to reconstruct a slug_cluster object from a
  // serialized buffer
  slug_cluster(const slug_cluster_buffer *buf,
	       const slug_PDF *imf_, 
	       const slug_tracks *tracks_, 
	       const slug_specsyn *specsyn_,
	       const slug_filter_set *filters_,
	       const slug_extinction *extinct_,
	       const slug_nebular *nebular_,
	       const slug_yields *yields_,
	       const slug_line_list *lines_,
	       slug_ostreams &ostreams_,
	       const slug_PDF *clf_ = nullptr,
	       const bool stoch_contrib_only_ = false);

  // Copy constructor
  slug_cluster(const slug_cluster &obj, const unsigned long id_ = 0);

  // Destructor
  ~slug_cluster() { }

  // Routines for manipulating serialized buffers

  // Return size in bytes needed for a buffer
  size_t buffer_size() const;

  // Allocate and fill a buffer
  slug_cluster_buffer *make_buffer() const;

  // Fill a buffer; user must allocate it and guarantee that it
  // contains sufficient memory
  void pack_buffer(slug_cluster_buffer *buf) const;

  // Free a buffer allocated by make_buffer
  void free_buffer(slug_cluster_buffer *buffer) const;

  // Routine to return the id
  unsigned long get_id() const { return id; }

  // Set the cluster back to its initial state; this involves
  // re-drawing all the stars. If keep_id is false, the new ID will be
  // set to the old ID plus 1.
  void reset(bool keep_id = false);

  // Routine to advance to a specified time
  void advance(double time);

  // Routine to report most mass still-living star
  double get_most_massive_star() const { 
    if (stars.size() > 0) return stars[stars.size()-1];
    else return 0.0;
  }

  // Routine to report various masses and similar data
  double get_target_mass() const { return targetMass; }
  double get_birth_mass() const { return birthMass; }
  double get_alive_mass() const { return aliveMass; }
  double get_non_stoch_alive_mass() const { return nonStochAliveMass; }
  double get_stellar_mass() const { return stellarMass; }
  double get_age() const { return curTime - formationTime; }
  double get_stellar_death_mass() const { return stellarDeathMass; }
  double get_remnant_mass() const {
    return stochRemnantMass + nonStochRemnantMass; }
  double get_non_stoch_remnant_mass() const { return nonStochRemnantMass; }

  // Routines to report the number of type II SN events
  double get_sn() const { return tot_sn; }
  int get_stoch_sn() const { return stoch_sn; }
  double get_non_stoch_sn() const
  { if (nonStochBirthMass > 0) return tot_sn - stoch_sn; else return 0.0; }

#ifdef WINDS_ON
  // Routines to report on winds
  double get_wind_mDot() const { set_winds(); return wind_mDot; }
  double get_wind_pDot() const { set_winds(); return wind_pDot; }
  double get_wind_lum() const { set_winds(); return wind_lum; }
#endif

  // Routine to return the element yield
  const std::vector<double>& get_yields() 
  { set_yield(); return all_yields; }

  // Pass through to routine that returns the mapping between isotope
  // number and weight and index in the yield table
  std::vector<double>::size_type
  isotope_index(unsigned int num, unsigned int wgt) const {
    return yields->isotope_index(num, wgt);
  }

  // Routine to report if this cluster is disrupted
  bool disrupted() const { return is_disrupted; }

  // Routine to return lifetime of cluster against disruption
  double get_lifetime() const { return lifetime; }

  // Routine to return the number of stochastic stars
  std::vector<double>::size_type get_nstars() const 
  { return stars.size(); }

  // Routines to return the stellar masses, other stellar data, the
  // bolometric luminosity, the spectrum, the wavelength data to
  // go with the spectrum, and the photometry. Where something must be
  // computed, these functions just invoke the corresponding set
  // method and then return.
  const std::vector<double> &get_stars() const { return stars; }

  const std::vector<slug_stardata> &get_isochrone()
  { set_isochrone(); return stardata; }

  // Routines to get bolometric luminosity, with or without extinction
  double get_Lbol() { set_Lbol(); return Lbol; }
  double get_Lbol_extinct() { set_Lbol(); return Lbol_ext; }

  // Methods to get spectrum, with or without extinction and nebular
  // contributions
  const std::vector<double> &get_spectrum();
  const std::vector<double> &get_spectrum_neb();
  const std::vector<double> &get_spectrum_extinct();
  const std::vector<double> &get_spectrum_neb_extinct();

  // Methods to get spectrum plus wavelength table
  void get_spectrum(std::vector<double> &lambda_out, 
		    std::vector<double> &L_lambda_out,
		    bool rest = false);
  void get_spectrum_neb(std::vector<double> &lambda_out, 
			std::vector<double> &L_lambda_out,
			bool rest = false);
  void get_spectrum_extinct(std::vector<double> &lambda_out, 
			    std::vector<double> &L_lambda_out,
			    bool rest = false);
  void get_spectrum_neb_extinct(std::vector<double> &lambda_out, 
				std::vector<double> &L_lambda_out,
				bool rest = false);

  // Methods to get photometry
  const std::vector<double> &get_photometry();
  const std::vector<double> &get_photometry_neb();
  const std::vector<double> &get_photometry_extinct();
  const std::vector<double> &get_photometry_neb_extinct();

  // Method to get the yields
  const std::vector<double> &get_yield();

  // Method to get equivalent widths & line names
  const std::vector<double> &get_ew();
  
  // Method to clear the spectrum and photometry
  void clear_spectrum();

  // Output routines. The optional argument cluster_only is used to
  // flag that only this cluster is being written, rather than a
  // collection of them, and that this cluster should therefore write
  // out its time and the number 1, indicating that it is the only
  // cluster. This only makes a difference if out_mode == BINARY.
  void write_prop(std::ofstream& outfile, const outputMode out_mode, 
		  const unsigned long trial,
		  const bool cluster_only = false,
		  const std::vector<double>& imfvp = {}) const;
  void write_spectrum(std::ofstream& outfile, const outputMode out_mode,
		      const unsigned long trial,
		      const bool cluster_only = false);
  void write_photometry(std::ofstream& outfile, const outputMode out_mode,
			const unsigned long trial,
			const bool cluster_only = false);
  void write_yield(std::ofstream& outfile, const outputMode out_mode,
		   const unsigned long trial,
		   const bool cluster_only = false);
  void write_sn(std::ofstream& outfile, const outputMode out_mode,
		const unsigned long trial,
		const bool cluster_only = false);

#ifdef ENABLE_FITS
  // These are identical to the previous three routines, but they
  // write to a FITS file instead of an ofstream
  void write_prop(fitsfile *out_fits, unsigned long trial, 
      const std::vector<double>& imfvp = {});
  void write_spectrum(fitsfile *out_fits, unsigned long trial);
  void write_photometry(fitsfile *out_fits, unsigned long trial);
  void write_yield(fitsfile *out_fits, unsigned long trial);
  void write_ew(fitsfile *out_fits, unsigned long trial);
  void write_sn(fitsfile *out_fits, unsigned long trial);
#endif

protected:
  slug_ostreams &ostreams;

private:

  // Routines to compute stellar data, bolometric luminosity, spectrum,
  // equivalent widths, yields, winds
  void set_isochrone();
  void set_Lbol();
  void set_spectrum();
  void set_photometry();
  void set_yield();
  void set_ew();
#ifdef WINDS_ON
  void set_winds();
#endif

  // Invariant data
  const double targetMass;            // Target mass
  const slug_PDF *imf;                // IMF
  const slug_PDF *clf;                // CLF
  const slug_tracks *tracks;          // Evolutionary track set
  const slug_specsyn *specsyn;        // Spectral synthesizer
  const slug_filter_set *filters;     // Photometric filters
  const slug_extinction *extinct;     // Extinction
  const slug_nebular *nebular;        // Nebular
  const slug_yields *yields;          // Yields
  const slug_line_list *lines;        // Line list  
  const slug_imf_integrator<double> integ; // Integration helper class

  // Malleable data
  unsigned long id;                   // Unique ID
  unsigned long stoch_sn;             // Number of stochastic SN
  double birthMass;                   // Actual mass at birth
  double aliveMass;                   // Mass of still-living stars
  double stochBirthMass;              // Birth mass of stochastic stars
  double stochAliveMass;              // Mass of still-living stochastic stars
  double nonStochBirthMass;           // Birth mass of non-stochastic stars
  double nonStochAliveMass;           // Mass of still-living non-stochastic stars
  double stochRemnantMass;            // Mass of stochastic stellar remnants
  double nonStochRemnantMass;         // Mass of non-stochastic stellar remnants
  double stellarMass;                 // Mass of all stars including remnants
  double stochStellarMass;            // Stochastic part of stellarMass
  double nonStochStellarMass;         // Non-stochastic part of stellarMass
  double formationTime;               // Cluster birth time
  double curTime;                     // Current time
  double clusterAge;                  // Current age of cluster
  double lifetime;                    // Time at which this cluster disrupts
  double stellarDeathMass;            // Mass of star dying at current age
  double A_V;                         // A_V value
  double A_Vneb;                      // Nebular A_V value
  double Lbol;                        // Bolometric luminosity
  double Lbol_ext;                    // Extincted bolometric luminosity
  double tot_sn;                      // Total number of SNe, stochastic and not
  double last_yield_time;             // Last time yields were updated
#ifdef WINDS_ON
  double wind_mDot;                   // Total wind mass loss rate
  double wind_pDot;                   // Total wind momentum flux
  double wind_lum;                    // Total wind mechanical luminosity
#endif
  std::vector<double> stars;          // List of stellar masses
  std::vector<double> dead_stars;     // Stars that died this time step
  std::vector<slug_stardata> stardata; // Data on stellar properties
  std::vector<double> L_lambda;       // Spectrum
  std::vector<double> phot;           // Photometry
  std::vector<double> L_lambda_ext;   // Spectrum after extinction
  std::vector<double> phot_ext;       // Photometry after extinction
  std::vector<double> L_lambda_neb;   // Specific star+nebular luminosity
  std::vector<double> phot_neb;       // Integrated star+nebular photometry
  std::vector<double> L_lambda_neb_ext; // Spec star+nebular lum w/extinction
  std::vector<double> phot_neb_ext;   // Integ star+nebular phot w/extinction
  std::vector<double> all_yields;     // Yield of all elements up to this point
  std::vector<double> stoch_yields;   // Yields from stochastic stars
  std::vector<double> non_stoch_yields; // Yields from non-stochastic stars
  std::vector<double> ew;             // Equivalent widths
  bool is_disrupted;                  // Is this cluster disrupted?
  bool data_set, Lbol_set, spec_set,
    phot_set, yield_set, ew_set;      // Status indicators
#ifdef WINDS_ON
  bool winds_set;
#endif
  bool stoch_contrib_only;            // Only include the stochastic
                                      // contribution to winds, yields,
                                      // and spectra
  std::vector<double> recspec_wl;     // WL for rectified spectrum
  std::vector<double> recspec;        // Rectified spectrum
};

#endif
// _slug_cluster_H_
