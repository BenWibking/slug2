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
// class slug_tracks
//
// The slug_tracks class and its children define how slug deals with
// evolutionary tracks, and provide the machinery to build isochrones,
// calculate the properties of individual stars, etc. The base
// slug_tracks class defined in this file is a pure virtual class the
// defines the common interface to all derived evolutionary track
// classes.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_tracks_H_
#define _slug_tracks_H_

#include <vector>
#include <string>
#include "../constants.H"
#include "../slug.H"
#include "../slug_IO.H"
#include "../interpolators/slug_interpolation.H"

namespace tracks {
  const double null_metallicity = -999.;
}

// Enum of known track sets
enum trackSet { NO_TRACK_SET, GENEVA_2013_VVCRIT_00, GENEVA_2013_VVCRIT_40,
		GENEVA_MDOT_STD, GENEVA_MDOT_ENHANCED,
		PADOVA_TPAGB_YES, PADOVA_TPAGB_NO, MIST_2016_VVCRIT_00,
		MIST_2016_VVCRIT_40, PARSEC_V1_2 };

// Methods to interpolate tracks in metallicity
enum ZInterpMethod { Z_UNSPECIFIED, Z_NEAR_NEIGHBOR, Z_LINEAR,
		     Z_AKIMA, Z_STEFFEN };

class slug_tracks {

public:
  // A trivial initialization routine; the real initialization is done
  // in the derived child classes
  slug_tracks(slug_ostreams& ostreams_) : ostreams(ostreams_) { }

  // Destructor
  virtual ~slug_tracks() { }

  // Maximum and minimum mass in the tracks
  virtual double min_mass() const = 0;
  virtual double max_mass() const = 0;

  // Return the mass(es) of star(s) dying at a particular time; this
  // is a polymorphic function with two forms. The first just takes a
  // time (and optionally a metallicity) as an argument, and returns a
  // single mass of star dying. This form is only suitable for tracks
  // that are monotonic, for which a unique death mass exists, and it
  // is an error to use this form for non-monotonic tracks. The second
  // form takes as an argument a vector of death masses, which on
  // return will contain the masses of all stars dying at a particular
  // time (of which there can be an arbitrary number). This form may
  // be called for either monotonic or non-monotonic tracks.
  virtual double
  death_mass(const double time,
	     const double Z = tracks::null_metallicity) const = 0;
  virtual void
  death_mass(const double time,
	     std::vector<double>& m_death,
	     const double Z = tracks::null_metallicity) const = 0;

  // Live mass range
  virtual std::vector<double>
  live_mass_range(const double time,
		  const double Z = tracks::null_metallicity) const = 0;

  // Check if the tracks are monotonic (meaning that death mass is a
  // always a non-increasing function of mass) or not
  virtual bool check_monotonic() const = 0;

  // Lifetime of a star of a specified mass
  virtual double
  star_lifetime(const double mass,
		const double Z = tracks::null_metallicity) const = 0;

  // Derivative of stellar lifetime with respect to stellar mass
  virtual double
  star_lifetime_deriv(const double mass,
		      const double Z = tracks::null_metallicity) const = 0;
  
  // Mass at death of a star of a specified starting mass
  virtual double
  star_mass_at_death(const double mass,
		     const double Z = tracks::null_metallicity) const = 0;

  // Remnant mass produced by a star of a given mass; if the age
  // argument is set to a non-negative value, the return value will be
  // 0 for stars that have not yet died as of the indicated age
  virtual double
  remnant_mass(const double mass,
	       const double age = -constants::big,
	       const double Z = tracks::null_metallicity) const = 0;
  virtual std::vector<double>
  remnant_mass(const std::vector<double>& mass,
	       const double age = -constants::big,
	       const double Z = tracks::null_metallicity) const = 0;

  // Method to evaluate the properties of a star at a particular mass,
  // age, metallicity
  virtual slug_stardata
  get_star(const double m, const double t,
	   const double Z = tracks::null_metallicity) const = 0;

  // Method to produce vector of stellar data from input vector of
  // mass and fixed time.
  virtual std::vector<slug_stardata>
  get_isochrone(const double t, const std::vector<double> &m,
		const double Z = tracks::null_metallicity) const = 0;

  // Method to return metallicity
  virtual double get_metallicity() const = 0;

protected:
  
  // IO handler
  slug_ostreams& ostreams;

  // Indices of variables
  static const
  boost::multi_array_types::size_type idx_log_cur_mass = 0,
    idx_log_L = 1,
    idx_log_Teff = 2,
    idx_log_mDot = 3;
};

#endif
// _slug_tracks_H_
