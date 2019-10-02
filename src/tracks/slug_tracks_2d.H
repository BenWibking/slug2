/*********************************************************************
Copyright (C) 2017 Robert da Silva, Michele Fumagalli, Mark Krumholz
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
// class slug_tracks_2d
//
// This class provides a general interface for evolutionary tracks
// that have two dimensions, meaning that a star is fully specified by
// its mass and age; metallicity is assumed to be the same for all
// stars. Note that many of the functions here take a metallicity
// argument, just for compatibility with the general slug_tracks
// interface, but this value is never used.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_tracks_2d_H_
#define _slug_tracks_2d_H_

#include <vector>
#include <string>
#include "slug_tracks.H"
#include "../slug.H"
#include "../interpolators/slug_mesh2d_interpolator.H"

class slug_tracks_2d : public slug_tracks {

public:
  // Here we just pass the ostreams pointer up to our parent, and
  // initialize some default values that the child class can override
  slug_tracks_2d(slug_ostreams& ostreams_,
		 double metallicity_ = tracks::null_metallicity,
		 ZInterpMethod Z_int_meth_ = Z_UNSPECIFIED,
		 bool force_mcur_monotonic_ = false) :
    slug_tracks(ostreams_),
    interp(nullptr),
    metallicity(metallicity_),
    monotonic(false),
    force_mcur_monotonic(force_mcur_monotonic_),
    Z_int_meth(Z_int_meth_),
    isochrone_age(-1.0)
  { }

  // Destructor
  virtual ~slug_tracks_2d();

  // Maximum and minimum mass in the tracks
  virtual double min_mass() const { return exp(interp->y_min()); }
  virtual double max_mass() const { return exp(interp->y_max()); }

  // Mass(es) of star(s) dying at a particular time; see slug_tracks.H
  // for an explanation of the polymorphism here.
  virtual double
  death_mass(const double time,
	     const double Z = tracks::null_metallicity) const;
  virtual void
  death_mass(const double time,
	     std::vector<double>& m_death,
	     const double Z = tracks::null_metallicity) const;
  
  // Live mass range
  virtual std::vector<double>
  live_mass_range(const double time,
		  const double Z = tracks::null_metallicity) const;

  // Check if the tracks are monotonic (meaning that death mass is a
  // non-increasing function of mass at all masses and ages) or not
  virtual bool check_monotonic() const { return monotonic; }

  // Lifetime of a star of a specified mass; note that this function
  // will return the lifetime of the least / most massive star in the
  // tracks if given a mass that is below / above the mass limits of
  // the tracks
  virtual double
  star_lifetime(const double mass,
		const double Z = tracks::null_metallicity) const {
    double logm = log(mass);
    if (logm < interp->y_min())
      return exp(interp->x_max(interp->y_min()));
    else if (logm > interp->y_max())
      return exp(interp->x_max(interp->y_max()));
    else					   
      return exp(interp->x_max(logm));
  };

  // Derivative of stellar lifetime with respect to stellar mass
  virtual double
  star_lifetime_deriv(const double mass,
		      const double Z = tracks::null_metallicity) const {
    double t = star_lifetime(mass, Z);
    double dlogm_dlogt = interp->dydx_hi(log(mass));
    return t / (mass * dlogm_dlogt);
  };
  
  // Mass at death of a star of a specified starting mass
  virtual double
  star_mass_at_death(const double mass,
		     const double Z = tracks::null_metallicity) const {
    return (*interp)(log(mass), mesh2d_xhi, idx_log_cur_mass);
  }

  // Remnant mass produced by a star of a given mass; if the age
  // argument is set to a non-negative value, the return value will be
  // 0 for stars that have not yet died as of the indicated age
  virtual double
  remnant_mass(const double mass,
	       const double age = -constants::big,
	       const double Z = tracks::null_metallicity) const;
  virtual std::vector<double>
  remnant_mass(const std::vector<double>& mass,
	       const double age = -constants::big,
	       const double Z = tracks::null_metallicity) const;

  // Method to evaluate the properties of a star at a particular mass,
  // age, metallicity
  virtual slug_stardata
  get_star(const double m, const double t,
	   const double Z = tracks::null_metallicity) const;

  // Method to produce vector of stellar data from input vector of
  // mass and fixed time. Input mass m is in Msun, and the time t is
  // in yr.
  virtual std::vector<slug_stardata>
  get_isochrone(const double t, const std::vector<double> &m,
		const double Z = tracks::null_metallicity) const;

  // Method to return whether stars are WR stars, and what type; this
  // is a pure virtual that must be implemented by the derived class,
  // because WR determinations are handled in different ways for
  // different track sets; the two forms differ in that one is
  // intended for use when an isochrone has already been constructed,
  // while the other is for use without one
  virtual void set_WR_type(const double m, const double t,
			   slug_stardata& star) const = 0;
  virtual void set_WR_type(const double m,
			   spl_arr_view_1d& isochrone_,
			   acc_arr_view_1d& isochrone_acc_,
			   slug_stardata& star) const = 0;

#ifdef WINDS_ON
  // Method to set stellar wind properties; this is a pure virtual
  // that must be implemented in a derived class
  virtual void set_wind(const double m, const double t,
			slug_stardata& star) const = 0;
  virtual void set_wind(const double m,
			spl_arr_view_1d& isochrone_,
			acc_arr_view_1d& isochrone_acc_,
			slug_stardata& star) const = 0;
#endif

  // Methods to get the metallicity (normalized to Solar)
  virtual double get_metallicity() const { return metallicity; }

  // Methods to return the list of files and metallicities for the
  // track set we're using; only non-empty if child class sets these
  ZInterpMethod trackset_Z_int_meth() const {
    return Z_int_meth;
  }
  const std::vector<std::string>& trackset_filenames() const {
    return filenames;
  }
  const std::vector<double>& trackset_metallicities() const {
    return Z_files;
  }
  
protected:
  
  // The interpolation machinery; must be set up by derived class
  slug_mesh2d_interpolator_vec *interp;
  
  // The metallicity of these tracks
  double metallicity;

  // Flag if tracks are monotonic
  bool monotonic;

  // Flag if we are going to use an interpolation scheme that
  // guarantees that the current stellar mass is monotonically
  // decreasing for all stars; this is in general not guaranteed for
  // higher-order and more accurate interpolation methods
  const bool force_mcur_monotonic;

  // File names and metallicities, and interpolation method, for track
  // sets
  std::vector<std::string> filenames;
  std::vector<double> Z_files;
  ZInterpMethod Z_int_meth;

private:

  // Machinery to free a cached isochrone; marked as const because it
  // only affects mutable objects
  void free_isochrone() const;
  
  // Storage to cache a copy of the interpolation machinery for an
  // isochrone; marked as mutable so that it can altered even by const
  // routines
  mutable double isochrone_age;
  mutable spl_arr_2d isochrone;
  mutable acc_arr_2d isochrone_acc;
  mutable array1d isochrone_logm_lim;
  
};

#endif
// _slug_tracks_2d_H_
