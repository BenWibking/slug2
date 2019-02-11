/*********************************************************************
Copyright (C) 2017 Mark Krumholz
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

#ifndef _slug_mesh2d_interpolator_H_
#define _slug_mesh2d_interpolator_H_

#include <vector>
#include "../constants.H"
#include "slug_interpolation.H"
#include "slug_mesh2d_grid.H"


////////////////////////////////////////////////////////////////////////
// class slug_mesh2d_interpolator
//
// This class performs interpolation of scalar-valued functions on
// semi-regular meshes. See the header to slug_mesh2d_grid.H for an
// explanation of the underlying geometry and logic for this class.
////////////////////////////////////////////////////////////////////////

enum mesh2d_edge_type { mesh2d_xlo, mesh2d_xhi, mesh2d_ylo,
			mesh2d_yhi };

class slug_mesh2d_interpolator {

public:
  
  // Constructor; needs arrays of x, y, and f points; an internal copy
  // of these data is made, so they can be modified or discarded after
  // construction. The option monotonicity-preserving tells the
  // interpolator to set up its coordinate system so that it is
  // possible to preserve monotonicity properties of the underlying
  // data. Monotonicity in a particular variable can be guaranteed
  // only if monotonicity-preserving is true AND the interpolation
  // type is set to gsl_interp_linear; if either of these conditions
  // is not met, monotonicity is not guaranteed
  slug_mesh2d_interpolator(const array2d& x, const array1d& y,
			   const array2d& f_, 
			   const gsl_interp_type *interp_type_ =
			   slug_default_interpolator,
			   const bool monotonicity_preserving = false);

  // Alternate constructor; completely identical except that it takes
  // views for x and f; this is to facilitate creating 2d meshes from
  // sub-arrays of 3d grids
  slug_mesh2d_interpolator(const view2d& x, const array1d& y,
			   const view2d& f_, 
			   const gsl_interp_type *interp_type_ =
			   slug_default_interpolator,
			   const bool monotonicity_preserving = false);

  // Destructor
  ~slug_mesh2d_interpolator();

  // Methods to query the mesh; these are just pass-throughs to the
  // methods of the same name that hang off the slug_mesh2d_grid
  // object
  double x_min(const double y = constants::big) const
  { return grid.x_min(y); }
  double x_max(const double y = constants::big) const
  { return grid.x_max(y); }
  double y_min() const { return grid.y_min(); }
  double y_max() const { return grid.y_max(); }
  std::vector<double> x_lim(const double y = constants::big) const
  { return grid.x_lim(y); }
  std::vector<double> y_lim(const double x = constants::big) const
  { return grid.y_lim(x); }
  double dydx_lo(const double x) const { return grid.dydx_lo(x); }
  double dydx_hi(const double x) const { return grid.dydx_hi(x); }
  bool convex() const { return grid.convex(); }
  const boost::multi_array_types::size_type* shape() const
  { return f.shape(); }
  
  // Methods to construct an interpolating function at fixed x or y;
  // for these functions, the user passes in either uninitialized
  // pointers to or empty arrays of gsl_spline * and gsl_interp_accel
  // * objects. These are allocated within the routine and returned
  // ready to use. The user is responsible for deleting these objects
  // via the appropriate calls for gsl_spline_free and
  // gsl_interp_accel_free. The user can specify that the interpolator
  // is to be build only over a limited range of x or y. Note the the
  // functions x_lim and y_lim can be used to determine the range of
  // validity of these spline functions. It is an error if the
  // supplied values of y or x do not intersect the mesh.
  void build_interp_const_x(const double x,
			    spl_arr_1d& spl,
			    acc_arr_1d& acc,
			    array1d& y_interp_lim,
			    const std::vector<double>& y_lim
			    = std::vector<double>()) const;
  void build_interp_const_y(const double x,
			    gsl_spline *& spl,
			    gsl_interp_accel *& acc,
			    array1d& x_interp_lim,
			    const std::vector<double>& x_lim 
			    = std::vector<double>()) const;
  
  // Method to interpolate to a specified (x, y) point in the mesh
  // interior; the fast_linear flag, if set, causes the interpolation
  // to be linear regardless of the overall interpolation type, and is
  // provided for speed
  double operator()(const double x, const double y,
  		    const bool fast_linear = false) const;

  // Method to interpolate to a specified point on the mesh edge; this
  // should be used to interpolate to points on the mesh edge instead
  // of the previous operator in order to avoid problems locating
  // edges with finite precision arithmetic; pos is the x or y
  // position along the mesh edge
  double operator()(const double pos, const mesh2d_edge_type edge) const;

  // Routine to initalize the interpolators along the tracks; if
  // called after object creation, it re-initializes the interpolation
  // to use the new data passed in, while keeping the mesh unchanged;
  // note that we declare this as a template method so that it can
  // accept either arrays or array views
  template<typename T> void init_interpolators(T& f);

private:

  // The grid
  slug_mesh2d_grid grid;
  const boost::multi_array_types::size_type nx, ny;

  // The data
  array2d f;

  // Interpolation machinery
  const gsl_interp_type *interp_type;
  const bool monotonic_grid;
  const unsigned int interp_npt;
  mutable spl_arr_1d spl_x, spl_s;
  mutable acc_arr_1d acc_x, acc_s;
};


////////////////////////////////////////////////////////////////////////
// class slug_mesh2d_interpolator_vec
//
// Same as slug_mesh2d_interpolator, but for vector-valued functions
////////////////////////////////////////////////////////////////////////

class slug_mesh2d_interpolator_vec {

public:
  
  
  // Constructor
  slug_mesh2d_interpolator_vec(const array2d& x, const array1d& y,
			       const array3d& f_, 
			       const std::vector<const gsl_interp_type *>
			       interp_type_ =
			       std::vector<const gsl_interp_type *>(),
			       const bool monotonicity_preserving = false);

  // Destructor
  ~slug_mesh2d_interpolator_vec();

  // Methods to query the mesh
  double x_min(const double y = constants::big) const
  { return grid.x_min(y); }
  double x_max(const double y = constants::big) const
  { return grid.x_max(y); }
  double y_min() const { return grid.y_min(); }
  double y_max() const { return grid.y_max(); }
  std::vector<double> x_lim(const double y = constants::big) const
  { return grid.x_lim(y); }
  std::vector<double> y_lim(const double x = constants::big) const
  { return grid.y_lim(x); }
  double dydx_lo(const double x) const { return grid.dydx_lo(x); }
  double dydx_hi(const double x) const { return grid.dydx_hi(x); }
  bool convex() const { return grid.convex(); }
  const boost::multi_array_types::size_type* shape() const
  { return f.shape(); }
  
  // Methods to construct an interpolating function at fixed x or y
  void build_interp_const_x(const double x,
			    spl_arr_2d& spl,
			    acc_arr_2d& acc,
			    array1d& y_interp_lim,
			    const std::vector<double>& y_lim
			    = std::vector<double>()) const;
  void build_interp_const_y(const double x,
			    spl_arr_1d& spl,
			    acc_arr_1d& acc,
			    array1d& x_interp_lim,
			    const std::vector<double>& x_lim 
			    = std::vector<double>()) const;
  
  // Methods to interpolate to a specified (x, y) point in mesh interior
  void operator()(const double x, const double y,
		  array1d& f_interp,
		  const bool fast_linear = false) const;
  double operator()(const double x, const double y,
		    const boost::multi_array_types::size_type f_idx,
		    const bool fast_linear = false) const;

  // Methods to interpolate to a point on mesh edge
  void operator()(const double pos, const mesh2d_edge_type edge,
		  array1d& f_interp) const;
  double operator()(const double pos, const mesh2d_edge_type edge,
		    const boost::multi_array_types::size_type f_idx) const;

  // Routine to initalize the interpolators along the tracks; if
  // called after object creation, it re-initializes the interpolation
  // to use the new data passed in, while keeping the mesh unchanged;
  // note that we declare this as a template method so that it can
  // accept either arrays or array views
  template<typename T> void init_interpolators(T& f_);

private:

  // The grid
  slug_mesh2d_grid grid;
  const boost::multi_array_types::size_type nx, ny;

  // The data
  const boost::multi_array_types::size_type nf;
  array3d f;

  // Interpolation machinery
  const std::vector<const gsl_interp_type *> interp_type;
  const bool monotonic_grid;
  mutable spl_arr_2d spl_x, spl_s;
  mutable acc_arr_2d acc_x, acc_s;
};

#endif
// _slug_mesh2d_interpolator_H_
