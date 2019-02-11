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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "slug_mesh2d_interpolator.H"
extern "C" {
#   include <gsl/gsl_errno.h>
}

using namespace std;
using namespace boost::multi_array_types;

////////////////////////////////////////////////////////////////////////
// class slug_mesh2d_interpolator, constructor and destructor
////////////////////////////////////////////////////////////////////////

slug_mesh2d_interpolator::
slug_mesh2d_interpolator(const array2d& x, const array1d& y,
			 const array2d& f_, 
			 const gsl_interp_type *interp_type_,
			 const bool monotonicity_preserving) :
  grid(x, y),
  nx(x.shape()[0]),
  ny(x.shape()[1]),
  interp_type(interp_type_),
  monotonic_grid(monotonicity_preserving),
  interp_npt(gsl_interp_type_min_size(interp_type_)) {

  // Allocate memory for the spine interpolation arrays; note that we
  // do not allocate the splines themselves here, because in some
  // directions we may have redundant points
  spl_x.resize(boost::extents[ny]);
  acc_x.resize(boost::extents[ny]);
  spl_s.resize(boost::extents[nx]);
  acc_s.resize(boost::extents[nx]);
  for (size_type j=0; j<ny; j++) {
    spl_x[j] = nullptr;
    acc_x[j] = gsl_interp_accel_alloc();
  }
  for (size_type i=0; i<nx; i++) {
    spl_s[i] = nullptr;
    acc_s[i] = gsl_interp_accel_alloc();
  }

  // Load the data into the interpolators
  init_interpolators(f_);
}

slug_mesh2d_interpolator::
slug_mesh2d_interpolator(const view2d& x, const array1d& y,
			 const view2d& f_, 
			 const gsl_interp_type *interp_type_,
			 const bool monotonicity_preserving) :
  grid(x, y),
  nx(x.shape()[0]),
  ny(x.shape()[1]),
  interp_type(interp_type_),
  monotonic_grid(monotonicity_preserving),
  interp_npt(gsl_interp_type_min_size(interp_type_)) {

  // Allocate memory for the spine interpolation arrays
  spl_x.resize(boost::extents[ny]);
  acc_x.resize(boost::extents[ny]);
  spl_s.resize(boost::extents[nx]);
  acc_s.resize(boost::extents[nx]);
  for (size_type j=0; j<ny; j++) {
    if (nx >= interp_npt)
      spl_x[j] = gsl_spline_alloc(interp_type, nx);
    else
      spl_x[j] = gsl_spline_alloc(gsl_interp_linear, nx);
    acc_x[j] = gsl_interp_accel_alloc();
  }
  for (size_type i=0; i<nx; i++) {
    if (ny >= interp_npt)
      spl_s[i] = gsl_spline_alloc(interp_type, ny);
    else
      spl_s[i] = gsl_spline_alloc(gsl_interp_linear, ny);
    acc_s[i] = gsl_interp_accel_alloc();
  }

  // Load the data into the interpolators
  init_interpolators(f_);
}

slug_mesh2d_interpolator::~slug_mesh2d_interpolator() {
  for (size_type j=0; j<ny; j++) {
    gsl_spline_free(spl_x[j]);
    gsl_interp_accel_free(acc_x[j]);
  }
  for (size_type i=0; i<nx; i++) {
    gsl_spline_free(spl_s[i]);
    gsl_interp_accel_free(acc_s[i]);
  }
}


////////////////////////////////////////////////////////////////////////
// Method to initialize the "spine" interpolators along the tracks
////////////////////////////////////////////////////////////////////////

template<typename T>
void slug_mesh2d_interpolator::
init_interpolators(T& f_) {

  // Safety check
  assert(f_.num_dimensions() == 2);
  assert(f_.shape()[0] == nx);
  assert(f_.shape()[1] == ny);

  // Copy the data
  f.resize(boost::extents[nx][ny]);
  f = f_;

  // Initialize the spines in the x direction; be careful to handle
  // case where grid points are overlapping
  double *x_tmp = new double[nx];
  double *f_tmp = new double[nx];
  for (size_type j=0; j<ny; j++) {
    size_type count = 0;
    for (size_type i=0; i<nx; i++) {
      if (i > 0) {
	if (grid.x_grid()[i-1][j] == grid.x_grid()[i][j]) continue;
      }
      x_tmp[count] = grid.x_grid()[i][j];
      f_tmp[count] = f[i][j];
      count++;
    }
    if (spl_x[j] == nullptr) {
      if (count >= interp_npt)
	spl_x[j] = gsl_spline_alloc(interp_type, count);
      else
	spl_x[j] = gsl_spline_alloc(gsl_interp_linear, count);
    }
    gsl_spline_init(spl_x[j], x_tmp, f_tmp, count);
  }
  delete[] x_tmp;
  delete[] f_tmp;

  // Initialize the spines in the y direction
  double *s_tmp = new double[ny];
  f_tmp = new double[ny];
  for (size_type i=0; i<nx; i++) {
    size_type count = 0;
    for (size_type j=0; j<ny; j++) {
      if (j > 0) {
	if (grid.s_grid()[i][j-1] == grid.s_grid()[i][j]) continue;
      }
      if (!monotonic_grid) s_tmp[count] = grid.s_grid()[i][j];
      else s_tmp[count] = grid.y_grid()[j];
      f_tmp[count] = f[i][j];
      count++;
    }
    if (spl_s[i] == nullptr) {
      if (count >= interp_npt)
	spl_s[i] = gsl_spline_alloc(interp_type, count);
      else
	spl_s[i] = gsl_spline_alloc(gsl_interp_linear, count);
    }
    gsl_spline_init(spl_s[i], s_tmp, f_tmp, count);
  }
  delete[] s_tmp;
  delete[] f_tmp;
}


////////////////////////////////////////////////////////////////////////
// Method to interpolate to get the value of a single point in the
// mesh interior
////////////////////////////////////////////////////////////////////////
double slug_mesh2d_interpolator::
operator()(const double x, const double y,
	   const bool fast_linear) const {

  // Safety assertion
  assert(grid.in_mesh(x,y));

  // Decide on method
  if (fast_linear) {

    // Find point in mesh
    size_type i, j;
    grid.ij_index(x, y, i, j);

    // Compute x and y weights
    const array2d& x_gr = grid.x_grid();
    const array1d& y_gr = grid.y_grid();
    double y_wgt = (y_gr[j+1] - y) / (y_gr[j+1] - y_gr[j]);
    double x_l = y_wgt * x_gr[i][j] + (1.0-y_wgt) * x_gr[i][j+1];
    double x_r = y_wgt * x_gr[i+1][j] + (1.0-y_wgt) * x_gr[i+1][j+1];
    double x_wgt = (x_r - x) / (x_r - x_l);

    // Return value
    return x_wgt * y_wgt * f[i][j] +
      (1.0-x_wgt) * y_wgt * f[i+1][j] +
      x_wgt * (1.0-y_wgt) * f[i][j+1] +
      (1.0-x_wgt) * (1.0-y_wgt) * f[i+1][j+1];

  } else {

    // This method gives back the same result as building an isochrone
    // at this x and then using it
  
    // Get back points required to build an interpolator in y around
    // this point
    vector<double> int_y, pos;
    vector<size_type> idx, edge;
    grid.intercept_const_x_n(x, y, interp_npt, int_y, pos,
			     idx, edge);

    // Turn off GSL error handling so that we can handle this on our own
    // and print out more informative error messages
    gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

    // Use the interpolators along the spines to generate the function
    // value at all intersection points
    vector<double> f_tmp(pos.size());
    for (vector<double>::size_type i=0; i<pos.size(); i++) {
      int gsl_errstat;
      gsl_spline *spl;
      gsl_interp_accel *acc;
      if (edge[i] == 0) {
        spl = spl_s[idx[i]];
        acc = acc_s[idx[i]];
	if (monotonic_grid)
	  pos[i] = int_y[i]; // If vertical position coordinate is y
			     // instead of s
      } else {
        spl = spl_x[idx[i]];
        acc = acc_x[idx[i]];
      }
      gsl_errstat = gsl_spline_eval_e(spl, pos[i], acc, f_tmp.data()+i);
      if (gsl_errstat) {
	stringstream ss;
	ss << "slug_mesh2d_interpolator::build_interp_const_x_n: "
	   << "bad interpolation evaluation at x = "
	   << setprecision(20) << x << ", pos = "
	   << setprecision(20) << pos[i]
	   << ", interpolation direction = " << edge[i]
	   << "; spline range limits are "
           << spl->x[0] << " to "
           << spl->x[spl->size-1]
           << "; input (x, y) are ("
           << x << ", " << y << ")"
           << "; gsl says: "
	   << gsl_strerror(gsl_errstat);
	throw runtime_error(ss.str());
      }
    }

    // Restore the GSL error handler
    gsl_set_error_handler(gsl_err);
    
    // Build interpolator
    gsl_interp *interp;
    if (int_y.size() == interp_npt)
      interp = gsl_interp_alloc(interp_type, int_y.size());
    else
      interp = gsl_interp_alloc(gsl_interp_linear, int_y.size());
    gsl_interp_init(interp, int_y.data(), f_tmp.data(), int_y.size());
    
    // Interpolate to desired point
    double f_interp
      = gsl_interp_eval(interp, int_y.data(), f_tmp.data(), y, nullptr);
    
    // Free interpolator
    gsl_interp_free(interp);
    
    // Return point
    return f_interp;
  }
}

////////////////////////////////////////////////////////////////////////
// Method to interpolate to get the value of a single point on the
// mesh edge
////////////////////////////////////////////////////////////////////////
double slug_mesh2d_interpolator::
operator()(const double pos, const mesh2d_edge_type edge) const {

  // Action depends on edge
  double f;
  switch (edge) {

  case mesh2d_xlo: // Fall through
  case mesh2d_xhi: {

    // Safety assertions
    double y = pos;
    assert(y >= grid.y_min());
    assert(y <= grid.y_max());
    
    // Get the s coordinate along this edge
    double s;
    size_type idx;
    double x;
    if (edge == mesh2d_xlo) {
      idx = 0;
      x = grid.x_min(y);
    } else {
      idx = nx-1;
      x = grid.x_max(y);
    }
    if (!monotonic_grid) {
      size_type j = grid.j_index(y);
      s = grid.s_grid()[idx][j] +
	sqrt((x - grid.x_grid()[idx][j]) * (x - grid.x_grid()[idx][j]) +
	     (y - grid.y_grid()[j]) * (y - grid.y_grid()[j]));
    } else {
      s = y;
    }
    
    // Interpolate to point
    f = gsl_spline_eval(spl_s[idx], s, acc_s[idx]);
    break;
  }
    
  case mesh2d_ylo: // Fall through
  case mesh2d_yhi: {

    // Safety assertions
    double x = pos;
    size_type idx;
    if (edge == mesh2d_ylo) {
      idx = 0;
    } else {
      idx = ny-1;
    }
    assert(x >= grid.x_grid()[0][idx]);
    assert(x <= grid.x_grid()[nx-1][idx]);

    // Interpolate to point
    f = gsl_spline_eval(spl_x[idx], x, acc_x[idx]);
    break;
  }
  };

  // Return
  return f;
}

////////////////////////////////////////////////////////////////////////
// Methods to build interpolators on the mesh
////////////////////////////////////////////////////////////////////////

// Build an interpolator at constant x
void slug_mesh2d_interpolator::
build_interp_const_x(const double x,
		     spl_arr_1d& spl,
		     acc_arr_1d& acc,
		     array1d& y_interp_lim,
		     const vector<double>& y_lim) const {

  // Safety assertion
  assert(x >= x_min() && x <= x_max());

  // Generate the intersection list
  vector<double> y, pos;
  vector<size_type> idx, edge, seq;
  grid.intercept_const_x(x, y, pos, idx, edge, seq, y_lim);

  // Safety check
  for (vector<double>::size_type i=0; i<y.size()-1; i++)
    assert(y[i] < y[i+1]);
  
  // Turn off GSL error handling so that we can handle this on our own
  // and print out more informative error messages
  gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

  // Use the interpolators along the spines to generate the function
  // value at all intersection points
  vector<double> f_tmp(pos.size());
  for (vector<double>::size_type i=0; i<pos.size(); i++) {
    int gsl_errstat;
    gsl_spline *spl;
    gsl_interp_accel *acc;
    if (edge[i] == 0) {
      spl = spl_s[idx[i]];
      acc = acc_s[idx[i]];
      if (monotonic_grid)
	pos[i] = y[i]; // If using y instead of s as vertical
		       // coordinate
    } else {
      spl = spl_x[idx[i]];
      acc = acc_x[idx[i]];
    }
    gsl_errstat = gsl_spline_eval_e(spl, pos[i], acc, f_tmp.data()+i);
    if (gsl_errstat) {
      stringstream ss;
      ss << "slug_mesh2d_interpolator::build_interp_const_x: "
	 << "bad interpolation evaluation at x = "
	 << setprecision(20) << x << ", pos = "
	 << setprecision(20) << pos[i]
	 << ", interpolation direction = " << edge[i]
         << "; spline range limits are "
         << spl->x[0] << " to "
         << spl->x[spl->size-1]
         << "; input x is "
         << x << ", full lits of y values is:";
      for (vector<int>::size_type j=0; j<y.size(); j++) ss << "  " << y[j];
      ss << "; gsl says: "
	 << gsl_strerror(gsl_errstat);
      throw runtime_error(ss.str());
    }
  }

  // Restore the GSL error handler
  gsl_set_error_handler(gsl_err);

  // Now build new interpolators; we need one for each distinct segment
  spl.resize(boost::extents[seq.size()/2]);
  acc.resize(boost::extents[seq.size()/2]);
  y_interp_lim.resize(boost::extents[seq.size()]);
  for (vector<double>::size_type i=0; i<seq.size(); i+=2) {
    const gsl_interp_type *itype;
    if (seq[i+1]-seq[i] >= interp_npt)
      itype = interp_type;
    else
      itype = gsl_interp_linear;
    spl[i/2] = gsl_spline_alloc(itype, seq[i+1]-seq[i]);
    acc[i/2] = gsl_interp_accel_alloc();
    gsl_spline_init(spl[i/2], y.data()+seq[i], f_tmp.data()+seq[i],
		    seq[i+1]-seq[i]);
    y_interp_lim[i] = y[seq[i]];
    y_interp_lim[i+1] = y[seq[i+1]-1];
  }
}

  
// Build an interpolator at constant y
void slug_mesh2d_interpolator::
build_interp_const_y(const double y,
		     gsl_spline *& spl,
		     gsl_interp_accel *& acc,
		     array1d& x_interp_lim,
		     const vector<double>& x_lim) const {

  // Safety assertion
  assert(y >= y_min() && y <= y_max());

  // Generate the intersection list
  vector<double> x, pos;
  vector<size_type> idx;
  grid.intercept_const_y(y, x, pos, idx, x_lim);

  // Safety check
  for (vector<double>::size_type i=0; i<x.size()-1; i++)
    assert(x[i] < x[i+1]);
  
  // Turn off GSL error handling so that we can handle this on our own
  // and print out more informative error messages
  gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

  // Use the interpolators along the spines to generate the function
  // value at all intersection points
  vector<double> f_tmp(pos.size());
  for (vector<double>::size_type i=0; i<pos.size(); i++) {
    int gsl_errstat;
    double f_eval;
    if (!monotonic_grid) {
      gsl_errstat =
	gsl_spline_eval_e(spl_s[idx[i]], pos[i], acc_s[idx[i]], &f_eval);
    } else {
      gsl_errstat =
	gsl_spline_eval_e(spl_s[idx[i]], y, acc_s[idx[i]], &f_eval);
    }
    f_tmp[i] = f_eval;
    if (gsl_errstat) {
      stringstream ss;
      ss << "slug_mesh2d_interpolator::build_interp_const_y: "
	 << "bad interpolation evaluation at y = "
	 << setprecision(20) << y << ", x = "
	 << setprecision(20) << x[i] << ", pos = "
	 << setprecision(20) << pos[i]
	 << "; gsl says: "
	 << gsl_strerror(gsl_errstat);
      throw runtime_error(ss.str());
    }
  }

  // Return error handler to prior state
  gsl_set_error_handler(gsl_err);

  // Now build new interpolator
  if (x.size() >= interp_npt)
    spl = gsl_spline_alloc(interp_type, x.size());
  else
    spl = gsl_spline_alloc(gsl_interp_linear, x.size());
  acc = gsl_interp_accel_alloc();
  gsl_spline_init(spl, x.data(), f_tmp.data(), x.size());
  x_interp_lim.resize(boost::extents[2]);
  x_interp_lim[0] = x[0];
  x_interp_lim[1] = x.back();
}





////////////////////////////////////////////////////////////////////////
// class slug_mesh2d_interpolator_vec, constructor and destructor
////////////////////////////////////////////////////////////////////////

slug_mesh2d_interpolator_vec::
slug_mesh2d_interpolator_vec(const array2d& x, const array1d& y,
			     const array3d& f_, 
			     const vector<const gsl_interp_type *>
			     interp_type_,
			     const bool monotonicity_preserving) :
  grid(x, y),
  nx(f_.shape()[0]),
  ny(f_.shape()[1]),
  nf(f_.shape()[2]),
  interp_type(interp_type_),
  monotonic_grid(monotonicity_preserving) {

  // Allocate memory for the spine interpolation arrays
  spl_x.resize(boost::extents[ny][nf]);
  acc_x.resize(boost::extents[ny][nf]);
  spl_s.resize(boost::extents[nx][nf]);
  acc_s.resize(boost::extents[nx][nf]);
  for (size_type n=0; n<nf; n++) {
    // Allocate memory
    for (size_type j=0; j<ny; j++) {
      spl_x[j][n] = nullptr;
      acc_x[j][n] = gsl_interp_accel_alloc();
    }
    for (size_type i=0; i<nx; i++) {
      spl_s[i][n] = nullptr;
      acc_s[i][n] = gsl_interp_accel_alloc();
    }
  }

  // Load the data into the interpolators
  init_interpolators(f_);
}

slug_mesh2d_interpolator_vec::~slug_mesh2d_interpolator_vec() {
  for (size_type n=0; n<nf; n++) {
    for (size_type j=0; j<ny; j++) {
      gsl_spline_free(spl_x[j][n]);
      gsl_interp_accel_free(acc_x[j][n]);
    }
    for (size_type i=0; i<nx; i++) {
      gsl_spline_free(spl_s[i][n]);
      gsl_interp_accel_free(acc_s[i][n]);
    }
  }
}

////////////////////////////////////////////////////////////////////////
// Method to initialize the "spine" interpolators along the tracks
////////////////////////////////////////////////////////////////////////

template<typename T>
void slug_mesh2d_interpolator_vec::
init_interpolators(T& f_) {

  // Safety check
  assert(f_.num_dimensions() == 3);
  assert(f_.shape()[0] == nx);
  assert(f_.shape()[1] == ny);
  assert(f_.shape()[2] == nf);

  // Copy the data
  f.resize(boost::extents[nx][ny][nf]);
  f = f_;

  // Initialize the spines in the x direction
  double *x_tmp = new double[nx];
  double *f_tmp = new double[nx];
  for (size_type j=0; j<ny; j++) {
    size_type count = 0;
    for (size_type i=0; i<nx; i++) {
      if (i > 0) {
	if (grid.x_grid()[i-1][j] == grid.x_grid()[i][j]) continue;
      }
      x_tmp[count] = grid.x_grid()[i][j];
      count++;
    }
    for (size_type n=0; n<nf; n++) {
      count = 0;
      for (size_type i=0; i<nx; i++) {
	if (i > 0) {
	  if (grid.x_grid()[i-1][j] == grid.x_grid()[i][j]) continue;
	}
	f_tmp[count] = f[i][j][n];
	count++;
      }
      if (spl_x[j][n] == nullptr) {
	// Set interpolation type for this quantity
	const gsl_interp_type *itype;
	if (interp_type.size() == 0) itype = slug_default_interpolator;
	else itype = interp_type[n];
	unsigned int interp_npt = gsl_interp_type_min_size(itype);
	if (count >= interp_npt)
	  spl_x[j][n] = gsl_spline_alloc(itype, count);
	else
	  spl_x[j][n] = gsl_spline_alloc(gsl_interp_linear, count);
      }
      gsl_spline_init(spl_x[j][n], x_tmp, f_tmp, count);
    }
  }
  delete[] x_tmp;
  delete[] f_tmp;

  // Initialize the spines in the y direction
  double *s_tmp = new double[ny];
  f_tmp = new double[ny];
  for (size_type i=0; i<nx; i++) {
    size_type count = 0;
    for (size_type j=0; j<ny; j++) {
      if (j > 0) {
	if (grid.s_grid()[i][j-1] == grid.s_grid()[i][j]) continue;
      }
      if (!monotonic_grid)
	s_tmp[count] = grid.s_grid()[i][j];
      else
	s_tmp[count] = grid.y_grid()[j];
      count++;
    }
    for (size_type n=0; n<nf; n++) {
      count = 0;
      for (size_type j=0; j<ny; j++) {
	if (j > 0) {
	  if (grid.s_grid()[i][j-1] == grid.s_grid()[i][j]) continue;
	}
	f_tmp[count] = f[i][j][n];
	count++;
      }
      if (spl_s[i][n] == nullptr) {
	// Set interpolation type for this quantity
	const gsl_interp_type *itype;
	if (interp_type.size() == 0) itype = slug_default_interpolator;
	else itype = interp_type[n];
	unsigned int interp_npt = gsl_interp_type_min_size(itype);
	if (count >= interp_npt)
	  spl_s[i][n] = gsl_spline_alloc(itype, count);
	else
	  spl_s[i][n] = gsl_spline_alloc(gsl_interp_linear, count);
      }
      gsl_spline_init(spl_s[i][n], s_tmp, f_tmp, ny);
    }
  }
  delete[] s_tmp;
  delete[] f_tmp;
}

////////////////////////////////////////////////////////////////////////
// Method to interpolate to get the value of a single point in the
// mesh interior 
////////////////////////////////////////////////////////////////////////
void slug_mesh2d_interpolator_vec::
operator()(const double x, const double y, array1d& f_interp,
	   const bool fast_linear) const {

  // Safety assertion
  assert(grid.in_mesh(x,y));

  // Resize output array
  f_interp.resize(boost::extents[nf]);

  // Decide on method
  if (fast_linear) {

    // Find point in mesh
    size_type i, j;
    grid.ij_index(x, y, i, j);

    // Compute x and y weights
    const array2d& x_gr = grid.x_grid();
    const array1d& y_gr = grid.y_grid();
    double y_wgt = (y_gr[j+1] - y) / (y_gr[j+1] - y_gr[j]);
    double x_l = y_wgt * x_gr[i][j] + (1.0-y_wgt) * x_gr[i][j+1];
    double x_r = y_wgt * x_gr[i+1][j] + (1.0-y_wgt) * x_gr[i+1][j+1];
    double x_wgt = (x_r - x) / (x_r - x_l);

    // Fill output array
    for (size_type n=0; n<nf; n++)
      f_interp[n] =  x_wgt * y_wgt * f[i][j][n] +
      (1.0-x_wgt) * y_wgt * f[i+1][j][n] +
      x_wgt * (1.0-y_wgt) * f[i][j+1][n] +
      (1.0-x_wgt) * (1.0-y_wgt) * f[i+1][j+1][n];

  } else {

    // This method gives back the same result as building an isochrone
    // at this x and then using it
  
    // Get back points required to build an interpolator in y around
    // this point
    vector<double> int_y, pos;
    vector<size_type> idx, edge;
    vector<unsigned int> interp_npt(nf);
    unsigned int interp_npt_max = 0;
    if (interp_type.size() > 0) {
      for (size_type n=0; n<nf; n++) {
	interp_npt[n] = gsl_interp_type_min_size(interp_type[n]);
	if (interp_npt[n] > interp_npt_max) interp_npt_max = interp_npt[n];
      }
    } else {
      for (size_type n=0; n<nf; n++)
	interp_npt[n] = gsl_interp_type_min_size(slug_default_interpolator);
      interp_npt_max = interp_npt[0];
    }
    grid.intercept_const_x_n(x, y, interp_npt_max, int_y, pos,
			     idx, edge);

    // Turn off GSL error handling so that we can handle this on our own
    // and print out more informative error messages
    gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

    // Loop over elements of the data
    for (size_type n=0; n<nf; n++) {
      
      // Use the interpolators along the spines to generate the function
      // value at all intersection points
      vector<double> f_tmp(pos.size());
      for (vector<double>::size_type i=0; i<pos.size(); i++) {
	int gsl_errstat;
	double f_eval;
	if (edge[i] == 0) {
	  if (!monotonic_grid)
	    gsl_errstat =
	      gsl_spline_eval_e(spl_s[idx[i]][n], pos[i], acc_s[idx[i]][n],
				&f_eval);
	  else
	    gsl_errstat =
	      gsl_spline_eval_e(spl_s[idx[i]][n], y, acc_s[idx[i]][n],
				&f_eval);
	} else
	  gsl_errstat =
	    gsl_spline_eval_e(spl_x[idx[i]][n], pos[i], acc_x[idx[i]][n],
			      &f_eval);
	f_tmp[i] = f_eval;
	if (gsl_errstat) {
	  stringstream ss;
	  ss << "slug_mesh2d_interpolator_vec::operator(): "
	     << "bad interpolation evaluation at x = "
	     << setprecision(20) << x << ", pos = "
	     << setprecision(20) << pos[i]
	     << ", interpolation direction = " << edge[i]
	     << "; gsl says: "
	     << gsl_strerror(gsl_errstat);
	  throw runtime_error(ss.str());
	}
      }

      // Set interpolator type for this quantity
      const gsl_interp_type *itype;
      if (interp_type.size() == 0) itype = slug_default_interpolator;
      else itype = interp_type[n];
      
      // Build interpolator
      gsl_interp *interp;
      if (int_y.size() == interp_npt[n])
	interp = gsl_interp_alloc(itype, int_y.size());
      else
	interp = gsl_interp_alloc(gsl_interp_linear, int_y.size());
      gsl_interp_init(interp, int_y.data(), f_tmp.data(), int_y.size());
    
      // Interpolate to desired point
      f_interp[n]
	= gsl_interp_eval(interp, int_y.data(), f_tmp.data(), y, nullptr);
    
      // Free interpolator
      gsl_interp_free(interp);
	
    }

    // Restore the GSL error handler
    gsl_set_error_handler(gsl_err);

  }
}

////////////////////////////////////////////////////////////////////////
// Method to interpolate to get the value of a single point for a
// single quantity in the mesh interior
////////////////////////////////////////////////////////////////////////
double slug_mesh2d_interpolator_vec::
operator()(const double x, const double y, const size_type f_idx,
	   const bool fast_linear) const {

  // Safety assertions
  assert(grid.in_mesh(x,y));
  assert(f_idx < nf);

  // Decide on method
  if (fast_linear) {

    // Find point in mesh
    size_type i, j;
    grid.ij_index(x, y, i, j);

    // Compute x and y weights
    const array2d& x_gr = grid.x_grid();
    const array1d& y_gr = grid.y_grid();
    double y_wgt = (y_gr[j+1] - y) / (y_gr[j+1] - y_gr[j]);
    double x_l = y_wgt * x_gr[i][j] + (1.0-y_wgt) * x_gr[i][j+1];
    double x_r = y_wgt * x_gr[i+1][j] + (1.0-y_wgt) * x_gr[i+1][j+1];
    double x_wgt = (x_r - x) / (x_r - x_l);

    // Fill output array
    double f_interp = x_wgt * y_wgt * f[i][j][f_idx] +
      (1.0-x_wgt) * y_wgt * f[i+1][j][f_idx] +
      x_wgt * (1.0-y_wgt) * f[i][j+1][f_idx] +
      (1.0-x_wgt) * (1.0-y_wgt) * f[i+1][j+1][f_idx];
    return f_interp;

  } else {

    // This method gives back the same result as building an isochrone
    // at this x and then using it
  
    // Get back points required to build an interpolator in y around
    // this point
    vector<double> int_y, pos;
    vector<size_type> idx, edge;
    unsigned int interp_npt = gsl_interp_type_min_size(interp_type[f_idx]);
    grid.intercept_const_x_n(x, y, interp_npt, int_y, pos,
			     idx, edge);

    // Turn off GSL error handling so that we can handle this on our own
    // and print out more informative error messages
    gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

    // Use the interpolators along the spines to generate the function
    // value at all intersection points
    vector<double> f_tmp(pos.size());
    for (vector<double>::size_type i=0; i<pos.size(); i++) {
      int gsl_errstat;
      double f_eval;
      gsl_spline *spl;
      gsl_interp_accel *acc;
      if (edge[i] == 0) {
	spl = spl_s[idx[i]][f_idx];
	acc = acc_s[idx[i]][f_idx];
	if (monotonic_grid)
	  pos[i] = int_y[i]; // If using y instead of s as coordinate
			     // along spines
      } else {
	spl = spl_x[idx[i]][f_idx];
	acc = acc_x[idx[i]][f_idx];
      }
      gsl_errstat =
	gsl_spline_eval_e(spl, pos[i], acc, &f_eval);
      f_tmp[i] = f_eval;
      if (gsl_errstat) {
	stringstream ss;
	ss << "slug_mesh2d_interpolator_vec::operator(): "
	   << "bad interpolation evaluation at x = "
	   << setprecision(20) << x << ", pos = "
	   << setprecision(20) << pos[i]
	   << ", interpolation direction = " << edge[i]
	   << "; gsl says: "
	   << gsl_strerror(gsl_errstat);
	throw runtime_error(ss.str());
      }
    }

    // Set interpolator type for this quantity
    const gsl_interp_type *itype;
    if (interp_type.size() == 0) itype = slug_default_interpolator;
    else itype = interp_type[f_idx];
      
    // Build interpolator
    gsl_interp *interp;
    if (int_y.size() == interp_npt)
      interp = gsl_interp_alloc(itype, int_y.size());
    else
      interp = gsl_interp_alloc(gsl_interp_linear, int_y.size());
    gsl_interp_init(interp, int_y.data(), f_tmp.data(), int_y.size());
    
    // Interpolate to desired point
    double f_interp
      = gsl_interp_eval(interp, int_y.data(), f_tmp.data(), y, nullptr);
    
    // Free interpolator
    gsl_interp_free(interp);

    // Restore the GSL error handler
    gsl_set_error_handler(gsl_err);

    // Return
    return f_interp;
  }
}

////////////////////////////////////////////////////////////////////////
// Method to interpolate to get the value of a single point on the
// mesh edge
////////////////////////////////////////////////////////////////////////
void slug_mesh2d_interpolator_vec::
operator()(const double pos, const mesh2d_edge_type edge,
	   array1d& f_interp) const {

  // Action depends on edge
  switch (edge) {
    
  case mesh2d_xlo: // Fall through
  case mesh2d_xhi: {

    // Safety assertions
    double y = pos;
    assert(y >= grid.y_min());
    assert(y <= grid.y_max());
    
    // Get the s coordinate along this edge
    size_type idx;
    double x;
    if (edge == mesh2d_xlo) {
      idx = 0;
      x = grid.x_min(y);
    } else {
      idx = nx-1;
      x = grid.x_max(y);
    }
    double s;
    if (!monotonic_grid) {
      size_type j = grid.j_index(y);
      s = grid.s_grid()[idx][j] +
        sqrt((x - grid.x_grid()[idx][j]) * (x - grid.x_grid()[idx][j]) +
	     (y - grid.y_grid()[j]) * (y - grid.y_grid()[j]));
    } else {
      s = y;
    }

    // Interpolate to point
    for (size_type n=0; n<nf; n++)
      f_interp[n] = gsl_spline_eval(spl_s[idx][n], s, acc_s[idx][n]);
    break;
  }

  case mesh2d_ylo: // Fall through
  case mesh2d_yhi: {

    // Safety assertions
    double x = pos;
    size_type idx;
    if (edge == mesh2d_ylo) {
      idx = 0;
    } else {
      idx = ny-1;
    }
    assert(x >= grid.x_grid()[0][idx]);
    assert(x <= grid.x_grid()[nx-1][idx]);

    // Interpolate to point
    for (size_type n=0; n<nf; n++)
      f_interp[n] = gsl_spline_eval(spl_x[idx][n], x, acc_x[idx][n]);
    break;
  }
  };
}


////////////////////////////////////////////////////////////////////////
// Method to interpolate to get the value of a single point on the
// mesh edge for a single quantity
////////////////////////////////////////////////////////////////////////
double slug_mesh2d_interpolator_vec::
operator()(const double pos, const mesh2d_edge_type edge,
	   size_type f_idx) const {
  
  // Action depends on edge
  double f_interp;
  switch (edge) {
    
  case mesh2d_xlo: // Fall through
  case mesh2d_xhi: {

    // Safety assertions
    double y = pos;
    assert(y >= grid.y_min());
    assert(y <= grid.y_max());
    
    // Get the s coordinate along this edge
    size_type idx;
    double x;
    if (edge == mesh2d_xlo) {
      idx = 0;
      x = grid.x_min(y);
    } else {
      idx = nx-1;
      x = grid.x_max(y);
    }
    double s;
    if (!monotonic_grid) {
      size_type j = grid.j_index(y);
      s = grid.s_grid()[idx][j] +
	sqrt((x - grid.x_grid()[idx][j]) * (x - grid.x_grid()[idx][j]) +
	     (y - grid.y_grid()[j]) * (y - grid.y_grid()[j]));
    } else {
      s = y;
    }

    // Interpolate to point
    f_interp = gsl_spline_eval(spl_s[idx][f_idx], s, acc_s[idx][f_idx]);
    break;
  }

  case mesh2d_ylo: // Fall through
  case mesh2d_yhi: {

    // Safety assertions
    double x = pos;
    size_type idx;
    if (edge == mesh2d_ylo) {
      idx = 0;
    } else {
      idx = ny-1;
    }
    assert(x >= grid.x_grid()[0][idx]);
    assert(x <= grid.x_grid()[nx-1][idx]);

    // Interpolate to point
    f_interp = gsl_spline_eval(spl_x[idx][f_idx], x, acc_x[idx][f_idx]);
    break;
  }
  };

  // Return
  return f_interp;
}

////////////////////////////////////////////////////////////////////////
// Methods to build interpolators on the mesh
////////////////////////////////////////////////////////////////////////

// Build an interpolator at constant x
void slug_mesh2d_interpolator_vec::
build_interp_const_x(const double x,
		     spl_arr_2d& spl,
		     acc_arr_2d& acc,
		     array1d& y_interp_lim,
		     const vector<double>& y_lim) const {

  // Safety assertion
  assert(x >= x_min() && x <= x_max());

  // Generate the intersection list
  vector<double> y, pos;
  vector<size_type> idx, edge, seq;
  grid.intercept_const_x(x, y, pos, idx, edge, seq, y_lim);

  // Safety check
  for (vector<double>::size_type i=0; i<y.size()-1; i++)
    assert(y[i] < y[i+1]);

  // Allocate output holders
  spl.resize(boost::extents[seq.size()/2][nf]);
  acc.resize(boost::extents[seq.size()/2][nf]);
  y_interp_lim.resize(boost::extents[seq.size()]);
  
  // Turn off GSL error handling so that we can handle this on our own
  // and print out more informative error messages
  gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

  // Loop over copmonents of f
  for (size_type n=0; n<nf; n++) {

    // Use the interpolators along the spines to generate the function
    // value at all intersection points
    vector<double> f_tmp(pos.size());
    for (vector<double>::size_type i=0; i<pos.size(); i++) {
      int gsl_errstat;
      gsl_spline *spl;
      gsl_interp_accel *acc;
      if (edge[i] == 0) {
        spl = spl_s[idx[i]][n];
        acc = acc_s[idx[i]][n];
	if (monotonic_grid)
	  pos[i] = y[i]; // If using y instead of s as coordinate
			 // along spines
      } else {
        spl = spl_x[idx[i]][n];
        acc = acc_x[idx[i]][n];
      }
      gsl_errstat = gsl_spline_eval_e(spl, pos[i], acc, f_tmp.data()+i);
      if (gsl_errstat
	  || fabs(f_tmp.data()[i]) > 1e100
	  ) {
        stringstream ss;
        ss << "slug_mesh2d_interpolator_vec::build_interp_const_x: "
           << "bad interpolation evaluation at x = "
           << setprecision(20) << x 
           << ", pos = " << pos[i]
           << ", y = " << y[i]
           << ", interpolation direction = " << edge[i]
           << "; spline range limits are "
           << spl->x[0] << " to "
           << spl->x[spl->size-1]
           << "; input x is "
           << x << ", full list of (y, pos, edge) values is:" << endl;
        for (vector<int>::size_type j=0; j<y.size(); j++) 
           ss << y[j] << "   " << pos[j] << "   " << edge[i] << endl;
        ss << "gsl says: "
           << gsl_strerror(gsl_errstat);
	throw runtime_error(ss.str());
      }
    }

    // Now build new interpolators; we need one for each distinct segment
    for (vector<double>::size_type i=0; i<seq.size(); i+=2) {
      // Set interpolator type for this quantity
      const gsl_interp_type *itype;
      if (interp_type.size() == 0) itype = slug_default_interpolator;
      else itype = interp_type[n];
      unsigned int interp_npt = gsl_interp_type_min_size(itype);
      if (seq[i+1]-seq[i] < interp_npt) itype = gsl_interp_linear;
      spl[i/2][n] = gsl_spline_alloc(itype, seq[i+1]-seq[i]);
      acc[i/2][n] = gsl_interp_accel_alloc();
      gsl_spline_init(spl[i/2][n], y.data()+seq[i], f_tmp.data()+seq[i],
		      seq[i+1]-seq[i]);
    }
  }

  // Restore the GSL error handler
  gsl_set_error_handler(gsl_err);
  
  // Set interpolator limits
  for (vector<double>::size_type i=0; i<seq.size(); i+=2) {
    y_interp_lim[i] = y[seq[i]];
    y_interp_lim[i+1] = y[seq[i+1]-1];
  }
}

  
// Build an interpolator at constant y
void slug_mesh2d_interpolator_vec::
build_interp_const_y(const double y,
		     spl_arr_1d& spl,
		     acc_arr_1d& acc,
		     array1d& x_interp_lim,
		     const vector<double>& x_lim) const {

  // Safety assertion
  assert(y >= y_min() && y <= y_max());

  // Generate the intersection list
  vector<double> x, pos;
  vector<size_type> idx;
  grid.intercept_const_y(y, x, pos, idx, x_lim);

  // Safety check
  for (vector<double>::size_type i=0; i<x.size()-1; i++)
    assert(x[i] < x[i+1]);
  
  // Turn off GSL error handling so that we can handle this on our own
  // and print out more informative error messages
  gsl_error_handler_t *gsl_err = gsl_set_error_handler_off();

  // Loop over components of f
  for (size_type n=0; n<nf; n++) {

    // Use the interpolators along the spines to generate the function
    // value at all intersection points
    vector<double> f_tmp(pos.size());
    for (vector<double>::size_type i=0; i<pos.size(); i++) {
      int gsl_errstat;
      double f_eval;
      if (!monotonic_grid)
	gsl_errstat =
	  gsl_spline_eval_e(spl_s[idx[i]][n], pos[i], acc_s[idx[i]][n],
			    &f_eval);
      else
	gsl_errstat =
	  gsl_spline_eval_e(spl_s[idx[i]][n], y, acc_s[idx[i]][n],
			    &f_eval);
      f_tmp[i] = f_eval;
      if (gsl_errstat) {
	stringstream ss;
	ss << "slug_mesh2d_interpolator_vec::build_interp_const_y: "
	   << "bad interpolation evaluation at y = "
	   << setprecision(20) << y << ", x = "
	   << setprecision(20) << x[i] << ", pos = "
	   << setprecision(20) << pos[i]
	   << "; gsl says: "
	   << gsl_strerror(gsl_errstat);
	throw runtime_error(ss.str());
      }
    }

    // Now build new interpolator
    const gsl_interp_type *itype;
    if (interp_type.size() == 0) itype = slug_default_interpolator;
    else itype = interp_type[n];
    unsigned int interp_npt = gsl_interp_type_min_size(itype);
    if (x.size() >= interp_npt)
      spl[n] = gsl_spline_alloc(itype, x.size());
    else
      spl[n] = gsl_spline_alloc(gsl_interp_linear, x.size());
    acc[n] = gsl_interp_accel_alloc();
    gsl_spline_init(spl[n], x.data(), f_tmp.data(), x.size());
  }

  // Restore the GSL error handler
  gsl_set_error_handler(gsl_err);

  // Set interpolator limits
  x_interp_lim.resize(boost::extents[2]);
  x_interp_lim[0] = x[0];
  x_interp_lim[1] = x.back();
}
