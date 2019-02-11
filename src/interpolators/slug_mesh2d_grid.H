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

////////////////////////////////////////////////////////////////////////
//
// The classes in this module and a few others solve a generic
// mathematical problem. We have a function whose value is specified
// at a set of data points that lie on a semi-regular grid of the form
//
// (x_{ij}, y_j), i = 0 ... N-1, j = 0 ... M-1
//
// The points are non-decreasing, so x_{i,j} <= x_{i+1,j} and y_j <=
// y_{j+1}, but need not be strictly increasing. There also need not
// be any particular ordering of x_{i,j} and x_{i,j+1}, which implies
// that the outer edge of the mesh need not be convex, and that edges
// can be degenerate, i.e., (x_{ij}, y_j) = (x_{i+1,j}, y_j) =
// (x_{i,j+1}, y_{j+1}), so the ij - i+1,j and ij - i,j+1 edges are
// identical.
//
// Given a set of function values f_{ij} defined at the grid points,
// we wish to define an interpolating function f(x,y) that will give
// us the value of f at any point in the interior of the semi-regular
// grid. We wish to extract this in two forms. First, we want to be
// able to specify a value (x,y) that lies in the interior of the grid
// and get back a value of f. Second, we want to be able to specify a
// value x and get back a function f_x(y) that interpolates in the y
// direction at fixed x, or specify a value of y and get back a
// function f_y(x) that interpolates as fixed y.
//
// In the context of stellar evolution, the grid consists of a set of
// evolutionary tracks. Each track is defined by an initial mass y_j =
// log m_j, and by a series of evolutionary points, each characterized
// by a time x_ij = log t_ij. The times t_ij are not the same for
// different initial masses, so in general x_{ij} != x_{i+1,j}. For
// each point in the evolutionary tracks, we have one or more stellar
// properties f_{ij}, and we want to be able to interpolate these
// properties, either by specifying a mass and an age directly, or by
// specifying an age and getting back an isochrone -- a function that
// returns stellar properties as a function of mass at fixed age.
//
// To solve this problem, we note that the points (x_ij, y_j) define a
// series of quadrilaterals in the (x, y) plane. These quadrilaterals
// have parallel edges in the x direction, but non-parallel edges in
// the y direction. Graphically, the structure is something like this:
//
//  (x_{i,j+1}, y_{j+1})   (x_{i+1}, y_{i+1,j+1})
//          o----------------------o
//         /                        \
//        /                          \
//       /                            \
//      o------------------------------o
//  (x_{ij}, y_j)              (x_{i+1,j}, y_j})
//
// For convenience we will refer to the edges at constant y as
// horizontal edges, and the edges that are not at constant y as
// vertical edges, though they in fact are generally not perfectly
// vertical. Let index (i,j) refer to the vertical edge that lies
// above point (x_{ij}, y_j). The slope of this edge is
//
// m_{ij} = (y_{j+1} - y_j) / (x_{i,j+1} - x_{ij}).
//
// We can also define a distance coordinate along the vertical cell
// edges. There are two possible ways to parmaeterize the coordinates
// in this direction. One is just by the y coordinate, but another
// choice is to define the distances to the grid points along the
// vertical cell edges recursively by
//
// s_{i0} = 0
// s_{i,j+1} = s_{i,j} + sqrt[ (x_{i,j+1} - x_{ij})^2 +
//                             (y_{j+1} - y_j)^2 ].
//
// This can be generalized naturally to give a coordinate s for an
// arbitrary point (x, y) that lies on one of the vertical edges:
//
// s(x,y) = s_{i,jmax} + sqrt[ (x - x_{i,jmax})^2 +
//                             (y - y_{jmax})^2 ]
//
// where jmax = max(j) s.t. y > y_j and the condition that (x, y) be
// on the vertical cell edge is equivalent to requiring that
//
// y - y_{jmax} = m_{i,jmax} (x - x_{i,jmax})
//
// for some i.
//
// Next, we define spline fits along both the horizontal and vertical
// edges. Specifically, consider the set of points along a horizontal
// edge
//
// p_j = (x_{0j}, f_{0j}), (x_{1,j}, f_{1,j}), ...
//       (x_{N-1,j}, f_{N-1,j})
//
// These points p_j can be used to define a spline function P_j(x)
// that accepts as an argument any x in [ x_{0j}, x_{N-1,j} ], and
// returns a function value at x. We refer to this spline as a "spine"
// along the x direction. The order and form of the spline is
// arbitrary, but we will use a true spline, i.e., one that is
// constrained to go through each data point, rather than a basis
// spline.
//
// We can similarly define the set of points in along the vertical
// edges. If we use the variable s defined above as our coordinate
// along the spines, this is:
//
// q_i = (s_{i0}, f_{i0}), (s_{i1}, f_{i1}), ...
//       (s_{i,M-1}, f_{i,M-1})
//
// and their corresponding "spine" spline functions Q_i(s(y)), where
//
// s(y) = s(x(y), y)
// x(y) = x_{i,jmax} + (y - y_{jmax}) / m_{i,jmax}
//
// Alternatively, if we use the y position as our coordinate across
// the spines the results are the same, except that we simply set
// s = y everywhere.
//
// With these definitions in place, we can now define the procedures
// for interpolating to a single point or alone a line of fixed
// x. First consider generating an interpolating function along a line
// of fixed x. We find all locations where a line of fixed x
// intersects the mesh; these locations can be alone either horizontal
// or vertical mesh edges. Let y_k be the set of intersection points,
// ordered so that y_{k+1} > y_k. For each intersection point, we
// define a corresponding interpolated function value f_k by
// evaluating either P(x) (if the point corresponds to intersecting a
// horizontal cell edge) or Q(s(y)) (if it corresponds to intersecting
// a vertical cell edge). Given the set (y_k, f_k), we can define a
// new interpolating function F(y) that allows us to evaluate the
// function at any y.
//
// Again, depending on our choice of coordinate along the spines, we
// can also just set s = y in the above. The reason for using s = y
// as an option is that this option, together with an overall linear
// interpolation scheme, makes it possible to enforce strict
// monotonicity, in the sense that if some quantity is monotonic in
// x along all tracks, and we wish all interpolants to be similarly
// monotonic, then we can guarantee that this will be so by setting
// s = y and using linear interpolation. Any other choice does not
// strictly guarantee monotonic behavior. However, the price of this
// choice is an overall reduction in the level of accuracy in the
// interpolation.
//
// The case of interpolating to a single point is analogous, except
// that we need not find all the intersection points y_k. We just find
// as many as needed on each side to have enough for the type of
// interpolating function we wish to use.
//
// Finally, the case of creating a function that interpolates in x at
// fixed y is the same as well, except that the intersection points
// are all vertical.
//
// The code to do all of this is broken up into two classes:
//
// slug_mesh2d_grid -- this handles all the geometry of defining a
//    grid, finding intersections of lines through it, etc. It does
//    not know anything about interpolation, and just handles
//    geometry.
// slug_mesh2d_interpolator -- this class wraps slug_mesh2d_grid, and
//    handles all the interpolation machinery.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_mesh2d_grid_H_
#define _slug_mesh2d_grid_H_

#include <vector>
#include "../constants.H"
#include "slug_interpolation.H"

class slug_mesh2d_grid {
  
public:
  
  // Constructor; takes arrays of x and y points. The x array can be
  // either a true 2d array or a 2d view. The class makes a copy of
  // the input data, so the data may be modified freely after being
  // passed.
  slug_mesh2d_grid(const array2d& x_, const array1d& y_);

  // Alternate constructor that takes a view rather than an array for x
  slug_mesh2d_grid(const view2d& x_, const array1d& y_);

  // Destructor
  ~slug_mesh2d_grid() { };
  
  // Methods to return the grid data
  const array2d& x_grid() const { return x; }
  const array1d& y_grid() const { return y; }
  const array2d& m_grid() const { return m; }
  const array2d& s_grid() const { return s; }

  // Method to return the shape
  const boost::multi_array_types::size_type* shape() const
  { return x.shape(); }

  // Methods to query the locations of the edge of the mesh; default
  // arguments are sentinel values, and for these values the minima /
  // maxima returned are the global min / max over the entire grid. If
  // the argument is not this sentinel value, it is an error if the
  // argument given misses the mesh entirely. Finally, note that y_min
  // and y_max return the global min and max only, while y_lim can
  // return an arbitrary number of intersections where the mesh begins
  // and ends.
  double x_min(const double y_ = constants::big) const;
  double x_max(const double y_ = constants::big) const;
  std::vector<double> x_lim(const double y_ = constants::big) const;
  double y_min() const { return ymin; };
  double y_max() const { return ymax; };
  std::vector<double> y_lim(const double x_ = constants::big) const;

  // Method to return the slopes of the mesh edges at a particular y;
  // formally, these functions return the slopes dy/dx on the low and
  // high x edges of the mesh
  double dydx_lo(const double x_) const;
  double dydx_hi(const double x_) const;

  // Method to ask of the mesh is convex
  bool convex() const { return convex_; }

  // Method to query if a point is in the mesh
  bool in_mesh(const double x_, const double y_) const;

  // Method to query if a point is exactly on a spine; return value is
  // true if yes, false otherwise; if the return value is true,
  // pos, idx, and edge are altered to give information about the
  // intersection: pos holds the position along the spine (x or s,
  // depending on whether it is a horizontal or vertical spine), idx
  // gives the index of the spine, and edge is 0 for a vertical spine,
  // 1 for a horizontal spine. If the routine returns false, the input
  // values of these quantities are not altered.
  bool on_spine(const double x_, const double y_,
		double &pos,
		boost::multi_array_types::size_type& idx,
		boost::multi_array_types::size_type& edge) const;
  
  // Methods to index points; note that bounds-checking is disabled if
  // compiled with NDEBUG, and the resuls if the point is outside the
  // mesh are unpredictable, so call in_mesh first if necessary. The
  // ij_index method returns the indices to the lower left of a
  // specified i and j. The j_index method returns the j index below a
  // specified y value. The i_index_j method value returns the i index
  // for a point along a particular j track.
  void ij_index(const double x_, const double y_,
		boost::multi_array_types::size_type &i,
		boost::multi_array_types::size_type &j) const;
  boost::multi_array_types::size_type j_index(const double y_) const;
  boost::multi_array_types::size_type
  i_index_j(const double x_,
	    const boost::multi_array_types::size_type j) const;

  // Routines to return the list of points where lines of constant y or
  // constant x intersect the mesh; in these functions:
  // x_ or y_ is the x or y coordinate of the intersection line
  // int_y or int_x is the y or x coordinate of the intersection point
  // int_pos is the intersection position; this is the x or s position
  //    along the intersection edge
  // int_index is the i or j index of the spine that is intersected
  // int_edge specifies which edge the intersection lies along; 0 =
  //    intersection along an x spine, 1 = intersection along a y
  //    spine
  // int_seq: if the mesh is not convex, a line of constant x may
  //    enter and exist the mesh multiple times; int_seq labels these
  //    entries and exits; even entries in int_seq indicate points in
  //    int_pos that are points of entry to the tracks, while odd
  //    entries indicate points of exit; for convex tracks, int_seq
  //    will therefore always have 2 elements, equal to 0 and to the
  //    number of elements in int_pos (unless the specified x or y
  //    miss the tracks completely, in which case it will be empty)
  // xlim / ylim is an optional 2-element array specifying that
  //    intersection should only be found over a limited range of x
  //    and y
  void intercept_const_x(const double x_,
			 std::vector<double>& int_y,
			 std::vector<double>& int_pos,
			 std::vector<boost::multi_array_types::size_type>&
			 int_index,
			 std::vector<boost::multi_array_types::size_type>&
			 int_edge,
			 std::vector<boost::multi_array_types::size_type>&
			 int_seq,
			 const std::vector<double>& ylim =
			 std::vector<double>()) const;
  void intercept_const_y(const double y_,
			 std::vector<double>& int_x,
			 std::vector<double>& int_pos,
			 std::vector<boost::multi_array_types::size_type>&
			 int_index,
			 const std::vector<double>& xlim =
			 std::vector<double>()) const;

  // This routine is like intercept_const_x, in that it returns points
  // where a line of constant x intersects the mesh. However, instead
  // of working on all y or over a specified range in y, it returns a
  // certain number of intercept points on either side of some
  // particular, specified y. The intent is to return enough points to
  // enable interpolation, without needing to traverse the entire
  // mesh. Arguments are:
  // x_ : constant x value
  // y_ : y value around which to search for intercepts
  // npt : total number of intercept points to return
  // int_pos, int_index, int_edge: same as for intercept_const_x
  // Note that the routine will attempt to find the n points uniformly
  // distributed around the input y (with one extra point below y if
  // npt is odd), but it will not go past a mesh edge. If it cannot
  // find enough points in one direction, it will add more points in
  // the other direction until it reaches the requested number, or
  // hits the mesh edge in both directions; if that happens, it will
  // return as many points as possible, which may be fewer than npt.
  void intercept_const_x_n(const double x_,
			   const double y_,
			   const boost::multi_array_types::size_type npt,
			   std::vector<double>& int_y,
			   std::vector<double>& int_pos,
			   std::vector<boost::multi_array_types::size_type>&
			   int_index,
			   std::vector<boost::multi_array_types::size_type>&
			   int_edge) const;

private:

  // Binary search routines; these are just lightly-modified versions
  // of the gsl_interp_bsearch routine from the gsl. The main
  // difference is that these are set up so that, if the search point
  // lies exactly on an index point, the routines will always return
  // an index pointer so that the input point lies on the left /
  // bottom side of the indexing interval. The routines do the
  // following:
  // i_bsearch_j: returns the i index that goes with a given x along
  //    a particular j track
  // i_bsearch_j2: returns the i index that goes with a given x for a
  //    certain offset between two j tracks
  // j_bsearch: returns the j index that goes with a given y
  boost::multi_array_types::size_type
  i_bsearch_j(const double x_,
	      const boost::multi_array_types::size_type j,
	      const boost::multi_array_types::size_type idx_lo,
	      const boost::multi_array_types::size_type idx_hi) const;
  boost::multi_array_types::size_type
  i_bsearch_j2(const double x_, const double dy,
	       const boost::multi_array_types::size_type j,
	       const boost::multi_array_types::size_type idx_lo,
	       const boost::multi_array_types::size_type idx_hi) const;
  boost::multi_array_types::size_type
  j_bsearch(const double y_,
	    const boost::multi_array_types::size_type idx_lo,
	    const boost::multi_array_types::size_type idx_hi) const;

  // Worker routines used by intercept_const_x
  void intercept_const_x_seg(const double x_,
			     const double y_start,
			     const double y_end,
			     const bool start_interior,
			     const bool end_interior,
			     std::vector<double>& int_y,
			     std::vector<double>& int_pos,
			     std::vector<boost::multi_array_types::size_type>&
			     int_index,
			     std::vector<boost::multi_array_types::size_type>&
    int_edge) const;
  bool find_next_y_intersect(const double x_,
			     double& y_,
			     const double y_end,
			     const bool end_interior,
			     const bool search_up,
			     bool& last_intersect_left,
			     bool& last_intersect_right,
			     std::vector<double>& int_y,
			     std::vector<double>& int_pos,
			     std::vector<boost::multi_array_types::size_type>&
			     int_index,
			     std::vector<boost::multi_array_types::size_type>&
			     int_edge) const;

  // Mesh data
  array2d x;
  array1d y;
  double xmin, xmax, ymin, ymax;
  const boost::multi_array_types::size_type nx, ny;

  // Data describing the mesh; m = slope, s = distances along vertical
  // spines of mesh, convex_ = flag indicating if mesh is convex.
  array2d m, s;
  bool convex_;

  // Cached indices to speed up searches
  mutable boost::multi_array_types::size_type i_cache, j_cache;

};

#endif
// slug_mesh2d_grid_H_
