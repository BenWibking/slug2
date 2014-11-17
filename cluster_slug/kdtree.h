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

/*********************************************************************/
/* This file contains routines to build and search KD trees          */
/*********************************************************************/

#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <stdlib.h>
#include "geometry.h"

/*********************************************************************/
/* Macros for traversing KD trees                                    */
/*********************************************************************/
#define ROOT 1
#define LEFT(i) (i<<1)
#define RIGHT(i) ((i<<1)+1)
#define PARENT(i) (i>>1)
#define SIBLING(i) ((i&1)?i-1:i+1)
#define SETNEXT(i) { while (i&1) i=i>>1; ++i; }

/*********************************************************************/
/* Structures for a node and a tree                                  */
/*********************************************************************/
typedef struct {
  /* Number of points in this node */
  unsigned int npt;
  /* Corners of box bounding region occupied by this node. */
  double *xlim[2];
  /* Corners of box tightly bounding the data points in this node */
  double *xbnd[2];
  /* Dimension along which this node splits. -1 for a leaf. */
  int splitdim;
  /* Pointer to the positions for this node */
  double *x;
  /* Extra data associated to the positions */
  void *dptr;
} KDnode;

typedef struct {
  /* Number of dimensions in the tree */
  unsigned int ndim;
  /* Number of levels in the tree */
  unsigned int levels;
  /* Number of leaves in the tree */
  unsigned int leaves;
  /* Number of nodes in the tree */
  unsigned int nodes;
  /* Size of extra data elements associated with positions */
  size_t dsize;
  /* Pointer to tree root */
  KDnode *tree;
} KDtree;


/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

KDtree* build_tree(double *x, unsigned int ndim, unsigned int npt, 
		   unsigned int leafsize, void *dptr, size_t dsize);
/* This routine builds a KD tree from the input data.

   Parameters:
      INPUT/OUTPUT x
         array of npt * ndim elements containing the positions,
         ordered so that element x[j + i*ndim] is the jth coordinate
         of point i; on return, this will be sorted into a tree
      INPUT ndim
         number of dimensions in the data set
      INPUT npt
         number of data points
      INPUT leafsize
         number of points to place in a leaf of the tree
      INPUT/OUTPUT dptr
         pointer to extra data associated to each position
      INPUT dsize
         size of each element of dptr

   Returns:
      OUTPUT tree
         a pointer to a KD tree decomposition of the data
*/

void free_tree(KDtree *tree);
/* Frees the memory associated with a KD tree.

   Parameters
      INPUT/OUTPUT tree
         The KDtree to be de-allocated

   Returns
      Nothing
*/

void neighbors(const KDtree *tree, const double *xpt, 
	       const unsigned int *dims, const unsigned int ndim, 
	       const unsigned int nneighbor,
	       const double *scale, double *pos,
	       void *dptr, double *d2);
/* Routine to find the N nearest neighbors to an input point; input
   points can have fewer dimensions than the search space, in which
   case the routine searches for the nearest neighbors to a line,
   plane, or higher-dimensional object.

   Parameters:
      INPUT tree
         the KD tree to be searched
      INPUT xpt
         array of ndim elements giving the position of the search
         point
      INPUT dims
         array specifying which dimensions in the kd tree are given in
         xpt; e.g., if the KD tree has tree dimensions, and dims = [0,
         2], then the two elements of x will be interpreted as
         specifying x and z coordinates, and the routine will search
         for the closest neighbors to a line as the specified x and
         z. If ndim is equal to the number of dimensions in the KD
         tree, this argument may be left as NULL.
      INPUT ndim
         number of elements in x and dims
      INPUT nneighbor
         the number of neighbors to find
      INPUT scale
         Array of ndim elements giving the scale factors for the
         Euclidean metric; if left as NULL, all scale factors are set to 1
      OUTPUT pos
         positions of the nearest neighbor points; on entry, this
         pointer must point to a block of at least
         tree->ndim*nneighbor elements, and on return element
	 x[i*tree->ndim+j] contains the jth coordinate for the ith
         neighbor found; points are sorted by distance from xpt
      OUTPUT dptr
         extra data associated with each of the nearest neighbor
         points in x; element dptr[i] is the extra data for the ith
         point; must point to a block of valid memory at least i
         elements long
      OUTPUT d2
         squared distances of all particles found from xpt; on entry,
         this pointer mut point to a block of nneighbor elements, and
         on return dist2[i] gives the distance from the ith point
         found to xpt
*/

unsigned int query_box(const KDtree *tree, const double *xbox[2], 
		       unsigned int ndim, const unsigned int *dim, 
		       const double *scale, double **x, void **dptr,
		       double **d2);
/* This routine searches the KD tree and finds all the points that lie
   within a specified box; the box may have fewer dimensions than the
   KD tree, in which case it is interpreted as a hyperslab.

   Parameters:
      INPUT tree
         The KD tree to be searched
      INPUT xbox
         2D array of 2 x ndim elements giving positions of the
         corners of the box
      INPUT ndim
         Number of coordinates in xbox; if ndim == tree->ndim, the
         search region is a closed box, while if ndim < tree->ndim it
         is hyperslab
      INPUT dim
         Array of ndim elements specifying which dimensions of the
         space are specified by the entries in xbox; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1]
      INPUT scale
         Array of tree->ndim elements giving the scale factors for the
         Euclidean metric; if left as NULL, all scale factors are set
         to 1
      OUTPUT x
         pointer to an array containing the locations of all the
         points found; on output (*x)[tree->ndim*i + j] is the jth
         coordinate of point i; (*x) should not point to valid memory
         on input, as the required memory is allocated within the
         routine.
      OUTPUT dptr
         pointer to an array of extra data associated with the points
         found; on output ((char *) (*dptr)) + i*tree->dsize is a
         pointer to the start of the data for point i; (*dptr) should
         not point to valid memory on input, as the required memory is
         allocated within the routine
      OUTPUT d2
         pointer to an array containing the square of the distance
         from each point returned by the query to the center of the
         hyperslab or the box; (*d2) should not point to valid
         memory on input, as the required memory is allocated within
         the routine

   Returns:
      OUTPUT npt
         the number of points returned by the search
*/


unsigned int query_sphere(const KDtree *tree, const double *xcen,
			  unsigned int ndim, const unsigned int *dim, 
			  const double radius, const double *scale,
			  double **x, void **dptr, double **d2);
/* This routine searches the KD tree and finds all the points that lie
   within a specified sphere; the sphere may have fewer dimensions than the
   KD tree, in which case it is interpreted as a hypercylinder that is
   infinite in the unspecified directions.

   Parameters:
      INPUT tree
         The KD tree to be searched
      INPUT xcen
         array of ndim elements giving the position of the center of
	 the sphere
      INPUT ndim
         Number of coordinates in xcen; if ndim < tree->ndim, the
         search region is a hypercylinder, while if ndim == tree->ndim
         the search region is a sphere
      INPUT dim
         Array of ndim elements specifying which dimensions of the
         space are specified by the entries in xbox; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1]
      INPUT radius
         radius of the sphere
      INPUT scale
         Array of tree->ndim elements giving the scale factors for the
         Euclidean metric; if left as NULL, all scale factors are set
         to 1
      OUTPUT x
         pointer to an array containing the locations of all the
         points found; on output (*x)[tree->ndim*i + j] is the jth
         coordinate of point i; (*x) should not point to valid memory
         on input, as the required memory is allocated within the
         routine.
      OUTPUT dptr
         pointer to an array of extra data associated with the points
         found; on output ((char *) (*dptr)) + i*tree->dsize is a
         pointer to the start of the data for point i; (*dptr) should
         not point to valid memory on input, as the required memory is
         allocated within the routine
      OUTPUT d2
         pointer to an array containing the square of the distance
         from each point returned by the query to the center of the
         sphere; (*d2) should not point to valid memory on input, as
         the required memory is allocated within the routine

   Returns:
      OUTPUT npt
         the number of points returned by the search
*/

#endif
/* _KDTREE_H_ */
