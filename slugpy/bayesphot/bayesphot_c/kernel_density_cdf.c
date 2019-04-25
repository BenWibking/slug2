/*********************************************************************
Copyright (C) 2019 Mark Krumholz
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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "geometry.h"
#include "kernel_density_cdf.h"

/*********************************************************************/
/* Useful macros                                                     */
/*********************************************************************/

#define SQR(x) ((x)*(x))
#define NODEBLOCKSIZE 16384

/*********************************************************************/
/* Static functions                                                  */
/*********************************************************************/

/* Functions to compute PDFs and integrals thereof on single tree
   nodes. We define these as static inline's for speed, since these
   are doing the bulk of our computation. */
static inline
void kd_cdf_node(const kernel_density *kd, const double *x, 
		 const unsigned long curnode, double *cdf,
		 double *cdferr);
static inline
void kd_cdf_node_int(const kernel_density *kd, const double *x,
		     const unsigned long *dims, const unsigned long ndim,
		     const unsigned long ndim_int, const double fac,
		     const unsigned long curnode, double *cdf,
		     double *cdferr);

/**********************************************************************/
/* Function to evaluate the PDF approximately using a kernel_density */
/* object                                                            */
/*********************************************************************/
double kd_cdf(const kernel_density *kd, const double *x,
	      const double reltol, const double abstol) {

  /* Make sure we have a Gaussian kernel, and return an impossible
     value as a flag if we don't */
  if (kd->ktype != gaussian) return -1.0;

  /* Analyze root node */
  double cdf, abserr;
  kd_cdf_node(kd, x, ROOT, &cdf, &abserr);

  /* Deal with special case where root node is a leaf */
  if (kd->tree->tree[ROOT].splitdim == -1) return(cdf);

  /* Allocate memory for node list */
  unsigned long *nodelist;
  if (!(nodelist = (unsigned long *) 
	calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf\n");
    exit(1);
  }
  double *nodecdf, *nodeerr;
  if (!(nodecdf = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf\n");
    exit(1);
  }
  if (!(nodeerr = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf\n");
    exit(1);
  }
  unsigned long nalloc = NODEBLOCKSIZE;  

  /* Push root node onto the node list */
  nodelist[0] = ROOT;
  nodecdf[0] = cdf;
  nodeerr[0] = abserr;
  unsigned long nnode = 1;
  double leafcdf = 0.0;

  /* Now work recursively through the tree, identifying the node in
     the tree that contributes most to the error and opening it until
     the error estimate is within our tolerance */
  while (1) {

    /* Compute the current estimate and total error, and find the node
       that is contributing the most to the error budget. Note that in
       principle we could do fewer arithmetic operations by keeping a
       running tally of abserr and cdf, and updating these each time
       we go through the loop, rather than recomputing them from
       scratch each time. However, in practice this involves a lot of
       cancelling additions and subtractions that can create problems
       with roundoff when the error tolerance is tight and a lot of
       node evaluations are needed to achieve it. Doing the
       calculation fresh each time through the loop avoids this
       problem, which is worth the small amount of extra work required. */
    cdf = leafcdf + nodecdf[0];
    unsigned long ptr = 0;
    double maxerr = nodeerr[0];
    abserr = nodeerr[0];
    for (unsigned long i=1; i<nnode; i++) {
      cdf += nodecdf[i];
      abserr += nodeerr[i];
      if (nodeerr[i] > maxerr) {
	ptr = i;
	maxerr = nodeerr[i];
      }
    }

    /* Check for convergence */
    double relerr = abserr / cdf;
    if ((abserr <= abstol) || (relerr <= reltol)) break;
    
    /* Compute estimates for that node's children */
    unsigned long lchild = LEFT(nodelist[ptr]);
    unsigned long rchild = RIGHT(nodelist[ptr]);
    double lcdf, lerr, rcdf, rerr;
    kd_cdf_node(kd, x, lchild, &lcdf, &lerr);
    kd_cdf_node(kd, x, rchild, &rcdf, &rerr);

    /* If children are leaves, add their contribution to the CDF from
       leaves */
    if (kd->tree->tree[lchild].splitdim == -1) leafcdf += lcdf;
    if (kd->tree->tree[rchild].splitdim == -1) leafcdf += rcdf;

    /* Remove the node we just analyzed from the list */
    for (unsigned long i=ptr; i<nnode-1; i++) {
      nodelist[i] = nodelist[i+1];
      nodecdf[i] = nodecdf[i+1];
      nodeerr[i] = nodeerr[i+1];
    }
    nnode--;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode+2 >= nalloc) {
      if (!(nodelist = (unsigned long *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	fprintf(stderr,
		"bayesphot: error: unable to allocate memory in kd_cdf\n");
	exit(1);
      }
      if (!(nodecdf = (double *) 
	    realloc(nodecdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr,
		"bayesphot: error: unable to allocate memory in kd_cdf\n");
	exit(1);
      }
      if (!(nodeerr = (double *) 
	    realloc(nodeerr, 2*nalloc*sizeof(double)))) {
	fprintf(stderr,
		"bayesphot: error: unable to allocate memory in kd_cdf\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* If children are not leaves, push them onto the node list */
    if (kd->tree->tree[lchild].splitdim != -1) {
      nodelist[nnode] = lchild;
      nodecdf[nnode] = lcdf;
      nodeerr[nnode] = lerr;
      nnode++;
    }
    if (kd->tree->tree[rchild].splitdim != -1) {
      nodelist[nnode] = rchild;
      nodecdf[nnode] = rcdf;
      nodeerr[nnode] = rerr;
      nnode++;
    }

    /* Safety check: bail out if no nodes left */
    if (nnode == 0) break;
  }

  /* Free memory */
  free(nodelist);
  free(nodecdf);
  free(nodeerr);

  /* Return */
  return(cdf);
}

/*********************************************************************/
/* Function to evaluate the CDF integrated over certain dimensions   */
/* approximately using a kernel_density object                       */
/*********************************************************************/
double kd_cdf_int(const kernel_density *kd, const double *x,
		  const unsigned long *dims, const unsigned long ndim,
		  const double reltol, const double abstol) {
  unsigned long i, ndim_int, nnode, nalloc, ptr, lchild, rchild;
  unsigned long *nodelist = NULL;
  double hprod, ds_n, fac;
  double *nodecdf = NULL, *nodeerr = NULL;
  double maxerr, relerr, abserr;
  double cdf, leafcdf, lcdf, rcdf, lerr, rerr;
  
  /* Make sure we have a Gaussian kernel, and return an impossible
     value as a flag if we don't */
  if (kd->ktype != gaussian) return -1.0;

  /* Pre-compute constant factor in the integrals we're evaluating;
     this is the part that depends only on h and the number of
     dimensions. */
  ndim_int = kd->tree->ndim - ndim;
  ds_n = ds(ndim_int);
  hprod = 1.0;
  for (i=0; i<kd->tree->ndim; i++) hprod *= kd->h[i];
  for (i=0; i<ndim; i++) hprod /= kd->h[dims[i]];
  fac = hprod * ds_n *
    pow(2.0, ndim_int/2.0 - 1) * gsl_sf_gamma(0.5*ndim_int);

  /* Process the root node */
  kd_cdf_node_int(kd, x, dims, ndim, ndim_int, fac, ROOT,
		  &cdf, &abserr);

  /* Initialize node list */
  if (kd->tree->tree[ROOT].splitdim == -1) {

    /* Special case: root node is a leaf, so just return the exact value */
    return(cdf);

  } else {

    /* The usual case: root node is not a leaf */

    /* Allocate memory for node list */
    if (!(nodelist = (unsigned long *) 
	  calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf_int\n");
      exit(1);
    }
    if (!(nodecdf = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf_int\n");
      exit(1);
    }
    if (!(nodeerr = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf_int\n");
      exit(1);
    }
    nalloc = NODEBLOCKSIZE;  

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodecdf[0] = cdf;
    nodeerr[0] = abserr;
    nnode = 1;
    leafcdf = 0.0;
  }

  /* Now work recursively through the tree, identifying the node in
     the tree that contributes most to the error and opening it until
     the error estimate is within our tolerance */
  while (1) {

    /* Compute the current estimate and total error, and find the node
       that makes the largest contribution to the total error */
    abserr = maxerr = nodeerr[0];
    cdf = leafcdf + nodecdf[0];
    ptr = 0;
    for (i=1; i<nnode; i++) {
      cdf += nodecdf[i];
      abserr += nodeerr[i];
      if (nodeerr[i] > maxerr) {
	ptr = i;
	maxerr = nodeerr[i];
      }
    }

    /* Check termination condition */
    relerr = abserr / cdf;
    if ((abserr <= abstol) || (relerr <= reltol)) break;

    /* Analyze the children of the node contributing the most error */
    lchild = LEFT(nodelist[ptr]);
    rchild = RIGHT(nodelist[ptr]);
    kd_cdf_node_int(kd, x, dims, ndim, ndim_int, fac, lchild, &lcdf, &lerr);
    kd_cdf_node_int(kd, x, dims, ndim, ndim_int, fac, rchild, &rcdf, &rerr);

    /* If children are leaves, add their contribution to the CDF from
       leaves */
    if (kd->tree->tree[lchild].splitdim == -1) leafcdf += lcdf;
    if (kd->tree->tree[rchild].splitdim == -1) leafcdf += rcdf;
    
    /* Remove the node we just analyzed from the list */
    for (i=ptr; i<nnode-1; i++) {
      nodelist[i] = nodelist[i+1];
      nodecdf[i] = nodecdf[i+1];
      nodeerr[i] = nodeerr[i+1];
    }
    nnode--;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode+2 >= nalloc) {
      if (!(nodelist = (unsigned long *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf_int\n");
	exit(1);
      }
      if (!(nodecdf = (double *) 
	    realloc(nodecdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf_int\n");
	exit(1);
      }
      if (!(nodeerr = (double *) 
	    realloc(nodeerr, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_cdf_int\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* If children are not leaves, push them onto the node list */
    if (kd->tree->tree[lchild].splitdim != -1) {
      nodelist[nnode] = lchild;
      nodecdf[nnode] = lcdf;
      nodeerr[nnode] = lerr;
      nnode++;
    }
    if (kd->tree->tree[rchild].splitdim != -1) {
      nodelist[nnode] = rchild;
      nodecdf[nnode] = rcdf;
      nodeerr[nnode] = rerr;
      nnode++;
    }
    
    /* Safety check: bail out if no nodes left */
    if (nnode == 0) break;
  }

  /* Free memory */
  free(nodelist);
  free(nodecdf);
  free(nodeerr);

  /* Return */
  return(cdf);
}



/*********************************************************************/
/* Function to estimate the contribution to a CDF from a node, and   */
/* the error on it; if the input node is a leaf, the routine         */
/* computes the exact CDF contribution, and sets the error to zero.  */
/* If it is not a leaf, the routine computes the minimum and maximum */
/* possible contributions from the node, returns their average as    */
/* the central estimate, and half their difference as the maximum    */
/* possible error.                                                   */
/*********************************************************************/
inline
void kd_cdf_node(const kernel_density *kd, const double *x, 
		 const unsigned long curnode, double *cdf,
		 double *cdferr) {

  /* Convenience variables */
  unsigned long ndim = kd->tree->ndim;
  KDnode node = kd->tree->tree[curnode];
  
  /* Is this node a leaf? If so, just sum over it and return that */
  if (node.splitdim == -1) {

    /* Loop over points, summing their contribution */
    *cdf = 0.0;
    for (unsigned long i=0; i<node.npt; i++) {

      /* Compute integral for this point */
      double cdf_pt = 1.0;
      for (unsigned long j=0; j<ndim; j++) {
	double xn = (x[j] - node.x[ndim*i+j]) / kd->h[j];
	cdf_pt *= kd->h[j] * (1.0 + erf(xn/sqrt(2.0)));
      }
      cdf_pt *= pow(M_PI/2.0, 0.5*ndim);

      /* If we have a weight for this point, apply it */
      if (node.dptr != NULL) cdf_pt *= ((double *) node.dptr)[i];

      /* Add to sum */
      *cdf += cdf_pt;
    }

    /* Apply overall normalization */
    *cdf *= kd->norm_tot;

    /* Set error to zero */
    *cdferr = 0.0;

  } else {

    /* Node is not a leaf. We therefore compute the minimum and
       maximum possible contributions, corresponding to assuming that
       all points in the node are in the upper right and lower left
       corners of the bounding box, respectively */
    double cdf_min = 1.0;
    double cdf_max = 1.0;

    /* Loop over dimensions, adding contribution */
    for (unsigned long j=0; j<ndim; j++) {
      double xn = (x[j] - node.xbnd[1][j]) / kd->h[j];
      cdf_min *= kd->h[j] * (1.0 + erf(xn/sqrt(2.0)));
      xn = (x[j] - node.xbnd[0][j]) / kd->h[j];
      cdf_max *= kd->h[j] * (1.0 + erf(xn/sqrt(2.0)));
    }

    /* Apply normalization factors */
    cdf_min *= pow(M_PI/2.0, 0.5*ndim) *
      kd->nodewgt[curnode] * kd->norm_tot;
    cdf_max *= pow(M_PI/2.0, 0.5*ndim) *
      kd->nodewgt[curnode] * kd->norm_tot;

    /* Set estimate to average of min and max, and error to half the
       difference between them. */
    *cdf = (cdf_min + cdf_max) / 2.0;
    *cdferr = (cdf_max - cdf_min) / 2.0;

  }
}


/*********************************************************************/
/* Function to estimate the contribution to an integrated CDF from a */
/* node. Functionality is identical to kdf_cdf_node, except that     */
/* certain dimensions are being integrated out.                      */
/*********************************************************************/
inline
void kd_cdf_node_int(const kernel_density *kd, const double *x,
		     const unsigned long *dims, const unsigned long ndim,
		     const unsigned long ndim_int, const double fac,
		     const unsigned long curnode, double *cdf,
		     double *cdferr) {

  /* Convenience variables */
  unsigned long ndim_tot = kd->tree->ndim;
  KDnode node = kd->tree->tree[curnode];

  /* Is this node a leaf? If so, just sum over it and return that */
  if (node.splitdim == -1) {

    /* Loop over points, summing their contribution */
    *cdf = 0.0;
    for (unsigned long i=0; i<node.npt; i++) {

      /* Compute integral for this point */
      double cdf_pt = 1.0;
      for (unsigned long j=0; j<ndim; j++) {
	double xn = (x[dims[j]] - node.x[ndim_tot*i+dims[j]])
	  / kd->h[dims[j]];
	cdf_pt *= kd->h[dims[j]] * (1.0 + erf(xn/sqrt(2.0)));
      }
      cdf_pt *= pow(M_PI/2.0, 0.5*ndim);

      /* If we have a weight for this point, apply it */
      if (node.dptr != NULL) cdf_pt *= ((double *) node.dptr)[i];

      /* Add to sum */
      *cdf += cdf_pt;
    }

    /* Normalize and return */
    *cdf *= fac * kd->norm_tot;
    *cdferr = 0.0;

  } else {

    /* This node is not a leaf, so compute the minimum and maximum
       possible contribution */
    double cdf_min = 1.0;
    double cdf_max = 1.0;

    /* Loop over dimensions, adding contribution */
    for (unsigned long j=0; j<ndim; j++) {
      double xn = (x[dims[j]] - node.xbnd[1][dims[j]]) / kd->h[dims[j]];
      cdf_min *= kd->h[dims[j]] * (1.0 + erf(xn/sqrt(2.0)));
      xn = (x[dims[j]] - node.xbnd[0][dims[j]]) / kd->h[dims[j]];
      cdf_max *= kd->h[dims[j]] * (1.0 + erf(xn/sqrt(2.0)));
    }

    /* Apply normalization factors */
    cdf_min *= pow(M_PI/2.0, 0.5*ndim) *
      fac * kd->nodewgt[curnode] * kd->norm_tot;
    cdf_max *= pow(M_PI/2.0, 0.5*ndim) *
      fac * kd->nodewgt[curnode] * kd->norm_tot;


    /* Set estimate to average of min and max, and error to half the
       difference between them. */
    *cdf = (cdf_min + cdf_max) / 2.0;
    *cdferr = (cdf_max - cdf_min) / 2.0;
  }
}

