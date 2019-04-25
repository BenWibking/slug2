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

/*********************************************************************/
/* This module contains routines to compute cumulative distributions */
/* from kernel density estimates stored in KD trees. The functions   */
/* in this module are closely analogous to those in kernel_density.h */
/* with the difference that instead of computing PDFs, in this       */
/* module the functions all compute CDFs.                            */
/*                                                                   */
/* IMPORTANT WARNING: for now all these functions ONLY work for      */
/* Gaussian kernels (because computing the required integrals for    */
/* the other kernel shapes in arbitrary dimensions is a pain in the  */
/* ass that I don't want to have to deal with right now. If you call */
/* any of these functions for any other kernel type, they return a   */
/* nonsense value to flag the error.                                 */
/*********************************************************************/

#ifndef _KERNEL_DENSITY_CDF_H_
#define _KERNEL_DENSITY_CDF_H_

#include "kernel_density.h"

/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

double kd_cdf(const kernel_density *kd, const double *x,
	      const double reltol, const double abstol);
/* This routine returns the value of the cumulative distribution
   function for a kernel_density object evaluated at a specified
   position x, i.e., if the PDF is given by 
      p(\vec{x}) = p(x_1, x_2, x_3, ...),
   then the value returned by this function is
      \int_{-\infty}^{x[0]} \int_{-\infty}^{x[1]} ... p(\vec{x}) d\vec{x}.
   This quantity is evaluated with some specified relative and absolute
   tolerances.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT x
         an ndim element array giving the position at which the PDF is
         to be evaluated
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / pdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUTPUT nodecheck
         Number of individual nodes examined during the evaluation;
         only if compiled with DIAGNOSTIC set
      OUTPUT leafcheck
         Number of leaves examined during the evaluation; only if
         compiled with DIAGNOSTIC set
      OUTPUT termcheck
         Number of nodes examined during the evaluation which were not
         futher sub-divided; only if compiled with DIAGNOSTIC set

   Returns:
      OUT cdf_approx
         an approximation to the CDF evaluated at x, satisfying the
         input error tolerances
*/


double kd_cdf_int(const kernel_density *kd, const double *x,
		  const unsigned long *dims, const unsigned long ndim,
		  const double reltol, const double abstol
#ifdef DIAGNOSTIC
		  , unsigned long *nodecheck, unsigned long *leafcheck,
		  unsigned long *termcheck
#endif
		  );
/* This routine returns the value of the cumulative distrbution
   function evaluated at a particular point x in certain dimensions,
   with all other dimensions integrated out (or equivalently, with the
   value of the input point for those dimensions set to infinity. For
   example, if the PDF depends on n variables, and is written out as

   p(x(0), x(1), x(2), ... x(n-1)),

   the input data point is

   x = [0.1, 0.3, 0.5],

   and the input list of dimensions is

   dims = [0, 1, 3],

   then the value returned will be

   \int_{-inf}^{0.1} dx(0) \int_{-inf}^{0.3} dx(1)
      \int_{-inf}^{inf} dx(2) \int_{-inf}^{0.5} dx(3)
      \int_{-inf}^{inf} dx(4) ... \int_{-inf}^{inf} dx(n-1) p(x)

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the CDF
      INPUT x
         an ndim element array giving the position at which the CDF is
         to be evaluated
      INPUT dims
         an ndim element array specifying the dimensions included in x
      INPUT ndim
         number of dimensions in x; must be less than the number in
         the kd tree
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value cdf_approx will be returned once the estimated error
	 | cdf_approx - cdf_true | / cdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value cdf_approx will be returned once the estimated error
	 | cdf_approx - cdf_true | < abstol.

   Returns:
      OUT cdf_approx
         an approximation to the output integral, accurate within the
         specified error tolerances
*/

#endif
/* _KERNEL_DENSITY_CDF_H_ */
