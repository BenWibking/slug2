"""
This defines a class that can be used to estimate the PDF of physical
quantities from a set of input photometry in various bands, together
with a training data set.
"""

import numpy as np
import scipy.interpolate as interp
import os
import os.path as osp
import ctypes
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_int
from ctypes import c_uint
from ctypes import c_ulong
from ctypes import c_double
from ctypes import c_bool
import numpy.ctypeslib as npct
import random
from copy import deepcopy
from warnings import warn
try:
    import emcee
    mc_avail = True
except:
    mc_avail = False
    pass
import multiprocessing as mp


##################################################################
# Define some types for use later                                #
##################################################################
array_1d_double = npct.ndpointer(dtype=c_double, ndim=1,
                                 flags="CONTIGUOUS")
array_1d_ulong = npct.ndpointer(dtype=c_ulong, ndim=1,
                               flags="CONTIGUOUS")
array_1d_int = npct.ndpointer(dtype=c_int, ndim=1,
                              flags="CONTIGUOUS")

##################################################################
# Define the cluster_slug class                                  #
##################################################################

class bp(object):
    """
    A class that can be used to estimate the PDF of the physical
    properties of stellar population from a training set plus a set of
    measured photometric values.

    Properties
       priors : array, shape (N) | callable | None
          prior probability on each data point; interpretation
          depends on the type passed; array, shape (N): values are
          interpreted as the prior probability of each data point;
          callable: the callable must take as an argument an array
          of shape (N, nphys), and return an array of shape (N)
          giving the prior probability at each data point; None:
          all data points have equal prior probability
       bandwidth : 'auto' | float | array, shape (M)
          bandwidth for kernel density estimation; if set to
          'auto', the bandwidth will be estimated automatically; if
          set to a scalar quantity, the same bandwidth is used for all
          dimensions
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, dataset, nphys, filters=None, bandwidth='auto',
                 ktype='gaussian', priors=None, sample_density=None,
                 reltol=1.0e-2, abstol=1.0e-6, leafsize=16,
                 nosort=None, thread_safe=True):
        """
        Initialize a bp object.

        Parameters
           dataset : array, shape (N, M)
              training data set; this is a set of N sample stellar
              populations, having M properties each; the first nphys
              represent physical properties (e.g., log mass, log age),
              while the next M - nphys represent photometric
              properties
           nphys : int
              number of physical properties in dataset
           filters : listlike of strings
              names of photometric filters; not used, but can be
              stored for convenience
           bandwidth : 'auto' | float | array, shape (M)
              bandwidth for kernel density estimation; if set to
              'auto', the bandwidth will be estimated automatically; if
              set to a scalar quantity, the same bandwidth is used for all
              dimensions
           ktype : string
              type of kernel to be used in densty estimation; allowed
              values are 'gaussian' (default), 'epanechnikov', and
              'tophat'; only Gaussian can be used with error bars
           priors : array, shape (N) | callable | None
              prior probability on each data point; interpretation
              depends on the type passed; array, shape (N): values are
              interpreted as the prior probability of each data point;
              callable: the callable must take as an argument an array
              of shape (N, nphys), and return an array of shape (N)
              giving the prior probability at each data point; None:
              all data points have equal prior probability
           sample_density : array, shape (N) | callable | 'auto' | None
              the density of the data samples at each data point; this
              need not match the prior density; interpretation depends
              on the type passed; array, shape (N): values are
              interpreted as the density of data sampling at each
              sample point; callable: the callable must take as an
              argument an array of shape (N, nphys), and return an
              array of shape (N) giving the sampling density at each
              point; 'auto': the sample density will be computed
              directly from the data set; note that this can be quite
              slow for large data sets, so it is preferable to specify
              this analytically if it is known; None: data are assumed
              to be uniformly sampled
           reltol : float
              relative error tolerance; errors on all returned
              probabilities p will satisfy either
              abs(p_est - p_true) <= reltol * p_est   OR
              abs(p_est - p_true) <= abstol,
              where p_est is the returned estimate and p_true is the
              true value
           abstol : float
              absolute error tolerance; see above
           leafsize : int
              number of data points in each leaf of the KD tree
           nosort : arraylike of bool, shape (N) | None
              if specified, this keyword causes the KD tree not to be
              sorted along the dimensions for which nosort is True
           thread_safe : bool
              if True, bayesphot will make extra copies of internals
              as needed to ensure thread safety when the computation
              routines (logL, mpdf, mcmc, bestmatch, make_approx_phot,
              make_approx_phys, mpdf_approx) are used with
              multiprocessing; this incurs a minor performance
              penalty, and can be disabled by setting to False if the
              code will not be run with the multiprocessing module

        Returns
           Nothing

        Raises
           IOError, if the bayesphot c library cannot be found
        """

        # Load the c library
        self.__clib = npct.load_library("bayesphot", 
                                        osp.realpath(__file__))

        # Check for diagnostic mode
        self.__clib.diagnostic_mode.restype = c_bool
        self.__clib.diagnostic_mode.argtypes = None
        self.__diag_mode = bool(self.__clib.diagnostic_mode())

        # Define interfaces to all the c library functions
        self.__clib.build_kd.restype = c_void_p
        self.__clib.build_kd.argtypes \
            = [ array_1d_double,   # x
                c_ulong,            # ndim
                c_ulong,           # npt
                ctypes.
                POINTER(c_double), # wgt
                c_ulong,           # leafsize
                array_1d_double,   # bandwidth
                c_int,             # ktype
                c_ulong ]          # minsplit

        self.__clib.build_kd_sortdims.restype = c_void_p
        self.__clib.build_kd_sortdims.argtypes \
            = [ array_1d_double,   # x
                c_ulong,           # ndim
                c_ulong,           # npt
                ctypes.
                POINTER(c_double), # wgt
                c_ulong,           # leafsize
                array_1d_double,   # bandwidth
                c_int,             # ktype
                array_1d_int ]     # nosort

        self.__clib.copy_kd.restype = c_void_p
        self.__clib.copy_kd.argtypes = [ c_void_p ] # kd

        self.__clib.free_kd.restype = None
        self.__clib.free_kd.argtypes = [ c_void_p ]

        self.__clib.free_kd_copy.restype = None
        self.__clib.free_kd_copy.argtypes = [ c_void_p ]

        self.__clib.kd_change_wgt.restype = None
        self.__clib.kd_change_wgt.argtypes \
            = [ POINTER(c_double), # wgt
                c_void_p ]         # kd

        self.__clib.kd_change_bandwidth.restype = None
        self.__clib.kd_change_bandwidth.argtypes \
            = [ array_1d_double,   # bandwidth
                c_void_p ]         # kd

        self.__clib.kd_neighbors.restype = None
        self.__clib.kd_neighbors.argtypes \
            = [ c_void_p,          # kd
                array_1d_double,   # xpt
                POINTER(c_ulong),  # dims
                c_ulong,           # ndim
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # d2
        self.__clib.kd_neighbors_vec.restype = None
        self.__clib.kd_neighbors_vec.argtypes \
            = [ c_void_p,          # kd
                array_1d_double,   # xpt
                POINTER(c_ulong),  # dims
                c_ulong,           # ndim
                c_ulong,           # npt
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # d2

        self.__clib.kd_neighbors_point.restype = None
        self.__clib.kd_neighbors_point.argtypes \
            = [ c_void_p,          # kd
                c_ulong,           # idxpt
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_ulong,     # idx
                array_1d_double ]  # d2
        self.__clib.kd_neighbors_point_vec.restype = None
        self.__clib.kd_neighbors_point_vec.argtypes \
            = [ c_void_p,          # kd
                array_1d_ulong,     # idxpt
                c_ulong,           # npt
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_ulong,     # idx
                array_1d_double ]  # d2

        self.__clib.kd_neighbors_all.restype = None
        self.__clib.kd_neighbors_all.argtypes \
            = [ c_void_p,          # kd
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_ulong,     # idx
                array_1d_double ]  # d2

        self.__clib.kd_pdf.restype = c_double
        if self.__diag_mode:
            self.__clib.kd_pdf.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_double,          # reltol
                    c_double,          # abstol
                    POINTER(c_ulong),  # nodecheck
                    POINTER(c_ulong),  # leafcheck
                    POINTER(c_ulong) ] # termcheck
        else:
            self.__clib.kd_pdf.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_double,          # reltol
                    c_double ]         # abstol

        self.__clib.kd_pdf_grid.restype = None
        self.__clib.kd_pdf_grid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,         # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgrid
                array_1d_ulong,         # dimgrid
                c_ulong,               # ndimgrid
                c_ulong,               # ngrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf
        self.__clib.kd_pdf_reggrid.restype = None
        self.__clib.kd_pdf_reggrid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,         # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgridlo
                array_1d_double,       # xgridhi
                array_1d_ulong,         # ngrid
                array_1d_ulong,         # dimgrid
                c_ulong,               # ndimgrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf

        self.__clib.kd_pdf_vec.restype = None
        if self.__diag_mode:
            self.__clib.kd_pdf_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double,   # pdf
                    array_1d_ulong,     # nodecheck
                    array_1d_ulong,     # leafcheck
                    array_1d_ulong ]    # termcheck
        else:
            self.__clib.kd_pdf_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double ]  # pdf
        self.__clib.kd_pdf_int.restype = c_double
        if self.__diag_mode:
            self.__clib.kd_pdf_int.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_ulong,     # dims
                    c_ulong,           # ndim
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_ulong,     # nodecheck
                    array_1d_ulong,     # leafcheck
                    array_1d_ulong ]    # termcheck
        else:
            self.__clib.kd_pdf_int.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_ulong,     # dims
                    c_ulong,           # ndim
                    c_double,          # reltol
                    c_double ]         # abstol
        self.__clib.kd_pdf_int_vec.restype = None
        if self.__diag_mode:
            self.__clib.kd_pdf_int_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_ulong,     # dims
                    c_ulong,           # ndim
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double,   # pdf
                    array_1d_ulong,     # nodecheck
                    array_1d_ulong,     # leafcheck
                    array_1d_ulong ]    # termcheck
        else:
            self.__clib.kd_pdf_int_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_ulong,     # dims
                    c_ulong,           # ndim
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double ]  # pdf
        self.__clib.kd_pdf_int_grid.restype = None
        self.__clib.kd_pdf_int_grid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,         # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgrid
                array_1d_ulong,         # dimgrid
                c_ulong,               # ndimgrid
                c_ulong,               # ngrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf
        self.__clib.kd_pdf_int_reggrid.restype = None
        self.__clib.kd_pdf_int_reggrid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,         # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgridlo
                array_1d_double,       # xgridhi
                array_1d_ulong,         # ngrid
                array_1d_ulong,         # dimgrid
                c_ulong,               # ndimgrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf                

        self.__clib.kd_rep.restype = c_ulong
        self.__clib.kd_rep.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # x
                array_1d_ulong,        # dims
                c_ulong,               # ndim
                c_double,              # reltol
                POINTER(c_ulong),      # dim_return
                c_ulong,               # ndim_return
                POINTER(POINTER(
                    c_double)),        # xpt
                POINTER(POINTER(
                    c_double)) ]       # wgts
        self.__clib.free_kd_rep.restype = None
        self.__clib.free_kd_rep.argtypes \
            = [ POINTER(POINTER(
                c_double)),            # xpt
               POINTER(POINTER(
                   c_double)) ]        # wgts
        self.__clib.squeeze_rep.restype = c_ulong
        self.__clib.squeeze_rep.argtypes \
            = [ c_ulong,               # npts
                c_ulong,               # ndim
                array_1d_double,       # h
                c_double,              # tol
                POINTER(POINTER(
                    c_double)),        # xpt
                POINTER(POINTER(
                    c_double)) ]       # wgts

        # Record some of the input parameters
        if (ktype == 'gaussian'):
            self.__ktype = 2
        elif (ktype == 'epanechnikov'):
            self.__ktype = 0
        elif (ktype == 'tophat'):
            self.__ktype = 1
        self.leafsize = leafsize
        self.abstol = abstol
        self.reltol = reltol
        self.__sden = sample_density
        self.__sample_density = None
        self.thread_safe = thread_safe

        # Store data set
        self.__dataset = np.ascontiguousarray(dataset)

        # Store list of available filters
        self.__filters = deepcopy(filters)

        # Initialize internal data
        self.__ndata = self.__dataset.shape[0]
        self.__nphys = nphys
        self.__nphot = self.__dataset.shape[1] - self.__nphys
        self.__auto_bw = None
        self.__auto_bw_set = False
        self.__priors = None
        self.__kd_phys = None

        # Build the initial kernel density estimation object, using a
        # dummy bandwidth
        self.__bandwidth = np.ones(self.__nphys + self.__nphot)
        if nosort is None:
            self.__kd = self.__clib.build_kd(
                np.ravel(self.__dataset), self.__dataset.shape[1],
                self.__ndata, None, self.leafsize, self.__bandwidth,
                self.__ktype, 0)
        else:
            nosort_c = np.zeros(nosort.shape, dtype=np.intc)
            nosort_c[:] = nosort == True
            self.__kd = self.__clib.build_kd_sortdims(
                np.ravel(self.__dataset), self.__dataset.shape[1],
                self.__ndata, None, self.leafsize, self.__bandwidth,
                self.__ktype, nosort_c)


        # Initialize the bandwidth
        self.bandwidth = bandwidth

        # Set priors
        self.priors = priors


    ##################################################################
    # De-allocation method
    ##################################################################
    def __del__(self):
        if self.__kd is not None:
            self.__clib.free_kd(self.__kd)
        if self.__kd_phys is not None:
            self.__clib.free_kd(self.__kd_phys)


    ##################################################################
    # Return a copy of the list of available filters
    ##################################################################
    def filters(self):
        return deepcopy(self.__filters)


    ##################################################################
    # Define the priors property
    ##################################################################

    @property
    def priors(self):
        """
        The current set of prior probabilities for every
        simulation in the library
        """
        return self.__priors

    @priors.setter
    def priors(self, pr):
        """
        This function sets the prior probabilities to use

        Parameters
           priors : array, shape (N) | callable | None
              prior probability on each data point; interpretation
              depends on the type passed:
                 array, shape (N) : 
                    values are interpreted as the prior probability of
                    each data point
                 callable : 
                    the callable must take as an argument an array of
                    shape (N, nphys), and return an array of shape (N)
                    giving the prior probability at each data point
                 None :
                    all data points have equal prior probability

        Returns
           Nothing
        """

        # If the prior is unchanged, do nothing
        if (type(pr) == np.ndarray) and \
           (type(self.__priors) == np.ndarray):
            if np.array_equal(pr, self.__priors):
                return
        elif pr == self.__priors:
            return

        # If priors is None, just remove all weighting
        if pr is None:
            self.__clib.kd_change_wgt(None, self.__kd)
            self.__priors = None
            return

        else:
            # If we're here, we have a non-trival prior

            # Evaluate the raw sample density at each point if we have
            # not previously done so
            if self.__sample_density is None:

                # Choose computation method
                if self.__sden is None:

                    # None means uniform sampling
                    self.__sample_density = np.ones(self.__ndata)

                elif hasattr(self.__sden, '__call__'):

                    # Callable, so pass the physical data to the
                    # callable and store the result
                    self.__sample_density \
                        = self.__sden(self.__dataset[:,:self.__nphys])

                elif type(self.__sden) is np.ndarray:

                    # Array, so treat treat this as the data
                    self.__sample_density = self.__sden

                elif self.__sden == 'auto':

                    # We've been asked to calculate the sample density
                    # ourselves, so do so

                    # Create unweighted kernel density object for just
                    # the physical parameters if we have not done so
                    # already
                    if self.__kd_phys is None:
                        self.__dataset_phys \
                            = np.copy(self.__dataset[:,:self.__nphys])
                        self.__kd_phys \
                            = self.__clib.build_kd(
                                np.ravel(self.__dataset_phys), 
                                self.__nphys, self.__ndata,
                                None, self.leafsize, self.__bandwidth,
                                self.__ktype, self.__nphys)

                        # Use the unweighted kernel density object to
                        # evaluate the raw sample density near each data
                        # point, or near a sub-sample which we can
                        # interpolate from
                        self.__sample_density = np.zeros(self.__ndata)
                        nsamp = 500
                        if self.__ndata < nsamp:
                            # Few data points, just use them all; note
                            # that we cannot pass self.__dataset_phys,
                            # because it is not in the same order as
                            # the full data set anymore
                            pts = np.ravel(self.__dataset[:,:self.__nphys])
                            if not self.__diag_mode:
                                self.__clib.kd_pdf_vec(
                                    self.__kd_phys, pts,
                                    self.__ndata, self.reltol, self.abstol,
                                    self.__sample_density)
                            else:
                                nodecheck = np.zeros(self.__ndata, dtype=c_ulong)
                                leafcheck = np.zeros(self.__ndata, dtype=c_ulong)
                                termcheck = np.zeros(self.__ndata, dtype=c_ulong)
                                self.__clib.kd_pdf_vec(
                                    self.__kd_phys, pts,
                                    self.__ndata, self.reltol, self.abstol,
                                    self.__sample_density, nodecheck, 
                                    leafcheck, termcheck)
                        else:
                            # Many data points, so choose a sample at random
                            idxpt = np.array(
                                random.sample(np.arange(self.__ndata), 
                                              nsamp), dtype=c_ulong)
                            pos = np.copy(self.__dataset_phys[idxpt,:])
                            # Add points at the edges of the dataset
                            # to ensure we enclose all the points
                            lowlim = np.amin(self.__dataset_phys,
                                             axis=0)
                            pos = np.append(pos, lowlim)
                            hilim = np.amax(self.__dataset_phys,
                                            axis=0)
                            pos = np.append(pos, hilim)
                            # Compute density at selected points
                            sample_density = np.zeros(nsamp+2)
                            if not self.__diag_mode:
                                self.__clib.kd_pdf_vec(
                                    self.__kd_phys, np.ravel(pos), 
                                    nsamp+2, self.reltol, self.abstol,
                                    sample_density)
                            else:
                                nodecheck = np.zeros(nsamp+2, dtype=c_ulong)
                                leafcheck = np.zeros(nsamp+2, dtype=c_ulong)
                                termcheck = np.zeros(nsamp+2, dtype=c_ulong)
                                self.__clib.kd_pdf_vec(
                                    self.__kd_phys, np.ravel(pos), 
                                    nsamp+2, self.reltol, self.abstol,
                                    sample_density, nodecheck, 
                                    leafcheck, termcheck)
                            # Now interpolate the sample points to all
                            # points in the data set
                            pts = np.ravel(self.__dataset[:,:self.__nphys])
                            self.__sample_density \
                                = np.exp(
                                    interp.griddata(pos, 
                                                    np.log(sample_density),
                                                    pts,
                                                    method='linear')).flatten()

            # We now have the sample density. Record the new prior,
            # and, if our prior is a callable, call it; otherwise
            # just record the input data
            self.__priors = pr
            if hasattr(self.__priors, '__call__'):
                prior_data \
                    = self.__priors(self.__dataset[:,:self.__nphys]).flatten()
            else:
                prior_data = self.__priors

            # Compute the weights from the ratio of the prior to
            # the sample density, then adjust the weights in the kd
            self.__wgt = prior_data / self.__sample_density
            self.__clib.kd_change_wgt(self.__wgt.ctypes.data_as(POINTER(c_double)),
                                      self.__kd)


    ##################################################################
    # Define the bandwidth property
    ##################################################################
 
    @property
    def bandwidth(self):
        """
        The current bandwidth
        """
        return deepcopy(self.__bandwidth)

    @bandwidth.setter
    def bandwidth(self, bw):

        if np.array_equal(self.__bandwidth, bw):
            # If new bandwidth equals old bandwidth, do nothing
            return

        elif bw != 'auto':
            # If we've been given a specified bandwidth, set to that
            if hasattr(bw, '__iter__'):
                self.__bandwidth = np.copy(bw)
            else:
                self.__bandwidth \
                    = np.zeros(self.__nphys + self.__nphot) + bw
            self.__auto_bw_set = False

        else:
            # Automatic bandwidth setting

            # Are we already set on auto? If so, just return
            if self.__auto_bw_set:
                return

            # Do we have a stored value for the automatic bandwidth?
            # If not, we need to compute it.
            if self.__auto_bw is None:

                # Find 10th nearest neighbors
                nneighbor=10
                if self.__ndata > 5000:
                    # For data sets with > 5000 samples, just use a
                    # sub-sample of 5,000, which is more than enough
                    # to get a reasonable estimate of the distribution
                    idxpt = np.array(
                        random.sample(np.arange(self.__ndata), 
                                      5000), dtype=c_ulong)
                    neighbors = np.zeros(nneighbor*5000, 
                                              dtype=c_ulong)
                    d2 = np.zeros(nneighbor*5000)
                    self.__clib.kd_neighbors_point_vec(
                        self.__kd, idxpt, 5000, nneighbor, False,
                        neighbors, d2)

                else:
                    # For smaller data sets, use it all
                    neighbors = np.zeros(nneighbor*self.__ndata, 
                                         dtype=c_ulong)
                    d2 = np.zeros(nneighbor*self.__ndata)
                    idxpt = np.arange(self.__ndata, dtype=c_ulong)
                    self.__clib.kd_neighbors_all(self.__kd, nneighbor, 
                                                 False, neighbors, d2)

                # Take the bandwidth in each dimension to be the 90th
                # percentile of the 10th nearest neighbor distance
                offset = np.abs(self.__dataset[idxpt,:] -
                                self.__dataset[
                                    neighbors[nneighbor-1::nneighbor],:])
                self.__auto_bw = np.zeros(self.__nphys+self.__nphot)
                for i in range(self.__nphys+self.__nphot):
                    self.__auto_bw[i] = np.percentile(offset[:,i], 95)

            # Set to the auto bandwidth
            self.__bandwidth = np.copy(self.__auto_bw)
            self.__auto_bw_set = True

        # Set the new bandwidth
        self.__clib.kd_change_bandwidth(self.__bandwidth, self.__kd)

        # If we have stored priors, and we're using automatic sample
        # density setting, we need to recompute them for the
        # new bandwidth; zero out the stored sample density, and
        # adjust the kernel density estimator for the physical
        # parameters to the new bandwidth before doing so
        if self.__sden == 'auto':
            self.__sample_density = None
            if self.__kd_phys is not None:
                self.__clib.kd_change_bandwidth(
                    self.__bandwidth[:self.__nphys], self.__kd_phys)
            pr = self.priors
            self.priors = None
            self.priors = pr

    ##################################################################
    # Utility methods to set a new bandwidth to account for
    # photometric errors, to reset to the original bandwidth, and to
    # free a temporary kd_tree holder, all in a thread-safe way. These
    # are intended for internal use.
    ##################################################################

    # Broaden the bandwidth by the photometric error, and return a KD
    # object with the new bandwidth; if we are thread-safe, this is a
    # newly-created KD object, which will need to be freed later; if a
    # value for kd_cur is passed in, the bandwidth will be changed for
    # it, and no allocation will be done regardless of thread safety
    def __change_bw_err(self, photerr, kd_cur=None, err_only=False):
        if err_only:
            bandwidth = np.copy(self.__bandwidth)
            bandwidth[self.__nphys:] = photerr
        else:
            err = np.zeros(self.__bandwidth.size)
            err[self.__nphys:] = photerr
            bandwidth = np.sqrt(self.__bandwidth**2+err**2)
        if kd_cur is not None:
            kd_tmp = kd_cur
        elif not self.thread_safe: 
            kd_tmp = self.__kd
        else: 
            kd_tmp = self.__clib.copy_kd(self.__kd)
        self.__clib.kd_change_bandwidth(bandwidth, kd_tmp)
        return kd_tmp

    # Free a temporary KD object
    def __free_bw_err(self, kd_tmp):
        if self.thread_safe:
            self.__clib.free_kd_copy(kd_tmp)

    # Restore bandwidth to the default, de-allocating the temporary KD
    # object if needed
    def __restore_bw_err(self, kd_tmp):
        if not self.thread_safe:
            self.__clib.kd_change_bandwidth(
                self.__bandwidth, self.__kd)
        else:
            self.__clib.free_kd_copy(kd_tmp)

    ##################################################################
    # Method to compute the log likelihood function for a particular
    # set of physical properties given a particular set of photometric
    # properties
    ##################################################################
    def logL(self, physprop, photprop, photerr=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters
           physprop : arraylike, shape (nphys) or (..., nphys)
              array giving values of the physical properties; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photprop : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions

        Returns
           logL : float or arraylike
              natural log of the likelihood function
        """

        # Safety check
        if (np.array(physprop).shape[-1] != self.__nphys) and \
           (self.__nphys > 1):
            raise ValueError("need " + str(self.__nphys) + 
                             " physical properties!")
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")

        # Reshape arrays if necessary
        if (np.array(physprop).shape[-1] != self.__nphys) and \
           (self.__nphys == 1):
            physprop = physprop.reshape(physprop.shape+(1,))
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            photprop = photprop.reshape(photprop.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Figure out number of distinct input sets of physical and
        # photometric properties
        nphys_in = np.array(physprop).size/self.__nphys
        nphot_in = np.array(photprop).size/self.__nphot

        # Figure out how many sets of photometric errors we have. We
        # unfortunately have to iterate over these, because each
        # distinct set requires changing the bandwidth of the kernel
        # density estimation.
        nphot_err = np.array(photerr).size/self.__nphot

        # Allocate an array to hold the results
        if photerr is not None:
            pdf = np.zeros(
                np.broadcast(np.array(physprop)[..., 0],
                             np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape)
        else:
            pdf = np.zeros(
                np.broadcast(np.array(physprop)[..., 0],
                             np.array(photprop)[..., 0]).shape)
        if self.__diag_mode:
            nodecheck = np.zeros(pdf.shape, dtype=c_ulong)
            leafcheck = np.zeros(pdf.shape, dtype=c_ulong)
            termcheck = np.zeros(pdf.shape, dtype=c_ulong)

        # Make an array suitable for passing data to c routines
        cdata = np.zeros(pdf.shape + (self.__nphys+self.__nphot,))
        cdata[..., :self.__nphys] \
            = np.vstack((physprop,) * (pdf.size/nphys_in)). \
            reshape(cdata[..., :self.__nphys].shape)
        cdata[..., self.__nphys:] \
            = np.vstack((photprop,) * (pdf.size/nphot_in)). \
            reshape(cdata[..., self.__nphys:].shape)

        # Separate cases with single / no photometric errors from
        # cases with multiple sets of photometric errors
        if nphot_err <= 1:

            # Case with at most one set of photometric errors

            # Set the bandwidth based on the photometric errors if we
            # were given some
            if photerr is None:
                kd_tmp = self.__kd
            else:
                kd_tmp = self.__change_bw_err(photerr)

            # Call the PDF computation routine
            if not self.__diag_mode:
                self.__clib.kd_pdf_vec(
                    kd_tmp, np.ravel(cdata), pdf.size, 
                    self.reltol, self.abstol, np.ravel(pdf))
            else:
                self.__clib.kd_pdf_vec(
                    kd_tmp, np.ravel(cdata), pdf.size, 
                    self.reltol, self.abstol, np.ravel(pdf),
                    np.ravel(nodecheck), np.ravel(leafcheck),
                    np.ravel(termcheck))

            # Set the bandwidth back to its default if necessary
            if photerr is not None:
                self.__restore_bw_err(kd_tmp)

            # Return
            if not self.__diag_mode:
                return np.log(pdf)
            else:
                return np.log(pdf), nodecheck, leafcheck, termcheck

        else:

            # Case with multiple sets of photometric errors

            # Loop over photometric errors
            kd_tmp = None
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                kd_tmp = self.__change_bw_err(photerr[i], kd_cur=kd_tmp)

                # Grab the corresponding portions of the arrays going
                # to and from the c code
                cdata_sub = cdata[i]
                pdf_sub = np.zeros(np.array(pdf[i]).shape)
                if self.__diag_mode:
                    nodecheck_sub = np.zeros(np.array(nodecheck[i]).shape)
                    leafcheck_sub = np.zeros(np.array(leafcheck[i]).shape)
                    termcheck_sub = np.zeros(np.array(termcheck[i]).shape)

                # Call kernel density estimate with this bandwidth
                if not self.__diag_mode:
                    self.__clib.kd_pdf_vec(
                        kd_tmp, np.ravel(cdata_sub), pdf_sub.size, 
                        self.reltol, self.abstol, np.ravel(pdf_sub))
                else:
                    self.__clib.kd_pdf_vec(
                        kd_tmp, np.ravel(cdata_sub), pdf_sub.size, 
                        self.reltol, self.abstol, np.ravel(pdf_sub),
                        np.ravel(nodecheck_sub), np.ravel(leafcheck_sub),
                        np.ravel(termcheck_sub))
                pdf[i] = pdf_sub
                if self.__diag_mode:
                    nodecheck[i] = nodecheck_sub
                    leafcheck[i] = leafcheck_sub
                    termcheck[i] = termcheck_sub

            # Restore the bandwidth if we changed it
            if photerr is not None:
                self.__restore_bw_err(kd_tmp)

            # Return
            if not self.__diag_mode:
                return np.log(pdf)
            else:
                return np.log(pdf), nodecheck, leafcheck, termcheck


    ##################################################################
    # Function to return the marginal distribution of one or more of
    # the physical properties for a specified set of photometric
    # properties
    ##################################################################
    def mpdf(self, idx, photprop, photerr=None, ngrid=128,
             qmin=None, qmax=None, grid=None, norm=True):
        """
        Returns the marginal probability for one or mode physical
        quantities for one or more input sets of photometric
        properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid.

        Parameters:
           idx : int or listlike containing ints
              index of the physical quantity whose PDF is to be
              computed; if this is an iterable, the joint distribution of
              the indicated quantities is returned
           photprop : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions
           ngrid : int or listlike containing ints
              number of points in each dimension of the output grid;
              if this is an iterable, it must have the same number of
              elements as idx
           qmin : float or listlike
              minimum value in the output grid in each quantity; if
              left as None, defaults to the minimum value in the
              library; if this is an iterable, it must contain the
              same number of elements as idx
           qmax : float or listlike
              maximum value in the output grid in each quantity; if
              left as None, defaults to the maximum value in the
              library; if this is an iterable, it must contain the
              same number of elements as idx
           grid : listlike of arrays
              set of values defining the grid on which the PDF is to
              be evaluated, in the same format used by meshgrid
           norm : bool
              if True, returned pdf's will be normalized to integrate
              to 1

        Returns:
           grid_out : array
              array of values at which the PDF is evaluated; contents
              are the same as returned by meshgrid
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input cluster; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """

        # Safety check
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")
        if (np.amax(idx) > self.__nphys) or (np.amin(idx) < 0) or \
           (not np.array_equal(np.squeeze(np.unique(np.array(idx))), 
                               np.squeeze(np.array([idx])))):
            raise ValueError("need non-repeating indices in " +
                             "the range 0 - {:d}!".
                             format(self.__nphys-1))

        # Reshape arrays if necessary
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            photprop = photprop.reshape(photprop.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Set up the grid of outputs, and the versions of it that we
        # will be passing to the c code
        if grid is not None:
            grid_out = np.array(grid)
            qmin_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            qmax_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            ngrid_tmp = np.array(grid_out.shape[1:], dtype=c_ulong)
        else:
            if qmin is None:
                qmin = np.amin(self.__dataset[:,idx], axis=0)
            if qmax is None:
                qmax = np.amax(self.__dataset[:,idx], axis=0)
            griddims = []
            if hasattr(idx, '__len__'):
                nidx = len(idx)
            else:
                nidx = 1
            if nidx > 1:
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = np.array(ngrid, dtype=c_ulong)
                else:
                    ngrid_tmp = np.array([ngrid]*len(idx), dtype=c_ulong)
                for i in range(len(idx)):
                    griddims.append(qmin[i] + np.arange(ngrid_tmp[i]) * 
                                    float(qmax[i]-qmin[i])/(ngrid_tmp[i]-1))
                grid_out = np.squeeze(np.array(np.meshgrid(*griddims,
                                                           indexing='ij')))
                out_shape = grid_out[0, ...].shape
                qmin_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)])
                qmax_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)])
            else:
                # Case for a single index
                ngrid_tmp = np.array([ngrid], dtype=c_ulong)
                qmin_tmp = np.array([qmin], dtype=np.double).reshape(1)
                qmax_tmp = np.array([qmax], dtype=np.double).reshape(1)
                grid_out = qmin + \
                           np.arange(ngrid) * \
                           float(qmax-qmin)/(ngrid-1)
                out_shape = grid_out.shape

        # Figure out how many distinct photometric values we've been
        # given, and how many sets of photometric errors
        nphot = np.array(photprop).size/self.__nphot
        nphot_err = np.array(photerr).size/self.__nphot

        # Set up a grid to hold the outputs
        if photerr is not None:
            pdf = np.zeros(
                np.broadcast(np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape +
                out_shape)
        else:
            pdf = np.zeros(np.array(photprop)[..., 0].shape +
                           out_shape)

        # Prepare data for c library
        if hasattr(idx, '__len__'):
            nidx = len(idx)
        else:
            nidx = 1
        dims = np.zeros(nidx+self.__nphot, dtype=c_ulong)
        dims[:nidx] = idx
        dims[nidx:] = self.__nphys+np.arange(self.__nphot, dtype=c_ulong)
        ndim = c_ulong(nidx + self.__nphot)
        phottmp = np.array(photprop, dtype=c_double)

        # Separate cases with single / no photometric errors from
        # cases with multiple sets of photometric errors
        if nphot_err <= 1:

            # Case with at most one set of photometric errors

            # Set the bandwidth based on the photometric errors if we
            # were given some
            if photerr is None:
                kd_tmp = self.__kd
            else:
                kd_tmp = self.__change_bw_err(photerr)

            # Call the PDF computation routine; note that we call
            # kd_pdf_int_vec if we're actually marginalizing over any
            # physical properties (which is the case if len(idx) <
            # nphys), but we invoke kd_pdf_vec if we're not actually
            # marginalizing (len(idx)==nphys) because then we don't
            # need to do any integration
            if nidx < self.__nphys:
                self.__clib.kd_pdf_int_reggrid(
                    kd_tmp, np.ravel(phottmp),
                    dims[nidx:], self.__nphot, nphot,
                    qmin_tmp, qmax_tmp, ngrid_tmp, dims[:nidx], nidx,
                    self.reltol, self.abstol, np.ravel(pdf))
            else:
                self.__clib.kd_pdf_reggrid(
                    kd_tmp, np.ravel(phottmp),
                    dims[nidx:], self.__nphot, nphot,
                    qmin_tmp, qmax_tmp, ngrid_tmp, dims[:nidx], nidx,
                    self.reltol, self.abstol, np.ravel(pdf))

        else:

            # Case with multiple sets of photometric errors

            # Loop over photometric errors
            kd_tmp = None
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                kd_tmp = self.__change_bw_err(photerr[i], kd_cur=kd_tmp)

                # Grab the corresponding portions of the arrays going
                # to and from the c code
                if photprop.size < photerr.size:
                    phot_sub = np.ravel(phottmp)
                else:
                    phot_sub = np.ravel(phottmp[i])
                pdf_sub = np.zeros(np.array(pdf[i]).shape)

                # Call kernel density estimate with this bandwidth
                if nidx < self.__nphys:
                    self.__clib.kd_pdf_int_reggrid(
                        kd_tmp, phot_sub,
                        dims[nidx:], self.__nphot, 
                        phot_sub.size/self.__nphot,
                        qmin_tmp, qmax_tmp, ngrid_tmp, 
                        dims[:nidx], nidx,
                        self.reltol, self.abstol, np.ravel(pdf_sub))
                else:
                    self.__clib.kd_pdf_reggrid(
                        kd_tmp, phot_sub,
                        dims[nidx:], self.__nphot, 
                        phot_sub.size/self.__nphot,
                        qmin_tmp, qmax_tmp, ngrid_tmp, dims[:nidx], nidx,
                        self.reltol, self.abstol, np.ravel(pdf_sub))
                pdf[i] = pdf_sub

        # Set the bandwidth back to its default if necessary
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # Normalize if requested
        if norm:

            # Compute the sizes of the output cells
            if nidx == 1:
                cellsize = np.zeros(grid_out.size)
                cellsize[1:-1] = 0.5*(grid_out[2:]-grid_out[:-2])
                cellsize[0] = grid_out[1] - grid_out[0]
                cellsize[-1] = grid_out[-1] - grid_out[-2]
            else:
                # Get the cell sizes in each dimension
                csize = []
                for i in range(nidx):
                    vec = grid_out[(i,)+i*(0,)+(slice(None),) + 
                                   (grid_out.shape[0]-i-1)*(0,)]
                    csize.append(np.zeros(vec.size))
                    csize[i][1:-1] = 0.5*(vec[2:]-vec[:-2])
                    csize[i][0] = vec[1] - vec[0]
                    csize[i][-1] = vec[-1] - vec[-2]
                # Take outer product to get grid of sizes
                cellsize = np.multiply.outer(csize[0], csize[1])
                for i in range(2, nidx):
                    cellsize = np.multiply.outer(cellsize, csize[i])

            # Compute integral
            normfac = np.sum(pdf*cellsize, axis = 
                             tuple(range(np.array(photprop).ndim-1, pdf.ndim)))

            # Normalize
            pdf = np.transpose(np.transpose(pdf)/normfac)

        # Return
        return grid_out, pdf


    ##################################################################
    # Method to return log likelihood function at a specified set of
    # physical properties for a particular set of photometric
    # variables; this is set up for use by emcee, and is not intended
    # for use by humans
    ##################################################################
    def __logL(self, *args):
        x = np.zeros(self.__nphys+self.__nphot)
        x[:self.__nphys] = args[0]
        x[self.__nphys:] = args[1:-1]
        kd_tmp = args[-1]
        if not self.__diag_mode:
            return np.log(self.__clib.kd_pdf(kd_tmp, x, self.reltol,
                                             self.abstol))
        else:
            nodecheck = c_ulong(0)
            leafcheck = c_ulong(0)
            termcheck = c_ulong(0)
            return np.log(self.__clib.kd_pdf(kd_tmp, x, self.reltol,
                                             self.abstol, nodecheck,
                                             leafcheck, termcheck))


    ##################################################################
    # Function to compute an MCMC sampler for a particular set of
    # photometric values
    ##################################################################
    def mcmc(self, photprop, photerr=None, mc_walkers=100,
             mc_steps=500, mc_burn_in=50):
        """
        This function returns a sample of MCMC walkers sampling the
        physical parameters at a specified set of photometric values.

        Parameters:
           photprop : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions
           mc_walkers : int
              number of walkers to use in the MCMC
           mc_steps : int
              number of steps in the MCMC
           mc_burn_in : int
              number of steps to consider "burn-in" and discard

        Returns
           samples : array
              array of sample points returned by the MCMC
        """

        # See if we have emcee
        if not mc_avail:
            raise ImportError("unable to import emcee")

        # Safety check
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")
 
        # Reshape arrays if necessary
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            photprop = photprop.reshape(photprop.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Prepare storage for output
        samples = []

        # Make dummy photometric errors if necessary
        if photerr is None:
            photerr = np.zeros(self.__nphot)

        # Loop over photometric errors
        kd_tmp = None
        for i in np.ndindex(*np.array(photerr).shape[:-1]):

            # Set bandwidth based on photometric error for this
            # iteration
            kd_tmp = self.__change_bw_err(photerr[i], kd_cur=kd_tmp)
 
            # Grab the clusters that go with this photometric error
            if photprop.ndim < photerr.ndim:
                ph = photprop
            else:
                ph = photprop[i]

            # Loop over clusters
            for j in np.ndindex(*np.array(ph).shape[:-1]):

                # Grab photometric values for this cluster
                ph_tmp = ph[j]

                # Search the data set for the sample cluster closest to
                # the observed luminosities; this will be the starting
                # point for the MCMC
                dims = np.arange(self.__nphys,
                                 self.__nphys+self.__nphot, 
                                 dtype=np.uint32)
                nearpt = np.zeros(self.__nphys+self.__nphot)
                wgt = np.zeros(1)
                dist2 = np.zeros(1)
                self.__clib. \
                    kd_neighbors(kd_tmp, ph_tmp,
                                 dims.ctypes.data_as(POINTER(c_ulong)),
                                 self.__nphot, 1, True, nearpt, 
                                 wgt.ctypes.data_as(POINTER(c_double)),
                                 dist2)

                # Generate a set of starting points by scattering walkers
                # around the starting position
                pos = [nearpt[:self.__nphys] + 
                       self.bandwidth[:self.__nphys] * 
                       np.random.randn(self.__nphys) 
                       for i in range(mc_walkers)]

                # Run the MCMC
                sampler=emcee.EnsembleSampler(mc_walkers, self.__nphys, 
                                              self.__logL, 
                                              args=[ph_tmp, kd_tmp])
                sampler.run_mcmc(pos, mc_steps)

                # Store the result
                samples.append(sampler.chain[:,mc_burn_in:,:].
                               reshape((-1,self.__nphys)))


        # Set bandwidth back to default if necessary
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # Reshape the samples
        samples = np.squeeze(
            np.array(samples).reshape(
                np.broadcast(np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape +
                samples[0].shape))

        # Return
        return samples


    ##################################################################
    # Function to return the N best matches in the library to an input
    # set of photometric properties
    ##################################################################
    def bestmatch(self, phot, photerr=None, nmatch=1, 
                  bandwidth_units=False):
        """
        Searches through the simulation library and returns the closest
        matches to an input set of photometry.

        Parameters:
           phot : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors, which must have the
              same shape as phot; if this is not None,
              then distances will be measured in units of the
              photometric error if bandwidth_units is False, or in
              units of the bandwidth added in quadrature with the
              errors if it is True
           nmatch : int
              number of matches to return; returned matches will be
              ordered by distance from the input
           bandwidth_units : bool
              if False, distances are computed based on the
              logarithmic difference in luminosity; if True, they are
              measured in units of the bandwidth

        Returns:
           matches : array, shape (..., nmatch, nphys + nfilter)
              best matches to the input photometry; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted
           dist : array, shape (..., nmatch)
              distances between the matches and the input photometry
        """

        # Safety check
        if (np.array(phot).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")
            if np.array(photerr).shape != np.array(phot).shape:
                raise ValueError("phot and photerr must have the same shape!")

        # Reshape arrays if necessary
        if (np.array(phot).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            phot = phot.reshape(phot.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Figure out how many distinct photometric values we've been
        # given, and how many sets of photometric errors
        nphot = np.array(phot).size/self.__nphot
        nphot_err = np.array(photerr).size/self.__nphot

        # Figure out what shape the output should have
        if nphot == 1 and nphot_err == 1:
            outshape = [self.__nphys + self.__nphot]
        elif photerr is None:
            outshape = list(phot.shape[:-1]) + \
                       [self.__nphys + self.__nphot]
        else:
            outshape = list(
                np.broadcast(np.array(phot)[..., 0],
                             np.array(photerr)[..., 0]).shape) + \
                [self.__nphys + self.__nphot]
        if nmatch > 1:
            outshape.insert(-1, nmatch)

        # Create array to hold output
        matches = np.zeros(outshape)
        if len(outshape) > 1:
            wgts = np.zeros(outshape[:-1])
            d2 = np.zeros(outshape[:-1])
        else:
            wgts = np.zeros([1])
            d2 = np.zeros([1])

        # Do we have errors?
        if photerr is None:

            # No, so just call the c neighbor-finding routine
            dims = np.arange(self.__nphys, self.__nphot+self.__nphys,
                             dtype=c_ulong)
            self.__clib.kd_neighbors_vec(
                self.__kd, np.ravel(phot), 
                dims.ctypes.data_as(POINTER(c_ulong)),
                self.__nphot, np.array(phot).size/self.__nphot,
                nmatch, bandwidth_units, np.ravel(matches),
                wgts.ctypes.data_as(POINTER(c_double)), 
                np.ravel(d2))

        elif nphot_err == 1:

            # Yes, we have errors, but only one set, so no need to
            # loop

            # Change the bandwidth to match the input errors
            kd_tmp = self.__change_bw_err(
                photerr, err_only = not bandwidth_units)

            # Now call c neighbor-finding routine
            dims = np.arange(self.__nphys, self.__nphot+self.__nphys,
                             dtype=c_ulong)
            self.__clib.kd_neighbors_vec(
                kd_tmp, np.ravel(phot), 
                dims.ctypes.data_as(POINTER(c_ulong)),
                self.__nphot, np.array(phot).size/self.__nphot,
                nmatch, True, np.ravel(matches),
                wgts.ctypes.data_as(POINTER(c_double)), 
                np.ravel(d2))

            # Restore bandwidth to previous value
            self.__restore_bw_err(kd_tmp)

        else:

            # We have multiple sets of errors, so loop
            ptr = 0
            kd_tmp = None
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                kd_tmp = self.__change_bw_err(
                    photerr[i], kd_cur = kd_tmp,
                    err_only = not bandwidth_units)

                # Call c neighbor-finding routine
                offset1 = ptr*nmatch
                offset2 = offset1*(self.__nphot+self.__nphys)
                dims = np.arange(self.__nphys, self.__nphot+self.__nphys,
                                 dtype=c_ulong)                
                self.__clib.kd_neighbors_vec(
                    self.__kd, np.ravel(phot[i]), 
                    dims.ctypes.data_as(POINTER(c_ulong)),
                    self.__nphot, np.array(phot[i]).size/self.__nphot,
                    nmatch, True, np.ravel(matches)[offset2:],
                    wgts.ctypes.data_as(POINTER(c_double)), 
                    np.ravel(d2)[offset1:])
                ptr = ptr+1

            # Restore bandwidth to previous value
            self.__restore_bw_err(kd_tmp)

        # Return
        return matches, np.sqrt(d2)

    ##################################################################
    # Functions to return an approximate kernel density representation
    # of a given point; that is, for an input physical point or
    # photometric point (there are two versions below), this routine
    # returns a set of points and weights that can be used to compute
    # an approximation to the PDF for it anywhere in space, much
    # faster than could be done using the full data set
    ##################################################################
    def make_approx_phot(self, phys, squeeze=True, filter_ignore=None):
        """
        Returns an object that can be used for a fast approximation of
        the PDF of photometric properties that corresponds to a set of
        physical properties. The PDF produced by summing over the
        points returned is guaranteed to account for at least 1-reltol
        of the marginal photometric probability, and to represent the
        shape of the PDF in photometric space within a local accuracy
        of reltol as well.

        Parameters:
           phys : arraylike, shape (nphys) or (N, nphys)
              the set or sets of physical properties for which the
              approximation is to be generated
           squeeze : bool
              if True, the representation returned will be squeezed to
              minimize the number of points included, using reltol as
              the error tolerance
           filter_ignore : None or listlike of bool
              if None, the kernel density representation returned
              covers all filters; otherwise this must be a listlike of
              bool, one entry per filter, with a value of False
              indicating that filter should be excluded from the
              values returned; suppressing filters can allow for more
              efficient representations

        Returns:
           x : array, shape (M, nphot), or a list of such arrays
              an array containing the list of points to be used for
              the approximation, where nphot is the number of
              photometric filters being returned
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points

        Notes:
           if the requested relative tolerance cannot be reached for
           numerical reasons (usually because the input point is too
           far from the library to allow accurate computation), x and
           wgts will be return as None, and a warning will be issued
        """

        # If given only one set of physical variables, make a
        # temporary list we can loop over
        if hasattr(phys[0], '__iter__'):
            phystmp = phys
        else:
            phystmp = [phys]

        # Specify dimensions to return, and set bandwidth for them
        if filter_ignore is None:
            dim_return_ptr = None
            ndim_return = self.__nphot
            bw = self.__bandwidth[self.__nphys:]
        else:
            dim_return = np.where(np.logical_not(
                np.array(filter_ignore)))[0] + self.__nphys
            dim_return_ptr = dim_return.ctypes.data_as(POINTER(c_ulong))
            ndim_return = len(dim_return)
            bw = self.__bandwidth[dim_return]

        # Loop over input physical variables
        x = []
        wgts = []
        for ph in phystmp:

            # Safety check
            if len(ph) != self.__nphys:
                raise ValueError("need " + str(self.__nphys) + 
                                 " physical properties!")

            # Call c library routine to make the representation
            xout = POINTER(c_double)()
            wgtsout = POINTER(c_double)()
            npts = self.__clib.kd_rep(self.__kd, np.array(ph), 
                                      np.arange(self.__nphys, dtype=c_ulong),
                                      self.__nphys, self.reltol, 
                                      dim_return_ptr, ndim_return,
                                      ctypes.byref(xout), 
                                      ctypes.byref(wgtsout))

            # Check if the c routine encountered an error due to
            # insufficient precision; if so, issue warning and return
            # None in place of xpt and wgts
            if npts == 0:
                warn("bp.make_approx_phot: requested precision cannot be reached")
                return None

            # Call c library routine to squeeze the representation
            if squeeze:
                npts = self.__clib.squeeze_rep(
                    npts, ndim_return, bw, self.reltol,
                    ctypes.byref(xout), ctypes.byref(wgtsout))

            # Convert the returned values into numpy arrays; copy the
            # data so that we can free the c buffers
            xsave = np.copy(npct.as_array(xout, shape=(npts, ndim_return)))
            wgtssave = np.copy(npct.as_array(wgtsout, shape=(npts,)))

            # Free the c buffers
            self.__clib.free_kd_rep(ctypes.byref(xout), 
                                    ctypes.byref(wgtsout))

            # Append to the lists we'll be returning
            x.append(xsave)
            wgts.append(wgtssave)

        # If we were given a single object, return something in the
        # same shape
        if not hasattr(phys[0], '__iter__'):
            x = x[0]
            wgts = wgts[0]

        # Return
        return x, wgts

    def make_approx_phys(self, phot, photerr=None, squeeze=True,
                         phys_ignore=None):
        """
        Returns an object that can be used for a fast approximation of
        the PDF of physical properties that corresponds to a set of
        photometric properties. The PDF produced by summing over the
        points returned is guaranteed to account for at least 1-reltol
        of the marginal photometric probability, and to represent the
        shape of the PDF in photometric space within a local accuracy
        of reltol as well.

        Parameters:
           phot : arraylike, shape (nfilter) or (N, nfilter)
              the set or sets of photometric properties for which the
              approximation is to be generated
           photerr : arraylike, shape (nfilter) or (N, nfilter)
              array giving photometric errors; the number of elements
              in the output lists will be the size that results from
              broadcasting together the leading dimensions of phot and
              photerr
           squeeze : bool
              if True, the representation returned will be squeezed to
              minimize the number of points included, using reltol as
              the error tolerance
           phys_ignore : None or listlike of bool
              if None, the kernel density representation returned
              covers all physical properties; otherwise this must be a
              listlike of bool, one entry per physical dimension, with
              a value of False indicating that dimension should be
              excluded from the values returned; suppressing
              dimensions can allow for more efficient representations

        Returns:
           x : array, shape (M, nphys), or a list of such arrays
              an array containing the list of points to be used for
              the approximation, where nphys is the number of
              physical dimensions being returned
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points

        Notes:
           if the requested relative tolerance cannot be reached for
           numerical reasons (usually because the input point is too
           far from the library to allow accurate computation), x and
           wgts will be return as None, and a warning will be issued
        """

        # If given only one set of photometric variables, make a
        # temporary list we can loop over; same for photometric errors
        if hasattr(phot[0], '__iter__'):
            phottmp = phot
        else:
            phottmp = [phot]
        if photerr is not None:
            if hasattr(photerr[0], '__iter__'):
                photerrtmp = photerr
            else:
                photerrtmp = [photerr]
        else:
            photerrtmp = [None]

        # Specify dimensions to return, and set bandwidth for them
        if phys_ignore is None:
            dim_return_ptr = None
            ndim_return = self.__nphys
            bw = self.__bandwidth[:self.__nphys]
        else:
            dim_return = np.where(np.logical_not(
                np.array(phys_ignore)))[0]
            dim_return_ptr = dim_return.ctypes.data_as(POINTER(c_ulong))
            ndim_return = len(dim_return)
            bw = self.__bandwidth[dim_return]

        # Loop over input photometric variables
        x = []
        wgts = []
        for ph in phottmp:

            # Safety check
            if len(ph) != self.__nphot:
                raise ValueError("need " + str(self.__nphot) + 
                                 " photometric properties!")

            # Loop over photometric errors
            kd_tmp = None
            for pherr in photerrtmp:

                # Change the bandwidth if necessary
                if pherr is not None:
                    if len(pherr) != self.__nphot:
                        raise ValueError("need " + str(self.__nphot) + 
                                         " photometric errors!")
                    kd_tmp = self.__change_bw_err(pherr, kd_cur=kd_tmp)

                # Call c library routine to make the representation
                xout = POINTER(c_double)()
                wgtsout = POINTER(c_double)()
                npts = self.__clib.kd_rep(
                    kd_tmp, np.array(ph), 
                    np.arange(self.__nphot, dtype=c_ulong)+self.__nphys,
                    self.__nphot, self.reltol, 
                    dim_return_ptr, ndim_return,
                    ctypes.byref(xout), 
                    ctypes.byref(wgtsout))

                # Check if the c routine encountered an error due to
                # insufficient precision; if so, issue warning and
                # return None in place of xpt and wgts
                if npts == 0:
                    warn("bp.make_approx_phot: requested precision cannot be reached")
                    if kd_tmp is not None:
                        self.__restore_bw_err(kd_tmp)
                    return None

                # Call c library routine to squeeze the representation
                if squeeze:
                    npts = self.__clib.squeeze_rep(
                        npts, ndim_return, bw, self.reltol,
                        ctypes.byref(xout), ctypes.byref(wgtsout))

                # Convert the returned values into numpy arrays; copy the
                # data so that we can free the c buffers
                xsave = np.copy(npct.as_array(xout, 
                                              shape=(npts, ndim_return)))
                wgtssave = np.copy(npct.as_array(wgtsout, shape=(npts,)))

                # Free the c buffers
                self.__clib.free_kd_rep(ctypes.byref(xout), 
                                        ctypes.byref(wgtsout))

                # Append to the lists we'll be returning
                x.append(xsave)
                wgts.append(wgtssave)

        # Reset the bandwidth
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # If we were given a single object, return something in the
        # same shape
        if not hasattr(phot[0], '__iter__'):
            if photerr is None:
                x = x[0]
                wgts = wgts[0]
            elif not hasattr(photerr[0], '__iter__'):
                x = x[0]
                wgts = wgts[0]

        # Return
        return x, wgts
            
    ##################################################################
    # Method to squeeze a representation that has already been created
    ##################################################################
    def squeeze_rep(self, x, wgts, dims=None):
        """
        Takes an input array of positions and weights that form a
        kernel density representation and approximates them using
        fewer points, using an error tolerance of reltol

        Parameters:
           x : array, shape (N, ndim)
              an array of points forming a kernel density
              representation; on exit, x will be resized to (M, ndim)
              with M <= N
           wgts : array, shape (N)
              an array of weights for the kernel density
              representation; on exit, wgts will be resized to (M),
              with M <= N
           dims : array, shape (ndim)
              array specifying which dimensions in the kernel density
              representation the coordinates in x correspond to; if
              left as None, they are assumed to correspond to the
              first ndim dimensions in the data set

        Returns:
           Nothing
        """

        # Get the bandwidth
        if dims is None:
            bw = self.__bandwidth[:x.shape[1]]
        else:
            bw = self.__bandwidth[dims]

        # Prepare data to pass to the c routine
        npt = x.shape[0]
        ndim = x.shape[1]
        xtmp = POINTER(c_double)(x.ctypes.data_as(POINTER(c_double)))
        wgtstmp = POINTER(wgts.ctypes.data_as(POINTER(c_double)))

        # Call the c squeeze routine
        npt_out = self.__clib.squeeze_rep(
            npt, ndim, bw, self.reltol,
            ctypes.byref(xtmp), ctypes.byref(wgtstmp))

        # Discard the old metadata for x and wgts, and rebuild them
        # using the new metadata
        x = npct.as_array(xtmp, shape=(npt_out, ndim))
        wgts = npct.as_array(wgtstmp, shape=(ndim,))


    ##################################################################
    # Routines to use the approximate representations returned by
    # make_approx_phys and make_approx_phot to compute marginal
    # posterior probabilities
    ##################################################################

    def mpdf_approx(self, x, wgts, dims='phys', dims_return=None,
                    ngrid=64, qmin=None, qmax=None, grid=None,
                    norm=True):
        """
        Returns the marginal posterior PDF computed from a kernel
        density approximation returned by make_approx_phys or
        make_approx_phot. Outputs are computed on a grid of values, in
        the same style as meshgrid.

        Parameters:
           x : array, shape (M, ndim), or a list of such arrays
              array of points retured by make_approx_phot or
              make_approx_phys
           wgts : array, shape (M) or a list of such arrays
              array of weights returned by make_approx_phot or
              make_approx_phys
           dims : 'phys' | 'phot' | arraylike of ints
              dimensions covered by x and wgts; the strings 'phys' or
              'phot' indicate that they cover all physical or
              photometric dimensions, and correspond to the defaults
              returned by make_approx_phys and make_approx_phot,
              respectively; if dims is an array of ints, these specify
              the dimensions covered by x and wgts, where the
              physical dimensions are numbered 0, 1, ... nphys-1, and
              the photometric ones are nphys, nphys+1,
              ... nphys+nphot-1
           dims_return : None or arraylike of ints
              if None, the output PDF has the same dimensions as
              specified in dms; if not, then dimreturn must be a
              subset of dim, and a marginal PDF in certain dimensions
              will be generated
           ngrid : int or listlike containing ints
              number of points in each dimension of the output grid;
              if this is an iterable, it must have the same number of
              elements as idx
           qmin : float or listlike
              minimum value in the output grid in each quantity; if
              left as None, the range is chosen automatically to zoom in
              on the maximum of the final PDF; if this is an iterable,
              it must contain the same number of elements as idx
           qmax : float or listlike
              maximum value in the output grid in each quantity; if
              left as None, the range is chosen automatically to zoom in
              on the maximum of the final PDF; if this is an iterable,
              it must contain the same number of elements as idx
           grid : listlike of arrays
              set of values defining the grid on which the PDF is to
              be evaluated, in the same format used by meshgrid
           norm : bool
              if True, returned pdf's will be normalized to integrate
              to 1

        Returns:
           grid_out : array
              array of values at which the PDF is evaluated; contents
              are the same as returned by meshgrid
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input cluster; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """

        # Set up input and output dimensions
        if dims == 'phys':
            dim_in = np.arange(self.__nphys)
        elif dims == 'phot':
            dim_in = np.arange(self.__nphot) + self.__nphys
        else:
            dim_in = dims
        if dims_return is None:
            dim_out = dim_in
        elif not hasattr(dims_return, '__iter__'):
            dim_out = [dims_return]
            if not dim_out in dim_in:
                raise ValueError("dims_return must be a subset of dims!")
        elif set(dims_return) <= set(dim_in):
            dim_out = dims_return
        else:
            raise ValueError("dims_return must be a subset of dims!")
        nidx = len(dim_out)

        # Set up the output grid
        if grid is not None:
            grid_out = np.meshgrid(grid)
            qmin \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            qmax \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            ngrid_tmp = np.squeeze(grid_out).shape
            for i in range(len(dim_out)):
                griddims.append(np.linspace(qmin[i], qmax[i], 
                                            ngrid_tmp[i]))
        else:
            if qmin is None or qmax is None:
                if qmin is None:
                    get_qmin = True
                    qmin = []
                else:
                    get_qmin = False
                if qmax is None:
                    qmax = []
                    get_qmax = True
                else:
                    get_qmax = False
                for d in dim_out:
                    xtmp = np.copy(x[:,d])
                    wgttmp = np.copy(wgts)
                    idx = np.argsort(xtmp)
                    xtmp = xtmp[idx]
                    wgttmp = wgttmp[idx]
                    wgtsum = np.cumsum(wgttmp)
                    if get_qmin:
                        qmin.append(xtmp[np.argmax(wgtsum > 0.001)])
                    if get_qmax:
                        if wgtsum[-1] > 0.999:
                            qmax.append(xtmp[np.argmax(wgtsum >
                                                       0.999)])
                        else:
                            qmax.append(xtmp[-1])
            griddims = []
            if nidx > 1:
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = np.array(ngrid)
                else:
                    ngrid_tmp = np.array([ngrid]*len(dim_out))
                for i in range(len(dim_out)):
                    griddims.append(np.linspace(qmin[i], qmax[i], 
                                                ngrid_tmp[i]))
                grid_out = np.squeeze(np.array(np.meshgrid(*griddims,
                                                           indexing='ij')))
            else:
                # Case for a single index
                grid_out = qmin + \
                           np.arange(ngrid) * \
                           float(qmax-qmin)/(ngrid-1)
                griddims = [grid_out]

        # Compute result; in the multi-dimensional case this involves
        # some tricky array indexing
        if nidx == 1:
            pdf = np.sum(wgts*np.exp(-(x[:,dim_out[0]]-grid_out)**2 / 
                                     (2*self.__bandwidth[dim_out[0]]**2)),
                         axis=1)
        else:
            compsum = np.exp(
                -np.subtract.outer(x[:,dim_out[0]], griddims[0])**2 /
                (2.0*self.__bandwidth[dim_out[0]]**2))
            for d, grd in zip(dim_out[1:], griddims[1:]):
                comp = np.exp(-np.subtract.outer(x[:,d], grd)**2 / \
                              (2.0*self.__bandwidth[d]**2))
                compsum = np.einsum('...j,...k->...jk', compsum, comp)
            pdf = np.einsum('i,i...', wgts, compsum)

        # Normalize if requested
        if norm:

            # Compute the sizes of the output cells
            if nidx == 1:
                cellsize = np.zeros(grid_out.size)
                cellsize[1:-1] = 0.5*(grid_out[2:]-grid_out[:-2])
                cellsize[0] = grid_out[1] - grid_out[0]
                cellsize[-1] = grid_out[-1] - grid_out[-2]
            else:
                # Get the cell sizes in each dimension
                csize = []
                for i in range(nidx):
                    vec = grid_out[(i,)+i*(0,)+(slice(None),) + 
                                   (grid_out.shape[0]-i-1)*(0,)]
                    csize.append(np.zeros(vec.size))
                    csize[i][1:-1] = 0.5*(vec[2:]-vec[:-2])
                    csize[i][0] = vec[1] - vec[0]
                    csize[i][-1] = vec[-1] - vec[-2]
                # Take outer product to get grid of sizes
                cellsize = np.multiply.outer(csize[0], csize[1])
                for i in range(2, nidx):
                    cellsize = np.multiply.outer(cellsize, csize[i])

            # Compute integral
            normfac = np.sum(pdf*cellsize)

            # Normalize
            pdf = pdf/normfac

        # Return
        return grid_out, pdf
