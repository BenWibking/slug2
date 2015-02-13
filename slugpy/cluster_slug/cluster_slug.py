"""
This defines a class that can be used to estimate the PDF of star
cluster properties (mass, age, extinction) from a set of input
photometry in various bands.
"""

import numpy as np
import copy
import os
import os.path as osp
import warnings
import urllib2

# Import the data reading and Bayesian inference stuff we need
from ..bayesphot import bp
from ..read_cluster_prop import read_cluster_prop
from ..read_cluster_phot import read_cluster_phot

##################################################################
# Define the cluster_slug class                                  #
##################################################################

class cluster_slug(object):
    """
    A class that can be used to estimate the PDF of star cluster
    properties (mass, age, extinction) from a set of input photometry
    in various bands.

    Properties
       priors : array, shape (N) | callable | None
          prior probability on each data point; interpretation
          depends on the type passed; array, shape (N): values are
          interpreted as the prior probability of each data point;
          callable: the callable must take as an argument an array
          of shape (N, nphys), and return an array of shape (N)
          giving the prior probability at each data point; None:
          all data points have equal prior probability
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, filters=None, photsystem=None,
                 bw_phys=0.1, bw_phot=None, ktype='gaussian', 
                 priors=None, sample_density=None, reltol=1.0e-2,
                 abstol=1.0e-6, leafsize=16, use_nebular=True):
        """
        Initialize a cluster_slug object.

        Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/cluster_slug/clusterslug_mw
           filters : iterable of stringlike
              list of filter names to be used for inferenence
           photsystem : None or string
              If photsystem is None, the library will be left in
              whatever photometric system was used to write
              it. Alternately, if it is a string, the data will be
              converted to the specified photometric system. Allowable
              values are 'L_nu', 'L_lambda', 'AB', 'STMAG', and
              'Vega', corresponding to the options defined in the SLUG
              code. Once this is set, any subsequent photometric data
              input are assumed to be in the same photometric system.
           bw_phys : 'auto' | float | array, shape (3)
              bandwidth for the physical quantities in the kernel
              density estimation; if set to 'auto', the bandwidth will
              be estimated automatically; if set to a scalar quantity,
              this will be used for all physical quantities
           bw_phot : None | 'auto' | float | array
              bandwidth for the photometric quantities; if set to
              None, defaults to 0.25 mag / 0.1 dex; if set to 'auto',
              bandwidth is estimated automatically; if set to a float,
              this bandwidth is used for all photometric dimensions;
              if set to an array, the array must have the same number
              of dimensions as len(filters)
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
           use_nebular : bool
              if True, photometry including nebular emission will be
              used if available; if not, nebular emission will be
              omitted

        Returns
           Nothing

        Raises
           IOError, if the library cannot be found
        """

        # If using the default library, assign the library name
        if libname is None:
            self.__libname = osp.join('cluster_slug', 'clusterslug_mw')
            if 'SLUG_DIR' in os.environ:
                self.__libname = osp.join(os.environ['SLUG_DIR'], 
                                          self.__libname)
        else:
            self.__libname = libname

        # Load the cluster physical properties
        try:
            prop = read_cluster_prop(self.__libname)
        except IOError:

            # If we're here, we failed to load the library. If we were
            # given a library file name explicitly, just raise an
            # error.
            if libname is not None:
                raise IOError("unable to open library {}".
                              format(self.__libname))

            # If we've made it to here, we were asked to open the
            # default library but failed to do so. Check if the
            # failure could be because we don't have astropy and thus
            # can't open fits files. If that's the cause, print out a
            # helpful error message.
            try:
                import astropy.io.fits as fits
            except ImportError:
                raise IOError("failed to read default cluster_slug " +
                              "library cluster_slug/clusterslug_mw " +
                              "due to missing " +
                              "astropy.io.fits; install astropy " +
                              "or specify a library in a non-FITS " +
                              "format")

            # If we're here, we couldn't open the default library, and
            # it's not because we don't have FITS capability. The file
            # must not exist, or must be damaged. Check if we're in
            # interactive mode. If not, just raise an error and
            # suggest the user to go get the library file.
            errstr = "Unable to open default sfr_slug " + \
                     "library file sfr_slug/SFR_SLUG."
            import __main__ as main
            if hasattr(main, '__file__'):
                # We're not interactive; just raise an error
                raise IOError(errstr + " " +
                              "Try downloading it from " +
                              "https://sites.google.com/site/runslug/data")

            # If we're here, we don't have hte library file, but we
            # are in interactive mode. Thus offer the user an option
            # to go download the file now.
            usr_response \
                = raw_input(errstr + " Would you like to download it "
                            "now (warning: 12 GB)? [y/n] ").\
                lower().strip()
            if not usr_response in ['yes', 'y', 'ye']:
                # User didn't say yes, so raise error
                raise IOError("Unable to proceeed")

            # If we're here, download the files
            print("Fetching clusterslug_mw_cluster_prop.fits " +
                  "(this may take a while)...")
            url = urllib2.urlopen(
                'https://dl.dropboxusercontent.com/s/oxuuxa1ci0zird4/clusterslug_mw_cluster_prop.fits?dl=0')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'clusterslug_mw_cluster_prop.fits'), 'wb')
            fp.write(rawdata)
            fp.close()
            print("Fetching clusterslug_mw_cluser_phot.fits " +
                  "(this make take a while)...")
            url = urllib2.urlopen(
                'https://dl.dropboxusercontent.com/s/arwy8u0xwt9dnz2/clusterslug_mw_cluster_phot.fits?dl=0')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'clusterslug_mw_cluster_phot.fits'), 'wb')
            fp.write(rawdata)
            fp.close()

            # Now try reading the data
            try:
                prop = read_integrated_prop(self.__libname)
            except IOError:
                raise IOError("still unable to open default library")

        # Store the physical properties
        self.__ds_phys = np.zeros((len(prop.id), 3))
        self.__ds_phys[:,0] = np.log10(prop.actual_mass)
        self.__ds_phys[:,1] = np.log10(prop.time - prop.form_time)
        self.__ds_phys[:,2] = prop.A_V

        # Record available filters
        filter_info = read_cluster_phot(self.__libname,
                                        filters_only=True, 
                                        nofilterdata=True)
        self.__allfilters = filter_info.filter_names
        self.__allunits = filter_info.filter_units

        # Record other stuff that we'll use later
        self.__photsystem = photsystem
        self.__use_nebular = use_nebular
        self.__ktype = ktype
        self.__priors = priors
        self.__sample_density = sample_density
        self.__reltol = reltol
        self.__abstol = abstol

        # Set the physical bandwidth
        self.__bw_phys = copy.deepcopy(bw_phys)

        # Initialize list of photometric data we've read to empty dict
        self.__photdata = {}
        self.__photbw = {}

        # Initialize an empty list of filter sets
        self.__filtersets = []

        # If we have been given a filter list, create the data set to
        # go with it
        if filters is not None:
            self.add_filters(filters, bandwidth=bw_phot)


    ##################################################################
    # Method to load the data off disk for a particular filter
    ##################################################################
    def load_data(self, filter_name):
        """
        Loads photometric data for the specified filter into memory

        Parameters:
           filter_name : string
              name of filter to load

        Returns:
           None

        Raises:
           ValueError, if filter_name is not one of the available
           filters
        """

        # Make sure we have this filter; if not, raise error
        if filter_name not in self.__allfilters:
            raise ValueError("no data available for filter {}".
                             format(filter_name))

        # Suppress obnoxious numpy warning messages here
        errstate = np.geterr()
        np.seterr(divide='ignore', invalid='ignore', over='ignore',
                  under='ignore')

        # Special case: for ionizing fluxes, always load the
        # non-nebular, non-extincted value, and don't do any
        # photometric system conversion
        if filter_name == 'QH0' or filter_name == 'QHe0' or \
           filter_name == 'QHe1':
            phot = read_cluster_phot(self.__libname, 
                                     read_filters=filter_name,
                                     read_nebular=False,
                                     read_extinct=False,
                                     phot_only=True)
        else:

            # Load data; first try for nebular, extincted data, unless
            # we've been told not to use nebular data
            phot = read_cluster_phot(self.__libname, 
                                     read_filters=filter_name,
                                     read_nebular=self.__use_nebular,
                                     read_extinct=True,
                                     phot_only=True,
                                     photsystem=self.__photsystem)

            # Make sure we got nebular data if we wanted it; if not, try
            # looking for non-nebular data
            if self.__use_nebular:
                if 'phot_neb_ex' in phot._fields:
                    # We want nebular data, and we have it
                    warn_nebular = False
                    nebular = True
                    phdata = phot.phot_neb_ex
                else:
                    # We want nebular data, but we don't have it
                    warn_nebular = True
                    nebular = False
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=False,
                        read_extinct=True,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    phdata = phot.phot_ex
            else:
                # We don't want nebular data
                warn_nebular = False
                nebular = False
                phot = read_cluster_phot(
                    self.__libname, 
                    read_filters=filter_name,
                    read_nebular=False,
                    read_extinct=True,
                    phot_only=True,
                    photsystem=self.__photsystem)
                phdata = phot.phot_ex
        
            # Make sure data aren't NaN's. If they are, that means we
            # don't have extinction for this filter, so we need to use
            # non-extincted data.
            if np.isnan(phdata[0]):
                warn_extinct = True
                if nebular:
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=True,
                        read_extinct=False,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    phdata = phot.phot_neb
                else:
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=False,
                        read_extinct=False,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    phdata = phot.phot
            else:
                warn_extinct = False

            # Issue warnings if necessary
            if warn_nebular:
                warnstr = ("cluster_slug: nebular data requested for "+ \
                           "filter {}, but is not available; using non-"+ \
                           "nebular data instead").format(filter_name)
                warnings.warn(warnstr)
            if warn_extinct:
                warnstr = ("cluster_slug: extincted data requested for "+ \
                           "filter {}, but is not available; using non-"+ \
                           "extinced data instead").format(filter_name)
                warnings.warn(warnstr)

        # Take the log of the photometric values if they're recorded
        # in a linear system; also set the bandwidth based on whether
        # we're in a magnitude system or not; we can override this
        # later if we want
        if 'mag' not in phot.filter_units[0]:
            phdata[phdata <= 0] = 1.0e-99
            phdata = np.log10(phdata)
            self.__photbw[filter_name] = 0.1
        else:
            self.__photbw[filter_name] = 0.25

        # Fix any inf's that were generated by photometric system
        # conversions or taking logs
        phdata[np.isinf(phdata)] = 99.0

        # Store the data
        self.__photdata[filter_name] = phdata

        # Restore the numpy error state
        np.seterr(divide=errstate['divide'], over=errstate['over'], 
                  under=errstate['under'], invalid=errstate['invalid'])


    ##################################################################
    # Method to return list of available filters and filter units
    ##################################################################
    def filters(self):
        """
        Returns list of all available filters

        Parameters:
           None

        Returns:
           filters : list of strings
              list of available filter names
        """

        return copy.deepcopy(self.__allfilters)

    def filter_units(self):
        """
        Returns list of all available filter units

        Parameters:
           None

        Returns:
           units : list of strings
              list of available filter units
        """

        return copy.deepcopy(self.__allunits)


    ##################################################################
    # Method to prepare to analyze a particular set of filters
    ##################################################################
    def add_filters(self, filters, bandwidth=None):
        """
        Add a set of filters to use for cluster property estimation

        Parameters
           filters : iterable of stringlike
              list of filter names to be used for inferenence
           bandwidth : None | 'auto' | float | array
              bandwidth for the photometric quantities; if set to
              None, defaults to 0.3 mag / 0.12 dex; if set to 'auto',
              bandwidth is estimated automatically; if set to a float,
              this bandwidth is used for all physical photometric
              dimensions; if set to an array, the array must have the
              same number of entries as 3+len(filters)

        Returns
           nothing
        """

        # If we already have this filter set in our dict, just set the
        # bandwidth and return
        for i, f in enumerate(self.__filtersets):
            if filters == f['filters']:
                if bandwidth is not None:
                    self.__filtersets[i]['bp'].bandwidth = bandwidth
                else:
                    bw = np.zeros(3+len(filters))
                    bw[:3] = self.__bw_phys
                    for j in range(len(filters)):
                        bw[3+j] = self.__photbw[filters[j]]
                    self.__filtersets[i]['bp'].bandwidth = bw
                return

        # We're adding a new filter set, so save its name
        newfilter = { 'filters' : copy.deepcopy(filters) }

        # Construct data set to use with this filter combination, and
        # add the physical property data
        newfilter['dataset'] = np.zeros((self.__ds_phys.shape[0], 
                                         3+len(filters)))
        newfilter['dataset'][:,:3] = self.__ds_phys

        # Loop over filters
        for i, f in enumerate(filters):

            # Do we have this filter loaded already? If not, read it.
            if f not in self.__photdata.keys():
                self.load_data(f)

            # Add data for this filter
            newfilter['dataset'][:,3+i] = np.squeeze(self.__photdata[f])

        # Set bandwidth
        if self.__bw_phys == 'auto' or bandwidth == 'auto':
            bw = 'auto'
        else:
            bw = np.zeros(3+len(filters))
            bw[:3] = self.__bw_phys
            if bandwidth is not None:
                bw[3:] = bandwidth
            else:
                for i in range(len(filters)):
                    bw[3+i] = self.__photbw[f]

        # Build the bp object
        newfilter['bp'] = bp(newfilter['dataset'], 3,
                             filters=filters,
                             bandwidth = bw,
                             ktype = self.__ktype,
                             priors = self.__priors,
                             sample_density = self.__sample_density,
                             reltol = self.__reltol,
                             abstol = self.__abstol)

        # Save to the master filter list
        self.__filtersets.append(newfilter)


    ##################################################################
    # Define the priors property. This just wraps around the
    # corresponding property defined for bp objects.
    ##################################################################
    @property
    def priors(self):
        return self.__priors

    @priors.setter
    def priors(self, pr):
        self.__priors = pr
        for f in self.__filtersets:
            f['bp'].priors = self.__priors

    ##################################################################
    # Define the abstol and reltol properties. These update the
    # current values of thes, and reset them for all child bp
    # objects.
    ##################################################################
    @property
    def abstol(self):
        return self.__abstol

    @abstol.setter
    def abstol(self, newtol):
        self.__abstol = newtol
        for f in self.__filtersets:
            f['bp'].abstol = self.__abstol

    @property
    def reltol(self):
        return self.__reltol

    @reltol.setter
    def reltol(self, newtol):
        self.__reltol = newtol
        for f in self.__filtersets:
            f['bp'].reltol = self.__reltol

    ##################################################################
    # Wrappers around the bp logL, mpdf, mcmc, and bestmatch functions
    ##################################################################
    def logL(self, physprop, photprop, photerr=None, filters=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters:
           physprop : arraylike, shape (3) or (..., 3)
              array giving values of the log M, log T, and A_V; for a
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
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           logL : float or arraylike
              natural log of the likelihood function
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    logL(physprop, photprop, photerr)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.logL(physprop, photprop, photerr)


    def mpdf(self, idx, photprop, photerr=None, ngrid=128,
             qmin=None, qmax=None, grid=None, norm=True,
             filters=None):
        """
        Returns the marginal probability for one or mode physical
        quantities for one or more input sets of photometric
        properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid

        Parameters:
           idx : int or listlike containing ints
              index of the physical quantity whose PDF is to be
              computed; 0 = log M, 1 = log T, 2 = A_V; if this is an
              iterable, the joint distribution of the indicated
              quantities is returned
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
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

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

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf(idx, photprop, photerr, ngrid,
                         qmin, qmax, grid, norm)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.mpdf(idx, photprop, photerr, ngrid,
                           qmin, qmax, grid, norm)


    def mcmc(self, photprop, photerr=None, mc_walkers=100,
             mc_steps=500, mc_burn_in=50, filters=None):
        """
        This function returns a sample of MCMC walkers for cluster
        mass, age, and extinction 

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
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns
           samples : array
              array of sample points returned by the MCMC
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mcmc(photprop, photerr, mc_walkers, mc_steps, 
                         mc_burn_in)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the mcmc method
            return bp.mcmc(photprop, photerr, mc_walkers, mc_steps, 
                           mc_burn_in)


    def bestmatch(self, phot, nmatch=1, bandwidth_units=False,
                  filters=None):
        """
        Searches through the simulation library and returns the closest
        matches to an input set of photometry.

        Parameters:
           phot : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           nmatch : int
              number of matches to return; returned matches will be
              ordered by distance from the input
           bandwidth_units : bool
              if False, distances are computed based on the
              logarithmic difference in luminosity; if True, they are
              measured in units of the bandwidth

        Returns:
           matches : array, shape (..., nmatch, 3 + nfilter)
              best matches to the input photometry; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted; in the
              final dimension, the first 3 elements give log M, log T,
              and A_V, while the last nfilter give the photometric
              values
           dist : array, shape (..., nmatch)
              distances between the matches and the input photometry
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    bestmatch(phot, nmatch, bandwidth_units)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the bestmatch method
            return bp.bestmatch(phot, nmatch, bandwidth_units)
