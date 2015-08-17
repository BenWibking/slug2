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
# Sample density of the default library                          #
##################################################################
def _default_sample_density(physprop):
    logm = physprop[:,0]
    logt = physprop[:,1]
    sden = np.ones(len(logm))
    sden[logm > 4] = sden[logm > 4] * 1.0/10.**(logm[logm > 4]-4)
    sden[logt > 8] = sden[logt > 8] * 1.0/10.**(logt[logt > 8]-8)
    return sden

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
       abstol : float
          absolute error tolerance for kernel density estimation
       reltol : float
          relative error tolerance for kernel density estimation

    Methods
       filters() : 
          returns list of filters available in the library
       filter_units() :
          returns units for available filters
       add_filters() : 
          adds a set of filters for use in parameter estimation
       logL() : 
          compute log likelihood at a particular set of physical and
          photometric parameters
       mpdf() : 
          computer marginal posterior probability distribution for a
          set of photometric measurements
       mcmc() :
          due MCMC estimation of the posterior PDF on a set of
          photometric measurments
       bestmatch() : 
          find the simulations in the library that are the closest
          matches to the input photometry
       make_approx_phot() :
          given a set of physical properties, return a set of points
          that can be used for fast approximation of the corresponding
          photometric properties
       make_approx_phys() :
          given a set of photometric properties, return a set of points
          that can be used for fast approximation of the corresponding
          physical properties
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, filters=None, photsystem=None,
                 bw_phys=0.1, bw_phot=None, ktype='gaussian', 
                 priors=None, sample_density=None, reltol=1.0e-2,
                 abstol=1.0e-8, leafsize=16, use_nebular=True,
                 use_extinction=True):
        """
        Initialize a cluster_slug object.

        Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/cluster_slug/modp020_chabrier_MW
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
           bw_phys : 'auto' | float | array, shape (2) | array, shape (3)
              bandwidth for the physical quantities in the kernel
              density estimation; if set to 'auto', the bandwidth will
              be estimated automatically; if set to a scalar quantity,
              this will be used for all physical quantities; if set to
              an array, the array must have 2 elements if
              use_extinction is False, or 3 if it is True
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
              to be uniformly sampled, or to be sampled as the default
              library is if libname is also None
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
           use_extinction : bool
              if True, photometry including extinction will be used;
              if not, it will be omitted, and in this case no results
              making use of the A_V dimension will be available

        Returns
           Nothing

        Raises
           IOError, if the library cannot be found
        """

        # If using the default library, assign the library name
        if libname is None:
            self.__libname = osp.join('cluster_slug', 'modp020_chabrier_MW')
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
            errstr = "Unable to open default cluster_slug " + \
                     "library file cluster_slug/clusterslug_mw."
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
                            "now (warning: 20 GB)? [y/n] ").\
                lower().strip()
            if not usr_response in ['yes', 'y', 'ye']:
                # User didn't say yes, so raise error
                raise IOError("Unable to proceeed")

            # If we're here, download the files
            print("Fetching modp020_chabrier_MW_cluster_prop " +
                  "(this may take a while)...")
            url = urllib2.urlopen(
                'https://www.dropbox.com/s/tg2ogad713nfux0/modp020_chabrier_MW_cluster_prop.fits?dl=0')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'modp020_chabrier_MW_cluster_phot.fits'), 'wb')
            fp.write(rawdata)
            fp.close()
            print("Fetching modp020_chabrier_MW_cluster_phot.fits " +
                  "(this make take a while)...")
            url = urllib2.urlopen(
                'https://www.dropbox.com/s/inmeeuldkxsazzs/modp020_chabrier_MW_cluster_phot.fits?dl=0')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'modp020_chabrier_MW_cluster_phot.fits'), 'wb')
            fp.write(rawdata)
            fp.close()

            # Now try reading the data
            try:
                prop = read_integrated_prop(self.__libname)
            except IOError:
                raise IOError("still unable to open default library")

        # Store the physical properties
        if use_extinction:
            self.__nphys = 3
            self.__ds_phys = np.zeros((len(prop.id), 3))
        else:
            self.__nphys = 2
            self.__ds_phys = np.zeros((len(prop.id), 2))
        self.__ds_phys[:,0] = np.log10(prop.actual_mass)
        self.__ds_phys[:,1] = np.log10(prop.time - prop.form_time)
        if use_extinction:
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
        self.__use_extinction = use_extinction
        self.__ktype = ktype
        self.__priors = priors
        if (sample_density is not None) or \
           (libname is not None):
            self.__sample_density = sample_density
        else:
            self.__sample_density = _default_sample_density
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
    def load_data(self, filter_name, bandwidth=None):
        """
        Loads photometric data for the specified filter into memory

        Parameters:
           filter_name : string
              name of filter to load
           bandwidth : float
              default bandwidth for this filter

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

            # Load data; first try reading the requested combination
            # of nebular and extincted values
            phot = read_cluster_phot(self.__libname, 
                                     read_filters=filter_name,
                                     read_nebular=self.__use_nebular,
                                     read_extinct=self.__use_extinction,
                                     phot_only=True,
                                     photsystem=self.__photsystem)

            # Make sure we have the data we want; if not, lots of ugly
            # special cases as fallbacks
            if (self.__use_nebular and self.__use_extinction):
                # Wanted nebular + extinction
                if 'phot_neb_ex' in phot._fields:
                    # Got it
                    phdata = phot.phot_neb_ex
                    warn_nebular = False
                    warn_extinct = False
                    nebular = True
                    extinct = True
                else:
                    # Didn't get it; try without nebular
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=False,
                        read_extinct=True,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    if 'phot_ex' in phot._fields:
                        phdata = phot.phot_ex
                        warn_nebular = True
                        warn_extinct = False
                        nebular = False
                        extinct = True
                    else:
                        # Couldn't get extinction only; try nebular only
                        phot = read_cluster_phot(
                            self.__libname, 
                            read_filters=filter_name,
                            read_nebular=True,
                            read_extinct=False,
                            phot_only=True,
                            photsystem=self.__photsystem)
                        if 'phot_neb' in phot._fields:
                            phdata = phot.phot_neb
                            warn_nebular = False
                            warn_extinct = True
                            nebular = True
                            extinct = False
                        else:
                            # Ultimate fallback: no nebular or extinction
                            phot = read_cluster_phot(
                                self.__libname, 
                                read_filters=filter_name,
                                read_nebular=False,
                                read_extinct=False,
                                phot_only=True,
                                photsystem=self.__photsystem)
                            phdata = phot.phot
                            warn_nebular = True
                            warn_extinct = True
                            nebular = False
                            extinct = False

            elif self.__use_nebular:
                # Want nebular, no extinction
                if 'phot_neb' in phot._fields:
                    # Got it
                    phdata = phot.phot_neb
                    warn_nebular = False
                    warn_extinct = False
                    nebular = True
                    extinct = False
                else:
                    # Didn't get it; go to no nebular or extinction
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=False,
                        read_extinct=False,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    phdata = phot.phot
                    warn_nebular = True
                    warn_extinct = False
                    nebular = False
                    extinct = False

            elif self._use_extinction:
                # Wanted extinction, no nebular
                if 'phot_ex' in phot._fields:
                    # Got it
                    phdata = phot.phot_ex
                    warn_nebular = False
                    warn_extinct = False
                    nebular = False
                    extinct = True
                else:
                    # Didn't get it; go to no nebular or extinction
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=False,
                        read_extinct=False,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    phdata = phot.phot
                    warn_nebular = False
                    warn_extinct = True
                    nebular = False
                    extinct = False

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
            if bandwidth is None:
                self.__photbw[filter_name] = 0.1
            else:
                self.__photbw[filter_name] = bandwidth
        else:
            if bandwidth is None:
                self.__photbw[filter_name] = 0.25
            else:
                self.__photbw[filter_name] = bandwidth

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
              None, defaults to 0.25 mag / 0.1 dex; if set to 'auto',
              bandwidth is estimated automatically; if set to a float,
              this bandwidth is used for all physical photometric
              dimensions; if set to an array, the array must have the
              same number of entries as nphys+len(filters)

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
                    bw = np.zeros(self.__nphys+len(filters))
                    bw[:self.__nphys] = self.__bw_phys
                    for j in range(len(filters)):
                        bw[self.__nphys+j] = self.__photbw[filters[j]]
                    self.__filtersets[i]['bp'].bandwidth = bw
                return

        # We're adding a new filter set, so save its name
        newfilter = { 'filters' : copy.deepcopy(filters) }

        # Construct data set to use with this filter combination, and
        # add the physical property data
        newfilter['dataset'] = np.zeros((self.__ds_phys.shape[0], 
                                         self.__nphys+len(filters)))
        newfilter['dataset'][:,:self.__nphys] = self.__ds_phys

        # Loop over filters
        for i, f in enumerate(filters):

            # Do we have this filter loaded already? If not, read it.
            if f not in self.__photdata.keys():
                if bandwidth is None or type(bandwidth) is str:
                    self.load_data(f)
                else:
                    if hasattr(bandwidth, '__iter__'):
                        self.load_data(f, bandwidth=bandwidth[i])
                    else:
                        self.load_data(f, bandwidth=bandwidth)

            # Add data for this filter
            newfilter['dataset'][:,self.__nphys+i] \
                = np.squeeze(self.__photdata[f])

        # Set bandwidth
        if self.__bw_phys == 'auto' or bandwidth == 'auto':
            bw = 'auto'
        else:
            bw = np.zeros(self.__nphys+len(filters))
            bw[:self.__nphys] = self.__bw_phys
            if bandwidth is not None:
                bw[self.__nphys:] = bandwidth
            else:
                for i in range(len(filters)):
                    bw[self.__nphys+i] = self.__photbw[f]

        # Build the bp object
        newfilter['bp'] = bp(newfilter['dataset'], self.__nphys,
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
    # The functions below just wrap around the bayesphot functions of
    # the same name. They all have in common that they accept an
    # additional keyword argument, filters, which specifies the filter
    # set to use.
    ##################################################################

    def logL(self, physprop, photprop, photerr=None, filters=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters:
           physprop : arraylike, shape (nhpys) or (..., nphys)
              array giving values of the log M, log T, and A_V; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions; if created with use_extinct =
              False, the A_V dimension should be omitted
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
                    mpdf(idx, photprop, photerr=photerr, 
                         ngrid=ngrid, qmin=qmin, qmax=qmax, 
                         grid=grid, norm=norm)
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
            return bp.mpdf(idx, photprop, photerr=photerr,
                           ngrid=ngrid, qmin=qmin, qmax=qmax,
                           grid=grid, norm=norm)


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


    def bestmatch(self, phot, photerr=None, nmatch=1, 
                  bandwidth_units=False, filters=None):
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
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           matches : array, shape (..., nmatch, nphys + nfilter)
              best matches to the input photometry; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted; in the
              final dimension, the first 3 elements give log M, log T,
              and A_V, while the last nfilter give the photometric
              values; if created with use_extinct = False, the A_V
              dimension is omitted
           dist : array, shape (..., nmatch)
              distances between the matches and the input photometry
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    bestmatch(phot, photerr, nmatch, bandwidth_units)
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
            return bp.bestmatch(phot, photerr, nmatch, bandwidth_units)

    def make_approx_phot(self, phys, squeeze=True, filter_ignore=None,
                         filters=None):
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
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           x : array, shape (M, nphot), or a list of such arrays
              an array containing the list of points to be used for
              the approximation
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    make_approx_phot(phys, squeeze=squeeze,
                                     filter_ignore=filter_ignore)
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

            # Call the method
            return bp.make_approx_phot(phys, squeeze=squeeze,
                                       filter_ignore=filter_ignore)


    def make_approx_phys(self, phot, photerr=None, squeeze=True, 
                         phys_ignore=None, filters=None):
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
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           x : array, shape (M, nphys), or a list of such arrays
              an array containing the list of points to be used for
              the approximation, where nphys is the number of
              physical dimensions being returned
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    make_approx_phys(phot, photerr=photerr,
                                     squeeze=squeeze,
                                     phys_ignore=phys_ignore)
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

            # Call the method
            return bp.make_approx_phys(phot, photerr=photerr,
                                       squeeze=squeeze,
                                       phys_ignore=phys_ignore)


    def mpdf_approx(self, x, wgts, dims='phys', dims_return=None,
                    ngrid=64, qmin=None, qmax=None, grid=None,
                    norm=True, filters=None):
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
              specified in dims; if not, then dimreturn must be a
              subset of dims, and a marginal PDF in certain dimensions
              will be generated
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
                    mpdf_approx(x, wgts, dims=dims, 
                                dims_return=dims_return,
                                ngrid=ngrid, qmin=qmin,
                                qmax=qmax, grid=grid, norm=norm)

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

            # Call the method
            return bp.mpdf_approx(x, wgts, dims=dims, 
                                  dims_return=dims_return,
                                  ngrid=ngrid, qmin=qmin,
                                  qmax=qmax, grid=grid, norm=norm)
