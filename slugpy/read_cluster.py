"""
Function to read all cluster data for a SLUG2 run.
"""

from collections import namedtuple
from read_cluster_prop import read_cluster_prop
from read_cluster_phot import read_cluster_phot
from read_cluster_spec import read_cluster_spec
from cloudy.read_cluster_cloudyphot import read_cluster_cloudyphot
from cloudy.read_cluster_cloudylines import read_cluster_cloudylines
from cloudy.read_cluster_cloudyspec import read_cluster_cloudyspec

def read_cluster(model_name, output_dir=None, fmt=None,
                 nofilterdata=False, photsystem=None, verbose=False,
                 read_filters=None, read_nebular=None, 
                 read_extinct=None, read_info=None):
    """
    Function to read all cluster data for a SLUG2 run.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set
       fmt : string
          Format for the file to be read. Allowed values are 'ascii',
          'bin' or 'binary, and 'fits'. If one of these is set, the code
          will only attempt to open ASCII-, binary-, or FITS-formatted
          output, ending in .txt., .bin, or .fits, respectively. If set
          to None, the code will try to open ASCII files first, then if
          it fails try binary files, and if it fails again try FITS
          files.
       nofilterdata : bool
          If True, the routine does not attempt to read the filter
          response data from the standard location
       photsystem : None or string
          If photsystem is None, the data will be returned in the same
          photometric system in which they were read. Alternately, if it
          is a string, the data will be converted to the specified
          photometric system. Allowable values are 'L_nu', 'L_lambda',
          'AB', 'STMAG', and 'Vega', corresponding to the options defined
          in the SLUG code. If this is set and the conversion requested
          involves a conversion from a wavelength-based system to a
          frequency-based one, nofilterdata must be False so that the
          central wavelength of the photometric filters is available.
       verbose : bool
          If True, verbose output is printed as code runs
       read_filters : None | string | listlike containing strings
          If this is None, photometric data on all filters is
          read. Otherwise only filters whose name(s) match the input
          filter names ar read.
       read_nebular : None | bool
          If True, only photometric data with the nebular contribution
          is read; if False, only data without it is read. Default
          behavior is to read all data.
       read_extinct : None | bool
          If True, only photometric data with extinction applied is
          read; if False, only data without it is read. Default
          behavior is to read all data.
       read_info : dict
          On return, this dict will contain the keys 'prop_name',
          'phot_name', 'spec_name', 'cloudyspec_name', 'cloudylines_name'
          and 'format', giving the names of the files read and the format
          they were in; 'format' will be one of 'ascii', 'binary', or
          'fits'. If one of the files is not present, the corresponding
          _name key will be omitted from the dict.

    Returns
       A namedtuple containing the following fields:

       (Always present)

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          time at which cluster's properties are being evaluated

       (Present if the run being read contains a cluster_prop file)

       form_time : array
          time when cluster formed
       lifetime : array
          time at which cluster will disrupt
       target_mass : array
          target cluster mass
       actual_mass : array
          actual mass at formation
       live_mass : array
          mass of currently living stars
       stellar_mass : array
          mass of all stars, living and stellar remnants
       num_star : array, dtype ulonglong
          number of living stars in cluster being treated stochastically
       max_star_mass : array
          mass of most massive living star in cluster

       (Present if the run being read contains a cluster_spec file)

       wl : array
          wavelength, in Angstrom
       spec : array, shape (N_cluster, N_wavelength)
          specific luminosity of each cluster at each wavelength, in erg/s/A
       wl_neb : array
          wavelength for the nebular spectrum, in Angstrom (present
          only if SLUG was run with nebular emission enabled)
       spec_neb : array, shape (N_cluster, N_wavelength)
          specific luminosity at each wavelength and each time for each
          trial, including emission and absorption by the HII region,
          in erg/s/A (present only if SLUG was run with nebular
          emission enabled)
       wl_ex : array
          wavelength for the extincted spectrum, in Angstrom (present
          only if SLUG was run with extinction enabled)
       spec_ex : array, shape (N_cluster, N_wavelength)
          specific luminosity at each wavelength in wl_ex and each
          time for each trial after extinction has been applied, in
          erg/s/A (present only if SLUG was run with extinction
          enabled)
       wl_neb_ex : array
          wavelength for the extincted spectrum with nebular emission,
          in Angstrom (present only if SLUG was run with both nebular
          emission and extinction enabled)
       spec_neb_ex : array, shape (N_cluster, N_wavelength)
          specific luminosity at each wavelength in wl_ex and each
          time for each trial including emission and absorption by the
          HII region, after extinction has been applied, in erg/s/A
          (present only if SLUG was run with nebular emission and
          extinction both enabled)

       (Present if the run being read contains a cluster_phot file)

       filter_names : list of string
          a list giving the name for each filter
       filter_units : list of string
          a list giving the units for each filter
       filter_wl_cen : list
          central wavelength of each filter; this is set to None for the
          filters Lbol, QH0, QHe0, and QHe1; omitted if nofilterdata is
          True
       filter_wl : list of arrays
          a list giving the wavelength table for each filter; this is
          None for the filters Lbol, QH0, QHe0, and QHe1; omitted if
          nofilterdata is True
       filter_response : list of arrays
          a list giving the photon response function for each filter;
          this is None for the filters Lbol, QH0, QHe0, and QHe1; omitted
          if nofilterdata is True 
       phot : array, shape (N_cluster, N_filter)
          photometric value in each filter for each cluster; units are as
          indicated in the units field

       If extinction is enabled, phot_ex will contain photometry  
       after extinction has been applied.     

       (Present if the run being read contains a cluster_cloudyspec file)

       cloudy_wl : array
          wavelength, in Angstrom
       cloudy_inc : array, shape (N_cluster, N_wavelength)
          specific luminosity of the cluster's stellar radiation field at
          each wavelength, in erg/s/A
       cloudy_trans :  array, shape (N_cluster, N_wavelength)
          specific luminosity of the stellar radiation field after it has
          passed through the HII region, at each wavelength, in erg/s/A
       cloudy_emit :  array, shape (N_cluster, N_wavelength)
          specific luminosity of the radiation field emitted by the HII
          region, at each wavelength, in erg/s/A
       cloudy_trans_emit :  array, shape (N_cluster, N_wavelength)
          the sum of the emitted and transmitted fields; this is what
          would be seen by an observer looking at both the star cluster
          and its nebula

       (Present if the run being read contains a cluster_cloudylines file)

       cloudy_linelabel : array, dtype='S4', shape (N_lines)
          labels for the lines, following cloudy's 4 character line label
          notation
       cloudy_linewl : array, shape (N_lines)
          rest wavelength for each line, in Angstrom
       cloudy_linelum : array, shape (N_cluster, N_lines)
          luminosity of each line at each time for each trial, in erg/s

       (Present if the run being read contains a cluster_cloudyphot file)

       cloudy_filter_names : list of string
          a list giving the name for each filter
       cloudy_filter_units : list of string
          a list giving the units for each filter
       cloudy_filter_wl_eff : list
          effective wavelength of each filter; this is set to None for the
          filters Lbol, QH0, QHe0, and QHe1; omitted if nofilterdata is
          True
       cloudy_filter_wl : list of arrays
          a list giving the wavelength table for each filter; this is
          None for the filters Lbol, QH0, QHe0, and QHe1; omitted if
          nofilterdata is True
       cloudy_filter_response : list of arrays
          a list giving the photon response function for each filter;
          this is None for the filters Lbol, QH0, QHe0, and QHe1; omitted
          if nofilterdata is True 
       cloudy_filter_beta : list
          powerlaw index beta for each filter; used to normalize the
          photometry
       cloudy_filter_wl_c : list
          pivot wavelength for each filter; used to normalize the photometry
       cloudy_phot_trans : array, shape (N_cluster, N_filter)
          photometric value for each cluster in each filter for the
          transmitted light (i.e., the starlight remaining after it has
          passed through the HII region); units are as indicated in
          the units field
       cloudy_phot_emit : array, shape (N_cluster, N_filter)
          photometric value for each cluster in each filter for the
          emitted light (i.e., the diffuse light emitted by the HII
          region); units are as indicated in the units field
       cloudy_phot_trans_emit : array, shape (N_cluster, N_filter)
          photometric value in each filter for each cluster for the
          transmitted plus emitted light (i.e., the light coming
          directly from the stars after absorption by the HII region,
          plus the diffuse light emitted by the HII region); units are as
          indicated in the units field 

    Raises
       IOError, if no photometry file can be opened
       ValueError, if photsystem is set to an unknown values
    """

    # Read properties
    try:
        prop = read_cluster_prop(model_name, output_dir, fmt=fmt, 
                                 verbose=verbose,
                                 read_info=read_info)
        if read_info is not None:
            read_info['prop_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        prop = None

    # Read spectra
    try:
        spec = read_cluster_spec(model_name, output_dir, fmt=fmt,
                                 verbose=verbose,
                                 read_info=read_info)
        if read_info is not None:
            read_info['spec_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        spec = None

    # Read photometry
    try:
        phot = read_cluster_phot(model_name, output_dir, fmt=fmt, 
                                 nofilterdata=nofilterdata,
                                 photsystem=photsystem, verbose=verbose,
                                 read_filters=read_filters,
                                 read_nebular=read_nebular,
                                 read_extinct=read_extinct,
                                 read_info=read_info)
        if read_info is not None:
            read_info['phot_name'] = read_info['fname']
            del read_info['fname']
    except IOError, e:
        if str(e.message)[:16] == 'requested filter' and \
           str(e.message)[-14:] == 'not available!':
            # If the IO error is that we didn't have the requested
            # filter, re-raise the error
            raise IOError(e.message)
        else:
            phot = None

    # Read cloudy spectra
    try:
        cloudyspec = read_cluster_cloudyspec(model_name, output_dir,
                                             fmt=fmt,
                                             verbose=verbose,
                                             read_info=read_info)
        if read_info is not None:
            read_info['cloudyspec_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudyspec = None

    # Read cloudy lines
    try:
        cloudylines \
            = read_cluster_cloudylines(model_name, output_dir, fmt=fmt,
                                       verbose=verbose, read_info=read_info)
        if read_info is not None:
            read_info['cloudylines_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudylines = None

    # Read cloudy photometry
    try:
        cloudyphot \
            = read_cluster_cloudyphot(model_name, output_dir, fmt=fmt,
                                      nofilterdata=nofilterdata,
                                      photsystem=photsystem,
                                      verbose=verbose, read_info=read_info)
        if read_info is not None:
            read_info['cloudyphot_name'] = read_info['fname']
            del read_info['fname']
    except IOError:
        cloudyphot = None

    # Build the output
    out_fields = ['id', 'trial', 'time']
    if prop is not None:
        out_data = [prop.id, prop.trial, prop.time]
    elif spec is not None:
        out_data = [spec.id, spec.trial, spec.time]
    elif phot is not None:
        out_data = [phot.id, phot.trial, phot.time]
    elif cloudyspec is not None:
        out_data = [cloudyspec.id, cloudyspec.trial, cloudyspec.time]
    elif cloudylines is not None:
        out_data = [cloudylines.id, cloudylines.trial, cloudylines.time]
    elif cloudyphot is not None:
        out_data = [cloudyphot.id, cloudyphot.trial, cloudyphot.time]
    else:
        raise IOError("unable to open any cluster files for run " +
                      model_name)
    if prop is not None:
        out_fields = out_fields + list(prop._fields[3:])
        out_data = out_data + list(prop[3:])
    if spec is not None:
        out_fields = out_fields + list(spec._fields[3:])
        out_data = out_data + list(spec[3:])
    if phot is not None:
        out_fields = out_fields + list(phot._fields[3:])
        out_data = out_data + list(phot[3:])
    if cloudyspec is not None:
        out_fields = out_fields + list(cloudyspec._fields[3:])
        out_data = out_data + list(cloudyspec[3:])
    if cloudylines is not None:
        out_fields = out_fields + list(cloudylines._fields[3:])
        out_data = out_data + list(cloudylines[3:])
    if cloudyphot is not None:
        out_fields = out_fields + list(cloudyphot._fields[3:])
        out_data = out_data + list(cloudyphot[3:])
    out_type = namedtuple('cluster_data', out_fields)
    out = out_type._make(out_data)

    # Return data
    return out

