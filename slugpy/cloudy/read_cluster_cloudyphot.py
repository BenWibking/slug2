"""
Function to read a SLUG2 cluster_cloudyphot file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open
from ..photometry_convert import photometry_convert
from ..read_filter import read_filter

def read_cluster_cloudyphot(model_name, output_dir=None, fmt=None,
                               nofilterdata=False, photsystem=None,
                               verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_cloudyphot file.

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
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
       A namedtuple containing the following fields:

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          times at which cluster spectra are output, in yr
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
       IOError, if no photometry file can be opened;
       ValueError, if photsystem is set to an unknown value
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_cloudyphot", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster cloudy photometry for model " +
              model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Read the list of filters
        line = fp.readline()
        filters = line.split()[2::3]
        nfilter = len(filters)

        # Burn a line
        line = fp.readline()

        # Read the list of units
        line = fp.readline()
        line = line.replace(')', '(').split('(') # split by ( and )
        units = []
        for l in line:
            if (not l.isspace()) and (len(l) > 0):
                units.append(l)
        units = units[1::3]   # Get rid of the units for time

        # Burn a line
        line = fp.readline()

        # Prepare holders for data
        cluster_id = []
        time = []
        trial = []
        phot_trans = []
        phot_emit = []
        phot_trans_emit = []

        # Read through data
        trialptr = 0
        for line in fp:
            if line[:3] == '---':
                trialptr = trialptr+1
                continue
            linesplit = line.split()
            cluster_id.append(int(linesplit[0]))
            time.append(float(linesplit[1]))
            phot_trans.append(linesplit[2::3])
            phot_emit.append(linesplit[3::3])
            phot_trans_emit.append(linesplit[4::3])
            trial.append(trialptr)

        # Convert to arrays
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time, dtype='float')
        trial = np.array(trial, dtype='uint')
        phot_trans = np.reshape(np.array(phot_trans, dtype='float'),
                                (len(time), len(filters)))
        phot_emit = np.reshape(np.array(phot_emit, dtype='float'),
                                (len(time), len(filters)))
        phot_trans_emit = np.reshape(np.array(phot_trans_emit, 
                                              dtype='float'),
                                     (len(time), len(filters)))

        # Close file
        fp.close()

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read number of filters
        nfilter = int(fp.readline())

        # Read filter names and units
        filters = []
        units = []
        for i in range(nfilter):
            line = fp.readline()
            filters.append(line.split()[0])
            units.append(line.split()[1])

        # Prepare holders for data
        cluster_id = []
        time = []
        trial = []
        phot_trans = []
        phot_emit = []
        phot_trans_emit = []

        # Go through the rest of the file
        while True:

            # Read number of clusters and time in next block, checking
            # if we've hit eof
            data = fp.read(struct.calcsize('LdL'))
            if len(data) < struct.calcsize('LdL'):
                break
            trialptr, t, ncluster = struct.unpack('LdL', data)

            # Skip if no clusters
            if ncluster == 0:
                continue

            # Add to time and trial arrays
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)

            # Read the next block of clusters
            data = fp.read(struct.calcsize('L')*ncluster + 
                           struct.calcsize('d')*ncluster*nfilter*3)
            data_list = struct.unpack(('L'+'d'*nfilter*3)*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[::3*nfilter+1])
            phot_trans.extend(
                [data_list[(nfilter+1)*i+1:(nfilter+1)*i+1+nfilter] 
                 for i in range(ncluster)])
            phot_emit.extend(
                [data_list[(nfilter+1)*i+1+nfilter:(nfilter+1)*i+1+2*nfilter] 
                 for i in range(ncluster)])
            phot_trans_emit.extend(
                [data_list[(nfilter+1)*i+1+2*nfilter:(nfilter+1)*i+1+3*nfilter] 
                 for i in range(ncluster)])

        # Convert to arrays
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time, dtype='float')
        trial = np.array(trial, dtype='uint')
        phot_trans = np.reshape(np.array(phot_trans, dtype='float'),
                                (len(time), len(filters)))
        phot_emit = np.reshape(np.array(phot_emit, dtype='float'),
                                (len(time), len(filters)))
        phot_trans_emit = np.reshape(np.array(phot_trans_emit, 
                                              dtype='float'),
                                     (len(time), len(filters)))

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Get cluster ID, trial, time
        cluster_id = fp[1].data.field('UniqueID')
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')

        # Get filter names and units
        filters = []
        units = []
        i = 4
        while 'TTYPE'+str(i) in fp[1].header.keys():
            filters.append(fp[1].header['TTYPE'+str(i)][:-12])
            units.append(fp[1].header['TUNIT'+str(i)])
            i = i+3

        # Get photometric data
        nfilter = len(filters)
        ntime = len(time)
        phot_trans = np.zeros((ntime, nfilter))
        phot_emit = np.zeros((ntime, nfilter))
        phot_trans_emit = np.zeros((ntime, nfilter))
        for i in range(len(filters)):
            phot_trans[:,i] = fp[1].data.field(filters[i]+'_Transmitted')
            phot_emit[:,i] = fp[1].data.field(filters[i]+'_Emitted')
            phot_trans_emit[:,i] \
                = fp[1].data.field(filters[i] + 
                                   '_Transmitted_plus_emitted')

    # Read filter data if requested
    if not nofilterdata:
        if verbose:
            print("Reading filter data")
        wl_eff, wavelength, response, beta, wl_c = read_filter(filters)

    # Do photometric system conversion if requested
    if photsystem is not None:
        if verbose:
            print("Converting photometric system")
        if nofilterdata:
            photometry_convert(photsystem, phot_trans, units,
                               filter_last=True)
            photometry_convert(photsystem, phot_emit, units,
                               filter_last=True)
            photometry_convert(photsystem, phot_trans_emit, units, 
                               filter_last=True)
        else:
            photometry_convert(photsystem, phot_trans, units, wl_eff,
                               filter_last=True)
            photometry_convert(photsystem, phot_emit, units, wl_eff,
                               filter_last=True)
            photometry_convert(photsystem, phot_trans_emit, units, wl_eff,
                               filter_last=True)

    # Construct return object
    if nofilterdata:
        out_type = namedtuple('cluster_phot',
                              ['id', 'trial', 'time',
                               'cloudy_filter_names', 
                               'cloudy_filter_units', 
                               'cloudy_phot_trans',
                               'cloudy_phot_emit',
                               'cloudy_phot_trans_emit'])
        out = out_type(cluster_id, trial, time, filters, units, 
                       phot_trans, phot_emit, phot_trans_emit)
    else:
        out_type = namedtuple('cluster_phot',
                              ['id', 'trial', 'time', 
                               'cloudy_filter_names', 
                               'cloudy_filter_units',
                               'cloudy_filter_wl_eff',
                               'cloudy_filter_wl',
                               'cloudy_filter_response',
                               'cloudy_filter_beta', 
                               'cloudy_filter_wl_c', 
                               'cloudy_phot_trans',
                               'cloudy_phot_emit',
                               'cloudy_phot_trans_emit'])
        out = out_type(cluster_id, trial, time, filters, units, wl_eff,
                       wavelength, response, beta, wl_c, 
                       phot_trans, phot_emit, phot_trans_emit)

    # Return
    return out
