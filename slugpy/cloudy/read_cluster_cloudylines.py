"""
Function to read a SLUG2 cluster_cloudylines file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open

def read_cluster_cloudylines(model_name, output_dir=None, fmt=None,
                             verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_cloudylines file.

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
       cloudy_linelabel : array, dtype='S4', shape (N_lines)
          labels for the lines, following cloudy's 4 character line label
          notation
       cloudy_linewl : array, shape (N_lines)
          rest wavelength for each line, in Angstrom
       cloudy_linelum : array, shape (N_cluster, N_lines)
          luminosity of each line at each time for each trial, in erg/s
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_cloudylines", 
                          output_dir=output_dir,
                          fmt=fmt)
    if read_info is not None:
        read_info['fname'] = fname

    # Print status
    if verbose:
        print("Reading cluster cloudy line luminosities for "
              "model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Prepare storage
        cluster_id = []
        time = []
        trial = []
        label = []
        wl = []
        lum = []

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read first line and store cluster data
        trialptr = 0
        entry = fp.readline()
        cluster_id.append(int(entry[0:11]))
        time.append(float(entry[14:25]))
        label.append(entry[35:39])
        wl.append(float(entry[42:53]))
        lum.append(float(entry[56:67]))
        trial.append(trialptr)

        # Read the rest of the data for first cluster
        while True:
            entry = fp.readline()

            # Check for EOF and separator lines
            if entry == '':
                break
            if entry[:3] == '---':
                trialptr = trialptr+1
                break

            # Get id, time, luminosity
            id_tmp = int(entry[0:11])
            time_tmp = float(entry[14:25])
            lum.append(float(entry[56:67]))

            # Stop when we find a different cluster or a different time
            if id_tmp != cluster_id[0] or time_tmp != time[0]:
                break

            # Still the same cluster, so append to line label and
            # wavelength lists
            label.append(entry[35:39])
            wl.append(float(entry[42:53]))

        # We have now read one full chunk, so we know how many
        # wavelength entries per cluster there are
        nl = len(wl)

        # Start of next chunk
        ptr = 1

        # Now read through rest of file
        while True:

            # Read a line
            entry = fp.readline()
            if entry == '':
                break
            if entry[:3] == '---':
                trialptr = trialptr+1
                continue
            lum.append(float(entry[56:67]))
            ptr = ptr+1

            # When we get to the end of a chunk, push cluster ID,
            # time, trial number list, then reset pointer
            if ptr == nl:
                cluster_id.append(int(entry[0:11]))
                time.append(float(entry[14:25]))
                trial.append(trialptr)
                ptr = 0

        # Convert to arrays
        wl = np.array(wl)
        label = np.array(label)
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time)
        trial = np.array(trial, dtype='uint')
        lum = np.reshape(np.array(lum), (len(cluster_id), nl))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read number of lines
        nl = int(fp.readline())

        # Read line labels
        label = []
        for i in range(nl):
            label.append(fp.readline().split()[0])

        # Read the line wavelengths
        data = fp.read(struct.calcsize('d')*nl)
        wl = np.array(struct.unpack('d'*nl, data))

        # Prepare storage
        cluster_id = []
        time = []
        trial = []
        lum = []

        # Go through the rest of the file
        trialptr = 0
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
                           struct.calcsize('d')*ncluster*nl)
            data_list = struct.unpack(('L'+'d'*nl)*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[::nl+1])
            lum.extend(
                [data_list[(nl+1)*i+1:(nl+1)*i+1+nl] 
                 for i in range(ncluster)])

        # Convert to arrays
        wl = np.array(wl)
        label = np.array(label)
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time)
        trial = np.array(trial, dtype='uint')
        lum = np.reshape(np.array(lum), (len(cluster_id), nl))

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Get the line labels and wavelengths from the first HDU
        label = fp[1].data.field('Line_Label').flatten()
        wl = fp[1].data.field('Wavelength').flatten()

        # Get time, trial, id, luminosity from second HDU
        cluster_id = fp[2].data.field('UniqueID')
        trial = fp[2].data.field('Trial')
        time = fp[2].data.field('Time')
        lum = fp[2].data.field('Line_luminosity')

    # Build namedtuple to hold output
    out_type = namedtuple('cluster_cloudylines',
                          ['id', 'trial', 'time', 'cloudy_linelist', 
                           'cloudy_linewl', 'cloudy_linelum'])
    out = out_type(cluster_id, trial, time, label, wl, lum)

    # Return
    return out
