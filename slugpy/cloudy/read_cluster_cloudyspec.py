"""
Function to read a SLUG2 cluster_cloudyspec file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open

def read_cluster_cloudyspec(model_name, output_dir=None, fmt=None,
                            verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_cloudyspec file.

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

    Raises
       IOError, if no spectrum file can be opened
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_cloudyspec", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading integrated cloudy spectra for model "+model_name)
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
        wavelength = []
        inc = []
        trans = []
        emit = []
        trans_emit = []

        # Burn the three header lines
        fp.readline()
        fp.readline()
        fp.readline()

        # Read first line and store cluster data
        trialptr = 0
        entry = fp.readline()
        data = entry.split()
        cluster_id.append(int(data[0]))
        time.append(float(data[1]))
        wavelength.append(float(data[2]))
        inc.append(float(data[3]))
        trans.append(float(data[4]))
        emit.append(float(data[5]))
        trans_emit.append(float(data[6]))
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

            # Split up data
            data = entry.split()
            id_tmp = int(data[0])
            time_tmp = float(data[1])
            inc.append(float(data[3]))
            trans.append(float(data[4]))
            emit.append(float(data[5]))
            trans_emit.append(float(data[6]))

            # Stop when we find a different cluster or a different time
            if id_tmp != cluster_id[0] or time_tmp != time[0]:
                break

            # Still the same cluster, so append to wavelength list
            wavelength.append(float(data[2]))

        # We have now read one full chunk, so we know how many
        # wavelength entries per cluster there are
        nl = len(wavelength)

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
            data = entry.split()
            inc.append(float(data[3]))
            trans.append(float(data[4]))
            emit.append(float(data[5]))
            trans_emit.append(float(data[6]))
            ptr = ptr+1

            # When we get to the end of a chunk, push cluster ID,
            # time, trial number list, then reset pointer
            if ptr == nl:
                cluster_id.append(int(data[0]))
                time.append(float(data[1]))
                trial.append(trialptr)
                ptr = 0

        # Convert to arrays
        wavelength = np.array(wavelength)
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time)
        trial = np.array(trial, dtype='uint')
        inc = np.reshape(np.array(inc), (len(cluster_id), len(wavelength)))
        trans = np.reshape(np.array(trans), (len(cluster_id), len(wavelength)))
        emit = np.reshape(np.array(emit), (len(cluster_id), len(wavelength)))
        trans_emit = np.reshape(np.array(trans_emit), 
                                (len(cluster_id), len(wavelength)))

        # Close file
        fp.close()

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # First read number of wavelengths and wavelength table
        data = fp.read(struct.calcsize('L'))
        nl, = struct.unpack('L', data)
        data = fp.read(struct.calcsize('d')*nl)
        wavelength = np.array(struct.unpack('d'*nl, data))

        # Prepare storage
        cluster_id = []
        time = []
        trial = []
        inc = []
        trans = []
        emit = []
        trans_emit = []

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
                           struct.calcsize('d')*4*ncluster*nl)
            data_list = struct.unpack(('L'+'d'*nl*4)*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[::4*nl+1])
            inc.extend(
                [data_list[(4*nl+1)*i+1:(4*nl+1)*i+1+nl] 
                 for i in range(ncluster)])
            trans.extend(
                [data_list[(4*nl+1)*i+1+nl:(4*nl+1)*i+1+2*nl] 
                 for i in range(ncluster)])
            emit.extend(
                [data_list[(4*nl+1)*i+1+2*nl:(4*nl+1)*i+1+3*nl] 
                 for i in range(ncluster)])
            trans_emit.extend(
                [data_list[(4*nl+1)*i+1+3*nl:(4*nl+1)*i+1+4*nl] 
                 for i in range(ncluster)])

        # Convert to arrays
        wavelength = np.array(wavelength)
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time)
        trial = np.array(trial, dtype='uint')
        inc = np.reshape(np.array(inc), (len(cluster_id), len(wavelength)))
        trans = np.reshape(np.array(trans), (len(cluster_id), len(wavelength)))
        emit = np.reshape(np.array(emit), (len(cluster_id), len(wavelength)))
        trans_emit = np.reshape(np.array(trans_emit), 
                                (len(cluster_id), len(wavelength)))

        # Close file
        fp.close()

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'
        wavelength = fp[1].data.field('Wavelength')
        wavelength = wavelength.flatten()
        cluster_id = fp[2].data.field('UniqueID')
        trial = fp[2].data.field('Trial')
        time = fp[2].data.field('Time')
        inc = fp[2].data.field('Incident_spectrum')
        trans = fp[2].data.field('Transmitted_spectrum')
        emit = fp[2].data.field('Emitted_spectrum')
        trans_emit = fp[2].data.field('Transmitted_plus_emitted_spectrum')

    # Build namedtuple to hold output
    out_type = namedtuple('cluster_cloudyspec',
                          ['id', 'trial', 'time', 'cloudy_wl', 
                           'cloudy_inc', 'cloudy_trans',
                           'cloudy_emit', 'cloudy_trans_emit'])
    out = out_type(cluster_id, trial, time, wavelength, inc,
                   trans, emit, trans_emit)

    # Return
    return out
