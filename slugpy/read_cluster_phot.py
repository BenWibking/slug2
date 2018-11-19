"""
Function to read a SLUG2 cluster_phot file.
"""

import numpy as np
from collections import namedtuple
from copy import deepcopy
import struct
import re
import errno
from .photometry_convert import photometry_convert
from .read_filter import read_filter
from .slug_open import slug_open

def read_cluster_phot(model_name, output_dir=None, fmt=None, 
                      nofilterdata=False, photsystem=None,
                      verbose=False, read_info=None,
                      filters_only=False, read_filters=None,
                      read_nebular=None, read_extinct=None,
                      phot_only=False):
    """
    Function to read a SLUG2 cluster_phot file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the SLUG2 output is located; if set to None,
          the current directory is searched, followed by the SLUG_DIR
          directory if that environment variable is set
       fmt : 'txt' | 'ascii' | 'bin' | 'binary' | 'fits' | 'fits2'
          Format for the file to be read. If one of these is set, the
          function will only attempt to open ASCII-('txt' or 'ascii'), 
          binary ('bin' or 'binary'), or FITS ('fits' or 'fits2')
          formatted output, ending in .txt., .bin, or .fits,
          respectively. If set to None, the code will try to open
          ASCII files first, then if it fails try binary files, and if
          it fails again try FITS files.
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
       filters_only : bool
          If True, the code only reads the data on the filters, not
          any of the actual photometry. If combined with nofilterdata,
          this can be used to return the list of available filters
          and nothing else.
       read_filters : None | string | listlike containing strings
          If this is None, data on all filters is read. Otherwise only
          filters whose name(s) match the input filter names ar
          read.
       read_nebular : None | bool
          If True, only data with the nebular contribution is read; if
          False, only data without it is read. Default behavior is to
          read all data.
       read_extinct : None | bool
          If True, only data with extinction applied is read; if
          False, only data without it is read. Default behavior is to
          read all data.
       phot_only : bool
          If true, id, trial, time, and filter information are not
          read, only photometry

    Returns
       A namedtuple, which can contain the following fields depending
       on the input options, and depending on which fields are present
       in the file being read:

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          times at which cluster spectra are output, in yr
       filter_names : list of string
          a list giving the name for each filter
       filter_units : list of string
          a list giving the units for each filter
       filter_wl_eff : list
          effective wavelength of each filter; this is set to None for the
          filters Lbol, QH0, QHe0, and QHe1; omitted if nofilterdata is
          True
       filter_wl : list of arrays
          a list giving the wavelength table for each filter; this is
          None for the filters Lbol, QH0, QHe0, and QHe1
       filter_response : list of arrays
          a list giving the photon response function for each filter;
          this is None for the filters Lbol, QH0, QHe0, and QHe1
       filter_beta : list
          powerlaw index beta for each filter; used to normalize the
          photometry
       filter_wl_c : list
          pivot wavelength for each filter; used to normalize the photometry
       phot : array, shape (N_cluster, N_filter)
          photometric value in each filter for each cluster; units are as
          indicated in the units field
       phot_neb : array, shape (N_filter, N_times, N_trials)
          same as phot, but for the light after it has passed through
          the HII region
       phot_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but after extinction has been applied
       phot_neb_ex : array, shape (N_filter, N_times, N_trials)
          same as phot, but for the light after it has passed through
          the HII region and then had extinction applied
       
    Raises
       IOError, if no photometry file can be opened
       ValueError, if photsystem is set to an unknown value
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_phot", 
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster photometry for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False

    # Prepare holders for data
    if not phot_only:
        cluster_id = []
        time = []
        trial = []
    phot = []

    # Read data
    if fname.endswith('.txt'):

        ########################################################
        # ASCII mode
        ########################################################

        if read_info is not None:
            read_info['format'] = 'ascii'

        # If this is a checkpoint file, skip the line stating how many
        # trials it contains; this line is not guaranteed to be
        # accurate, and is intended for the C++ code, not for us
        if checkpoint:
            fp.readline()
            
        # Read the list of filters
        line = fp.readline()
        filters = line.split()[2:]
        nfilter = len(filters)

        # Read the list of units
        line = fp.readline()
        line = line.replace(')', '(').split('(') # split by ( and )
        units = []
        for l in line:
            if (not l.isspace()) and (len(l) > 0):
                units.append(l)
        units = units[1:]    # Get rid of the units for time

        # Search for filters with names that end in _n, _ex, or _nex,
        # indicating that they include the effects of the nebula,
        # extinction, or both
        neb = []
        ex = []
        for i in range(nfilter):
            if len(filters[i]) > 2:
                if filters[i][-2:] == '_n':
                    neb.append(i)
            if len(filters[i]) > 3:
                if filters[i][-3:] == '_ex':
                    ex.append(i)
        if len(neb) > 0:
            nebular = True
            phot_neb = []
        else:
            nebular = False
        if len(ex) > 0:
            extinct = True
            phot_ex = []
        else:
            extinct = False
        if nebular and extinct:
            phot_neb_ex = []

        # If we have nebular emission or extinction, reshape filter
        # and units lists
        nuniq = nfilter // ((1+nebular)*(1+extinct))
        nfilter = nuniq
        filters = filters[:nfilter]
        units = units[:nfilter]

        # If given a list of filters to read, make sure that we
        # haven't been given ones that are not available; if we have,
        # raise an error now
        if read_filters is not None:
            if type(read_filters) is not str and \
               type(read_filters) is not np.str_:
                for f in read_filters:
                    if not f in filters:
                        raise IOError(errno.EIO,
                                      "requested filter {:s} not available!".
                                      format(f))
            else:
                if not read_filters in filters:
                    raise IOError(errno.EIO, 
                                  "requested filter {:s} not available!".
                                  format(read_filters))

        # If only reading filter data, don't read any further
        if not filters_only:

            # If requested to read only certain filters, figure out
            # which ones are in our list to read
            if read_filters is not None:
                fread = []
                if type(read_filters) is not str and \
                   type(read_filters) is not np.str_:
                    for f in filters:
                        if f in read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)
                else:
                    for f in filters:
                        if f == read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)

            # Burn a line
            line = fp.readline()

            # Read through data
            trialptr = 0
            for line in fp:
                if line[:3] == '---':
                    trialptr = trialptr+1
                    continue
                linesplit = line.split()
                if not phot_only:
                    cluster_id.append(int(linesplit[0]))
                    time.append(float(linesplit[1]))
                if read_nebular is not True and read_extinct is not True:
                    if read_filters is None:
                        phot.append(np.array(linesplit[2:2+nfilter],
                                             dtype='float'))
                    else:
                        tmp = []
                        for i, j in enumerate(range(2,2+nfilter)):
                            if fread[i]:
                                tmp.append(linesplit[j])
                        phot.append(np.array(tmp, dtype='float'))
                if nebular and read_nebular is not False and \
                   read_extinct is not True:
                    if read_filters is None:
                        phot_neb.append(np.array(
                            linesplit[nfilter+2:2*nfilter+2],
                            dtype='float'))
                    else:
                        tmp = []
                        for i, j in enumerate(range(nfilter+2,2+2*nfilter)):
                            if fread[i]:
                                tmp.append(linesplit[j])
                        phot_neb.append(np.array(tmp, dtype='float'))
                if extinct and read_nebular is not True and \
                   read_extinct is not False:
                    tmp_ex = []
                    for i in range(nfilter):
                        substr = line[21*(2+(1+nebular)*nfilter+i):
                                      21*(2+(1+nebular)*nfilter+i+1)]
                        if read_filters is None:
                            if substr.isspace():
                                tmp_ex.append(np.nan)
                            else:
                                tmp_ex.append(float(substr))
                        else:
                            if fread[i]:
                                if substr.isspace():
                                    tmp_ex.append(np.nan)
                                else:
                                    tmp_ex.append(float(substr))
                    phot_ex.append(np.array(tmp_ex, dtype='float'))
                if nebular and extinct and read_nebular is not False \
                   and read_extinct is not False:
                    tmp_neb_ex = []
                    for i in range(nfilter):
                        substr = line[21*(2+3*nfilter+i):
                                      21*(2+3*nfilter+i+1)]
                        if read_filters is None:
                            if substr.isspace():
                                tmp_neb_ex.append(np.nan)
                            else:
                                tmp_neb_ex.append(float(substr))
                        else:
                            if fread[i]:
                                if substr.isspace():
                                    tmp_neb_ex.append(np.nan)
                                else:
                                    tmp_neb_ex.append(float(substr))
                    phot_neb_ex.append(np.array(tmp_neb_ex, dtype='float'))
                if not phot_only:
                    trial.append(trialptr)

    elif fname.endswith('.bin'):

        ########################################################
        # Binary mode
        ########################################################

        if read_info is not None:
            read_info['format'] = 'binary'

        # If this is a checkpoint, skip the bytes specifying how many
        # trials we have; this is inteded for the C++ code, not for
        # us, since we will determine that on our own
        if checkpoint:
            data = fp.read(struct.calcsize('i'))
        
        # Read number of filters
        nfilter = int(fp.readline())

        # Read filter names and units
        filters = []
        units = []
        for i in range(nfilter):
            line = fp.readline()
            filters.append(line.split()[0].decode("utf-8"))
            units.append(line.split()[1].decode("utf-8"))

        # If given a list of filters to read, make sure that we
        # haven't been given ones that are not available; if we have,
        # raise an error now
        if read_filters is not None:
            if type(read_filters) is not str and \
               type(read_filters) is not np.str_:
                for f in read_filters:
                    if not f in filters:
                        raise IOError(errno.EIO,
                                      "requested filter {:s} not available!".
                                      format(f))
            else:
                if not read_filters in filters:
                    raise IOError(errno.EIO,
                                  "requested filter {:s} not available!".
                                  format(read_filters))

        # Read the bits that tells us if we're using nebular emission
        # and extinction
        data = fp.read(struct.calcsize('b'))
        nebular = struct.unpack('b', data)[0] != 0
        data = fp.read(struct.calcsize('b'))
        extinct = struct.unpack('b', data)[0] != 0
        if nebular:
            phot_neb = []
        if extinct:
            phot_ex = []
        if nebular and extinct:
            phot_neb_ex = []
        nftot = nfilter*(1+nebular)*(1+extinct)

        # If only reading the filters, skip the rest of this
        if not filters_only:

            # If requested to read only certain filters, figure out
            # which ones are in our list to read
            if read_filters is not None:
                fread = []
                if type(read_filters) is not str and \
                   type(read_filters) is not np.str_:
                    for f in filters:
                        if f in read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)
                else:
                    for f in filters:
                        if f == read_filters:
                            fread.append(True)
                        else:
                            fread.append(False)

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
                if not phot_only:
                    time.extend([t]*ncluster)
                    trial.extend([trialptr]*ncluster)

                # Handle things differently if we're reading
                # everything versus if we're reading just some filters
                if read_filters is None:

                    # Read the next block of clusters
                    chunkstr = ('L'+'d'*nftot)*ncluster
                    data = fp.read(struct.calcsize(chunkstr))
                    data_list = struct.unpack(chunkstr, data)

                    # Pack clusters into data list
                    if not phot_only:
                        cluster_id.extend(data_list[::nftot+1])
                    if read_nebular is not True and \
                       read_extinct is not True:
                        phot.extend(
                            [data_list[(nftot+1)*i+1:
                                       (nftot+1)*i+1+nfilter] 
                             for i in range(ncluster)])
                    ptr = 1
                    if nebular:
                        if read_nebular is not False and \
                           read_extinct is not True:
                            phot_neb.extend(
                                [data_list[(nftot+1)*i+1+ptr*nfilter:
                                           (nftot+1)*i+1+(ptr+1)*nfilter] 
                                 for i in range(ncluster)])
                        ptr = ptr+1
                    if extinct:
                        if read_nebular is not True and \
                           read_extinct is not False:
                            phot_ex.extend(
                                [data_list[(nftot+1)*i+1+ptr*nfilter:
                                           (nftot+1)*i+1+(ptr+1)*nfilter] 
                                 for i in range(ncluster)])
                        ptr = ptr+1
                    if nebular and extinct:
                        if read_nebular is not False and \
                           read_extinct is not False:
                            phot_neb_ex.extend(
                                [data_list[(nftot+1)*i+1+ptr*nfilter:
                                           (nftot+1)*i+1+(ptr+1)*nfilter] 
                                 for i in range(ncluster)])
                        ptr = ptr+1

                else:

                    # Loop over clusters
                    for i in range(ncluster):

                        # Read cluster ID
                        if not phot_only:
                            data = fp.read(struct.calcsize('L'))
                            data_list = struct.unpack('L', data)
                            cluster_id.extend(data_list)
                        else:
                            fp.seek(struct.calcsize('L'), 1)

                        # Loop over filters
                        for j in range(nfilter):
                            if fread[j] and \
                               read_nebular is not True and \
                               read_extinct is not True:
                                data = fp.read(struct.calcsize('d'))
                                data_list = struct.unpack('d', data)
                                phot.extend(data_list)
                            else:
                                fp.seek(struct.calcsize('d'), 1)

                        # Repeat for nebular emission
                        if nebular:
                            for j in range(nfilter):
                                if fread[j] and \
                                   read_nebular is not False and \
                                   read_extinct is not True:
                                    data = fp.read(struct.calcsize('d'))
                                    data_list = struct.unpack('d', data)
                                    phot_neb.extend(data_list)
                                else:
                                    fp.seek(struct.calcsize('d'), 1)

                        # Repeat for extincted emission
                        if extinct:
                            for j in range(nfilter):
                                if fread[j] and \
                                   read_nebular is not True and \
                                   read_extinct is not False:
                                    data = fp.read(struct.calcsize('d'))
                                    data_list = struct.unpack('d', data)
                                    phot_ex.extend(data_list)
                                else:
                                    fp.seek(struct.calcsize('d'), 1)

                        # Repeat for nebular plus extincted emission
                        if nebular and extinct:
                            for j in range(nfilter):
                                if fread[j] and \
                                   read_nebular is not False and \
                                   read_extinct is not False:
                                    data = fp.read(struct.calcsize('d'))
                                    data_list = struct.unpack('d', data)
                                    phot_neb_ex.extend(data_list)
                                else:
                                    fp.seek(struct.calcsize('d'), 1)

    elif fname.endswith('.fits'):

        ########################################################
        # FITS mode
        ########################################################

        if read_info is not None:
            read_info['format'] = 'fits'

        # Figure out which of the two fits formats we're using; if
        # there's > 2 HDUs, we're using fits2, otherwise we're using
        # fits1
        if len(fp) == 2:

            ########################################################
            # FITS1 format: all the data is stored in a single HDU
            ########################################################

            # Get filter names and units
            filters = []
            units = []
            i = 4
            while 'TTYPE'+str(i) in fp[1].header.keys():
                filters.append(fp[1].header['TTYPE'+str(i)])
                units.append(fp[1].header['TUNIT'+str(i)])
                i = i+1
            nfilter = len(filters)

            # If given a list of filters to read, make sure that we
            # haven't been given ones that are not available; if we have,
            # raise an error now
            if read_filters is not None:
                if type(read_filters) is not str and \
                   type(read_filters) is not np.str_:
                    for f in read_filters:
                        if not f in filters:
                            raise IOError(errno.EIO,
                                          "requested filter {:s} not available!".
                                          format(f))
                else:
                    if not read_filters in filters:
                        raise IOError(errno.EIO, 
                                      "requested filter {:s} not available!".
                                      format(read_filters))

            # Search for filters with names that end in _neb, _ex, or
            # _neb_ex, indicating that they include the effects of the
            # nebula, extinction, or both
            neb = []
            ex = []
            for i in range(nfilter):
                if len(filters[i]) > 4:
                    if filters[i][-4:] == '_neb':
                        neb.append(i)
                if len(filters[i]) > 3:
                    if filters[i][-3:] == '_ex':
                        ex.append(i)
            if len(neb) > 0:
                nebular = True
            else:
                nebular = False
            if len(ex) > 0:
                extinct = True
            else:
                extinct = False

            # If we have nebular emission or extinction, reshape
            # filter and units lists
            nuniq = nfilter // ((1+nebular)*(1+extinct))
            nfilter = nuniq
            filters = filters[:nfilter]
            units = units[:nfilter]

            # If only reading filters, skip the rest
            if not filters_only:

                # If requested to read only certain filters, figure
                # out which ones are in our list to read
                if read_filters is not None:
                    fread = []
                    nf_final = 0
                    if type(read_filters) is not str and \
                       type(read_filters) is not np.str_:
                        for f in filters:
                            if f in read_filters:
                                fread.append(True)
                                nf_final = nf_final+1
                            else:
                                fread.append(False)
                    else:
                        for f in filters:
                            if f == read_filters:
                                fread.append(True)
                                nf_final = nf_final+1
                            else:
                                fread.append(False)
                else:
                    fread = [True] * len(filters)
                    nf_final = len(filters)

                # Get cluster ID, trial, time
                if not phot_only:
                    cluster_id = fp[1].data.field('UniqueID')
                    trial = fp[1].data.field('Trial')
                    time = fp[1].data.field('Time')

                # Get photometric data
                if read_nebular is not True and \
                   read_extinct is not True:
                    phot = np.zeros((fp[1].header['NAXIS2'], nf_final))
                if nebular and \
                   read_nebular is not False and \
                   read_extinct is not True:
                    phot_neb = np.zeros((fp[1].header['NAXIS2'], nf_final))
                if extinct and \
                   read_nebular is not True and \
                   read_extinct is not False:
                    phot_ex = np.zeros((fp[1].header['NAXIS2'], nf_final))
                if nebular and extinct and \
                   read_nebular is not False and \
                   read_extinct is not False:
                    phot_neb_ex = np.zeros((fp[1].header['NAXIS2'], nf_final))
                ptr = 0
                for i in range(len(filters)):
                    if fread[i]:
                        if read_nebular is not True and \
                           read_extinct is not True:
                            phot[:,ptr] = fp[1].data.field(filters[i])
                        if nebular and \
                           read_nebular is not False and \
                           read_extinct is not True:
                            phot_neb[:,ptr] = fp[1].data.field(
                                filters[i]+"_neb")
                        if extinct and \
                           read_nebular is not True and \
                           read_extinct is not False:
                            phot_ex[:,ptr] = fp[1].data.field(
                                filters[i]+"_ex")
                        if nebular and extinct and \
                           read_nebular is not False and \
                           read_extinct is not False:
                            phot_neb_ex[:,ptr] = fp[1].data. \
                                                 field(filters[i]+"_neb_ex")
                        ptr = ptr+1

        else:

            ########################################################
            # FITS2 format: each filter in its own HDU
            ########################################################

            # Get filter names and units
            filters = [fp[i].header['TTYPE1'] for i in 
                       range(2,len(fp))]
            units = [fp[i].header['TUNIT1'] for i in 
                       range(2,len(fp))]
            nfilter = len(filters)

            # If given a list of filters to read, make sure that we
            # haven't been given ones that are not available; if we have,
            # raise an error now
            if read_filters is not None:
                if type(read_filters) is not str and \
                   type(read_filters) is not np.str_:
                    for f in read_filters:
                        if not f in filters:
                            raise IOError(errno.EIO,
                                          "requested filter {:s} not available!".
                                          format(f))
                else:
                    if not read_filters in filters:
                        raise IOError(errno.EIO,
                                      "requested filter {:s} not available!".
                                      format(read_filters))

            # Search for filters with names that end in _neb, _ex, or
            # _neb_ex, indicating that they include the effects of the
            # nebula, extinction, or both
            neb = []
            ex = []
            for i in range(nfilter):
                if len(filters[i]) > 4:
                    if filters[i][-4:] == '_neb':
                        neb.append(i)
                if len(filters[i]) > 3:
                    if filters[i][-3:] == '_ex':
                        ex.append(i)
            if len(neb) > 0:
                nebular = True
            else:
                nebular = False
            if len(ex) > 0:
                extinct = True
            else:
                extinct = False

            # If we have nebular emission or extinction, reshape
            # filter and units lists
            nuniq = nfilter // ((1+nebular)*(1+extinct))
            nfilter = nuniq
            filters = filters[:nfilter]
            units = units[:nfilter]

            # If only reading filters, skip the rest
            if not filters_only:

               # If requested to read only certain filters, figure
                # out which ones are in our list to read
                if read_filters is not None:
                    fread = []
                    nf_final = 0
                    if type(read_filters) is not str and \
                       type(read_filters) is not np.str_:
                        for f in filters:
                            if f in read_filters:
                                fread.append(True)
                                nf_final = nf_final+1
                            else:
                                fread.append(False)
                    else:
                        for f in filters:
                            if f == read_filters:
                                fread.append(True)
                                nf_final = nf_final+1
                            else:
                                fread.append(False)
                else:
                    fread = [True] * len(filters)
                    nf_final = len(filters)

                # Get cluster ID, trial, time
                if not phot_only:
                    cluster_id = fp[1].data.field('UniqueID')
                    trial = fp[1].data.field('Trial')
                    time = fp[1].data.field('Time')

                # Get photometric data
                if read_nebular is not True and \
                   read_extinct is not True:
                    phot = np.zeros((fp[1].header['NAXIS2'], nf_final))
                if nebular and \
                   read_nebular is not False and \
                   read_extinct is not True:
                    phot_neb = np.zeros((fp[1].header['NAXIS2'], nf_final))
                if extinct and \
                   read_nebular is not True and \
                   read_extinct is not False:
                    phot_ex = np.zeros((fp[1].header['NAXIS2'], nf_final))
                if nebular and extinct and \
                   read_nebular is not False and \
                   read_extinct is not False:
                    phot_neb_ex = np.zeros((fp[1].header['NAXIS2'], nf_final))
                ptr1 = 0
                ptr2 = 2
                for i in range(len(filters)):
                    if fread[i] and \
                       read_nebular is not True and \
                       read_extinct is not True:
                        phot[:,ptr1] \
                            = fp[ptr2].data.field(filters[i])
                        ptr1 = ptr1 + 1
                    ptr2 = ptr2+1
                ptr1 = 0
                if nebular:
                    for i in range(len(filters)):
                        if fread[i] and \
                           read_nebular is not False and \
                           read_extinct is not True:
                            phot_neb[:,ptr1] \
                                = fp[ptr2].data.field(filters[i]+'_neb')
                            ptr1 = ptr1 + 1
                        ptr2 = ptr2+1
                ptr1 = 0
                if extinct:
                    for i in range(len(filters)):
                        if fread[i] and \
                           read_nebular is not True and \
                           read_extinct is not False:
                            phot_ex[:,ptr1] \
                                = fp[ptr2].data.field(filters[i]+'_ex')
                            ptr1 = ptr1 + 1
                        ptr2 = ptr2+1
                ptr1 = 0
                if nebular and extinct:
                    for i in range(len(filters)):
                        if fread[i] and \
                           read_nebular is not False and \
                           read_extinct is not False:
                            phot_neb_ex[:,ptr1] \
                                = fp[ptr2].data.field(filters[i]+'_neb_ex')
                            ptr1 = ptr1 + 1
                        ptr2 = ptr2+1
                ptr1 = 0

    ########################################################
    # End of file reading block
    ########################################################

    # Close file
    fp.close()

    # If using only a subset of filters, truncate the filter list now
    if read_filters is not None:
        filters_tmp = []
        units_tmp = []
        for i in range(len(filters)):
            if fread[i]:
                filters_tmp.append(filters[i])
                units_tmp.append(units[i])
        filters = filters_tmp
        units = units_tmp

    # Convert to arrays
    if not filters_only:
        if not phot_only:
            cluster_id = np.array(cluster_id, dtype='uint')
            time = np.array(time, dtype='float')
            trial = np.array(trial, dtype='uint')
        if read_nebular is not True and \
           read_extinct is not True:
            phot = np.array(phot, dtype='float')
            phot = np.reshape(phot, (phot.size//len(filters), 
                                     len(filters)))
        if nebular and \
           read_nebular is not False and \
           read_extinct is not True:
            phot_neb = np.array(phot_neb, dtype='float')
            phot_neb = np.reshape(phot_neb, (phot_neb.size//len(filters),
                                             len(filters)))
        if extinct and \
           read_nebular is not True and \
           read_extinct is not False:
            phot_ex = np.array(phot_ex, dtype='float')
            phot_ex = np.reshape(phot_ex, (phot_ex.size//len(filters), 
                                           len(filters)))
        if nebular and extinct and \
           read_nebular is not False and \
           read_extinct is not False:
            phot_neb_ex = np.array(phot_neb_ex, dtype='float')
            phot_neb_ex = np.reshape(phot_neb_ex, 
                                     (phot_neb_ex.size//len(filters), 
                                      len(filters)))

    # Read filter data if requested
    if not nofilterdata:
        if verbose:
            print("Reading filter data")
        wl_eff, wavelength, response, beta, wl_c = read_filter(filters)

    # Do photometric system conversion if requested
    if photsystem is not None and not filters_only:
        if verbose:
            print("Converting photometric system")
        if nofilterdata:
            units_save = deepcopy(units)
            if read_nebular is not True and \
               read_extinct is not True:
                photometry_convert(photsystem, phot, units, filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and \
               read_nebular is not False and \
               read_extinct is not True:
                photometry_convert(photsystem, phot_neb, units, 
                                   filter_last=True)            
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if extinct and \
               read_nebular is not True and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_ex, units, 
                                   filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and extinct and \
               read_nebular is not False and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_neb_ex, units, 
                                   filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            units = units_out
        else:
            units_save = deepcopy(units)
            if read_nebular is not True and \
               read_extinct is not True:
                photometry_convert(photsystem, phot, units, wl_eff, 
                                   filter_names=filters, filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and \
               read_nebular is not False and \
               read_extinct is not True:
                photometry_convert(photsystem, phot_neb, units, 
                                   wl_eff, filter_names=filters,
                                   filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if extinct and \
               read_nebular is not True and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_ex, units, wl_eff, 
                                   filter_names=filters,
                                   filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            if nebular and extinct and \
               read_nebular is not False and \
               read_extinct is not False:
                photometry_convert(photsystem, phot_neb_ex, units, 
                                   wl_eff, filter_names=filters,
                                   filter_last=True)
                units_out = deepcopy(units)
                units = deepcopy(units_save)
            units = units_out

    # Construct return object
    if filters_only or phot_only:
        fieldnames = ['filter_names', 'filter_units']
        fields = [filters, units]
    else:
        fieldnames = ['id', 'trial', 'time', 'filter_names', 'filter_units']
        fields = [cluster_id, trial, time, filters, units]
    if not nofilterdata:
        fieldnames = fieldnames + ['filter_wl_eff', 'filter_wl',
                                   'filter_response', 'filter_beta',
                                   'filter_wl_c']
        fields = fields + [wl_eff, wavelength, response, beta, wl_c]
    if not filters_only and read_nebular is not True and \
       read_extinct is not True:
        fieldnames = fieldnames + ['phot']
        fields = fields + [phot]
    if nebular and not filters_only and \
       read_nebular is not False and \
       read_extinct is not True:
        fieldnames = fieldnames + ['phot_neb']
        fields = fields + [phot_neb]
    if extinct and not filters_only and \
       read_nebular is not True and \
       read_extinct is not False:
        fieldnames = fieldnames + ['phot_ex']
        fields = fields + [phot_ex]
    if nebular and extinct and not filters_only and \
       read_nebular is not False and \
       read_extinct is not False:
            fieldnames = fieldnames + ['phot_neb_ex']
            fields = fields + [phot_neb_ex]
    out_type = namedtuple('cluster_phot', fieldnames)
    out = out_type(*fields)

    # Return
    return out
