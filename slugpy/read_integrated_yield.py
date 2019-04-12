"""
Function to read a SLUG2 integrated_yield file.
"""

import numpy as np
from collections import namedtuple
import struct
import re
from .slug_open import slug_open

def read_integrated_yield(model_name, output_dir=None, fmt=None,
                          read_info=None, verbose=False):

    """
    Function to read a SLUG2 integrated_yield file.

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
       verbose : bool
          If True, verbose output is printed as code runs
       read_info : dict
          On return, this dict will contain the keys 'fname' and
          'format', giving the name of the file read and the format it
          was in; 'format' will be one of 'ascii', 'binary', or 'fits'

    Returns
       A namedtuple containing the following fields:

       time : array, shape (N_times) or shape (N_trials)
          Times at which data are output; shape is either N_times (if
          the run was done with fixed output times) or N_trials (if
          the run was done with random output times)
       isotope_name : array of strings, shape (N_iso)
          Atomic symbols of the isotopes included in the yield table
       isotope_Z : array of int, shape (N_iso)
          Atomic numbers of the isotopes included in the yield table
       isotope_A : array of int, shape (N_iso)
          Atomic mass number of the isotopes included in the yield table
       yld : array, shape (N_iso, N_times) or (N_iso, N_trials)
          Yield of each isotope, defined as the instantaneous amount
          produced up to that time; for unstable isotopes, this
          includes the effects of decay since production
    """
    
    # Open file
    fp, fname = slug_open(model_name+"_integrated_yield", 
                          output_dir=output_dir,
                          fmt=fmt)

    # See if this file is a checkpoint file
    if len(re.findall('_chk\d\d\d\d', model_name)) != 0:
        checkpoint = True
    else:
        checkpoint = False
        
    # Print status and record
    if verbose:
        print("Reading integrated yield for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # If this is a checkpoint file, skip the line stating how many
        # trials it contains; this line is not guaranteed to be
        # accurate, and is intended for the C++ code, not for us
        if checkpoint:
            fp.readline()
            
        # Burn 3 header lines
        hdr = fp.readline()
        hdr = fp.readline()
        hdr = fp.readline()

        # Read data
        time = []
        isotope_name = []
        isotope_Z = []
        isotope_A = []
        yld = []
        yldtmp = []
        first_trial = True
        for entry in fp:
            # See if this is a new trial
            if entry[:3] == '---':
                yld.append(yldtmp)
                yldtmp = []
                first_trial = False
                continue       # Skip separator lines
            # Read data
            data = entry.split()
            time.append(float(data[0]))
            if first_trial:
                isotope_name.append(data[1].title())
                isotope_Z.append(float(data[2]))
                isotope_A.append(float(data[3]))
            yldtmp.append(float(data[4]))

        # Append last set of yields
        yld.append(yldtmp)

        # Truncate repeats in the isotope list that correspond to the
        # same set of isotopes at different times, and get unique times
        if len(isotope_name) > 1:
            isotopes = zip(isotope_name, isotope_Z, isotope_A)
            if isotopes[0] in isotopes[1:]:
                niso = isotopes[1:].index(isotopes[0])+1
                isotope_name = isotope_name[:niso]
                isotope_Z = isotope_Z[:niso]
                isotope_A = isotope_A[:niso]
                time = time[::niso]

        # Convert to arrays
        isotope_name = np.array(isotope_name)
        isotope_Z = np.array(isotope_Z, dtype=int)
        isotope_A = np.array(isotope_A, dtype=int)
        time = np.array(time)

        # If times are repeats, truncate
        if time.size > len(yld):
            ntime = time.size//len(yld)
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]

        # Make yield array
        yld = np.array(yld)
        yld = yld.reshape((yld.shape[0], time.size,
                           isotope_name.size))

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # If this is a checkpoint, skip the bytes specifying how many
        # trials we have; this is inteded for the C++ code, not for
        # us, since we will determine that on our own
        if checkpoint:
            data = fp.read(struct.calcsize('i'))
        
        # Read number of isotopes
        niso = struct.unpack('L', fp.read(struct.calcsize('L')))[0]

        # Read isotope data
        data = fp.read(struct.calcsize(('c'*4+'II')*niso))
        data_list = struct.unpack(('c'*4+'II')*niso, data)
        isotope_name = np.array(
            [ (data_list[6*i]+data_list[6*i+1]+
               data_list[6*i+2]+data_list[6*i+3]).strip().
              title() for i in range(niso) ])
        isotope_Z = np.array(data_list[4::6], dtype=int)
        isotope_A = np.array(data_list[5::6], dtype=int)

        # Now read remainder of file
        buf = fp.read()

        # Parse
        blockstr = 'Ld'+niso*'d'
        nout = len(buf) // struct.calcsize(blockstr)
        data_list = struct.unpack(blockstr*nout, buf)
        trial = np.array(data_list[::2+niso], dtype=int)
        time = np.array(data_list[1::2+niso], dtype=float)
        yld = np.array([np.array(data_list[2+(2+niso)*i:
                                           (2+niso)*(i+1)]) for
                        i in range(nout)])

        # Reformat
        idx = np.argmax(trial != trial[0])
        if idx > 1:
            if np.amin(time[:idx] == time[idx:2*idx]):
                time = time[:idx]
            yld = yld.reshape((len(trial)//idx, time.size,
                               isotope_name.size))
        else:
            yld = yld.reshape((1, time.size,
                               isotope_name.size))

    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Read data
        isotope_name = fp[1].data['Name']
        isotope_name = np.array([iso.strip().title() for iso in isotope_name])
        isotope_Z = fp[1].data['Z']
        isotope_A = fp[1].data['A']
        trial = fp[2].data['Trial']
        time = fp[2].data['Time']
        yld = fp[2].data['Yield']

        # Re-arrange data into desired shape
        ntrial = len(np.unique(trial))
        ntime = len(time)//ntrial
        if ntime != len(time):
            if np.amin(time[:ntime] == time[ntime:2*ntime]):
                time = time[:ntime]
        yld = yld.reshape((ntrial, ntime, isotope_name.size))

    # Close file
    fp.close()

    # Build output holder
    fieldnames = ['time', 'isotope_name', 'isotope_Z', 'isotope_A', 'yld']
    fields = [ time, isotope_name, isotope_Z, isotope_A, np.transpose(yld)]
    out_type = namedtuple('integrated_yield', fieldnames)
    out = out_type(*fields)

    # Return
    return out


