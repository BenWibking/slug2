"""
Function to read a SLUG2 cluster_cloudyparams file.
"""

import numpy as np
from collections import namedtuple
import struct
from ..slug_open import slug_open

def read_cluster_cloudyparams(model_name, output_dir=None, fmt=None,
                              verbose=False, read_info=None):
    """
    Function to read a SLUG2 cluster_cloudyparams file.

    Parameters
       model_name : string
          The name of the model to be read
       output_dir : string
          The directory where the output is located; if set to None,
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

       id : array, dtype uint
          unique ID of cluster
       trial: array, dtype uint
          which trial was this cluster part of
       time : array
          time at which cluster's properties are being evaluated
       cloudy_hden : array
          number density of H nuclei at the inner edge of the ionized
          region simulated by cloudy
       cloudy_r0 : array
          inner radius of the ionized region simulated by cloudy
       cloudy_r1 : array
          outer radius of the ionized region simulated by cloudy (approximate!)
       cloudy_QH0 : array
          ionizing luminosity used in the cloudy computation
       cloudy_covFac : array
          covering factor assumed in the cloudy computation; only a
          fraction covFac of the ionizing photons are assumed to
          produce emission within the HII region, while the remainder
          are assumed to escape
       cloudy_U : array
          volume-averaged ionization parameter of the HII region
          simulated by cloudy
       cloudy_U0 : array
          ionization parameter of the HII reegion at its inner edge
       cloudy_Omega : array
          Yeh & Matzner (2012) wind parameter for the HII region
          simulated by cloudy
       cloudy_zeta : array
          Krumholz & Matzner (2009) radiation pressure parameter for
          the HII region, again approximate; values of zeta >~1
          indicate that radiation pressure is dynamically important

       The following fields may or may not be present, depending on
       what is found in the output file:
       
       cloudy_hden_out : array
          volume-averaged number density produced by the cloudy
          calculation
       cloudy_r1_out : array
          HII region outer radius produced by cloudy
       cloudy_Omega_out : array
          value of Omega computed using cloudy_r1_out
       cloudy_zeta_out : array
          value of zeta computed using cloudy_r1_out

    Notes

       The relationships between U, U0, Omega, r0, r1, hden, and
       QH0 used in deriving various parameters are valid only in
       the limit of negligible radiation pressure. They may be
       significantly off if radiation pressure is significant,
       i.e., if zeta >~ 1.
    """

    # Open file
    fp, fname = slug_open(model_name+"_cluster_cloudyparams",
                          output_dir=output_dir,
                          fmt=fmt)

    # Print status
    if verbose:
        print("Reading cluster cloudy parameters for model "+model_name)
    if read_info is not None:
        read_info['fname'] = fname

    # Read data
    if fname.endswith('.txt'):

        # ASCII mode
        if read_info is not None:
            read_info['format'] = 'ascii'

        # Read the first header line and see what fields are present
        line = fp.readline()
        if 'Hden_out' in line: hden_out_set = True
        else: hden_out_set = False
        if 'R1_out' in line: r1_out_set = True
        else: r1_out_set = False
        if 'Omega_out' in line: Omega_out_set = True
        else: Omega_out_set = False
        if 'Zeta_out' in line: zeta_out_set = True
        else: zeta_out_set = False
                
        # Burn two header lines
        line = fp.readline()
        line = fp.readline()

        # Prepare data holders
        cluster_id = []
        time = []
        trial = []
        hden = []
        r0 = []
        r1 = []
        QH0 = []
        covFac = []
        U = []
        U0 = []
        Omega = []
        zeta = []
        if hden_out_set: hden_out = []
        if r1_out_set: r1_out = []
        if Omega_out_set: Omega_out = []
        if zeta_out_set: zeta_out = []

        # Read data
        trialptr = 0
        for line in fp:
            if line[:3] == '---':
                trialptr = trialptr+1
                continue
            linesplit = line.split()
            cluster_id.append(long(linesplit[0]))
            time.append(float(linesplit[1]))
            hden.append(float(linesplit[2]))
            r0.append(float(linesplit[3]))
            r1.append(float(linesplit[4]))
            QH0.append(float(linesplit[5]))
            covFac.append(float(linesplit[6]))
            U.append(float(linesplit[7]))
            U0.append(float(linesplit[8]))
            Omega.append(float(linesplit[9]))
            zeta.append(float(linesplit[10]))
            ptr = 11
            if hden_out_set:
                hden_out.append(float(linesplit[ptr]))
                ptr += 1
            if r1_out_set:
                r1_out.append(float(linesplit[ptr]))
                ptr += 1
            if Omega_out_set:
                Omega_out.append(float(linesplit[ptr]))
                ptr += 1
            if zeta_out_set:
                zeta_out.append(float(linesplit[ptr]))
                ptr += 1
            trial.append(trialptr)

        # Convert to arrays
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time, dtype='float')
        trial = np.array(trial, dtype='uint')
        hden = np.array(hden, dtype='float')
        r0 = np.array(r0, dtype='float')
        r1 = np.array(r1, dtype='float')
        QH0 = np.array(QH0, dtype='float')
        covFac = np.array(covFac, dtype='float')
        U = np.array(U, dtype='float')
        U0 = np.array(U0, dtype='float')
        Omega = np.array(Omega, dtype='float')
        if hden_out_set: hden_out = np.array(hden_out, dtype='float')
        if r1_out_set: r1_out = np.array(r1_out, dtype='float')
        if Omega_out_set: Omega_out = np.array(Omega_out, dtype='float')
        if zeta_out_set: zeta_out = np.array(zeta_out, dtype='float')

    elif fname.endswith('.bin'):

        # Binary mode
        if read_info is not None:
            read_info['format'] = 'binary'

        # Read leading bits to see which fields we have
        data = fp.read(struct.calcsize('bbbb'))
        field_bits = struct.unpack('bbbb', data)
        hden_out_set = field_bits[0] != 0
        r1_out_set = field_bits[1] != 0
        Omega_out_set = field_bits[2] != 0
        zeta_out_set = field_bits[3] != 0
        
        # Prepare storage
        datastr = 'ddddddddd'
        nfield = len(datastr)
        cluster_id = []
        time = []
        trial = []
        hden = []
        r0 = []
        r1 = []
        QH0 = []
        covFac = []
        U = []
        U0 = []
        Omega = []
        zeta = []
        if hden_out_set:
            hden_out = []
            nfield += 1
            datastr += 'd'
        if r1_out_set:
            r1_out = []
            nfield += 1
            datastr += 'd'
        if Omega_out_set:
            Omega_out = []
            nfield += 1
            datastr += 'd'
        if zeta_out_set:
            zeta_out = []
            nfield += 1
            datastr += 'd'

        # Go through file
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

            # Read the next block of clusters
            data = fp.read(struct.calcsize(datastr*ncluster))
            data_list = struct.unpack(datastr*ncluster, data)

            # Pack clusters into data list
            cluster_id.extend(data_list[0::nfield])
            time.extend([t]*ncluster)
            trial.extend([trialptr]*ncluster)
            hden.extend(data_list[0::nfield])
            r0.extend(data_list[1::nfield])
            r1.extend(data_list[2::nfield])
            QH0.extend(data_list[3::nfield])
            covFac.extend(data_list[4::nfield])
            U.extend(data_list[5::nfield])
            U0.extend(data_list[6::nfield])
            Omega.extend(data_list[7::nfield])
            zeta.extend(data_list[8::nfield])
            ptr = 9
            if hden_out_set:
                hden_out.extend(data_list[ptr::nfield])
                ptr += 1
            if r1_out_set:
                r1_out.extend(data_list[ptr::nfield])
                ptr += 1
            if Omega_out_set:
                Omega_out.extend(data_list[ptr::nfield])
                ptr += 1
            if zeta_out_set:
                zeta_out.extend(data_list[ptr::nfield])
                ptr += 1

        # Convert to arrays
        cluster_id = np.array(cluster_id, dtype='uint')
        time = np.array(time, dtype='float')
        trial = np.array(trial, dtype='uint')
        hden = np.array(hden, dtype='float')
        r0 = np.array(r0, dtype='float')
        r1 = np.array(r1, dtype='float')
        QH0 = np.array(QH0, dtype='float')
        covFac = np.array(covFac, dtype='float')
        U = np.array(U, dtype='float')
        U0 = np.array(U0, dtype='float')
        Omega = np.array(Omega, dtype='float')
        zeta = np.array(zeta, dtype='float')
        if hden_out_set: hden_out = np.array(hden_out, dtype='float')
        if r1_out_set: r1_out = np.array(r1_out, dtype='float')
        if Omega_out_set: Omega_out = np.array(Omega_out, dtype='float')
        if zeta_out_set: zeta_out = np.array(zeta_out, dtype='float')
        
    elif fname.endswith('.fits'):

        # FITS mode
        if read_info is not None:
            read_info['format'] = 'fits'

        # Read data
        cluster_id = fp[1].data.field('UniqueID')
        trial = fp[1].data.field('Trial')
        time = fp[1].data.field('Time')
        hden = fp[1].data.field('Hden')
        r0 = fp[1].data.field('R0')
        try:
            # Compatibility for some older outputs, which didn't write R1
            r1 = fp[1].data.field('R1')
        except KeyError:
            r1 = None
        QH0 = fp[1].data.field('QH0')
        covFac = fp[1].data.field('covFac')
        U = fp[1].data.field('U')
        try:
            # Compatibility for some older outputs, which didn't write U0
            U0 = fp[1].data.field('U0')
        except KeyError:
            U0 = None
        Omega = fp[1].data.field('Omega')
        try:
            # Compatibility for some older outputs, which didn't write U0
            zeta = fp[1].data.field('zeta')
        except KeyError:
            zeta = None
        try:
            hden_out = fp[1].data.field('Hden_out')
            hden_out_set = True
        except KeyError:
            hden_out_set = False
        try:
            r1_out = fp[1].data.field('R1_out')
            r1_out_set = True
        except KeyError:
            r1_out_set = False
        try:
            Omega_out = fp[1].data.field('Omega_out')
            Omega_out_set = True
        except KeyError:
            Omega_out_set = False
        try:
            zeta_out = fp[1].data.field('zeta_out')
            zeta_out_set = True
        except KeyError:
            zeta_out_set = False
                         
    # Close file
    fp.close()

    # Build namedtuple to hold output
    fields = ['id', 'trial', 'time',
              'cloudy_hden', 'cloudy_r0', 'cloudy_r1',
              'cloudy_QH0', 'cloudy_covFac',
              'cloudy_U', 'cloudy_U0', 'cloudy_Omega', 'cloudy_zeta']
    data = [cluster_id, trial, time, hden, r0, r1, QH0, covFac, U,
            U0, Omega, zeta]
    if hden_out_set:
        fields.append('cloudy_hden_out')
        data.append(hden_out)
    if r1_out_set:
        fields.append('cloudy_r1_out')
        data.append(r1_out)
    if Omega_out_set:
        fields.append('cloudy_Omega_out')
        data.append(Omega_out)
    if zeta_out_set:
        fields.append('cloudy_zeta_out')
        data.append(zeta_out)
    cloudyparams_type = namedtuple('cluster_cloudyparams', fields)
    out = cloudyparams_type(*data)

    # Return
    return out
