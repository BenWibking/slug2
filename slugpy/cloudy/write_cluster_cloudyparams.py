"""
This fucntion writes out the parameters that were used for a cloudy_slug run.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_cluster_cloudyparams(data, model_name, fmt):
    """
    Write out photometry for nebular emission computed by cloudy on a
    slug spectrum for a series of clusters

    Parameters
       data : namedtuple
          Cluster cloudy parameter data; a namedtuple containing the
          fields id, trial, time, cloudy_hden, cloudy_r0, 
          cloudy_r1, cloudy_QH0, cloudy_covFac, cloudy_U, cloudy_U0,
          cloudy_Omega, and cloudy_zeta; may also optionally contain
          the fields cloudy_r1_out, cloudy_hden_out, cloudy_Omega_out,
          and cloudy_zeta_out
       model_name : string
          Base file name to give the model to be written. Can include a
          directory specification if desired.
       fmt : string
          Format for the output file. Allowed values are 'ascii', 'bin'
          or 'binary, and 'fits'.

    Returns
       Nothing
    """

    # Make sure fmt is valid
    if fmt != 'ascii' and fmt != 'bin' and fmt != 'binary' and \
       fmt != 'fits':
        raise ValueError("fmt must be ascii, bin, binary, or fits")

    # Make sure we're not trying to do fits if we don't have astropy
    if fmt == 'fits' and fits is None:
        raise ValueError("Couldn't import astropy, so fits format "+
                         "is unavailable.")

    if fmt == 'ascii':

        # ASCII mode
        fp = open(model_name+'_cluster_cloudyparams.txt', 'w')

        # Write header lines
        fields = ['UniqueID', 'Time', 'Hden', 'R0', 'R1',
                  'QH0', 'CovFac', 'U', 'U0', 'Omega', 'Zeta']
        units = ['', '(yr)', '(cm^-3)', '(cm)', '(cm)',
                 '(s^-1)', '', '', '', '', '']
        if 'cloudy_hden_out' in data._fields:
            fields.append('Hden_out')
            units.append('(cm^-3)')
        if 'cloudy_r1_out' in data._fields:
            fields.append('R1_out')
            units.append('(cm)')
        if 'cloudy_Omega_out' in data._fields:
            fields.append('Omega_out')
            units.append('')
        if 'cloudy_zeta_out' in data._fields:
            fields.append('Zeta_out')
            units.append('')
        nfields = len(fields)
        for f in fields:
            fp.write("{:<14s}".format(f))
        fp.write("\n")
        for u in units:
            fp.write("{:<14s}".format(u))
        fp.write("\n")
        for i in range(nfields):
            fp.write("{:<14s}".format('-----------'))
        fp.write("\n")

        # Write data
        for i in range(data.trial.size):
            # If this is a new trial, write a separator
            if i != 0:
                if data.trial[i] != data.trial[i-1]:
                    fp.write("-"*(nfields*14-3)+"\n")
            fp.write(("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   " +
                      "{:11.5e}   {:11.5e}   {:11.5e}   {:11.5e}   " +
                      "{:11.5e}   {:11.5e}   {:11.5e}").format(
                          data.id[i], data.time[i], 
                          data.cloudy_hden[i], data.cloudy_r0[i],
                          data.cloudy_r1[i],
                          data.cloudy_QH0[i], data.cloudy_covFac[i],
                          data.cloudy_U[i], data.cloudy_U0[i],
                          data.cloudy_Omega[i],
                          data.cloudy_zeta[i]))
            if 'cloudy_hden_out' in data._fields:
                fp.write("   {:11.5e}".
                         format(data.cloudy_hden_out[i]))
            if 'cloudy_r1_out' in data._fields:
                fp.write("   {:11.5e}".
                         format(data.cloudy_r1_out[i]))
            if 'cloudy_Omega_out' in data._fields:
                fp.write("   {:11.5e}".
                         format(data.cloudy_Omega_out[i]))
            if 'cloudy_zeta_out' in data._fields:
                fp.write("   {:11.5e}".
                         format(data.cloudy_zeta_out[i]))
            fp.write("\n")
                
        # Close
        fp.close()

        
    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_cluster_cloudyparams.bin', 'wb')

        # Write out bytes indicating if we have various optional
        # fields
        if 'cloudy_hden_out' in data._fields:
            fp.write(str(bytearray([1])))
        else:
            fp.write(str(bytearray([0])))
        if 'cloudy_r1_out' in data._fields:
            fp.write(str(bytearray([1])))
        else:
            fp.write(str(bytearray([0])))
            fp.write(str(bytearray([0])))
        if 'cloudy_Omega_out' in data._fields:
            fp.write(str(bytearray([1])))
        else:
            fp.write(str(bytearray([0])))
        if 'cloudy_zeta_out' in data._fields:
            fp.write(str(bytearray([1])))
        else:
            fp.write(str(bytearray([0])))

        # Break data into blocks of clusters with the same time
        # and trial number
        ptr = 0
        while ptr < data.trial.size:

            # Find the next cluster that differs from this one in
            # either time or trial number
            diff = np.where(
                np.logical_or(data.trial[ptr+1:] != data.trial[ptr],
                              data.time[ptr+1:] != data.time[ptr]))[0]
            if diff.size == 0:
                block_end = data.trial.size
            else:
                block_end = ptr + diff[0] + 1

            # Write out time and number of clusters
            ncluster = block_end - ptr
            fp.write(np.uint(data.trial[ptr]))
            fp.write(data.time[ptr])
            fp.write(ncluster)

            # Loop over clusters and write them
            for k in range(ptr, block_end):
                fp.write(data.cloudy_hden[k])
                fp.write(data.cloudy_r0[k])
                fp.write(data.cloudy_r1[k])
                fp.write(data.cloudy_QH0[k])
                fp.write(data.cloudy_covFac[k])
                fp.write(data.cloudy_U[k])
                fp.write(data.cloudy_U0[k])
                fp.write(data.cloudy_Omega[k])
                fp.write(data.cloudy_zeta[k])
                if 'cloudy_hden_out' in data._fields:
                    fp.write(data.cloudy_hden_out[k])
                if 'cloudy_r1_out' in data._fields:
                    fp.write(data.cloudy_r1_out[k])
                if 'cloudy_Omega_out' in data._fields:
                    fp.write(data.cloudy_Omega_out[k])
                if 'cloudy_zeta_out' in data._fields:
                    fp.write(data.cloudy_zeta_out[k])

            # Move pointer
            ptr = block_end
                
        # Close file
        fp.close()
   
    elif fmt == 'fits' or fmt == 'fits2':

        # Convert data to FITS columns
        cols = []
        cols.append(fits.Column(name="Trial", format="1K",
                                unit="", array=data.trial))
        cols.append(fits.Column(name="UniqueID", format="1K",
                                unit="", array=data.id))
        cols.append(fits.Column(name="Time", format="1D",
                                unit="yr", array=data.time))
        cols.append(fits.Column(name="Hden", format="1D",
                                unit="cm^-3", array=data.cloudy_hden))
        cols.append(fits.Column(name="R0", format="1D",
                                unit="cm", array=data.cloudy_r0))
        cols.append(fits.Column(name="R1", format="1D",
                                unit="cm", array=data.cloudy_r1))
        cols.append(fits.Column(name="QH0",
                                format="1D",
                                unit="s^-1", array=data.cloudy_QH0))
        cols.append(fits.Column(name="covFac", format="1D",
                                unit="", array=data.cloudy_covFac))
        cols.append(fits.Column(name="U", format="1D",
                                unit="", array=data.cloudy_U))
        cols.append(fits.Column(name="U0", format="1D",
                                unit="", array=data.cloudy_U0))
        cols.append(fits.Column(name="Omega", format="1D",
                                unit="", array=data.cloudy_Omega))
        cols.append(fits.Column(name="zeta", format="1D",
                                unit="", array=data.cloudy_zeta))
        if 'cloudy_hden_out' in data._fields:
            cols.append(fits.Column(
                name="Hden_out", format="1D",
                unit="cm^-3", array=data.cloudy_hden_out))
        if 'cloudy_r1_out' in data._fields:
            cols.append(fits.Column(
                name="R1_out", format="1D",
                unit="cm", array=data.cloudy_r1_out))
        if 'cloudy_Omega_out' in data._fields:
            cols.append(fits.Column(
                name="Omega_out", format="1D",
                unit="", array=data.cloudy_Omega_out))
        if 'cloudy_zeta_out' in data._fields:
            cols.append(fits.Column(
                name="zeta_out", format="1D",
                unit="", array=data.cloudy_zeta_out))
        fitscols = fits.ColDefs(cols)

        # Create the binary table HDU
        tbhdu = fits.BinTableHDU.from_columns(fitscols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(model_name+'_cluster_cloudyparams.fits',
                        overwrite=True)

