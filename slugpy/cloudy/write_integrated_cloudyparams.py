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

def write_integrated_cloudyparams(data, model_name, fmt):
    """
    Write out photometry for nebular emission computed by cloudy on a
    slug spectrum for a series of clusters

    Parameters
       data : namedtuple
          Cluster cloudy parameter data; a namedtuple containing the
          fields time, cloudy_hden, cloudy_r0, cloudy_r1, 
          cloudy_QH0, cloudy_covFac, cloudy_U, cloudy_U0, cloudy_Omega, and
          cloudy_zeta; may also optionally contain the fields
          cloudy_r1_out, cloudy_hden_out, cloudy_Omega_out, and
          cloudy_zeta_out
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
        fp = open(model_name+'_integrated_cloudyparams.txt', 'w')

        # Write header lines
        fields = ['Time', 'Hden', 'R0', 'R1',
                  'QH0', 'CovFac', 'U', 'U0', 'Omega', 'Zeta']
        units = ['(yr)', '(cm^-3)', '(cm)', '(cm)',
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
        ntime = data.cloudy_hden.shape[-2]
        ntrial = data.cloudy_hden.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        for i in range(ntrial):
            if i != 0:
                fp.write("-"*(nfields*14-3)+"\n")
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                fp.write(("{:11.5e}   {:11.5e}   {:11.5e}   " +
                          "{:11.5e}   {:11.5e}   {:11.5e}   " +
                          "{:11.5e}   {:11.5e}   {:11.5e}   " +
                          "{:11.5e}")
                         .format(t_out, data.cloudy_hden[j,i],
                                 data.cloudy_r0[j,i],
                                 data.cloudy_r1[j,i],
                                 data.cloudy_QH0[j,i],
                                 data.cloudy_covFac[j,i],
                                 data.cloudy_U[j,i],
                                 data.cloudy_U0[j,i],
                                 data.cloudy_Omega[j,i],
                                 data.cloudy_zeta[j,i]))
                if 'cloudy_hden_out' in data._fields:
                    fp.write("   {:11.5e}".
                             format(data.cloudy_hden_out[j,i]))
                if 'cloudy_r1_out' in data._fields:
                    fp.write("   {:11.5e}".
                             format(data.cloudy_r1_out[j,i]))
                if 'cloudy_Omega_out' in data._fields:
                    fp.write("   {:11.5e}".
                             format(data.cloudy_Omega_out[j,i]))
                if 'cloudy_zeta_out' in data._fields:
                    fp.write("   {:11.5e}".
                             format(data.cloudy_zeta_out[j,i]))
            fp.write("\n")

        # Close
        fp.close()
        
    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_integrated_cloudyparams.bin', 'wb')

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

        # Write data
        ntime = data.cloudy_hden.shape[-2]
        ntrial = data.cloudy_hden.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        for i in range(ntrial):
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                fp.write(np.uint(i))
                fp.write(t_out)
                fp.write(data.cloudy_hden[j,i])
                fp.write(data.cloudy_r0[j,i])
                fp.write(data.cloudy_r1[j,i])
                fp.write(data.cloudy_QH0[j,i])
                fp.write(data.cloudy_covFac[j,i])
                fp.write(data.cloudy_U[j,i])
                fp.write(data.cloudy_U0[j,i])
                fp.write(data.cloudy_Omega[j,i])
                fp.write(data.cloudy_zeta[j,i])
                if 'cloudy_hden_out' in data._fields:
                    fp.write(data.cloudy_hden_out[j,i])
                if 'cloudy_r1_out' in data._fields:
                    fp.write(data.cloudy_r1_out[j,i])
                if 'cloudy_Omega_out' in data._fields:
                    fp.write(data.cloudy_Omega_out[j,i])
                if 'cloudy_zeta_out' in data._fields:
                    fp.write(data.cloudy_zeta_out[j,i])

        # Close
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Figure out number of trials, and tile arrays
        ntrial = data.cloudy_hden.shape[-1]
        ntimes = data.cloudy_hden.shape[-2]
        trial = np.transpose(np.tile(
            np.arange(ntrial, dtype='int64'), (ntimes,1))).\
            flatten()
        if len(data.time) > ntimes:
            times = data.time
        else:
            times = np.tile(data.time, ntrial)

        # Convert data to FITS columns
        cols = []
        cols.append(fits.Column(name="Trial", format="1K",
                                unit="", array=trial))
        cols.append(fits.Column(name="Time", format="1D",
                                unit="yr", array=times))
        cols.append(fits.Column(name="HDen", format="1D",
                                unit="cm^-3",
                                array=np.transpose(data.cloudy_hden).
                                flatten()))
        cols.append(fits.Column(name="R0", format="1D",
                                unit="cm",
                                array=np.transpose(data.cloudy_r0).
                                flatten()))
        cols.append(fits.Column(name="R1", format="1D",
                                unit="cm",
                                array=np.transpose(data.cloudy_r1).
                                flatten()))
        cols.append(fits.Column(name="QH0", format="1D",
                                unit="s^-1",
                                array=np.transpose(data.cloudy_QH0).
                                flatten()))
        cols.append(fits.Column(name="covFac", format="1D",
                                array=np.transpose(data.cloudy_covFac).
                                flatten()))
        cols.append(fits.Column(name="U", format="1D",
                                array=np.transpose(data.cloudy_U).
                                flatten()))
        cols.append(fits.Column(name="U0", format="1D",
                                array=np.transpose(data.cloudy_U0).
                                flatten()))
        cols.append(fits.Column(name="Omega", format="1D",
                                array=np.transpose(data.cloudy_Omega).
                                flatten()))
        cols.append(fits.Column(name="zeta", format="1D",
                                array=np.transpose(data.cloudy_zeta).
                                flatten()))
        if 'cloudy_hden_out' in data._fields:
            cols.append(fits.Column(
                name="Hden_out", format="1D",
                unit="cm^-3", array=np.transpose(
                    data.cloudy_hden_out).flatten()))
        if 'cloudy_r1_out' in data._fields:
            cols.append(fits.Column(
                name="R1_out", format="1D",
                unit="cm", array=np.transpose(
                    data.cloudy_r1_out).flatten()))
        if 'cloudy_Omega_out' in data._fields:
            cols.append(fits.Column(
                name="Omega_out", format="1D",
                unit="", array=np.transpose(
                    data.cloudy_Omega_out).flatten()))
        if 'cloudy_zeta_out' in data._fields:
            cols.append(fits.Column(
                name="zeta_out", format="1D",
                unit="", array=
                np.transpose(data.cloudy_zeta_out).flatten()))
        fitscols = fits.ColDefs(cols)

        # Create the binary table HDU
        tbhdu = fits.BinTableHDU.from_columns(fitscols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(model_name+'_integrated_cloudyparams.fits',
                        overwrite=True)
