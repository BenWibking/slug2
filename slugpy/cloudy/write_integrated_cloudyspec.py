"""
This function writes out integrated spectra computed by cloudy from a slug
run.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_integrated_cloudyspec(data, model_name, fmt):
    """
    Write out data computed by cloudy on a slug spectrum

    Parameters
       data : namedtuple
          Integrated cloudy spectral data to be written; a namedtuple
          containing the field time, cloudy_wl, cloudy_inc, cloudy_trans,
          cloudy_emit, and cloudy_trans_emit
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
        fp = open(model_name+'_integrated_cloudyspec.txt', 'w')

        # Write header lines
        fp.write(("{:<14s}"*6).
                 format('Time', 'Wavelength', 'Incident', 'Transmitted',
                        'Emitted', 'Trans+Emit') + "\n")
        fp.write(("{:<14s}"*6).
                 format('(yr)', '(Angstrom)', '(erg/s/A)',
                        '(erg/s/A)', '(erg/s/A)', '(erg/s/A)') + "\n")
        fp.write(("{:<14s}"*6).
                 format('-----------', '-----------', '-----------',
                        '-----------', '-----------', '-----------')
                 + "\n")

        # Write data
        ntime = data.cloudy_inc.shape[-2]
        ntrial = data.cloudy_inc.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        nl = len(data.cloudy_wl)
        for i in range(ntrial):
            if i != 0:
                fp.write("-"*(6*14-3)+"\n")
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                for k in range(nl):
                    fp.write(("{:11.5e}   {:11.5e}   {:11.5e}   " +
                              "{:11.5e}   {:11.5e}   {:11.5e}\n")
                             .format(t_out, data.cloudy_wl[k],
                                     data.cloudy_inc[k,j,i],
                                     data.cloudy_trans[k,j,i],
                                     data.cloudy_emit[k,j,i],
                                     data.cloudy_trans_emit[k,j,i]))

        # Close
        fp.close()

    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_integrated_cloudyspec.bin', 'wb')

        # Write out wavelength data
        fp.write(np.int64(len(data.cloudy_wl)))
        fp.write(data.cloudy_wl)

        # Write out times and spectra
        ntime = data.cloudy_inc.shape[-2]
        ntrial = data.cloudy_inc.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        for i in range(ntrial):
            for j in range(ntime):
                fp.write(np.uint(i))
                if random_time:
                    fp.write(data.time[i])
                else:
                    fp.write(data.time[j])
                # Note: need to put data in contiguous block before
                # writing
                tmp = np.copy(data.cloudy_inc[:,j,i])
                fp.write(tmp)
                tmp = np.copy(data.cloudy_trans[:,j,i])
                fp.write(tmp)
                tmp = np.copy(data.cloudy_emit[:,j,i])
                fp.write(tmp)
                tmp = np.copy(data.cloudy_trans_emit[:,j,i])
                fp.write(tmp)

        # Close file
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Convert wavelength data to FITS columns and make an HDU
        # from it; complication: astropy expects the dimensions of
        # the array to be (n_entries, n_wavelengths)
        nl = data.cloudy_wl.shape[0]
        fmtstring = str(nl)+"D"
        wlcols = [fits.Column(name="Wavelength",
                              format=fmtstring,
                              unit="Angstrom", 
                              array=data.cloudy_wl.reshape(1,nl))]
        wlfits = fits.ColDefs(wlcols)
        wlhdu = fits.BinTableHDU.from_columns(wlcols)

        # Figure out number of trials, and tile arrays
        ntrial = data.cloudy_inc.shape[-1]
        ntimes = data.cloudy_inc.shape[-2]
        trial = np.transpose(np.tile(
            np.arange(ntrial, dtype='int64'), (ntimes,1))).\
            flatten()
        if len(data.time) > ntimes:
            times = data.time
        else:
            times = np.tile(data.time, ntrial)

        # Convert spectra to FITS columns, and make an HDU from them
        speccols = []
        speccols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=trial))
        speccols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=times))
        speccols.append(fits.Column(name="Incident_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=np.transpose(data.cloudy_inc).
                                    reshape(ntimes*ntrial, nl)))
        speccols.append(fits.Column(name="Transmitted_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=np.transpose(data.cloudy_trans).
                                    reshape(ntimes*ntrial, nl)))
        speccols.append(fits.Column(name="Emitted_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=np.transpose(data.cloudy_emit).
                                    reshape(ntimes*ntrial, nl)))
        speccols.append(fits.Column(name="Transmitted_plus_emitted_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=np.transpose(data.cloudy_trans_emit).
                                    reshape(ntimes*ntrial, nl)))
        specfits = fits.ColDefs(speccols)
        spechdu = fits.BinTableHDU.from_columns(specfits)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
        hdulist.writeto(model_name+'_integrated_cloudyspec.fits', 
                        overwrite=True)
