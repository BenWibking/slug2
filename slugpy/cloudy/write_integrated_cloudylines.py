"""
This function writes out a list of line luminosities computed by
cloudy from a slug run.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_integrated_cloudylines(data, model_name, fmt):
    """
    Write out line luminosities computed by cloudy on a slug spectrum

    Parameters
       data : namedtuple
          Integrated cloudy line data to be written; a namedtuple
          containing the fields time, cloudy_linelist, cloudy_linewl, 
          cloudy_linelum
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
        fp = open(model_name+'_integrated_cloudylines.txt', 'w')

        # Write header lines
        fp.write(("{:<14s}"*4).
                 format('Time', 'LineLabel', 'Wavelength', 
                        'Luminosity') + "\n")
        fp.write(("{:<14s}"*4).
                 format('(yr)', '', '(Angstrom)', '(erg/s)') + "\n")
        fp.write(("{:<14s}"*4).
                 format('-----------', '-----------', '-----------',
                        '-----------')
                 + "\n")

        # Write data
        ntime = data.cloudy_linelum.shape[-2]
        ntrial = data.cloudy_linelum.shape[-1]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False
        nline = len(data.cloudy_linewl)
        for i in range(ntrial):
            if i != 0:
                fp.write("-"*(4*14-3)+"\n")
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                for k in range(nline):
                    fp.write(("{:11.5e}   {:>11s}   {:11.5e}   " +
                              "{:11.5e}\n")
                             .format(t_out, 
                                     data.cloudy_linelabel[k],
                                     data.cloudy_linewl[k],
                                     data.cloudy_linelum[k,j,i]))

        # Close
        fp.close()

    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_integrated_cloudylines.bin', 'wb')

        # Write out number of lines and line labels as ASCII, one per
        # line
        nlabel = len(data.cloudy_linelabel)
        fp.write(str(nlabel)+"\n")
        for i in range(nlabel):
            fp.write(data.cloudy_linelabel[i]+"\n")

        # Write line wavelengths
        fp.write(data.cloudy_linewl)

        # Write line luminosities
        ntime = data.cloudy_linelum.shape[1]
        ntrial = data.cloudy_linelum.shape[2]
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
                # This next line is needed to put the data into a
                # contiguous block before writing
                tmp = np.copy(data.cloudy_linelum[:,j,i])
                fp.write(tmp)

        # Close file
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Create a first HDU containing the line wavelengths and labels
        nl = len(data.cloudy_linewl)
        fmtstring = "A4"
        wlcols = [fits.Column(name="Line_Label",
                              format=fmtstring,
                              array=data.cloudy_linelabel)]
        fmtstring = "D"
        wlcols.append(fits.Column(name="Wavelength",
                                  format=fmtstring,
                                  unit="Angstrom", 
                                  array=data.cloudy_linewl))
        wlfits = fits.ColDefs(wlcols)
        wlhdu = fits.BinTableHDU.from_columns(wlcols)

        # Figure out number of trials, and tile arrays
        ntimes = data.cloudy_linelum.shape[1]
        ntrial = data.cloudy_linelum.shape[2]
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
        cols.append(fits.Column(name="Line_Luminosity",
                                format=str(nl)+"D",
                                unit="erg/s", 
                                array=np.transpose(data.cloudy_linelum).
                                reshape(ntimes*ntrial, nl)))
        linelum_cols = fits.ColDefs(cols)
        linelumhdu = fits.BinTableHDU.from_columns(linelum_cols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, wlhdu, linelumhdu])
        hdulist.writeto(model_name+'_integrated_cloudylines.fits', 
                        overwrite=True)
 
