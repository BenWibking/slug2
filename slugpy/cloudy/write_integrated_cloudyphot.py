"""
This function writes out photometry computed on a cloudy output
nebular spectrum on the integrated light from a galaxy.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_integrated_cloudyphot(data, model_name, fmt):
    """
    Write out photometry for nebular emission computed by cloudy on a
    slug spectrum

    Parameters
       data : namedtuple
          Integrated cloudy photometry data to be written; a namedtuple
          containing the fields time, cloudy_filter_names, 
          cloudy_filter_units, cloudy_phot_trans, cloudy_phot_emit,
          and cloudy_phot_trans_emit
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
        fp = open(model_name+'_integrated_cloudyphot.txt', 'w')

        # Get shape of output
        ntime = data.cloudy_phot_trans.shape[1]
        ntrial = data.cloudy_phot_trans.shape[2]
        if len(data.time) > ntime:
            random_time = True
        else:
            random_time = False

        # Write header lines
        fp.write("{:<18s}".format('Time'))
        for f in data.cloudy_filter_names:
            fp.write("{:<18s}{:<18s}{:<18s}".format(f,f,f))
        fp.write("\n")
        fp.write("{:<18s}".format(""))
        for f in data.cloudy_filter_names:
            fp.write("{:<18s}{:<18s}{:<18s}".
                     format("Trans", "Emit", "Trans+Emit"))
        fp.write("\n")
        fp.write("{:<18s}".format('(yr)'))
        for f in data.cloudy_filter_units:
            for i in range(3):
                fp.write("({:s}".format(f)+")"+" "*(16-len(f)))
        fp.write("\n")
        fp.write("{:<18s}".format('---------------'))
        nf = len(data.cloudy_filter_names)
        for j in range(3):
            for i in range(nf):
                fp.write("{:<18s}".format('---------------'))
        fp.write("\n")

        # Write data
        for i in range(ntrial):
            # Write separator between trials
            if i != 0:
                fp.write("-"*((1+3*nf)*18-3)+"\n")
            for j in range(ntime):
                if random_time:
                    t_out = data.time[i]
                else:
                    t_out = data.time[j]
                fp.write("    {:11.5e}".format(t_out))
                for k in range(nf):
                    fp.write("       {:11.5e}".
                             format(data.cloudy_phot_trans[k,j,i]))
                    fp.write("       {:11.5e}".
                             format(data.cloudy_phot_emit[k,j,i]))
                    fp.write("       {:11.5e}".
                             format(data.cloudy_phot_trans_emit[k,j,i]))
                fp.write("\n")

        # Close
        fp.close()

    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_integrated_cloudyphot.bin', 'wb')

        # Write number of filters and filter names as ASCII
        nf = len(data.cloudy_filter_names)
        fp.write(str(nf)+"\n")
        for i in range(nf):
            fp.write(data.cloudy_filter_names[i] + " " + 
                     data.cloudy_filter_units[i] + "\n")

        # Write data
        ntime = data.cloudy_phot_trans.shape[1]
        ntrial = data.cloudy_phot_trans.shape[2]
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
                tmp = np.copy(data.cloudy_phot_trans[:,j,i])
                fp.write(tmp)
                tmp = np.copy(data.cloudy_phot_emit[:,j,i])
                fp.write(tmp)
                tmp = np.copy(data.cloudy_phot_trans_emit[:,j,i])
                fp.write(tmp)

        # Close file
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Figure out number of trials, and tile arrays
        ntimes = data.cloudy_phot_trans.shape[1]
        ntrial = data.cloudy_phot_trans.shape[2]
        trial = np.transpose(np.tile(
            np.arange(ntrial, dtype='int64'), (ntimes,1))).\
            flatten()
        if len(data.time) > ntimes:
            times = data.time
        else:
            times = np.tile(data.time, ntrial)
        nf = len(data.cloudy_filter_names)

        # Convert data to FITS columns
        cols = []
        cols.append(fits.Column(name="Trial", format="1K",
                                unit="", array=trial))
        cols.append(fits.Column(name="Time", format="1D",
                                unit="yr", array=times))
        for i in range(nf):
            cols.append(
                fits.Column(name=data.cloudy_filter_names[i]+'_Transmitted',
                            unit=data.cloudy_filter_units[i],
                            format="1D",
                            array=np.transpose(data.cloudy_phot_trans[i,:,:]).
                            flatten()))
            cols.append(
                fits.Column(name=data.cloudy_filter_names[i]+'_Emitted',
                            unit=data.cloudy_filter_units[i],
                            format="1D",
                            array=np.transpose(data.cloudy_phot_emit[i,:,:]).
                            flatten()))
            cols.append(
                fits.Column(name=data.cloudy_filter_names[i]+
                            '_Transmitted_plus_emitted',
                            unit=data.cloudy_filter_units[i],
                            format="1D",
                            array=np.transpose(data.
                                               cloudy_phot_trans_emit[i,:,:]).
                            flatten()))
        fitscols = fits.ColDefs(cols)

        # Create the binary table HDU
        tbhdu = fits.BinTableHDU.from_columns(fitscols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(model_name+'_integrated_cloudyphot.fits',
                        overwrite=True)
