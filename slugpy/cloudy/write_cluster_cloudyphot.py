"""
This function writes out photometry computed on a cloudy output
nebular spectrum for a series of star clusters.
"""

import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_cluster_cloudyphot(data, model_name, fmt):
    """
    Write out photometry for nebular emission computed by cloudy on a
    slug spectrum for a series of clusters

    Parameters
       data : namedtuple
          Cluster cloudy photometry data to be written; a namedtuple
          containing the fields id, time, cloudy_filter_names, 
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
        fp = open(model_name+'_cluster_cloudyphot.txt', 'w')

        # Write header lines
        fp.write("{:<18s}{:<18s}".format('UniqueID', 'Time'))
        for f in data.cloudy_filter_names:
            fp.write("{:<18s}{:<18s}{:<18s}".format(f,f,f))
        fp.write("\n")
        fp.write("{:<18s}{:<18s}".format("", ""))
        for f in data.cloudy_filter_names:
            fp.write("{:<18s}{:<18s}{:<18s}".
                     format("Trans", "Emit", "Trans+Emit"))
        fp.write("\n")
        fp.write("{:<18s}{:<18s}".format('', '(yr)'))
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
        for i in range(data.cloudy_phot_trans.shape[0]):
            # If this is a new trial, write a separator
            if i != 0:
                if data.trial[i] != data.trial[i-1]:
                    fp.write("-"*((2+3*nf)*18-3)+"\n")
            fp.write("    {:11d}       {:11.5e}"
                     .format(data.id[i], data.time[i]))
            for j in range(nf):
                fp.write("       {:11.5e}".
                         format(data.cloudy_phot_trans[i,j]))
                fp.write("       {:11.5e}".
                         format(data.cloudy_phot_emit[i,j]))
                fp.write("       {:11.5e}".
                         format(data.cloudy_phot_trans_emit[i,j]))
            fp.write("\n")

        # Close
        fp.close()

    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_cluster_cloudyphot.bin', 'wb')

        # Write number of filters and filter names as ASCII
        nf = len(data.cloudy_filter_names)
        fp.write(str(nf)+"\n")
        for i in range(nf):
            fp.write(data.cloudy_filter_names[i] + " " + 
                     data.cloudy_filter_units[i] + "\n")

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
                fp.write(data.id[k])
                fp.write(data.cloudy_phot_trans[k,:])
                fp.write(data.cloudy_phot_emit[k,:])
                fp.write(data.cloudy_phot_trans_emit[k,:])

            # Move pointer
            ptr = block_end

        # Close file
        fp.close()

    elif fmt == 'fits':

        # FITS mode

        # Convert data to FITS columns
        cols = []
        cols.append(fits.Column(name="Trial", format="1K",
                                unit="", array=data.trial))
        cols.append(fits.Column(name="UniqueID", format="1K",
                                unit="", array=data.id))
        cols.append(fits.Column(name="Time", format="1D",
                                unit="yr", array=data.time))
        for i in range(len(data.cloudy_filter_names)):
            cols.append(
                fits.Column(name=data.cloudy_filter_names[i]+'_Transmitted',
                            unit=data.cloudy_filter_units[i],
                            format="1D",
                            array=data.cloudy_phot_trans[:,i]))
            cols.append(
                fits.Column(name=data.cloudy_filter_names[i]+'_Emitted',
                            unit=data.cloudy_filter_units[i],
                            format="1D",
                            array=data.cloudy_phot_emit[:,i]))
            cols.append(
                fits.Column(name=data.cloudy_filter_names[i]+
                            '_Transmitted_plus_emitted',
                            unit=data.cloudy_filter_units[i],
                            format="1D",
                            array=data.cloudy_phot_trans_emit[:,i]))
        fitscols = fits.ColDefs(cols)

        # Create the binary table HDU
        tbhdu = fits.BinTableHDU.from_columns(fitscols)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, tbhdu])
        hdulist.writeto(model_name+'_cluster_cloudyphot.fits',
                        overwrite=True)
