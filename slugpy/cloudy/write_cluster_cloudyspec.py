"""
This function writes out cluster spectra computed by cloudy from a slug
run.
"""

from collections import namedtuple
import numpy as np
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")

def write_cluster_cloudyspec(data, model_name, fmt):
    """
    Write out data computed by cloudy on a slug spectrum

    Parameters
       data : namedtuple
          Cloudy spectral data for clusters to be written; a namedtuple
          containing the fields id, time, cloudy_wl, cloudy_inc, cloudy_trans,
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
        fp = open(model_name+'_cluster_cloudyspec.txt', 'w')

        # Write header lines
        fp.write(("{:<14s}"*7).
                 format('UniqueID', 'Time', 'Wavelength', 'Incident',
                        'Transmitted', 'Emitted', 'Trans+Emit') + "\n")
        fp.write(("{:<14s}"*7).
                 format('', '(yr)', '(Angstrom)', '(erg/s/A)',
                        '(erg/s/A)', '(erg/s/A)', '(erg/s/A)') + "\n")
        fp.write(("{:<14s}"*7).
                 format('-----------', '-----------', '-----------',
                        '-----------', '-----------', '-----------',
                        '-----------')
                 + "\n")

        # Write data
        for i in range(data.cloudy_inc.shape[0]):
            # If this is a new trial, write a separator
            if i != 0:
                if data.trial[i] != data.trial[i-1]:
                    fp.write("-"*(7*14-3)+"\n")
            for j in range(data.cloudy_wl.shape[0]):
                fp.write(("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   " +
                          "{:11.5e}   {:11.5e}   {:11.5e}\n")
                         .format(data.id[i], data.time[i], 
                                 data.cloudy_wl[j],
                                 data.cloudy_inc[i,j],
                                 data.cloudy_trans[i,j],
                                 data.cloudy_emit[i,j],
                                 data.cloudy_trans_emit[i,j]))

        # Close
        fp.close()

    elif fmt == 'bin' or fmt == 'binary':

        # Binary mode
        fp = open(model_name+'_cluster_cloudyspec.bin', 'wb')

        # Write out wavelength data
        fp.write(np.int64(len(data.cloudy_wl)))
        fp.write(data.cloudy_wl)

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
                fp.write(data.cloudy_inc[k,:])
                fp.write(data.cloudy_trans[k,:])
                fp.write(data.cloudy_emit[k,:])
                fp.write(data.cloudy_trans_emit[k,:])

            # Move pointer
            ptr = block_end

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

        # Convert spectra to FITS columns, and make an HDU from
        # them
        speccols = []
        speccols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=data.trial))
        speccols.append(fits.Column(name="UniqueID", format="1K",
                                    unit="", array=data.id))
        speccols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=data.time))
        speccols.append(fits.Column(name="Incident_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=data.cloudy_inc))
        speccols.append(fits.Column(name="Transmitted_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=data.cloudy_trans))
        speccols.append(fits.Column(name="Emitted_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=data.cloudy_emit))
        speccols.append(fits.Column(name="Transmitted_plus_emitted_spectrum",
                                    format=fmtstring,
                                    unit="erg/s/A",
                                    array=data.cloudy_trans_emit))
        specfits = fits.ColDefs(speccols)
        spechdu = fits.BinTableHDU.from_columns(specfits)

        # Create dummy primary HDU
        prihdu = fits.PrimaryHDU()

        # Create HDU list and write to file
        hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
        hdulist.writeto(model_name+'_cluster_cloudyspec.fits', 
                        overwrite=True)
