"""
Function to write a set of output cluster files in SLUG2 format,
starting from a cluster data set as returned by read_cluster. This can
be used to translate data from one format to another (e.g., bin to
fits), or to consolidate multiple runs into a single output file.
"""

import numpy as np
from scipy.interpolate import interp1d
from cloudy import write_cluster_cloudyphot
from cloudy import write_cluster_cloudylines
from cloudy import write_cluster_cloudyspec
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")


def write_cluster(data, model_name, fmt):
    """
    Function to write a set of output cluster files in SLUG2 format,
    starting from a cluster data set as returned by read_cluster.

    Parameters
       data : namedtuple
          Cluster data to be written, in the namedtuple format returned
          by read_cluster
       model_name : string
          Base file name to give the model to be written. Can include a
          directory specification if desired.
       fmt : 'txt' | 'ascii' | 'bin' | 'binary' | 'fits' | 'fits2'
          Format for the output file; 'txt' and 'ascii' produce ASCII
          text files, 'bin' or 'binary' produce binary files, and
          'fits' or 'fits2' product FITS files; 'fits2' uses an
          ordering that allows for more efficient querying of outputs
          too large to fit in memory

    Returns
       Nothing
    """

    # Make sure fmt is valid
    if fmt != 'ascii' and fmt != 'txt' and fmt != 'bin' and \
       fmt != 'binary' and fmt != 'fits' and fmt != 'fits2':
        raise ValueError("unrecognized format {}".format(fmt))

    # Make sure we're not trying to do fits if we don't have astropy
    if (fmt == 'fits' or fmt == 'fits2') and fits is None:
        raise ValueError("Couldn't import astropy, so fits format "+
                         "is unavailable.")

    ################################################################
    # Write the properties file if we have the data for it
    ################################################################
    if 'form_time' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_prop.txt', 'w')

            # Write header lines
            if 'A_V' in data._fields:
                fp.write(("{:<14s}"*11).
                         format('UniqueID', 'Time', 'FormTime',
                                'Lifetime', 'TargetMass',
                                'BirthMass', 'LiveMass', 'StellarMass',
                                'NumStar', 'MaxStarMass', 'A_V') + "\n")
                fp.write(("{:<14s}"*11).
                         format('', '(yr)', '(yr)',
                                '(yr)', '(Msun)', '(Msun)'
                                '(Msun)', '(Msun)',
                                '', '(Msun)', '(mag)') + "\n")
                fp.write(("{:<14s}"*11).
                         format('-----------', '-----------', '-----------',
                                '-----------', '-----------', '-----------',
                                '-----------', '-----------',
                                '-----------', '-----------',
                                '-----------') + "\n")
            else:
                fp.write(("{:<14s}"*10).
                         format('UniqueID', 'Time', 'FormTime',
                                'Lifetime', 'TargetMass',
                                'BirthMass', 'LiveMass', 'StellarMass',
                                'NumStar', 'MaxStarMass') + "\n")
                fp.write(("{:<14s}"*10).
                         format('', '(yr)', '(yr)',
                                '(yr)', '(Msun)', '(Msun)',
                                '(Msun)', '(Msun)',
                                '', '(Msun)') + "\n")
                fp.write(("{:<14s}"*10).
                         format('-----------', '-----------', '-----------',
                                '-----------', '-----------', '-----------',
                                '-----------', '-----------',
                                '-----------', '-----------') + "\n")

            # Write data
            for i in range(len(data.id)):
                if 'A_V' in data._fields:
                    # If this is a new trial, write a separator
                    if i != 0:
                        if data.trial[i] != data.trial[i-1]:
                            fp.write("-"*(11*14-3)+"\n")
                    fp.write("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   "
                             "{:11.5e}   {:11.5e}   {:11.5e}   {:11.5e}   "
                             "{:11d}   {:11.5e}   {:11.5e}\n".
                             format(data.id[i], data.time[i], 
                                    data.form_time[i], data.lifetime[i],
                                    data.target_mass[i],
                                    data.actual_mass[i],
                                    data.live_mass[i],
                                    data.stellar_mass[i],
                                    data.num_star[i],
                                    data.max_star_mass[i],
                                    data.A_V[i]))
                else:
                    if i != 0:
                        if data.trial[i] != data.trial[i-1]:
                            fp.write("-"*(10*14-3)+"\n")
                    fp.write("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   "
                             "{:11.5e}   {:11.5e}   {:11.5e}   {:11.5e}   "
                             "{:11d}   {:11.5e}\n".
                             format(data.id[i], data.time[i], 
                                    data.form_time[i], data.lifetime[i],
                                    data.target_mass[i],
                                    data.actual_mass[i],
                                    data.live_mass[i],
                                    data.stellar_mass[i],
                                    data.num_star[i],
                                    data.max_star_mass[i]))

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_prop.bin', 'wb')

            # Write out a byte indicating extinction or no extinction
            if 'A_V' in data._fields:
                fp.write(str(bytearray([1])))
            else:
                fp.write(str(bytearray([0])))

            # Construct lists of times and trials
            trials = np.unique(data.trial)
            times = np.unique(data.time)

            # Break data into blocks of clusters with the same time
            # and trial number
            ptr = 0
            for i in range(len(trials)):
                for j in range(len(times)):

                    # Find block of clusters with this time and trial
                    block_match \
                        = np.logical_and(data.trial[ptr:] == trials[i],
                                         data.time[ptr:] == times[j])
                    block_end = ptr + np.argmin(block_match)

                    # Special case: if we found no clusters in this
                    # block, make sure that's not because all the
                    # remaining clusters belong in it, and the search
                    # ran off the end of the data. If that is the
                    # case, adjust block_end appropriately.
                    if block_end == ptr:
                        if np.sum(block_match) > 0:
                            block_end = len(data.trial)

                    # Special case: if block_end is the last entry in
                    # the data, check if the last entry is the same as
                    # the previous one. If so, move block_end one
                    # space, to off the edge of the data.
                    if block_end == len(data.trial)-1 and \
                       data.trial[-1] == trials[i] and \
                       data.time[-1] == times[j]:
                        block_end = block_end+1

                    # Write out time and number of clusters
                    ncluster = block_end - ptr
                    fp.write(np.uint(i))
                    fp.write(times[j])
                    fp.write(ncluster)

                    # Loop over clusters and write them
                    for k in range(ptr, block_end):
                        fp.write(data.id[k])
                        fp.write(data.form_time[k])
                        fp.write(data.lifetime[k])
                        fp.write(data.target_mass[k])
                        fp.write(data.actual_mass[k])
                        fp.write(data.live_mass[k])
                        fp.write(data.stellar_mass[k])
                        fp.write(data.num_star[k])
                        fp.write(data.max_star_mass[k])
                        if 'A_V' in data._fields:
                            fp.write(data.A_V[k])

                    # Move pointer
                    ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits' or fmt == 'fits2':

            ########################################################
            # FITS mode
            ########################################################

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=data.trial))
            cols.append(fits.Column(name="UniqueID", format="1K",
                                    unit="", array=data.id))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=data.time))
            cols.append(fits.Column(name="FormTime", format="1D",
                                    unit="yr", array=data.form_time))
            cols.append(fits.Column(name="Lifetime", format="1D",
                                    unit="yr", array=data.lifetime))
            cols.append(fits.Column(name="TargetMass", format="1D",
                                    unit="Msun", array=data.target_mass))
            cols.append(fits.Column(name="BirthMass", format="1D",
                                    unit="Msun", array=data.actual_mass))
            cols.append(fits.Column(name="LiveMass", format="1D",
                                    unit="Msun", array=data.live_mass))
            cols.append(fits.Column(name="StellarMass", format="1D",
                                    unit="Msun", array=data.stellar_mass))
            cols.append(fits.Column(name="NumStar", format="1K",
                                    unit="", array=data.num_star))
            cols.append(fits.Column(name="MaxStarMass", format="1D",
                                    unit="Msun", array=data.max_star_mass))
            if 'A_V' in data._fields:
                cols.append(fits.Column(name="A_V", format="1D",
                                        unit="mag", array=data.A_V))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_cluster_prop.fits',
                            clobber=True)


    ################################################################
    # Write spectra file if we have the data for it
    ################################################################
    if 'spec' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_spec.txt', 'w')

            # If we have nebular data, and this data wasn't originally
            # read in ASCII mode, then the stellar and nebular spectra
            # won't be on the same grids. Since in ASCII mode they'll
            # be written out in the same grid, we need to interpolate
            # the stellar spectra onto the nebular grid
            if ('wl_neb' in data._fields) and \
               (len(data.wl_neb) > len(data.wl)):
                wl = data.wl_neb
                # Suppress the numpy warnings we're going to generate
                # if any of the entries in spec are 0
                save_err = np.seterr(divide='ignore', invalid='ignore')
                ifunc = interp1d(np.log(data.wl), np.log(data.spec))
                spec = np.exp(ifunc(np.log(data.wl_neb)))
                # Restore original error settings
                np.seterr(divide=save_err['divide'], 
                          invalid=save_err['invalid'])
                # Fix NaN's
                spec[np.isnan(spec)] = 0.0
                if 'wl_neb_ex' in data._fields:
                    # Same for extincted nebular data
                    wl_ex = data.wl_neb_ex
                    save_err = np.seterr(divide='ignore', invalid='ignore')
                    ifunc = interp1d(np.log(data.wl_ex),
                                     np.log(data.spec_ex))
                    spec_ex = np.exp(ifunc(np.log(data.wl_neb_ex)))
                    np.seterr(divide=save_err['divide'],
                              invalid=save_err['invalid'])
                    spec_ex[np.isnan(spec_ex)] = 0.0
            else:
                # If no nebular data, just replicate the original
                # stellar grid
                wl = data.wl
                spec = data.spec
                if 'wl_ex' in data._fields:
                    wl_ex = data.wl_ex
                    spec_ex = data.spec_ex

            # Construct header lines
            line1 = ("{:<14s}"*4).format('UniqueID', 'Time', 
                                         'Wavelength', 'L_lambda')
            line2 = ("{:<14s}"*4).format('', '(yr)', '(Angstrom)', 
                                         '(erg/s/A)')
            line3 = ("{:<14s}"*4).format('-----------', '-----------', 
                                         '-----------', '-----------')
            sep_length = 4*14-3
            out_line = "{:11d}   {:11.5e}   {:11.5e}   {:11.5e}"
            if 'spec_neb' in data._fields:
                line1 = line1 + "{:<14s}".format("L_l_neb")
                line2 = line2 + "{:<14s}".format("(erg/s/A)")
                line3 = line3 + "{:<14s}".format("-----------")
                sep_length = sep_length + 14
                out_line = out_line + "   {:11.5e}"
            if 'spec_ex' in data._fields:
                line1 = line1 + "{:<14s}".format("L_lambda_ex")
                line2 = line2 + "{:<14s}".format("(erg/s/A)")
                line3 = line3 + "{:<14s}".format("-----------")
                sep_length = sep_length + 14
                out_line1 = out_line + "   {:11.5e}"
            if 'spec_neb_ex' in data._fields:
                line1 = line1 + "{:<14s}".format("L_l_neb_ex")
                line2 = line2 + "{:<14s}".format("(erg/s/A)")
                line3 = line3 + "{:<14s}".format("-----------")
                sep_length = sep_length + 14
                out_line1 = out_line1 + "   {:11.5e}"

            # Write header lines
            fp.write(line1+"\n")
            fp.write(line2+"\n")
            fp.write(line3+"\n")

            # Write data
            nl = len(data.wl)
            if 'spec_ex' in data._fields:
                offset = np.where(wl_ex[0] == wl)[0][0]
                nl_ex = len(wl_ex)
            else:
                offset = 0
                nl_ex = 0
            for i in range(len(data.id)):

                # If this is a new trial, write a separator
                if i != 0:
                    if data.trial[i] != data.trial[i-1]:
                        fp.write("-"*sep_length+"\n")

                # Write data for this trial
                for j in range(nl):
                    out_data = [data.id[i], data.time[i],
                                wl[j], spec[i,j]]
                    if 'spec_neb' in data._fields:
                        out_data = out_data + [data.spec_neb[i,j]]
                    if j >= offset and j < offset + nl_ex:
                        out_data = out_data + [spec_ex[i,j-offset]]
                        if 'spec_neb_ex' in data._fields:
                            out_data = out_data + \
                                       [data.spec_neb_ex[i,j-offset]]
                        out_fmt = out_line1
                    else:
                        out_fmt = out_line
                    fp.write(out_fmt.format(*out_data)+"\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_spec.bin', 'wb')

            # Write out bytes indicating nebular or no nebular, and
            # extinction or no extinction
            if 'spec_neb' in data._fields:
                fp.write(str(bytearray([1])))
            else:
                fp.write(str(bytearray([0])))
            if 'spec_ex' in data._fields:
                fp.write(str(bytearray([1])))
            else:
                fp.write(str(bytearray([0])))

            # Write out wavelength data
            fp.write(np.int64(len(data.wl)))
            fp.write(data.wl)
            if 'spec_neb' in data._fields:
                fp.write(np.int64(len(data.wl_neb)))
                fp.write(data.wl_neb)
            if 'spec_ex' in data._fields:
                fp.write(np.int64(len(data.wl_ex)))
                fp.write(data.wl_ex)
            if 'spec_neb_ex' in data._fields:
                fp.write(np.int64(len(data.wl_neb_ex)))
                fp.write(data.wl_neb_ex)

            # Construct lists of times and trials
            trials = np.unique(data.trial)
            times = np.unique(data.time)

            # Break data into blocks of clusters with the same time
            # and trial number
            ptr = 0
            for i in range(len(trials)):
                for j in range(len(times)):

                    # Find block of clusters with this time and trial
                    block_match \
                        = np.logical_and(data.trial[ptr:] == trials[i],
                                         data.time[ptr:] == times[j])
                    block_end = ptr + np.argmin(block_match)

                    # Special case: if we found no clusters in this
                    # block, make sure that's not because all the
                    # remaining clusters belong in it, and the search
                    # ran off the end of the data. If that is the
                    # case, adjust block_end appropriately.
                    if block_end == ptr:
                        if np.sum(block_match) > 0:
                            block_end = len(data.trial)

                    # Special case: if block_end is the last entry in
                    # the data, check if the last entry is the same as
                    # the previous one. If so, move block_end one
                    # space, to off the edge of the data.
                    if block_end == len(data.trial)-1 and \
                       data.trial[-1] == trials[i] and \
                       data.time[-1] == times[j]:
                        block_end = block_end+1

                    # Write out time and number of clusters
                    ncluster = block_end - ptr
                    fp.write(np.uint(i))
                    fp.write(times[j])
                    fp.write(ncluster)

                    # Loop over clusters and write them
                    for k in range(ptr, block_end):
                        fp.write(data.id[k])
                        fp.write(data.spec[k,:])
                        if 'spec_neb' in data._fields:
                            fp.write(data.spec_neb[k,:])
                        if 'spec_ex' in data._fields:
                            fp.write(data.spec_ex[k,:])
                        if 'spec_neb_ex' in data._fields:
                            fp.write(data.spec_neb_ex[k,:])

                    # Move pointer
                    ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits' or fmt == 'fits2':

            ########################################################
            # FITS mode
            ########################################################

            # Convert wavelength data to FITS columns and make an HDU
            # from it; complication: astropy expects the dimensions of
            # the array to be (n_entries, n_wavelengths)
            nl = data.wl.shape[0]
            fmtstring = str(nl)+"D"
            wlcols = [fits.Column(name="Wavelength",
                                  format=fmtstring,
                                  unit="Angstrom", 
                                  array=data.wl.reshape(1,nl))]
            if 'spec_neb' in data._fields:
                nl_neb = data.wl_neb.shape[0]
                fmtstring_neb = str(nl_neb)+"D"
                wlcols.append(
                    fits.Column(name="Wavelength_neb",
                                format=fmtstring_neb,
                                unit="Angstrom", 
                                array=data.wl_neb.reshape(1,nl_neb)))
            if 'spec_ex' in data._fields:
                nl_ex = data.wl_ex.shape[0]
                fmtstring_ex = str(nl_ex)+"D"
                wlcols.append(
                    fits.Column(name="Wavelength_ex",
                                format=fmtstring_ex,
                                unit="Angstrom", 
                                array=data.wl_ex.reshape(1,nl_ex)))
            if 'spec_neb_ex' in data._fields:
                nl_neb_ex = data.wl_neb_ex.shape[0]
                fmtstring_neb_ex = str(nl_neb_ex)+"D"
                wlcols.append(
                    fits.Column(name="Wavelength_neb_ex",
                                format=fmtstring_neb_ex,
                                unit="Angstrom", 
                                array=data.wl_neb_ex.reshape(1,nl_neb_ex)))
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
            speccols.append(fits.Column(name="L_lambda",
                                        format=fmtstring,
                                        unit="erg/s/A",
                                        array=data.spec))
            if 'spec_neb' in data._fields:
                speccols.append(fits.Column(name="L_lambda_neb",
                                            format=fmtstring_neb,
                                            unit="erg/s/A",
                                            array=data.spec_neb))
            if 'spec_ex' in data._fields:
                speccols.append(fits.Column(name="L_lambda_ex",
                                            format=fmtstring_ex,
                                            unit="erg/s/A",
                                            array=data.spec_ex))
            if 'spec_neb_ex' in data._fields:
                speccols.append(fits.Column(name="L_lambda_neb_ex",
                                            format=fmtstring_neb_ex,
                                            unit="erg/s/A",
                                            array=data.spec_neb_ex))
            specfits = fits.ColDefs(speccols)
            spechdu = fits.BinTableHDU.from_columns(specfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
            hdulist.writeto(model_name+'_cluster_spec.fits', 
                            clobber=True)

    ################################################################
    # Write photometry file if we have the data for it
    ################################################################
    if 'phot' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_phot.txt', 'w')

            # Write header lines
            fp.write(("{:<21s}"*2).format('UniqueID', 'Time'))
            fac = 1
            for f in data.filter_names:
                fp.write("{:<21s}".format(f))
            if 'phot_neb' in data._fields:
                fac = fac + 1
                for f in data.filter_names:
                    fp.write("{:<21s}".format(f+'_n'))
            if 'phot_ex' in data._fields:
                fac = fac + 1
                for f in data.filter_names:
                    fp.write("{:<21s}".format(f+'_ex'))
            if 'phot_neb_ex' in data._fields:
                fac = fac + 1
                for f in data.filter_names:
                    fp.write("{:<21s}".format(f+'_nex'))
            fp.write("\n")
            fp.write(("{:<21s}"*2).format('', '(yr)'))
            for i in range(fac):
                for f in data.filter_units:
                    fp.write("({:s}".format(f)+")"+" "*(19-len(f)))
            fp.write("\n")
            nf = len(data.filter_names)
            fp.write(("{:<21s}"*2).
                     format('------------------', '------------------'))
            for j in range(fac):
                for i in range(nf):
                    fp.write("{:<21s}".format('------------------'))
            fp.write("\n")

            # Write data
            for i in range(len(data.id)):
                # If this is a new trial, write a separator
                if i != 0:
                    if data.trial[i] != data.trial[i-1]:
                        fp.write("-"*((2+fac*nf)*21-3)+"\n")
                fp.write("       {:11d}          {:11.5e}"
                         .format(data.id[i], data.time[i]))
                for j in range(nf):
                    fp.write("          {:11.5e}".format(data.phot[i,j]))
                if 'phot_neb' in data._fields:
                    for j in range(nf):
                        fp.write("          {:11.5e}".
                                 format(data.phot_neb[i,j]))
                if 'phot_ex' in data._fields:
                    for j in range(nf):
                        if np.isnan(data.phot_ex[i,j]):
                            fp.write("          {:11s}".
                                     format(""))
                        else:
                            fp.write("          {:11.5e}".
                                     format(data.phot_ex[i,j]))
                if 'phot_neb_ex' in data._fields:
                    for j in range(nf):
                        if np.isnan(data.phot_ex[i,j]):
                            fp.write("          {:11s}".
                                     format(""))
                        else:
                            fp.write("          {:11.5e}".
                                     format(data.phot_neb_ex[i,j]))
                fp.write("\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_phot.bin', 'wb')

            # Write number of filters and filter names as ASCII
            nf = len(data.filter_names)
            fp.write(str(nf)+"\n")
            for i in range(nf):
                fp.write(data.filter_names[i] + " " + 
                         data.filter_units[i] + "\n")

            # Write out bytes indicating nebular or no nebular, and
            # extinction or no extinction
            if 'spec_neb' in data._fields:
                fp.write(str(bytearray([1])))
            else:
                fp.write(str(bytearray([0])))
            if 'spec_ex' in data._fields:
                fp.write(str(bytearray([1])))
            else:
                fp.write(str(bytearray([0])))

            # Construct lists of times and trials
            trials = np.unique(data.trial)
            times = np.unique(data.time)

            # Break data into blocks of clusters with the same time
            # and trial number
            ptr = 0
            for i in range(len(trials)):
                for j in range(len(times)):

                    # Find block of clusters with this time and trial
                    block_match \
                        = np.logical_and(data.trial[ptr:] == trials[i],
                                         data.time[ptr:] == times[j])
                    block_end = ptr + np.argmin(block_match)

                    # Special case: if we found no clusters in this
                    # block, make sure that's not because all the
                    # remaining clusters belong in it, and the search
                    # ran off the end of the data. If that is the
                    # case, adjust block_end appropriately.
                    if block_end == ptr:
                        if np.sum(block_match) > 0:
                            block_end = len(data.trial)

                    # Special case: if block_end is the last entry in
                    # the data, check if the last entry is the same as
                    # the previous one. If so, move block_end one
                    # space, to off the edge of the data.
                    if block_end == len(data.trial)-1 and \
                       data.trial[-1] == trials[i] and \
                       data.time[-1] == times[j]:
                        block_end = block_end+1

                    # Write out time and number of clusters
                    ncluster = block_end - ptr
                    fp.write(np.uint(i))
                    fp.write(times[j])
                    fp.write(ncluster)

                    # Loop over clusters and write them
                    for k in range(ptr, block_end):
                        fp.write(data.id[k])
                        fp.write(data.phot[k,:])
                        if 'phot_neb' in data._fields:
                            fp.write(data.phot_neb[k,:])
                        if 'phot_ex' in data._fields:
                            fp.write(data.phot_ex[k,:])
                        if 'phot_neb_ex' in data._fields:
                            fp.write(data.phot_neb_ex[k,:])

                    # Move pointer
                    ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=data.trial))
            cols.append(fits.Column(name="UniqueID", format="1K",
                                    unit="", array=data.id))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=data.time))
            for i in range(len(data.filter_names)):
                cols.append(fits.Column(name=data.filter_names[i],
                                        unit=data.filter_units[i],
                                        format="1D",
                                        array=data.phot[:,i]))
            if 'phot_neb' in data._fields:
                for i in range(len(data.filter_names)):
                    cols.append(
                        fits.Column(name=data.filter_names[i]+'_neb',
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=data.phot_neb[:,i]))
            if 'phot_ex' in data._fields:
                for i in range(len(data.filter_names)):
                    cols.append(
                        fits.Column(name=data.filter_names[i]+'_ex',
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=data.phot_ex[:,i]))
            if 'phot_neb_ex' in data._fields:
                for i in range(len(data.filter_names)):
                    cols.append(
                        fits.Column(name=data.filter_names[i]+'_neb_ex',
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=data.phot_neb_ex[:,i]))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_cluster_phot.fits',
                            clobber=True)

        elif fmt == 'fits2':

            ########################################################
            # FITS2 mode
            ########################################################

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Convert trial, time, unique ID to FITS columns
            hdu1cols = []
            hdu1cols.append(fits.Column(name="Trial", format="1K",
                                        unit="", array=data.trial))
            hdu1cols.append(fits.Column(name="UniqueID", format="1K",
                                        unit="", array=data.id))
            hdu1cols.append(fits.Column(name="Time", format="1D",
                                        unit="yr", array=data.time))

            # Make HDU containing trial, time, unique ID
            hdu1fitscol = fits.ColDefs(hdu1cols)
            hdu1 = fits.BinTableHDU.from_columns(hdu1fitscol)

            # Create master HDU list
            hdulist = [prihdu, hdu1]

            # Now loop over filters; each will be stored as a single
            # column in its own HDU
            for i in range(len(data.filter_names)):

                # Create column
                col = [fits.Column(name=data.filter_names[i],
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=data.phot[:,i])]

                # Create HDU from column
                hdu = fits.BinTableHDU.from_columns(
                    fits.ColDefs(col))

                # Add to master HDU list
                hdulist.append(hdu)

            # Repeat loop over filters for nebular, extincted, and
            # extincted+nebular photometry if we have those
            if 'phot_neb' in data._fields:
                for i in range(len(data.filter_names)):
                    col = [fits.Column(name=data.filter_names[i]+'_neb',
                                       unit=data.filter_units[i],
                                       format="1D",
                                       array=data.phot_neb[:,i])]
                    hdu = fits.BinTableHDU.from_columns(
                        fits.ColDefs(col))
                    hdulist.append(hdu)
            if 'phot_ex' in data._fields:
                for i in range(len(data.filter_names)):
                    col = [fits.Column(name=data.filter_names[i]+'_ex',
                                       unit=data.filter_units[i],
                                       format="1D",
                                       array=data.phot_ex[:,i])]
                    hdu = fits.BinTableHDU.from_columns(
                        fits.ColDefs(col))
                    hdulist.append(hdu)
            if 'phot_neb_ex' in data._fields:
                for i in range(len(data.filter_names)):
                    col = [fits.Column(name=data.filter_names[i]+'_neb_ex',
                                       unit=data.filter_units[i],
                                       format="1D",
                                       array=data.phot_neb_ex[:,i])]
                    hdu = fits.BinTableHDU.from_columns(
                        fits.ColDefs(col))
                    hdulist.append(hdu)

            # Create final HDU list and write to file
            hdulist = fits.HDUList(hdulist)
            hdulist.writeto(model_name+'_cluster_phot.fits',
                            clobber=True)

    ################################################################
    # Write cloudy files if we have the data for them
    ################################################################
    if 'cloudy_inc' in data._fields:
        write_cluster_cloudyspec(data, model_name, fmt=fmt)
    if 'cloudy_linelum' in data._fields:
        write_cluster_cloudylines(data, model_name, fmt=fmt)
    if 'cloudy_phot_trans' in data._fields:
        write_cluster_cloudyphot(data, model_name, fmt=fmt)
