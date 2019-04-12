"""
Function to write a set of output cluster files in SLUG2 format,
starting from a cluster data set as returned by read_cluster. This can
be used to translate data from one format to another (e.g., bin to
fits), or to consolidate multiple runs into a single output file.
"""

import numpy as np
import struct
import os
import sys
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if not on_rtd:
    from scipy.interpolate import interp1d
else:
    # Dummy interp1d function for RTD
    def interp1d(dummy1, dummy2, axis=None):
        pass
from .cloudy import write_cluster_cloudyparams
from .cloudy import write_cluster_cloudyphot
from .cloudy import write_cluster_cloudylines
from .cloudy import write_cluster_cloudyspec
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

            # Figure out which fields we have
            fields = ['UniqueID', 'Time', 'FormTime',
                      'Lifetime', 'TargetMass',
                      'BirthMass', 'LiveMass', 'StellarMass',
                      'NumStar', 'MaxStarMass']
            units = ['', '(yr)', '(yr)',
                     '(yr)', '(Msun)', '(Msun)',
                     '(Msun)', '(Msun)',
                     '', '(Msun)']
            if 'A_V' in data._fields:
                fields.append('A_V')
                units.append('(mag)')
            if 'A_Vneb' in data._fields:
                fields.append('A_Vneb')
                units.append('(mag)')
            vp_i = 0
            while 'VP'+repr(vp_i) in data._fields:
                fields.append('VP'+repr(vp_i))
                units.append('')
                vp_i += 1

            # Write headers
            fp.write(("{:<14s}"*len(fields)).format(*fields))
            fp.write("\n")
            fp.write(("{:<14s}"*len(units)).format(*units))
            fp.write("\n")
            fp.write(("{:<14s}"*len(fields)).
                     format(*(['-----------']*len(fields))))
            fp.write("\n")

            # Write data
            for i in range(len(data.id)):
                # If this is a new trial, write a separator
                if i != 0:
                    if data.trial[i] != data.trial[i-1]:
                        fp.write("-"*(len(fields)*14-3)+"\n")
                # Write fields that are always present
                fp.write("{:11d}   {:11.5e}   {:11.5e}   {:11.5e}   "
                         "{:11.5e}   {:11.5e}   {:11.5e}   {:11.5e}   "
                         "{:11d}   {:11.5e}".format(
                             data.id[i], data.time[i], 
                             data.form_time[i], data.lifetime[i],
                             data.target_mass[i],
                             data.actual_mass[i],
                             data.live_mass[i],
                             data.stellar_mass[i],
                             data.num_star[i],
                             data.max_star_mass[i]))
                # Write optional fields
                if 'A_V' in data._fields:
                    fp.write("   {:11.5e}".format(data.A_V[i]))
                if 'A_Vneb' in data._fields:
                    fp.write("   {:11.5e}".format(data.A_Vneb[i]))
                vp_i = 0
                while 'VP'+repr(vp_i) in data._fields:
                    current_vp=getattr(data, "VP"+repr(vp_i))
                    fp.write("   {:11.5e}".format(current_vp[i]))
                    vp_i+=1

                # End line
                fp.write("\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_prop.bin', 'wb')

            # Write out a bytes indicating extinction or no
            # extinction, and nebular extinction or no nebular
            # extinction
            if 'A_V' in data._fields:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([1])))
                else:
                    fp.write(b'\x01')
            else:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([0])))
                else:
                    fp.write(b'\x00')
            if 'A_Vneb' in data._fields:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([1])))
                else:
                    fp.write(b'\x01')
            else:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([0])))
                else:
                    fp.write(b'\x00')

            # Write out number of variable parameters
            nvp = 0
            for field in data._fields:
                if field.startswith("VP"):
                    nvp+=1
            fp.write( struct.pack('i', nvp) )   
           
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
                    if 'A_Vneb' in data._fields:
                        fp.write(data.A_Vneb[k])
                    for field in sorted(data._fields):
                        if field.startswith("VP"):
                            fp.write(getattr(data, field)[k])

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
            if 'A_Vneb' in data._fields:
                cols.append(fits.Column(name="A_Vneb", format="1D",
                                        unit="mag", array=data.A_Vneb))
            vp_i = 0
            while 'VP'+repr(vp_i) in data._fields:
                cols.append(fits.Column(name="VP"+repr(vp_i), format="1D",
                                        unit="",
                                        array=getattr(data,
                                                      "VP"+repr(vp_i))))
                vp_i += 1
            fitscols = fits.ColDefs(cols)
            
            
            
            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_cluster_prop.fits',
                            overwrite=True)


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
                if sys.version_info < (3,):
                    fp.write(str(bytearray([1])))
                else:
                    fp.write(b'\x01')
            else:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([0])))
                else:
                    fp.write(b'\x00')
            if 'spec_ex' in data._fields:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([1])))
                else:
                    fp.write(b'\x01')
            else:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([0])))
                else:
                    fp.write(b'\x00')

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

            if 'spec_rec' in data._fields:
                nl_rec = data.wl_r.shape[0]
                fmtstring_rec = str(nl_rec)+"D"
                wlcols.append(
                    fits.Column(name="Wavelength_Rectified",
                                format=fmtstring_rec,
                                unit="Angstrom",
                                array=data.wl_r.reshape(1,nl_rec)))
                                
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
                                            
            if 'spec_rec' in data._fields:
                speccols.append(fits.Column(name="Rectified_Spec",
                                            format=fmtstring_rec,
                                            unit="",
                                            array=data.spec_rec))
            specfits = fits.ColDefs(speccols)
            spechdu = fits.BinTableHDU.from_columns(specfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
            hdulist.writeto(model_name+'_cluster_spec.fits', 
                            overwrite=True)

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
            if sys.version_info < (3,):
                fp.write(str(nf)+"\n")
                for i in range(nf):
                    fp.write(data.filter_names[i] + " " + 
                             data.filter_units[i] + "\n")
            else:
                fp.write(bytes(str(nf)+"\n", "ascii"))
                for i in range(nf):
                    fp.write(bytes(
                        data.filter_names[i] + " " + 
                        data.filter_units[i] + "\n", "ascii"))

            # Write out bytes indicating nebular or no nebular, and
            # extinction or no extinction
            if 'phot_neb' in data._fields:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([1])))
                else:
                    fp.write(b'\x01')
            else:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([0])))
                else:
                    fp.write(b'\x00')
            if 'phot_ex' in data._fields:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([1])))
                else:
                    fp.write(b'\x01')
            else:
                if sys.version_info < (3,):
                    fp.write(str(bytearray([0])))
                else:
                    fp.write(b'\x00')

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
                            overwrite=True)

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
                            overwrite=True)


    ################################################################
    # Write equivalent width file if we have the data for it
    ################################################################
    if 'ew' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################
            
            raise NotImplementedError(
                "ERROR: EW available for fits format output only")


        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            raise NotImplementedError(
                "ERROR: EW available for fits format output only")
            
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
            for i in range(len(data.line_names)):
                #print data.line_names
                cols.append(fits.Column(name=data.line_names[i],
                                        unit=data.line_units[i],
                                        format="1D",
                                        array=data.ew[:,i]))
            
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_cluster_ew.fits',
                            overwrite=True)

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

            # Now loop over lines; each will be stored as a single
            # column in its own HDU
            for i in range(len(data.line_names)):

                # Create column
                col = [fits.Column(name=data.line_names[i],
                                    unit=data.line_units[i],
                                    format="1D",
                                    array=data.ew[:,i])]
                # Create HDU from column
                hdu = fits.BinTableHDU.from_columns(
                    fits.ColDefs(col))

                # Add to master HDU list
                hdulist.append(hdu)

            # Create final HDU list and write to file
            hdulist = fits.HDUList(hdulist)
            hdulist.writeto(model_name+'_cluster_ew.fits',
                            overwrite=True)


    ################################################################
    # Write the supernova file if we have the data for it
    ################################################################
    if 'tot_sn' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_sn.txt', 'w')
            
            # Write header
            fp.write(("{:<14s}"*4+"\n").\
                     format('UniqueID', 'Time', 'TotSN',
                            'StochSN'))
            fp.write(("{:<14s}"*4+"\n").\
                     format('', '(yr)', '', '', '', '(Msun)'))
            fp.write(("{:<14s}"*4+"\n").\
                     format('-----------', '-----------', 
                            '-----------', '-----------'))
            sep_length = 4*14-3
            out_line = "{:>11d}   {:11.5e}   {:>11.5e}   " + \
                       "{:>11d}\n"

            # Write data
            for i in range(len(data.id)):
                fp.write(out_line.
                         format(data.id[i], data.time[i],
                                data.tot_sn[i], data.stoch_sn[i]))

            # Close file
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':
            
            ########################################################
            # Binary mode
            ########################################################
            
            fp = open(model_name+'_cluster_sn.bin', 'wb')

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
                    fp.write(data.tot_sn[k])
                    fp.write(data.stoch_sn[k])
                   
                # Move pointer
                ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################
            
            # Create a HDU containing SN information
            sncols = []
            sncols.append(fits.Column(name="Trial", format="1K",
                                      unit="", array=data.trial))
            sncols.append(fits.Column(name="UniqueID", format="1K",
                                      unit="", array=data.id))
            sncols.append(fits.Column(name="Time", format="1D",
                                      unit="yr", array=data.time))
            sncols.append(fits.Column(name="TotSN", format="1D",
                                      unit="", array=data.tot_sn))
            sncols.append(fits.Column(name="StochSN", format="1K",
                                      unit="", array=data.stoch_sn))
            snfits = fits.ColDefs(sncols)
            snhdu = fits.BinTableHDU.from_columns(snfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, snhdu])
            hdulist.writeto(model_name+'_cluster_sn.fits', 
                            overwrite=True)
            
            
    ################################################################
    # Write the yield file if we have the data for it
    ################################################################
    if 'yld' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_cluster_yield.txt', 'w')

            # Write header
            fp.write(("{:<14s}"*6+"\n").\
                     format('UniqueID', 'Time', 'Symbol',
                            'Z', 'A', 'Yield'))
            fp.write(("{:<14s}"*6+"\n").\
                     format('', '(yr)', '', '', '', '(Msun)'))
            fp.write(("{:<14s}"*6+"\n").\
                     format('-----------', '-----------', 
                            '-----------', '-----------', 
                            '-----------', '-----------'))
            sep_length = 6*14-3
            out_line = "{:>11d}   {:11.5e}   {:>11s}   " + \
                       "{:>11d}   {:>11d}   {:11.5e}\n"

            # Write data
            for i in range(len(data.id)):
                for j in range(len(data.isotope_name)):
                    fp.write(out_line.
                             format(data.id[i], data.time[i],
                                    data.isotope_name[j],
                                    data.isotope_Z[j],
                                    data.isotope_A[j],
                                    data.yld[i,j]))

            # Close file
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_cluster_yield.bin', 'wb')

            # Write isotope data; note that we need to use struct to
            # force things to match up byte by byte, so that alignment
            # doesn't screw things
            fp.write(np.uint64(data.isotope_name.size))
            for i in range(data.isotope_name.size):
                tempstr = "{:<4s}".format(data.isotope_name[i])
                if sys.version_info < (3,):
                    fp.write(struct.pack('ccccII',
                                         tempstr[0],
                                         tempstr[1],
                                         tempstr[2],
                                         tempstr[3],
                                         data.isotope_Z[i],
                                         data.isotope_A[i]))
                else:
                    fp.write(struct.pack('ccccII',
                                         bytes(tempstr[0], "ascii"),
                                         bytes(tempstr[1], "ascii"),
                                         bytes(tempstr[2], "ascii"),
                                         bytes(tempstr[3], "ascii"),
                                         data.isotope_Z[i],
                                         data.isotope_A[i]))

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
                    fp.write(data.yld[k,:])
                    
                # Move pointer
                ptr = block_end

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Store isotope information in the first HDU
            niso = data.isotope_name.size
            isocols = []
            isocols.append(fits.Column(name="Name", format="3A", unit="",
                                       array=data.isotope_name))
            isocols.append(fits.Column(name="Z", format="1K", unit="",
                                       array=data.isotope_Z))
            isocols.append(fits.Column(name="A", format="1K", unit="",
                                       array=data.isotope_A))
            isofits = fits.ColDefs(isocols)
            isohdu = fits.BinTableHDU.from_columns(isofits)

            # Create a second HDU containing yield information
            yldcols = []
            yldcols.append(fits.Column(name="Trial", format="1K",
                                       unit="", array=data.trial))
            yldcols.append(fits.Column(name="UniqueID", format="1K",
                                       unit="", array=data.id))
            yldcols.append(fits.Column(name="Time", format="1D",
                                       unit="yr", array=data.time))
            yldcols.append(fits.Column(name="Yield",
                                       format="{:d}D".format(niso),
                                       unit="Msun",
                                       array=data.yld))
            yldfits = fits.ColDefs(yldcols)
            yldhdu = fits.BinTableHDU.from_columns(yldfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, isohdu, yldhdu])
            hdulist.writeto(model_name+'_cluster_yield.fits', 
                            overwrite=True)


    ################################################################
    # Write cloudy files if we have the data for them
    ################################################################
    if 'cloudy_hden' in data._fields:
        write_cluster_cloudyparams(data, model_name, fmt=fmt)
    if 'cloudy_inc' in data._fields:
        write_cluster_cloudyspec(data, model_name, fmt=fmt)
    if 'cloudy_linelum' in data._fields:
        write_cluster_cloudylines(data, model_name, fmt=fmt)
    if 'cloudy_phot_trans' in data._fields:
        write_cluster_cloudyphot(data, model_name, fmt=fmt)
