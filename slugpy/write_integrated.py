"""
Function to write a set of output integrated files in SLUG2 format,
starting from an integrated data set as returned by
read_integrated. This can be used to translate data from one format to
another (e.g., bin to fits), or to consolidate multiple runs into a
single output file.
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
from .cloudy import write_integrated_cloudyparams
from .cloudy import write_integrated_cloudyphot
from .cloudy import write_integrated_cloudylines
from .cloudy import write_integrated_cloudyspec
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
    import warnings
    warnings.warn("Unable to import astropy. FITS funtionality" +
                  " will not be available.")


def write_integrated(data, model_name, fmt):
    """
    Function to write a set of output integrated files in SLUG2 format,
    starting from an integrated data set as returned by
    read_integrated.

    Parameters
       data : namedtuple
          Integrated data to be written, in the namedtuple format returned
          by read_integrated
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

    ################################################################
    # Write the properties file if we have the data for it
    ################################################################
    if 'target_mass' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_prop.txt', 'w')

            # Write header lines
            fp.write(("{:<14s}"*9).
                     format('Time', 'TargetMass', 'ActualMass',
                            'LiveMass', 'StellarMass', 'ClusterMass', 
                            'NumClusters', 'NumDisClust', 'NumFldStar',))
                            
            vp_i=0   
            adding_vcols = True                                     
            while adding_vcols == True:
                                        
                if 'VP'+repr(vp_i) in data._fields:            
                    fp.write(("{:<14s}"*1).format('VP'+repr(vp_i),))
                    vp_i+=1
                else:
                    nvp = vp_i
                    adding_vcols = False                   
                    fp.write("\n")            
                            
            fp.write(("{:<14s}"*9).
                     format('(yr)', '(Msun)', '(Msun)', '(Msun)',
                            '(Msun)', '(Msun)', '', '', '',))
            vp_i=0   
            adding_vcols = True  
            while adding_vcols == True:
                                        
                if 'VP'+repr(vp_i) in data._fields:            
                    fp.write(("{:<14s}"*1).format('',))
                    vp_i+=1
                else:
                    adding_vcols = False                   
                    fp.write("\n")
                            
                            
                            
            fp.write(("{:<14s}"*9).
                     format('-----------', '-----------', '-----------',
                            '-----------', '-----------', '-----------',
                            '-----------', '-----------',
                            '-----------',))
            vp_i=0   
            adding_vcols = True  
            while adding_vcols == True:
                                        
                if 'VP'+repr(vp_i) in data._fields:            
                    fp.write(("{:<14s}"*1).format('-----------',))
                    vp_i+=1
                else:
                    adding_vcols = False                   
                    fp.write("\n")         
            # Write data
            ntime = data.actual_mass.shape[-2]
            ntrial = data.actual_mass.shape[-1]
            if len(data.time) > ntime:
                random_time = True
            else:
                random_time = False
            for i in range(ntrial):
                if i != 0:
                    if nvp == 0:
                        fp.write("-"*(9*14-3)+"\n")
                    else:
                        fp.write("-"*((9+nvp)*14-3)+"\n")
                for j in range(ntime):
                    if random_time:
                        t_out = data.time[i]
                    else:
                        t_out = data.time[j]
                    fp.write(("{:11.5e}   {:11.5e}   {:11.5e}   " +
                              "{:11.5e}   {:11.5e}   {:11.5e}   " +
                              "{:11d}   {:11d}   {:11d}")
                             .format(t_out, 
                                     data.target_mass[j,i],
                                     data.actual_mass[j,i],
                                     data.live_mass[j,i],
                                     data.stellar_mass[j,i],
                                     data.cluster_mass[j,i],
                                     data.num_clusters[j,i],
                                     data.num_dis_clusters[j,i],
                                     data.num_fld_stars[j,i]))


                    vp_i=0   
                    adding_vcols = True  
                    while adding_vcols == True:                                                
                        if 'VP'+repr(vp_i) in data._fields:
                            current_vp=getattr(data, "VP"+repr(vp_i))
                            fp.write("   {:11.5e}".format(current_vp[j,i]))
                            vp_i+=1
                        else:
                            adding_vcols = False                   
                            fp.write("\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_prop.bin', 'wb')

            # Write out number of variable parameters
            nvp = 0
            for field in data._fields:
                if field.startswith("VP"):
                    nvp+=1
            fp.write( struct.pack('i', nvp) )

            # Write data
            ntime = data.actual_mass.shape[-2]
            ntrial = data.actual_mass.shape[-1]
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
                    fp.write(data.target_mass[j,i])
                    fp.write(data.actual_mass[j,i])
                    fp.write(data.live_mass[j,i])
                    fp.write(data.stellar_mass[j,i])
                    fp.write(data.cluster_mass[j,i])
                    fp.write(data.num_clusters[j,i])
                    fp.write(data.num_dis_clusters[j,i])
                    fp.write(data.num_fld_stars[j,i])
                    for field in sorted(data._fields):
                        if field.startswith("VP"):
                            fp.write(getattr(data, field)[j,i])

            # Close
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Figure out number of trials, and tile arrays
            ntrial = data.actual_mass.shape[-1]
            ntimes = data.actual_mass.shape[-2]
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
            cols.append(
                fits.Column(name="TargetMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.target_mass).flatten()))
            cols.append(
                fits.Column(name="ActualMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.actual_mass).flatten()))
            cols.append(
                fits.Column(name="LiveMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.live_mass).flatten()))
            cols.append(
                fits.Column(name="StellarMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.stellar_mass).flatten()))
            cols.append(
                fits.Column(name="ClusterMass", format="1D",
                            unit="Msun", 
                            array=np.transpose(data.cluster_mass).flatten()))
            cols.append(
                fits.Column(name="NumClusters", format="1K",
                            unit="", 
                            array=np.transpose(data.num_clusters).flatten()))
            cols.append(
                fits.Column(name="NumDisClust", format="1K",
                            unit="", 
                            array=np.transpose(data.num_dis_clusters).
                            flatten()))
            cols.append(
                fits.Column(name="NumFldStar", format="1K",
                            unit="", 
                            array=np.transpose(data.num_fld_stars).flatten()))
            
                                        
            vp_i=0   
            adding_vcols = True                                     
            while adding_vcols == True:
                                        
                if 'VP'+repr(vp_i) in data._fields:
                    cols.append(fits.Column(name="VP"+repr(vp_i), format="1D",
                                        unit="", array=np.transpose(getattr(data, "VP"+repr(vp_i)) ).flatten() ))
                    vp_i+=1
                else:

                    adding_vcols = False  
            
            
            fitscols = fits.ColDefs(cols)
            
            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_integrated_prop.fits',
                            clobber=True)

    ################################################################
    # Write the spectra file if we have the data for it
    ################################################################
    if 'spec' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_spec.txt', 'w')

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
                ifunc = interp1d(np.log(data.wl), np.log(data.spec), axis=0)
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
                                     np.log(data.spec_ex), 
                                     axis=0)
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
            line1 = ("{:<14s}"*3).format('Time', 'Wavelength',
                                         'L_lambda')
            line2 = ("{:<14s}"*3).format('(yr)', '(Angstrom)', 
                                         '(erg/s/A)')
            line3 = ("{:<14s}"*3).format('-----------', '-----------', 
                                         '-----------')
            sep_length = 3*14-3
            out_line = "{:11.5e}   {:11.5e}   {:11.5e}"

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
            ntime = spec.shape[-2]
            ntrial = spec.shape[-1]
            if len(data.time) > ntime:
                random_time = True
            else:
                random_time = False
            nl = len(wl)
            if 'spec_ex' in data._fields:
                offset = np.where(wl_ex[0] == wl)[0][0]
                nl_ex = len(wl_ex)
            else:
                offset = 0
                nl_ex = 0
            for i in range(ntrial):
                # Write trial separator
                if i != 0:
                    fp.write("-"*sep_length+"\n")
                # Write data for this time
                for j in range(ntime):
                    if random_time:
                        t_out = data.time[i]
                    else:
                        t_out = data.time[j]
                    for k in range(nl):
                        out_data = [t_out, wl[k], spec[k,j,i]]
                        if 'spec_neb' in data._fields:
                            out_data = out_data + [data.spec_neb[k,j,i]]
                        if k >= offset and k < offset + nl_ex:
                            out_data = out_data + [spec_ex[k-offset,j,i]]
                            if 'spec_neb_ex' in data._fields:
                                out_data = out_data + \
                                           [data.spec_neb_ex[k-offset,j,i]]
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

            fp = open(model_name+'_integrated_spec.bin', 'wb')

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

            # Write out times and spectra
            ntime = data.spec.shape[-2]
            ntrial = data.spec.shape[-1]
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
                    tmp = np.copy(data.spec[:,j,i])
                    fp.write(tmp)
                    if 'spec_neb' in data._fields:
                        tmp = np.copy(data.spec_neb[:,j,i])
                        fp.write(tmp)
                    if 'spec_ex' in data._fields:
                        tmp = np.copy(data.spec_ex[:,j,i])
                        fp.write(tmp)
                    if 'spec_neb_ex' in data._fields:
                        tmp = np.copy(data.spec_neb_ex[:,j,i])
                        fp.write(tmp)

            # Close file
            fp.close()

        elif fmt == 'fits':

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
            wlhdu = fits.BinTableHDU.from_columns(wlfits)

            # Figure out number of trials, and tile arrays
            ntrial = data.spec.shape[-1]
            ntimes = data.spec.shape[-2]
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            if len(data.time) > ntimes:
                times = data.time
            else:
                times = np.tile(data.time, ntrial)

            # Convert spectra to FITS columns, and make an HDU from
            # them
            speccols = []
            speccols.append(fits.Column(name="Trial", format="1K",
                                        unit="", array=trial))
            speccols.append(fits.Column(name="Time", format="1D",
                                        unit="yr", array=times))
            speccols.append(fits.Column(name="L_lambda",
                                        format=fmtstring,
                                        unit="erg/s/A",
                                        array=np.transpose(data.spec).
                                        reshape(ntimes*ntrial, nl)))
            if 'spec_neb' in data._fields:
                speccols.append(
                    fits.Column(name="L_lambda_neb",
                                format=fmtstring_neb,
                                unit="erg/s/A",
                                array=np.transpose(data.spec_neb).
                                reshape(ntimes*ntrial, nl_neb)))
            if 'spec_ex' in data._fields:
                speccols.append(
                    fits.Column(name="L_lambda_ex",
                                format=fmtstring_ex,
                                unit="erg/s/A",
                                array=np.transpose(data.spec_ex).
                                reshape(ntimes*ntrial, nl_ex)))
            if 'spec_neb_ex' in data._fields:
                speccols.append(
                    fits.Column(name="L_lambda_neb_ex",
                                format=fmtstring_neb_ex,
                                unit="erg/s/A",
                                array=np.transpose(data.spec_neb_ex).
                                reshape(ntimes*ntrial, nl_neb_ex)))
            specfits = fits.ColDefs(speccols)
            spechdu = fits.BinTableHDU.from_columns(specfits)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, wlhdu, spechdu])
            hdulist.writeto(model_name+'_integrated_spec.fits', 
                            clobber=True)
                
    ################################################################
    # Write photometry file if we have the data for it
    ################################################################
    if 'phot' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_phot.txt', 'w')

            # Write header lines
            fp.write("{:<21s}".format('Time'))
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
            fp.write("{:<21s}".format('(yr)'))
            for i in range(fac):
                for f in data.filter_units:
                    fp.write("({:s}".format(f)+")"+" "*(19-len(f)))
            fp.write("\n")
            fp.write("{:<21s}".format('------------------'))
            nf = len(data.filter_names)
            for j in range(fac):
                for i in range(nf):
                    fp.write("{:<21s}".format('------------------'))
            fp.write("\n")

            # Write data
            ntime = data.phot.shape[1]
            ntrial = data.phot.shape[2]
            if len(data.time) > ntime:
                random_time = True
            else:
                random_time = False
            for i in range(ntrial):
                # Write separator between trials
                if i != 0:
                    fp.write("-"*((1+fac*nf)*21-3)+"\n")
                for j in range(ntime):
                    if random_time:
                        t_out = data.time[i]
                    else:
                        t_out = data.time[j]
                    fp.write("       {:11.5e}".format(t_out))
                    for k in range(nf):
                        fp.write("          {:11.5e}".
                                 format(data.phot[k,j,i]))
                    if 'phot_neb' in data._fields:
                        for k in range(nf):
                            fp.write("          {:11.5e}".
                                     format(data.phot_neb[k,j,i]))
                    if 'phot_ex' in data._fields:
                        for k in range(nf):
                            if np.isnan(data.phot_ex[k,j,i]):
                                fp.write("          {:11s}".format(""))
                            else:
                                fp.write("          {:11.5e}".
                                         format(data.phot_ex[k,j,i]))
                    if 'phot_neb_ex' in data._fields:
                        for k in range(nf):
                            if np.isnan(data.phot_neb_ex[k,j,i]):
                                fp.write("          {:11s}".format(""))
                            else:
                                fp.write("          {:11.5e}".
                                         format(data.phot_neb_ex[k,j,i]))
                    fp.write("\n")

            # Close
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_phot.bin', 'wb')

            # Write number of filters and filter names as ASCII
            nf = len(data.filter_names)
            if sys.version_info < (3,):
                # Python 2 version
                fp.write(str(nf)+"\n")
                for i in range(nf):
                    fp.write(data.filter_names[i] + " " + 
                             data.filter_units[i] + "\n")
            else:
                # Python 3 version
                fp.write((str(nf)+"\n").encode())
                for i in range(nf):
                    fp.write((data.filter_names[i] + " " + 
                              data.filter_units[i] + "\n").encode())

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

            # Write data
            ntime = data.phot.shape[1]
            ntrial = data.phot.shape[2]
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
                    tmp = np.copy(data.phot[:,j,i])
                    fp.write(tmp)
                    if 'phot_neb' in data._fields:
                        tmp = np.copy(data.phot_neb[:,j,i])
                        fp.write(tmp)
                    if 'phot_ex' in data._fields:
                        tmp = np.copy(data.phot_ex[:,j,i])
                        fp.write(tmp)
                    if 'phot_neb_ex' in data._fields:
                        tmp = np.copy(data.phot_neb_ex[:,j,i])
                        fp.write(tmp)

            # Close file
            fp.close()

        elif fmt == 'fits':

            ########################################################
            # FITS mode
            ########################################################

            # Figure out number of trials, and tile arrays
            ntimes = data.phot.shape[1]
            ntrial = data.phot.shape[2]
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            if len(data.time) > ntimes:
                times = data.time
            else:
                times = np.tile(data.time, ntrial)
            nf = len(data.filter_names)

            # Convert data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=trial))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=times))
            for i in range(len(data.filter_names)):
                cols.append(
                    fits.Column(name=data.filter_names[i],
                                unit=data.filter_units[i],
                                format="1D",
                                array=np.transpose(data.phot[i,:,:]).
                                flatten()))
            if 'phot_neb' in data._fields:
                for i in range(len(data.filter_names)):
                    cols.append(
                        fits.Column(name=data.filter_names[i]+'_neb',
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=np.transpose(data.phot_neb[i,:,:]).
                                    flatten()))
            if 'phot_ex' in data._fields:
                for i in range(len(data.filter_names)):
                    cols.append(
                        fits.Column(name=data.filter_names[i]+'_ex',
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=np.transpose(data.phot_ex[i,:,:]).
                                    flatten()))
            if 'phot_neb_ex' in data._fields:
                for i in range(len(data.filter_names)):
                    cols.append(
                        fits.Column(name=data.filter_names[i]+'_neb_ex',
                                    unit=data.filter_units[i],
                                    format="1D",
                                    array=np.transpose(data.phot_neb_ex[i,:,:]).
                                    flatten()))
            fitscols = fits.ColDefs(cols)

            # Create the binary table HDU
            tbhdu = fits.BinTableHDU.from_columns(fitscols)

            # Create dummy primary HDU
            prihdu = fits.PrimaryHDU()

            # Create HDU list and write to file
            hdulist = fits.HDUList([prihdu, tbhdu])
            hdulist.writeto(model_name+'_integrated_phot.fits',
                            clobber=True)

    ################################################################
    # Write the supernova file if we have the data for it
    ################################################################
    if 'tot_sn' in data._fields:

        if fmt == 'ascii' or fmt == 'txt':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_sn.txt', 'w')

            # Write header
            fp.write(("{:<14s}"*3+"\n").
                     format('Time', 'TotSN', 'StochSN'))
            fp.write(("{:<14s}"*3+"\n").format('(yr)', '', ''))
            fp.write(("{:<14s}"*3+"\n").format('-----------', '-----------', 
                                               '-----------'))
            sep_length = 3*14-3
            out_line = "{:11.5e}   {:>11.5e}   {:>11d}\n"
             
            # Write data
            ntime = data.tot_sn.shape[0]
            ntrial = data.tot_sn.shape[1]
            if len(data.time) > ntime:
                random_time = True
            else:
                random_time = False
            for i in range(ntrial):
                # Write separator between trials
                if i != 0:
                    fp.write('-'*sep_length+'\n')
                # Write data for this time
                for j in range(ntime):
                    if random_time:
                        t_out = data.time[i]
                    else:
                        t_out = data.time[j]
                    fp.write(out_line.format(
                        t_out, data.tot_sn[j,i], data.stoch_sn[j,i]))

            # Close file
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_sn.bin', 'wb')
            
            # Write data
            ntime = data.tot_sn.shape[-2]
            ntrial = data.tot_sn.shape[-1]
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
                    fp.write(data.tot_sn[j,i])
                    fp.write(data.stoch_sn[j,i])

            # Close file
            fp.close()
                    
        elif fmt == 'fits':
            
            ########################################################
            # FITS mode
            ########################################################
           
            # Figure out number of times and trials, and tile the
            # arrays
            ntimes = data.tot_sn.shape[0]
            ntrial = data.tot_sn.shape[1]
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            if len(data.time) > ntimes:
                times = data.time
            else:
                times = np.tile(data.time, ntrial)

            # Convert yield data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=trial))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=times))
            cols.append(fits.Column(name="TotSN", format="1D",
                                    unit="",
                                    array=np.transpose(data.tot_sn).flatten()))
            cols.append(fits.Column(name="StochSN", format="1K",
                                    unit="",
                                    array=np.transpose(data.tot_sn).flatten()))
            snfits = fits.ColDefs(cols)
            snhdu = fits.BinTableHDU.from_columns(snfits)

            # Create HDU list and write file
            prihdu = fits.PrimaryHDU()
            hdulist = fits.HDUList([prihdu, snhdu])
            hdulist.writeto(model_name+'_integrated_sn.fits',
                            clobber=True)
            
    ################################################################
    # Write yield file if we have the data for it
    ################################################################
    if 'yld' in data._fields:

        if fmt == 'ascii':

            ########################################################
            # ASCII mode
            ########################################################

            fp = open(model_name+'_integrated_yield.txt', 'w')

            # Write header
            fp.write(("{:<14s}"*5+"\n").format('Time', 'Symbol',
                                               'Z', 'A', 'Yield'))
            fp.write(("{:<14s}"*5+"\n").format('(yr)', '', '', '', '(Msun)'))
            fp.write(("{:<14s}"*5+"\n").format('-----------', '-----------', 
                                               '-----------', '-----------', 
                                               '-----------'))
            sep_length = 5*14-3
            out_line = "{:11.5e}   {:>11s}   {:>11d}   {:>11d}   {:11.5e}\n"
           
            # Write data
            ntime = data.yld.shape[1]
            ntrial = data.yld.shape[2]
            niso = data.yld.shape[0]
            if len(data.time) > ntime:
                random_time = True
            else:
                random_time = False
            for i in range(ntrial):
                # Write separator between trials
                if i != 0:
                    fp.write('-'*sep_length+'\n')
                # Write data for this time
                for j in range(ntime):
                    if random_time:
                        t_out = data.time[i]
                    else:
                        t_out = data.time[j]
                    for k in range(niso):
                        fp.write(out_line.format(
                            t_out, data.isotope_name[k],
                            data.isotope_Z[k],
                            data.isotope_A[k],
                            data.yld[k,j,i]))

            # Close file
            fp.close()

        elif fmt == 'bin' or fmt == 'binary':

            ########################################################
            # Binary mode
            ########################################################

            fp = open(model_name+'_integrated_yield.bin', 'wb')

            # Write isotope data; note that we need to use struct to
            # force things to match up byte by byte, so that alignment
            # doesn't screw things
            fp.write(np.uint64(data.isotope_name.size))
            for i in range(data.isotope_name.size):
                tempstr = "{:<4s}".format(data.isotope_name[i])
                if sys.version_info < (3,):
                    # Python 2 version
                    fp.write(struct.pack('ccccII',
                                         tempstr[0],
                                         tempstr[1],
                                         tempstr[2],
                                         tempstr[3],
                                         data.isotope_Z[i],
                                         data.isotope_A[i]))
                else:
                    # Python 3 version
                    fp.write(struct.pack('ccccII',
                                         bytes(tempstr[0].encode()),
                                         bytes(tempstr[1].encode()),
                                         bytes(tempstr[2].encode()),
                                         bytes(tempstr[3].encode()),
                                         data.isotope_Z[i],
                                         data.isotope_A[i]))

            # Write remainder of data
            ntime = data.yld.shape[1]
            ntrial = data.yld.shape[2]
            if len(data.time) > ntime:
                random_time = True
            else:
                random_time = False
            for i in range(ntrial):
                for j in range(ntime):
                    # Write trial number and time
                    fp.write(np.uint(i))
                    if random_time:
                        fp.write(data.time[i])
                    else:
                        fp.write(data.time[j])
                    # Write yields
                    tmp = np.copy(data.yld[:,j,i])
                    fp.write(tmp)

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

            # Figure out number of times and trials, and tile the
            # arrays
            ntimes = data.yld.shape[1]
            ntrial = data.yld.shape[2]
            trial = np.transpose(np.tile(
                np.arange(ntrial, dtype='int64'), (ntimes,1))).\
                flatten()
            if len(data.time) > ntimes:
                times = data.time
            else:
                times = np.tile(data.time, ntrial)

            # Convert yield data to FITS columns
            cols = []
            cols.append(fits.Column(name="Trial", format="1K",
                                    unit="", array=trial))
            cols.append(fits.Column(name="Time", format="1D",
                                    unit="yr", array=times))
            yldtmp = np.transpose(data.yld)
            cols.append(fits.Column(name="Yield",
                                    format="{:d}D".format(niso),
                                    unit="Msun", 
                                    array=yldtmp.reshape((ntrial*ntimes,
                                                          niso))))
            yldfits = fits.ColDefs(cols)
            yldhdu = fits.BinTableHDU.from_columns(yldfits)

            # Create HDU list and write file
            prihdu = fits.PrimaryHDU()
            hdulist = fits.HDUList([prihdu, isohdu, yldhdu])
            hdulist.writeto(model_name+'_integrated_yield.fits',
                            clobber=True)

    ################################################################
    # Write cloudy files if we have the data for them
    ################################################################
    if 'cloudy_hden' in data._fields:
        write_integrated_cloudyparams(data, model_name, fmt=fmt)
    if 'cloudy_inc' in data._fields:
        write_integrated_cloudyspec(data, model_name, fmt=fmt)
    if 'cloudy_linelum' in data._fields:
        write_integrated_cloudylines(data, model_name, fmt=fmt)
    if 'cloudy_phot_trans' in data._fields:
        write_integrated_cloudyphot(data, model_name, fmt=fmt)
