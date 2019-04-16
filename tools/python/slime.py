"""
SLug Ingestion of Mesa Evolution (SLIME) -- this script processes MIST
track sets taken directly from the MIST database and exports minimal
information from them that is sufficient to run SPS models. It writes
them out in gzip'ed FITS format, which vastly reduces file size.
"""

import argparse
import os
import os.path as osp
import glob
import numpy as np
from astropy.io import fits

# Read arguments
parser = argparse.ArgumentParser(
    description="SLug Ingestion of Mesa Evolution (SLIME); "
    " this script reads MIST database directories, extracts"
    " the minimal information required for a stellar population"
    " synthesis calculation, and writes it out in gzip'ed FITS"
    " file, thereby reducing the data size by a large factor")
parser.add_argument("directory", nargs="*",
                    help="names of MIST track directories "
                    "to process; default is to process all "
                    "subdirectories of current directory")
parser.add_argument("-o", "--outdir", help="name of directory into "
                    "which the processed outputs should be written; "
                    "default is current directory",
                    default=None)
parser.add_argument("-nc", "--noclobber", default=False,
                    action="store_true",
                    help="do not overwrite existing output files; "
                    "default behavior is to overwrite")
parser.add_argument("-v", "--verbose", help="produce verbose output",
                    default=False, action="store_true")
args = parser.parse_args()
cwd = os.getcwd()

# Get list of MIST directories to process
if len(args.directory) > 0:
    dirs = args.directory
else:
    dirs_tmp = glob.glob(osp.join(cwd, "MIST_*"))
    dirs = []
    for d in dirs_tmp:
        if not d.endswith('.gz') and not d.endswith('.zip') \
           and not d.endswith('.fits') and not d.endswith('.txz'):
            dirs.append(d)
if args.verbose:
    print("Directories to be processed:")
    for d in dirs:
        print(d)

# Loop through directories
for dirname in dirs:

    # Construct name of output file
    outfile = osp.basename(dirname)+'.fits.gz'
    if args.outdir is not None:
        outfile = osp.join(args.outdir, outfile)

    # If noclobber is set, skip if file exists
    if args.noclobber:
        if osp.isfile(outfile):
            if args.verbose:
                print("File "+outfile+" exists, skipping "
                      "directory "+dirname)
            continue

    # Print if verbose
    if args.verbose:
        print("Compressing directory "+dirname+" to "+
              outfile+"...")
    
    # Grab list of eep files and interpolated eep files
    fnames = glob.glob(osp.join(dirname, "*.track.eep*"))
    fnames.sort()

    # Replace incomplete tracks with the interpolated ones
    files = []
    for f in fnames:
        if f+"_INTERP" not in fnames:
            files.append(f)

    # Loop through files
    data = []
    for f in files:

        # Open file
        fp = open(f, "r")

        # Ingest the file header
        for i in range(5):
            line = fp.readline()
        splt = line[1:].split()
        y_init = float(splt[0])
        z_init = float(splt[1])
        fe_h = float(splt[2])
        a_fe = float(splt[3])
        v_vcrit = float(splt[4])
        for i in range(3):
            line = fp.readline()
        splt = line[1:].split()
        m_init = float(splt[0])
        n_pts = int(splt[1])
        n_eep = int(splt[2])
        line = fp.readline()
        splt = line[1:].split()
        eeps = np.array([int(s) for s in splt[1:]])
        for i in range(3):
            line = fp.readline()
        cols = line[1:].split()

        # Read remainder of file
        fdat = np.loadtxt(fp)

        # Close file
        fp.close()

        # Extract parts of data that we need to save
        t = fdat[:,cols.index('star_age')]
        m = fdat[:,cols.index('star_mass')]
        mdot = -fdat[:,cols.index('star_mdot')]  # Note sign convention change
        log_L = fdat[:,cols.index('log_L')]
        log_Teff = fdat[:,cols.index('log_Teff')]
        h_surf = fdat[:,cols.index('surface_h1')]
        he_surf = fdat[:,cols.index('surface_he3')] + \
                  fdat[:,cols.index('surface_he4')]
        c_surf = fdat[:,cols.index('surface_c12')] + \
                 fdat[:,cols.index('surface_c13')]
        n_surf = fdat[:,cols.index('surface_n14')]
        o_surf = fdat[:,cols.index('surface_o16')]
        phase = np.array(fdat[:,cols.index('phase')], dtype=int)

        # Store file data
        file_data = {
            'm_init' : m_init,
            'y_init' : y_init,
            'z_init' : z_init,
            'v_vcrit' : v_vcrit,
            'fe_h' : fe_h,
            'a_fe' : a_fe,
            'n_pts' : n_pts,
            'n_eep' : n_eep,
            'eeps' : eeps,
            'age' : t,
            'mass' : m,
            'mdot' : mdot,
            'log_L' : log_L,
            'log_Teff' : log_Teff,
            'h_surf' : h_surf,
            'he_surf' : he_surf,
            'c_surf' : c_surf,
            'n_surf' : n_surf,
            'o_surf' : o_surf,
            'phase' : phase,
            'interp' : f.endswith('_INTERP')
            }
        data.append(file_data)

    # Find maximum number of points in any file, and pad shorter data
    # to that length
    npts_max = np.amax([d['n_pts'] for d in data])
    flds = ['age', 'mass', 'mdot', 'log_L', 'log_Teff',
            'h_surf', 'he_surf', 'c_surf', 'n_surf', 'o_surf']
    for d in data:
        if d['n_pts'] < npts_max:
            for f in flds:
                d[f] = np.concatenate(
                    (d[f][:-1],
                     np.ones(npts_max-d['n_pts']+1)*d[f][-1]))
            d['phase'] = np.concatenate(
                (d['phase'][:-1],
                 np.ones(npts_max-d['n_pts']+1, dtype='int')*d['phase'][-1]))
        d['n_pts'] = npts_max

    # Create primary HDU to store metadata
    prihdr = fits.Header()
    prihdr['Y_init'] = data[0]['y_init']
    prihdr['Z_init'] = data[0]['z_init']
    prihdr['V_Vcrit'] = data[0]['v_vcrit']
    prihdr['Fe_H'] = data[0]['fe_h']
    prihdr['Alpha_Fe'] = data[0]['a_fe']
    prihdr['n_time'] = npts_max
    prihdr['n_mass'] = len(data)
    prihdr['m_min'] = data[0]['m_init']
    prihdr['m_max'] = data[-1]['m_init']
    prihdu = fits.PrimaryHDU(header=prihdr)

    # Create HDU that stores masses and information on which tracks
    # are interpolated
    mass_hdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='m_init', format='E',
                     array = np.array([d['m_init'] for d in data])),
         fits.Column(name='interpolated', format='I',
                     array = np.array([d['interp'] for d in data],
                                      dtype='int'))
         ])

    # Create table HDUs that store data for each star mass
    hdus = [prihdu, mass_hdu]
    for d in data:
        cols = []
        for f in flds:
            cols.append(fits.Column(name=f, format='E',
                                    array=d[f]))
        cols.append(fits.Column(name='phase', format='I',
                                array=d['phase']))
        hdus.append(fits.BinTableHDU.from_columns(cols))

    # Dump data to file
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfile, overwrite=True)

