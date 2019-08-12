"""
SLug Ingestion of Mesa Evolution (SLIME) -- this script processes MIST
track sets taken directly from the MIST database and exports minimal
information from them that is sufficient to run SPS models. It writes
them out in gzip'ed FITS format, which vastly reduces file size. A
later version extends this to the PARSEC tracks with EEPS produced by
Phil Rosenfeld (https://philrosenfield.github.io/padova_tracks/)
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
    " file, thereby reducing the data size by a large factor."
    " This script can also perform similar processing on the"
    " PARSEC tracks with EEPS")
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
parser.add_argument("-pc", "--parsec", default=False,
                    action="store_true",
                    help="process PARSEC tracks instead of MIST "
                    "tracks")
parser.add_argument("-v", "--verbose", help="produce verbose output",
                    default=False, action="store_true")
args = parser.parse_args()
cwd = os.getcwd()

# Get list of directories to process
if len(args.directory) > 0:
    dir_tmp = args.directory
    # Trim final separators
    dirs = []
    for d in dir_tmp:
        if d[-1] == os.sep: dirs.append(d[:-1])
        else: dirs.append(d)
else:
    if not args.parsec:
        dirs_tmp = glob.glob(osp.join(cwd, "MIST_*"))
        dirs = []
        for d in dirs_tmp:
            if not d.endswith('.gz') and not d.endswith('.zip') \
               and not d.endswith('.fits') and not d.endswith('.txz'):
                dirs.append(d)
            else:
                dirs_tmp = glob.glob(osp.join(cwd, "Z*Y*"))
                dirs = [ d for d in dirs_tmp if osp.isdirectory(d) ]
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

    # MIST procedure
    if not args.parsec:
        
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
            
    else:

        # PARSEC tracks procedure

        # Extract Y and Z from name of directory
        dname = osp.basename(dirname)
        y_init = float(dname[1:].split("Y")[1])
        z_init = float(dname[1:].split("Y")[0])

        # Grab list of all files
        fnames = glob.glob(osp.join(dirname, "match*.dat"))

        # Grab list of HB files
        hbnames = glob.glob(osp.join(dirname, "match*.HB.dat"))

        # Get list of non-HB files
        files = [ f for f in fnames if f not in hbnames ]

        # Get list of masses to go with file names
        m_f = [ float(f.split("M")[-1][:-4]) for f in files ]
        m_hb = [ float(f.split("M")[-1][:-7]) for f in hbnames ]

        # Sort by mass
        idx = np.argsort(m_f)
        m_f = list(np.array(m_f)[idx])
        files = list(np.array(files)[idx])

        # Loop through files
        data = []
        for m_init, f in zip(m_f, files):

            # Open file, burn header line, then import remainder
            fp = open(f, "r")
            fp.readline()
            fdat = np.loadtxt(fp)
            fp.close()

            # Extract the data we need
            t = 10.**fdat[:,0]   # Age in file is stored as log age
            m = fdat[:,1]
            log_Teff = fdat[:,2]
            log_L = (4.77 - fdat[:,3])/2.5 # File contains M_bol, we
                                           # want L
            # PLACEHOLDERS
            mdot = np.zeros(t.shape)
            h_surf = np.zeros(t.shape)
            he_surf = np.zeros(t.shape)
            c_surf = np.zeros(t.shape)
            n_surf = np.zeros(t.shape)
            o_surf = np.zeros(t.shape)

            # Assign phases
            phase = np.zeros(t.shape, dtype='int')
            phase[:200] = -1          # Pre-main sequence
            phase[200:600] = 0        # Main sequence
            if t.size > 600:
                phase[600:1100] = 2   # Red giant branch
            if t.size > 1100:
                phase[1100:1630] = 3  # Core He fuction branch
                phase[1630:1730] = 4  # EAGB, post-core He exhaustion
            if t.size > 1730:
                phase[1730:1930] = 5  # TP-AGB

            # See if this track has a horizontal branch file to go
            # with it; note that the track list includes a many more
            # masses at the HB than for the main files, so that if one
            # wanted to assume some mass loss on the RGB or during the
            # He flash, one could do so. However, this has never been
            # implemented, so we will assume no mass loss here, and
            # just match the masses.
            if m_init in m_hb:

                # Open file, burn header line, then import remainder
                fp = open(hbnames[m_hb.index(m_init)], "r")
                fp.readline()
                fdat = np.loadtxt(fp)
                fp.close()
                
                # Extract the data we need; note that time in the HB
                # file is stored as time since start of HB, which we
                # will assume is coincident with last time in the
                # corresponding file running up to the He flash
                t1 = t[-1] + 10.**fdat[:,0]
                m1 = fdat[:,1]
                log_Teff1 = fdat[:,2]
                log_L1 = (4.77 - fdat[:,3])/2.5
                # PLACEHOLDERS
                mdot1 = np.zeros(t1.shape)
                h_surf1 = np.zeros(t1.shape)
                he_surf1 = np.zeros(t1.shape)
                c_surf1 = np.zeros(t1.shape)
                n_surf1 = np.zeros(t1.shape)
                o_surf1 = np.zeros(t1.shape)
            
                # Contatenate on to existing data
                t = np.concatenate([t, t1])
                m = np.concatenate([m, m1])
                log_Teff = np.concatenate([log_Teff, log_Teff1])
                log_L = np.concatenate([log_L, log_L1])
                mdot = np.concatenate([mdot, mdot1])
                h_surf = np.concatenate([h_surf, h_surf1])
                he_surf = np.concatenate([he_surf, he_surf1])
                c_surf = np.concatenate([c_surf, c_surf1])
                n_surf = np.concatenate([n_surf, n_surf1])
                o_surf = np.concatenate([o_surf, o_surf1])
                phase = np.concatenate([phase,
                                        np.zeros(t1.size, dtype='int')+3])
        
            # Store file data
            file_data = {
                'm_init' : m_init,
                'y_init' : y_init,
                'z_init' : z_init,
                'n_pts' : t.size,
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
                'phase' : phase
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
                 np.ones(npts_max-d['n_pts']+1, dtype='int') * \
                 d['phase'][-1]))
        d['n_pts'] = npts_max

    # Create primary HDU to store metadata
    prihdr = fits.Header()
    prihdr['Y_init'] = data[0]['y_init']
    prihdr['Z_init'] = data[0]['z_init']
    if not args.parsec:
        prihdr['V_Vcrit'] = data[0]['v_vcrit']
        prihdr['Fe_H'] = data[0]['fe_h']
        prihdr['Alpha_Fe'] = data[0]['a_fe']
    prihdr['n_time'] = npts_max
    prihdr['n_mass'] = len(data)
    prihdr['m_min'] = data[0]['m_init']
    prihdr['m_max'] = data[-1]['m_init']
    prihdu = fits.PrimaryHDU(header=prihdr)

    # Create HDU that stores masses and information on which tracks
    # are interpolated (for MIST)
    if not args.parsec:
        mass_hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='m_init', format='E',
                         array = np.array([d['m_init'] for d in data])),
             fits.Column(name='interpolated', format='I',
                         array = np.array([d['interp'] for d in data],
                                          dtype='int'))
            ])
    else:
        mass_hdu = fits.BinTableHDU.from_columns(
            [fits.Column(name='m_init', format='E',
                         array = np.array([d['m_init'] for d in data]))])

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
