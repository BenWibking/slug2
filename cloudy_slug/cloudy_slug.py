"""
This is a script that takes spectra output by SLUG and calls cloudy on
them in order to calculate the resulting nebular emission.
"""

######################################################################
# Import libraries
######################################################################

import argparse
try:
    import astropy.io.fits as fits
except ImportError:
    fits = None
from collections import namedtuple
import copy
import multiprocessing
import numpy as np
import os
import os.path as osp
try:
    from Queue import Queue    # python 2.x
except ImportError:
    from queue import Queue    # python 3.x
import subprocess
import sys
from threading import Thread
from time import sleep
import warnings
import errno
try:
    from slugpy import *    # If slugpy is already in our path
    from slugpy.cloudy import *
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        from slugpy.cloudy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

######################################################################
# Set some constants; change to cgs units
######################################################################

from scipy.constants import c
from scipy.constants import k as kB
from scipy.constants import m_e
from scipy.constants import m_p
from scipy.constants import physical_constants as physcons
c = c*1e2
kB = kB*1e7
m_e = m_e * 1e-3
m_p = m_p * 1e-3
eps0 = physcons['Rydberg constant times hc in J'][0] * 1e7
mH = m_e + m_p       # Hydrogen atom mass
alphaB = 2.59e-13    # Case B recombination coefficient
muH = 1.4            # Mean mass per H nucleus for standard cosmic composition
fe = 1.1             # Electrons per H nuclues


######################################################################
# Step 1: set up and read command line arguments
######################################################################
parser = argparse. \
         ArgumentParser(
             description="Script to run cloudy on slug outputs")

# Positional arguments
parser.add_argument("slug_model_name", 
                    help="name of the SLUG model output to be " +
                    "processed")
parser.add_argument("start_spec", nargs="?", type=int, default=0,
                    help="starting cluster or trial number (default: 0)")
parser.add_argument("end_spec", nargs="?", type=int, default=-1,
                    help="ending cluster or trial number " + 
                    "(default: last spectrum)")

# Optional arguments
parser.add_argument('-a', '--agemax', default=10, type=float,
                    help="maximum cluster age in Myr for which to " +
                    "compute nebular emission; only used in clustermode " +
                    "(default: 10 Myr)")
parser.add_argument("--cloudypath", default=None, type=str,
                    help="path to the cloudy executable (default: "+
                    "$CLOUDY_DIR/cloudy.exe)")
parser.add_argument("--cloudytemplate", default=None, type=str,
                    help="template cloudy input file (default: "+
                    "$SLUG_DIR/cloudy_slug/cloudy.in_template)")
parser.add_argument("-cm", "--clustermode", action='store_true',
                    default=False, help="process cluster by "+
                    "by cluster spectra through cloudy "+
                    "(default: integrated mode, process integrated "
                    "spectra only)")
parser.add_argument('-cf', '--coveringfac', default=None, type=str,
                    help="covering factor; can be a PDF")
parser.add_argument('-d', '--dynamic', action='store_true', default=False,
                    help='use dynamic mode to compute HII region ' +
                    'radii; only allowed in combination with with clustermode')
parser.add_argument('-hd', '--hden', default=None, type=str,
                    help='hydrogen number density in cm^-3; either '+
                    'mean internal density for calculations not in '+
                    'dynamic mode, or ambient density if dynamic '+
                    'mode is set; can be a PDF')
parser.add_argument('-ip', '--ionparam', default=None, type=str,
                    help="volume-averaged ionization parameter; can be a PDF")
parser.add_argument('-ip0', '--ionparam0', default=None, type=str,
                    help="inner radius ionization parameter; "+
                    "can be a PDF")
parser.add_argument('--ionparammax', default=None, type=str,
                    help="maximum inner radius ionization parameter"+
                    "; if value exceeds this, response is set "+
                    "by the paramsafety option selected")
parser.add_argument('--ionparammin', default=None, type=str,
                    help="minimum inner radius ionization parameter"+
                    "; if value exceeds this, response is set "+
                    "by the paramsafety option selected")
parser.add_argument('-nl', '--nicelevel', default=0, type=int,
                    help="nice level of the cloudy processes " +
                    "(default: 0)")
parser.add_argument('-n', '--nproc', default=None, type=int,
                    help="number of cloudy processes (default: "+
                    "number of logical cores)")
parser.add_argument('-ps', '--paramsafety', default='warn',
                    help="parameter safety setting; options are: " +
                    "'warn' -- forbidden parameter combinations " +
                    "are corrected and a warning is issued; " +
                    "'skip' -- forbidden parameter combinations " +
                    "are skipped, exactly as for cases where " +
                    "Q(H0) < qH0min; 'halt' -- code halts if a " +
                    "forbidden parameter combination is entered; " +
                    "'redraw' -- if values were drawn from a PDF, " +
                    "draw again (up to 100 tries), then skip if "+
                    "still unsuccessful; " +
                    "default behavior is 'warn'")
parser.add_argument('-fix', '--fix_quantity', default=None, type=str,
                    help="when fixing out of range parameters, "+
                    "which quantity should be corrected; see help "+
                    "for hiiregparam class for details; only uses "+
                    "if paramsafety is set to warn")
parser.add_argument('-qm', '--qH0min', default=0.0, type=float,
                    help="minimum ionizing luminosity for which to "+
                    "compute nebular emission (default = 0)")
parser.add_argument('-r0', default=None, type=str,
                    help='inner radius of HII region in cm; can be a PDF')
parser.add_argument('-r1', default=None, type=str,
                    help='outer radius of HII region in cm; can be a PDF')
parser.add_argument('-s', '--save', default=False,
                    action='store_true', help='save full cloudy ' +
                    'output (default: delete after extracting data)')
parser.add_argument('--slugformat', default=None, type=str,
                    help="format of SLUG output data to use (valid "+
                    "values = 'ascii', 'bin' / 'binary', or 'fits'; "+
                    "default behavior = detect automatically")
parser.add_argument("--slugpath", default=None, type=str,
                    help="path to the SLUG output data (default: "+
                    "check cwd, then $SLUG_DIR/output)")
parser.add_argument('-t', '--tmpdir', default=None, type=str,
                    help="location of directory in which to store "+
                    "temporary files; default: " +
                    "./cloudy_tmp_MODEL_NAME")
parser.add_argument('-v', '--verbose', action='store_true',
                    default=False, help="produce verbose output")
parser.add_argument('-wp', '--windparam', default=None, type=str,
                    help='wind parameter Omega; can be a PDF')
parser.add_argument('-wr', '--writeparams', default=False,
                    action='store_true',
                    help='write run parameters as text in the cloudy'+
                    'output directory; only applied if --save is set')
parser.add_argument('-wf', '--warnfail', default=False,
                    action='store_true',
                    help='on failure to read the cloudy output '+
                    '(usually a result of cloudy crashing), '+
                    'issue a warning and exit normally rather '+
                    'than raising an error; useful for scripting')
args = parser.parse_args()


######################################################################
# Step 2: set paths
######################################################################

# Get current directory
cwd = os.getcwd()

# Set path to cloudy
if args.cloudypath is None:
    if 'CLOUDY_DIR' in os.environ:
        cloudypath = osp.join(os.environ['CLOUDY_DIR'], 'cloudy.exe')
    else:
        cloudypath = 'cloudy.exe'
else:
    cloudypath = args.cloudypath
if args.nproc is None:
    nproc = multiprocessing.cpu_count()
else:
    nproc = args.nproc

# Set location of cloudy template
if args.cloudytemplate is None:
    if 'SLUG_DIR' in os.environ:
        cloudytemplate = osp.join(os.environ['SLUG_DIR'],
                                  'cloudy_slug',
                                  'cloudy.in_template')
    else:
        cloudytemplate = osp.join('cloudy_slug',
                                  'cloudy.in_template')
else:
    cloudytemplate = args.cloudytemplate

# Set input directory
if args.slugpath is None:
    if 'SLUG_DIR' in os.environ:
        slugdir=osp.join(os.environ['SLUG_DIR'], 'output')
    else:
        slugdir='output'
else:
    slugdir=args.slugpath

# Set temporary directory
basename = osp.basename(args.slug_model_name)
if args.tmpdir is None:
    tmpdirname = osp.join(cwd, 'cloudy_tmp_'+basename)
else:
    tmpdirname = args.tmpdir

    
######################################################################
# Step 3: validate inputs
#
# Here we have to make sure we have enough information to proceed; we
# need one of these conditions to be met
# (1) a density and radius command are both in the template;
# (2) values for 2/6 of the following command line options: hden,
#     ionparam, ionparam0, r0, r1, windparam (but not the disallowed
#     combinations r0, ionparam and r1, ionparam0); or
# (3) dynamic must be true and hden and windparam must be set
#
# We issue warnings if parameter are double-specified.
######################################################################

# Check that the parameter safety setting is valid
if args.paramsafety != 'warn' and \
   args.paramsafety != 'skip' and \
   args.paramsafety != 'halt' and \
   args.paramsafety != 'redraw':
    raise ValueError("cloudy_slug: paramsafety must be 'warn', "+
                     "'skip', 'halt', or 'redraw'")

# Check if the template contains hden and radius commands
fp = open(cloudytemplate, 'r')
tempfile = fp.read().split('\n')
fp.close()
hden_template = None
radius_template = None
coveringfac_template = None
for line in tempfile:
    # Trim comments
    ltrim = line.partition(';')[0]
    ltrim = line.partition('//')[0]
    ltrim = line.partition('#')[0]
    ltrim = line.partition('%')[0]
    # Skip blank lines and lines that start with c, indicating a
    # comment
    if len(line) == 0:
        continue
    if line[0] == 'c':
        continue
    # Check for contents
    if line.split()[0] == 'hden':
        hden_template = 10.**float(line.split()[1])
    elif line.split()[0] == 'radius':
        radius_template = 10.**float(line.split()[1])
    elif line.split()[0] == 'covering' and \
         line.split()[1] == 'factor':
        coveringfac_template = float(line.split()[2])
        if args.coveringfac is not None:
            warnings.warn(
                "cloudy_slug: found covering factor "
                "command in " + cloudytemplate +
                " but also got command line argument; "
                "value in template file will be ignored")
template_complete = (hden_template is not None) and \
                    (radius_template is not None)

# Now check command line values for physical parameters; if we're not
# in dynamic mode, we need there to be exactly 0 or 2 of these. Issue
# errors and warnings as appropriate.
cmd_complete = False
if not args.dynamic:

    # See which parameters we got
    hden_set = args.hden is not None
    U_set = args.ionparam is not None
    U0_set = args.ionparam0 is not None
    r0_set = args.r0 is not None
    r1_set = args.r1 is not None
    Omega_set = args.windparam is not None
    nset = hden_set + U_set + r0_set + r1_set + Omega_set + U0_set

    # Flag on incorrect number or combination
    if nset != 0 and nset != 2:
        raise ValueError(
            "cloudy_slug: for parameters hden, ionparam, "+
            "ionparam0, r0, r1, windparam," +
            "can set either 0 or 2, not any other number")
    if r0_set and U_set:
        raise ValueError(
            "cloudy_slug: cannot use the parameter combination r0, "+
            "ionparam; this combination does not define unique parameters")
    if r1_set and U0_set:
        raise ValueError(
            "cloudy_slug: cannot use the parameter combination r1, "+
            "ionparam0; this combination does not define unique parameters")

    # Warn about double-setting
    if hden_template is not None and nset != 0:
        warnings.warn(
            "cloudy_slug: found hden command in " +
            cloudytemplate + " but also got 2 nebular " +
            "parameters from command line; calculation will " +
            "proceed, but hden command in template will be ignored")
    if radius_template is not None and nset != 0:
        warnings.warn(
            "cloudy_slug: found radius command in " +
            cloudytemplate + " but also got 2 nebular " +
            "parameters from command line; calculation will " +
            "proceed, but hden command in template will be ignored")

    # Flag if command line parameters are properly set
    if nset == 2:
        cmd_complete = True

else:
    # Dynamic mode; make sure we are in cluster mode, and that we have
    # an ambient density and wind paramater
    if not args.clustermode:
        raise ValueError("cloudy_slug: can only use --dynamic in "+
                         "conjunction with --clustermode")
    if args.hden is None or args.windparam is None:
        raise ValueError("cloudy_slug: in dynamic mode, need to "
                         "set both hden and windparam")
    cmd_complete = True

# Make sure we got what we needed at the command line or in the
# template; if not, bail out
if not (cmd_complete or template_complete):
    raise ValueError("cloudy_slug: insufficient input information "+
                     "to proceed; either set hden and radius in "+
                     "template, or set command line parameters to "+
                     "specify physical conditions")
    

######################################################################
# Step 4: read the SLUG output to be processed, and check that it
# contains all the required data
######################################################################
file_info = {}
if (args.clustermode):
    if args.slugpath is None:
        # If no directory is set, check cwd, then default directory
        try:
            data = read_cluster(args.slug_model_name,
                                fmt=args.slugformat,
                                read_info=file_info,
                                output_dir=cwd)
        except IOError:
            data = read_cluster(args.slug_model_name,
                                fmt=args.slugformat,
                                read_info=file_info,
                                output_dir=slugdir)
    else:
        data = read_cluster(args.slug_model_name,
                            fmt=args.slugformat,
                            read_info=file_info,
                            output_dir=slugdir)
else:
    if args.slugpath is None:
        # If no directory is set, check cwd, then default directory
        try:
            data = read_integrated(args.slug_model_name,
                                   fmt=args.slugformat,
                                   read_info=file_info,
                                   output_dir=cwd)
        except IOError:
            data = read_integrated(args.slug_model_name,
                                   fmt=args.slugformat,
                                   read_info=file_info,
                                   output_dir=slugdir)
    else:
        data = read_integrated(args.slug_model_name,
                               fmt=args.slugformat,
                               read_info=file_info,
                               output_dir=slugdir)
outpath = osp.dirname(file_info['spec_name'])
valid = True
if 'spec' not in data._fields:
    valid = False
if 'filter_names' not in data._fields:
    valid = False
else:
    if 'QH0' not in data.filter_names:
        valid = False
    elif 'QH0' in data.filter_names:
        qH0idx = data.filter_names.index('QH0')
if not valid:
    raise IOError("cloudy_slug: error: input slug data must " +
                  "contain spectra and ionizing luminosity")
if args.clustermode and 'form_time' not in data._fields:
    raise IOError("cloudy_slug: error: for cluster mode, "+
                  "input slug data must " +
                  "contain cluster physical properties")
freq = c/(data.wl*1e-8)           # Frequency in Hz
logfreq = np.log10(freq)          # Log frequency in Hz
if args.clustermode:
    if args.end_spec != -1:
        end_spec = min(args.end_spec, len(data.id))
    else:
        end_spec = len(data.id)
else:
    if args.end_spec != -1:
        end_spec = min(args.end_spec, data.spec.shape[-1])
    else:
        end_spec = data.spec.shape[-1]

# Figure out the photometric system we're using
if 'erg/s/Hz' in data.filter_units:
    photsystem = 'L_nu'
elif 'erg/s/A' in data.filter_units:
    photsystem = 'L_lambda'
elif 'AB mag' in data.filter_units:
    photsystem = 'AB'
elif 'ST mag' in data.filter_units:
    photsystem = 'STMAG'
elif 'Vega mag' in data.filter_units:
    photsystem = 'Vega'
else:
    photsystem = 'L_nu'

# Get number of filters
nfilt = len(data.filter_names)


######################################################################
# Step 5: read the template cloudy input file template, and set up
# storage for what we'll be computing
######################################################################

# Parse the template to see what will be written out
fp = open(cloudytemplate, 'r')
tempfile = fp.read().split('\n')
fp.close()
compute_hcon = False
compute_hion = False
compute_continuum = False
compute_lines = False
for line in tempfile:
    # Trim comments
    ltrim = line.partition(';')[0]
    ltrim = line.partition('//')[0]
    ltrim = line.partition('#')[0]
    ltrim = line.partition('%')[0]
    # Skip lines that start with c, indicating a comment
    if len(line) > 0:
        if line[0] == 'c':
            continue
    if 'save last continuum' in line:
        compute_continuum = True
    elif 'save last line list emergent absolute column' in line:
        compute_lines = True
    elif 'save last hydrogen conditions' in line:
        compute_hcon = True
    elif 'save last hydrogen ionization' in line:
        compute_hion = True

# Prepare storage
if compute_continuum:
    cloudywl = []
    cloudyspec = []
    cloudyphot = []
    # Create dummy holders for cloudy spectra
    if args.clustermode:
        for i in range(end_spec-args.start_spec):
            cloudywl.append(None)
            cloudyspec.append(None)
            cloudyphot.append(None)
    else:
        for j in range(data.spec.shape[-2]):
            cloudywl.append([])
            cloudyspec.append([])
            cloudyphot.append([])
            for i in range(args.start_spec, end_spec):
                cloudywl[j].append(None)
                cloudyspec[j].append(None)
                cloudyphot[j].append(None)
if compute_lines:
    linelist = None
    linewl = None
    linelum = []
    # Create dummy holders for cloudy lines
    if args.clustermode:
        for i in range(end_spec-args.start_spec):
            linelum.append(None)
    else:
        for j in range(data.spec.shape[-2]):
            linelum.append([])
            for i in range(args.start_spec, end_spec):
                linelum[j].append(None)
cloudy_params = []
if args.clustermode:
    for i in range(end_spec-args.start_spec):
        cloudy_params.append(None)
else:
    for j in range(data.spec.shape[-2]):
        cloudy_params.append([])
        for i in range(args.start_spec, end_spec):
            cloudy_params[j].append(None)

            
######################################################################
# Step 6: queue up the SLUG runs
######################################################################
slug_queue = Queue()
if args.clustermode:
    for i in range(args.start_spec, end_spec):
        slug_queue.put(i)
else:
    for i in range(args.start_spec, end_spec):
        for j in range(data.spec.shape[-2]):
            slug_queue.put((i, j))


######################################################################
# Step 7: define the worker function to do a cloudy run
######################################################################
def do_cloudy_run(thread_num, thread_return, q):

    # Declare global variables
    global compute_continuum
    global compute_lines
    global compute_hcon
    global compute_hion
    global data
    global basename
    global tmpdirname
    global cloudy_params
    if compute_continuum:
        global cloudywl
        global cloudyspec
        global cloudyphot
    if compute_lines:
        global linelist
        global linewl
        global linelum

    # Initialize the return value to None
    thread_return[thread_num] = None
        
    # Terminate when queue empties
    while not q.empty():

        # Do entire thread inside a try... except block, so that if
        # any unanticipated error occurs we can exit the thread
        # gracefully and the caller won't hang waiting for it
        try:

            ##########################################################
            # Step 7a: fetch a task from the queue, get the data we
            # need, and generate a file name extension to go with it
            ##########################################################
            if args.clustermode:

                # Cluster mode
                cluster_num = q.get()

                # Get data and construct output file name
                spec = data.spec[cluster_num,:]
                ext = "_n{0:09d}".format(cluster_num)
                qH0 = data.phot[cluster_num,qH0idx]
                outstr = "launching cloudy on cluster {0:d} of {1:d}" \
                         .format(cluster_num+1, len(data.id))

                # Check if this cluster is above our maximum age or below
                # our minimum ionizing luminosity; if so, flag it
                skip = False
                if data.time[cluster_num]-data.form_time[cluster_num] > \
                   args.agemax*1e6*365.25*24.*3600.:
                    skip = True
                    outstr = "skipping cluster {:d} of {:d} (age > {:f} Myr)" \
                             .format(cluster_num+1, len(data.id), args.agemax)
                if qH0 < args.qH0min:
                    skip = True
                    outstr = "skipping cluster {:d} of {:d} (QH0 < {:e} s^-1)" \
                             .format(cluster_num+1, len(data.id), args.qH0min)
            
            else:

                # Integrated mode
                trial, time = q.get()
                spec = data.spec[:, time, trial]
                ext = "_tr{0:05d}_ti{1:05d}".format(trial,time)
                qH0 = data.phot[qH0idx,time,trial]
                outstr = ("launching cloudy on trial {:d} of {:d}, " +
                          "time {:d} of {:d}").\
                          format(trial+1, data.spec.shape[-1],
                                 time+1, data.spec.shape[-2])

                # Check if this time is below our minimum ionizing
                # luminosity; if so, flag it
                skip = False
                if qH0 < args.qH0min:
                    skip = True
                    outstr = ("skipping trial {:d} of {:d}, " +
                              "time {:d} of {:d} (QH0 < {:e})").\
                              format(trial+1, data.spec.shape[-1],
                                     time+1,
                                     data.spec.shape[-2], args.qH0min)

            ##########################################################
            # Step 7b: if we have command line arguments that set the
            # physical conditions, use them to derive the inner radius
            # and density we need to give cloudy. Conversely, if we do
            # not have command line parameters, and we're getting
            # density and radius from the template, create an object
            # we will use below to figure out what the other
            # parameters are.
            ##########################################################
            if cmd_complete:

                # Data are coming from command line
            
                # For arguments that are PDFs, draw from them;
                # otherwise convert to float
                if args.hden is not None:
                    try:
                        hden = float(args.hden)
                    except ValueError:
                        pdf = slug_pdf(args.hden)
                        hden = pdf.draw()
                else:
                    hden = None
                if args.r0 is not None:
                    try:
                        r0 = float(args.r0)
                    except ValueError:
                        pdf = slug_pdf(args.r0)
                        r0 = pdf.draw()
                else:
                    r0 = None
                if args.r1 is not None:
                    try:
                        r1 = float(args.r1)
                    except ValueError:
                        pdf = slug_pdf(args.r1)
                        r1 = pdf.draw()
                else:
                    r1 = None
                if args.ionparam is not None:
                    try:
                        U = float(args.ionparam)
                    except ValueError:
                        pdf = slug_pdf(args.ionparam)
                        U = pdf.draw()
                else:
                    U = None
                if args.ionparam0 is not None:
                    try:
                        U0 = float(args.ionparam0)
                    except ValueError:
                        pdf = slug_pdf(args.ionparam0)
                        U0 = pdf.draw()
                else:
                    U0 = None
                if args.windparam is not None:
                    try:
                        Omega = float(args.windparam)
                    except ValueError:
                        pdf = slug_pdf(args.windparam)
                        Omega = pdf.draw()
                else:
                    Omega = None
                if args.coveringfac is not None:
                    try:
                        coveringfac = float(args.coveringfac)
                    except ValueError:
                        pdf = slug_pdf(args.coveringfac)
                        coveringfac = pdf.draw()
                else:
                    coveringfac = 1.0

                # Compute the inner radius and density we'll be giving
                # to cloudy
                if not args.dynamic:
                    
                    # Apply parameter safety check as instructed
                    try:
                        hp = hiiregparam(
                            qH0, nII=hden, r0=r0, r1=r1,
                            U=U, U0=U0, Omega=Omega,
                            warn=(args.paramsafety == 'warn'),
                            fix_quantity=args.fix_quantity)
                        if args.ionparammax is not None:
                            if hp.U0 > args.ionparammax:
                                if args.paramsafety == 'warn':
                                    hp = hiiregparam(
                                        qH0, nII=hden, r0=r0, r1=r1,
                                        U=U, U0=args.ionparammax,
                                        Omega=Omega)
                                else:
                                    raise ValueError
                        if args.ionparammin is not None:
                            if hp.U0 < args.ionparammin:
                                if args.paramsafety == 'warn':
                                    hp = hiiregparam(
                                        qH0, nII=hden, r0=r0, r1=r1,
                                        U=U, U0=args.ionparamin,
                                        Omega=Omega,
                                        fix_quantity=args.fix_quantity)
                                else:
                                    raise ValueError
                    except ValueError as ve:
                        # If we're here, the parameter combination is
                        # physically not allowed, and the safety setting
                        # is skip, halt, or redraw; do whichever we
                        # were told to do
                        if args.paramsafety == 'redraw':
                            ntries = 1
                            while ntries < 100:
                                if args.hden is not None:
                                    try:
                                        hden = float(args.hden)
                                    except ValueError:
                                        pdf = slug_pdf(args.hden)
                                        hden = pdf.draw()
                                else:
                                    hden = None
                                if args.r0 is not None:
                                    try:
                                        r0 = float(args.r0)
                                    except ValueError:
                                        pdf = slug_pdf(args.r0)
                                        r0 = pdf.draw()
                                else:
                                    r0 = None
                                if args.r1 is not None:
                                    try:
                                        r1 = float(args.r1)
                                    except ValueError:
                                        pdf = slug_pdf(args.r1)
                                        r1 = pdf.draw()
                                else:
                                    r1 = None
                                if args.ionparam is not None:
                                    try:
                                        U = float(args.ionparam)
                                    except ValueError:
                                        pdf = slug_pdf(args.ionparam)
                                        U = pdf.draw()
                                else:
                                    U = None
                                if args.ionparam0 is not None:
                                    try:
                                        U0 = float(args.ionparam0)
                                    except ValueError:
                                        pdf = slug_pdf(args.ionparam0)
                                        U0 = pdf.draw()
                                else:
                                    U = None
                                if args.windparam is not None:
                                    try:
                                        Omega = float(args.windparam)
                                    except ValueError:
                                        pdf = slug_pdf(args.windparam)
                                        Omega = pdf.draw()
                                else:
                                    Omega = None
                                try:
                                    hp = hiiregparam(
                                        qH0, nII=hden, r0=r0, r1=r1,
                                        U=U, Omega=Omega,
                                        warn=(args.paramsafety == 'warn'),
                                        fix_quantity=args.fix_quantity)
                                    break
                                except ValueError as ve:
                                    ntries += 1
                        elif args.paramsafety == 'skip' or ntries == 100:
                            skip = True
                            hp = None
                            outstr \
                                = "skipping cluster {:d} of {:d}; "\
                                .format(cluster_num+1, len(data.id))
                            outstr += "(unphysical parameter "+ \
                                      "combination "
                            nset = 0
                            if hden is not None:
                                outstr += "hden = {:e}".format(hden)
                                if nset == 0: outstr += ", "
                                nset += 1
                            if r0 is not None:
                                outstr += "r0 = {:e}".format(r0)
                                if nset == 0: outstr += ", "
                                nset += 1
                            if r1 is not None:
                                outstr += "r1 = {:e}".format(r1)
                                if nset == 0: outstr += ", "
                                nset += 1
                            if U is not None:
                                outstr += "U = {:e}".format(U)
                                if nset == 0: outstr += ", "
                                nset += 1
                            if U0 is not None:
                                outstr += "U0 = {:e}".format(U0)
                                if nset == 0: outstr += ", "
                                nset += 1
                            if Omega is not None:
                                outstr += "Omega = {:e}".format(Omega)
                                if nset == 0: outstr += ", "
                                nset += 1
                            outstr += ", Q(H0) = {:e})".format(qH0)
                        else:
                            # Re-raise in this case
                            raise ValueError(str(ve))
                else:
                    # Note that the combination Omega and n0 is always
                    # valid, so need to check for safety
                    hp = hiiregparam(qH0, Omega=Omega, n0=hden,
                                     t=data.time[cluster_num] -
                                     data.form_time[cluster_num])

            else:

                # Data are coming from template
                hp = hiiregparam(qH0, r0=radius_template,
                                 nII=hden_template)

                
            ##########################################################
            # Step 7c. Write out the cloudy input file header,
            # substituting a custom name for OUTPUT_FILENAME; omit
            # hden and radius lines if we're overriding these from the
            # command line; extract the names of the output files
            # we'll be reading back later
            ##########################################################
            if not skip:
                cloudy_in_fname = osp.join(tmpdirname, 'cloudy.in'+ext)
                fpout = open(cloudy_in_fname, 'w')
            continuum_file = None
            lines_file = None
            hcon_file = None
            for line in tempfile:
                linesplit = line.split()
                if len(linesplit) > 0:
                    if linesplit[0] == 'hden' and cmd_complete:
                        continue
                    elif linesplit[0] == 'radius' and cmd_complete:
                        continue
                    elif linesplit[1] == 'covering' and \
                         args.converingfac is not None:
                        continue
                    elif '$OUTPUT_FILENAME' in line:
                        if not skip:
                            newline \
                                = line.replace('$OUTPUT_FILENAME',
                                               osp.join(tmpdirname,
                                                        basename+ext))
                            fpout.write(newline + '\n')
                            if 'continuum' in newline:
                                lquote = newline.find('"')
                                rquote = newline.rfind('"')
                                continuum_file = newline[lquote+1:rquote]
                            elif 'line list' in newline or \
                                 'linelist' in newline:
                                lquote = newline.find('"')
                                rquote = newline[lquote+1:].find('"')
                                lines_file = newline[lquote+1:lquote+1+rquote]
                            elif 'hydrogen conditions' in newline:
                                lquote = newline.find('"')
                                rquote = newline[lquote+1:].find('"')
                                hcon_file = newline[lquote+1:lquote+1+rquote]
                    else:
                        if not skip:
                            fpout.write(line+'\n')

            # Write inner radius and density to cloudy input file if
            # necessary
            if not skip and cmd_complete:
                fpout.write("hden {0:f}\n".format(np.log10(hp.nII)))
                fpout.write("radius {0:f}\n".format(np.log10(hp.r0)))

            # Write the ionizing luminosity and covering factor to the
            # cloudy input file
            if not skip:
                fpout.write("Q(H) = {0:f}\n".format(np.log10(qH0)))
                if args.coveringfac is not None:
                    fpout.write("covering factor {:f}\n".format(coveringfac))

            ##########################################################
            # Step 7d: write the spectral shape into the cloudy input
            # file. Cloudy is very particular about the format of
            # these data, so this has to be done with a bit of care to
            # correctly handle the following issues:
            # 1. We need to prepend and append low fluxes for
            # frequencies outside the range covered by SLUG.
            # 2. We need to make sure to handle zero fluxes
            # gracefully; cloudy inputs are logarithmic, and we can't
            # give it NaN's.
            # 3. We need to check to make sure that we don't put the
            # same frequency twice in the cloudy input. This can occur
            # if the SLUG output is in ASCII, and does not contain
            # enough digits of precision for two adjacent wavelengths
            # to be separated. Should never occur for SLUG output in
            # FITS or binary format.
            # 4. We have to write out the frequencies to an absurd
            # number of digits of precision, to prevent cases where
            # two adjacent frequencies are close enough that they
            # would appear identical when printed out as ASCII in the
            # cloudy input file.
            ##########################################################
            if not skip:
                specclean = np.copy(spec)
                specclean[spec == 0.0] = np.amin(spec[spec > 0])*1e-4
                logL_nu = np.log10(specclean*c/freq**2)
                fpout.write("interpolate")
                fpout.write(" ({0:15.12f} {1:f})".
                            format(7.51, np.amin(logL_nu)-4))
                fpout.write(" ({0:15.12f} {1:f})".
                            format(logfreq[-1]-0.01, np.amin(logL_nu)-4))
                for i in range(len(logL_nu)):
                    if i % 4 == 0:
                        fpout.write("\ncontinue")
                    if logfreq[-i-1] == logfreq[-(i-1)-1]:
                        continue
                    fpout.write(" ({0:15.12f} {1:f})".
                                format(logfreq[-i-1], logL_nu[-i-1]))
                fpout.write("\ncontinue ({0:15.12f} {1:f})".
                            format(logfreq[0]+0.01, 
                                   np.amin(logL_nu)-4))
                fpout.write(" ({0:f} {1:f})\n".
                            format(22.4, np.amin(logL_nu)-4))

                # Close cloudy input file
                fpout.close()

            ##########################################################
            # Step 7e: record the run parameters
            ##########################################################
            if hp is not None:
                params = { 'hden' : hp.nII, 'r0' : hp.r0,
                           'r1' : hp.r1, 'QH0' : qH0, 'U' : hp.U,
                           'U0' : hp.U0, 'Omega' : hp.Omega,
                           'zeta' : hp.zeta() }
            else:
                # Special case where the parameter values we got are
                # invalid, so we can't construct an hp object, and
                # we're skipping this run. We'll enter 0 for all
                # values except the ones we were given
                params = { 'QH0' : qH0, 'zeta' : 0.0 }
                if hden is not None: params['hden'] = hden
                else: params['hden'] = 0.0
                if r0 is not None: params['r0'] = r0
                else: params['r0'] = 0.0
                if r1 is not None: params['r1'] = r1
                else: params['r1'] = 0.0
                if U is not None: params['U'] = U
                else: params['U'] = 0.0
                if U0 is not None: params['U0'] = U0
                else: params['U0'] = 0.0
                if Omega is not None: params['Omega'] = Omega
                else: params['Omega'] = 0.0
            if args.coveringfac is not None:
                params['coveringfac'] = coveringfac
            elif coveringfac_template is not None:
                params['coveringfac'] = coveringfac_template
            else:
                params['coveringfac'] = 1.0
            if args.clustermode:
                params['id'] = data.id[cluster_num]
                params['trial'] = data.trial[cluster_num]
                params['time'] = data.time[cluster_num]
            else:
                params['time'] = data.time[time]
            
            ##########################################################
            # Step 7f:launch the cloudy process and wait for it to
            # complete
            ##########################################################
            if args.verbose:
                print("thread {0:d}: ".format(thread_num+1) + outstr)
            if not skip:
                cloudy_out_fname = osp.join(tmpdirname, 'cloudy.out'+ext)
                cmd = cloudypath + " < " + cloudy_in_fname + \
                      " > " + cloudy_out_fname
                if args.nicelevel > 0:
                    cmd = "nice -n " + str(args.nicelevel) + " " + cmd
                proc = subprocess.call(cmd, shell=True)

            ##########################################################
            # Step 7g: read and store cloudy outputs
            ##########################################################

            # Read and store the cloudy continuum output
            if continuum_file is not None:
                while os.stat(continuum_file).st_size == 0:
                    sleep(2)
                try:
                    cdata = read_cloudy_continuum(continuum_file, r0=hp.r0)
                except:
                    errstr =  "cloudy_slug: failed to read continuum " \
                        "file " + continuum_file + "; suspect " \
                        "cloudy has crashed"
                    if args.warnfail:
                        warnings.warn(errstr+ "; cloudy_slug exiting")
                    raise IOError(errno.EIO, errstr)
                if args.clustermode:
                    cloudywl[cluster_num] = cdata.wl
                    cloudyspec[cluster_num] = cdata.L_lambda
                else:
                    cloudywl[time][trial] = cdata.wl
                    cloudyspec[time][trial] = cdata.L_lambda

            # Read and store the cloudy line luminosity output
            if lines_file is not None and compute_lines:
                while os.stat(lines_file).st_size == 0:
                    sleep(2)
                try:
                    ldata = read_cloudy_linelist(lines_file)
                except:
                    errstr =  "cloudy_slug: failed to read lines " \
                        "file " + lines_file + "; suspect " \
                        "cloudy has crashed"
                    if args.warnfail:
                        warnings.warn(errstr+ "; cloudy_slug exiting")
                    raise IOError(errno.EIO, errstr)
                if linelist is None:
                    linelist = ldata.label
                if linewl is None:
                    linewl = ldata.wl
                if args.clustermode:
                    linelum[cluster_num] = ldata.lum
                else:
                    linelum[time][trial] = ldata.lum

            # Read and store cloudy physical conditions output
            if hcon_file is not None:
                while os.stat(hcon_file).st_size == 0:
                    sleep(2)
                try:
                    r1_out, nII_out, Omega_out \
                        = read_cloudy_hcon(hcon_file,
                                           r0 = hp.r0)
                except:
                    errstr =  "cloudy_slug: failed to read hcon " \
                        "file " + hcon_file + "; suspect " \
                        "cloudy has crashed"
                    if args.warnfail:
                        warnings.warn(errstr+ "; cloudy_slug exiting")
                    raise IOError(errno.EIO, errstr)
                params['r1_out'] = r1_out
                params['hden_out'] = nII_out
                params['Omega_out'] = Omega_out
                params['zeta_out'] = hp.rch() / r1_out

            ##########################################################
            # Step 7h: compute photometry from the cloudy data
            ##########################################################
            if continuum_file is not None:
                trans_phot \
                    = compute_photometry(cdata.wl, cdata.L_lambda[1,:],
                                         data.filter_names, 
                                         photsystem=photsystem,
                                         filter_wl=data.filter_wl,
                                         filter_response=data.filter_response,
                                         filter_beta=data.filter_beta,
                                         filter_wl_c=data.filter_wl_c)
                emit_phot \
                    = compute_photometry(cdata.wl, cdata.L_lambda[2,:],
                                         data.filter_names, 
                                         photsystem=photsystem,
                                         filter_wl=data.filter_wl,
                                         filter_response=data.filter_response,
                                         filter_beta=data.filter_beta,
                                         filter_wl_c=data.filter_wl_c)
                trans_emit_phot \
                    = compute_photometry(cdata.wl, cdata.L_lambda[3,:],
                                         data.filter_names, 
                                         photsystem=photsystem,
                                         filter_wl=data.filter_wl,
                                         filter_response=data.filter_response,
                                         filter_beta=data.filter_beta,
                                         filter_wl_c=data.filter_wl_c)
                phot = np.array([trans_phot, emit_phot,
                                 trans_emit_phot])
                if args.clustermode:
                    cloudyphot[cluster_num] = phot
                else:
                    cloudyphot[time][trial] = phot
                
            ##########################################################
            # Step 7i: clean up the cloudy output unless requested to keep
            # it
            ##########################################################
            if not args.save:
                if continuum_file is not None:
                    os.remove(continuum_file)
                if lines_file is not None:
                    os.remove(lines_file)
                if hcon_file is not None:
                    os.remove(hcon_file)
                if not skip:
                    os.remove(cloudy_in_fname)
                    os.remove(cloudy_out_fname)

            ##########################################################
            # Step 7j: store the run parameters, and write them out if
            # requested
            ##########################################################
            if args.clustermode:
                cloudy_params[cluster_num] = params
            else:
                cloudy_params[time][trial] = params
            if args.writeparams and args.save:
                param_fname = osp.join(tmpdirname,
                                       'cloudy_slug.param'+ext)
                fpp = open(param_fname, 'w')
                for k in params.keys():
                    fpp.write(str(k) + '    ' + str(params[k]) + '\n')
                fpp.close()
            
            # Declare that we're done
            q.task_done()

        except Exception as e:
            
            # If we're here, something has gone wrong; record the
            # error, mark this thread as done, then re-raise to bail
            # out
            q.task_done()
            errstr \
                = (("thread {:d}: caught exception {:s} on line {:d}".
                    format(thread_num+1, str(e),
                           sys.exc_info()[2].tb_lineno)) + 
                   "; problem parameters were: " + 
                   ("Q(H) = {:e}, hden = {:e}, " +
                    "r0 = {:e}, r1 = {:e}, " +
                    "U = {:e}, U0 = {:e}, Omega = {:e}, " +
                    "zeta = {:e}").format(
                        params['QH0'], params['hden'],
                        params['r0'], params['r1'],
                        params['U'], params['U0'],
                        params['Omega'], params['zeta']))
            thread_return[thread_num] = sys.exc_info()[0](errstr)
            if not args.warnfail: 
                raise sys.exc_info()[0](errstr)
            else:
                warnings.warn(errstr)


######################################################################
# Step 8: start a set of threads to do the job; make sure all threads
# exit succesfully, and bail out if not
######################################################################
try: 
    os.mkdir(osp.join(tmpdirname))  # Temporary working directory
except OSError: pass           # Probably failed because dir exists
if __name__ == '__main__' :
    thread_return = [None]*nproc
    for i in range(nproc):
        p = Thread(target=do_cloudy_run,
                   args=(i, thread_return, slug_queue))
        p.start()
slug_queue.join()
for r in thread_return:
    if r is not None:
        if not args.warnfail:
            raise r
        else:
            sys.exit(0)

######################################################################
# Step 9: take the cloudy output spectra and put them into the SLUG
# format. This involves padding the spectra we get from cloudy to make
# them all the same size, so that the can be converted to arrays; the
# continuum spectra are assumed to follow standard cloudy output
# format, so that they all have the same wavelength spacing and the
# same maximum wavelength, but different minimum wavelengths. Thus we
# are padding the beginnings of the arrays.
######################################################################
if args.verbose:
    print("cloudy computation done, consolidating output...")
if compute_continuum:

    # Get the maximum length wavelength array
    cloudywl_max = np.zeros(0)
    if args.clustermode:
        for i in range(len(cloudywl)):
            if cloudywl[i] is not None:
                if cloudywl[i].shape[0] > cloudywl_max.shape[0]:
                    cloudywl_max = cloudywl[i]
    else:
        for i in range(len(cloudywl)):
            for j in range(len(cloudywl[0])):
                if cloudywl[i][j] is None:
                    continue # skipped trials
                elif cloudywl[i][j].shape[0] > cloudywl_max.shape[0]:
                    cloudywl_max = cloudywl[i][j]

    # Now loop over stored spectra, padding array beginnings
    if args.clustermode:
        for i in range(len(cloudywl)):
            if cloudywl[i] is not None:
                offset = len(cloudywl_max) - len(cloudywl[i])
                if offset > 0:
                    cloudyspec[i] \
                        = np.insert(cloudyspec[i], 0,
                                    np.zeros((offset,4)), axis=1)
            else:
                cloudyspec[i] = np.zeros((4,len(cloudywl_max)))
    else:
        for i in range(len(cloudywl)):
            for j in range(len(cloudywl[0])):
                if cloudywl[i][j] is not None:
                    offset = len(cloudywl_max) - len(cloudywl[i][j])
                    if offset > 0:
                        cloudyspec[i][j] \
                            = np.insert(cloudyspec[i][j], 0,
                                        np.zeros((offset,4)), axis=1)
                else:
                    cloudyspec[i][j] = np.zeros((4,len(cloudywl_max)))

    # Now that we've made the spectra the same size, turn them into an
    # array
    cloudyspec = np.array(cloudyspec)

    # Final step: make a namedtuple to hold the data
    if args.clustermode:
        cloudyspec_type = namedtuple('cluster_cloudyspec',
                                     ['id', 'trial', 'time', 
                                      'cloudy_wl', 'cloudy_inc', 
                                      'cloudy_trans', 'cloudy_emit', 
                                      'cloudy_trans_emit'])
        cloudyspec_data \
            = cloudyspec_type(data.id, data.trial, data.time, cloudywl_max,
                              cloudyspec[:,0,:], cloudyspec[:,1,:],
                              cloudyspec[:,2,:], cloudyspec[:,3,:])
    else:
        cloudyspec_type = namedtuple('integrated_cloudyspec',
                                     ['time', 'cloudy_wl', 'cloudy_inc', 
                                      'cloudy_trans', 'cloudy_emit', 
                                      'cloudy_trans_emit'])
        cloudyspec_data \
            = cloudyspec_type(data.time, cloudywl_max,
                              np.transpose(cloudyspec[:,:,0,:], (2,0,1)),
                              np.transpose(cloudyspec[:,:,1,:], (2,0,1)),
                              np.transpose(cloudyspec[:,:,2,:], (2,0,1)),
                              np.transpose(cloudyspec[:,:,3,:], (2,0,1)))

######################################################################
# Step 10: write the cloudy spectra to file
######################################################################
if compute_continuum:
    if args.clustermode:
        write_cluster_cloudyspec(
            cloudyspec_data,
            osp.join(outpath,
                     osp.basename(args.slug_model_name)),
            file_info['format'])
    else:
        write_integrated_cloudyspec(
            cloudyspec_data,
            osp.join(outpath,
                     osp.basename(args.slug_model_name)),
            file_info['format'])

######################################################################
# Step 11: write the line data to file
######################################################################
if compute_lines:
    if args.clustermode:
 
        # Set clusters we skipped to have line luminosities of zero
        nline = 1
        for i in range(len(linelum)):
            if linelum[i] is not None:
                nline = linelum[i].shape[0]
                break
        for i in range(len(linelum)):
            if linelum[i] is None:
                linelum[i] = np.zeros(nline)

        # Write line data
        linelum = np.array(linelum)
        cloudylines_type = namedtuple('cluster_cloudylines',
                                      ['id', 'trial', 'time', 
                                       'cloudy_linelist',
                                       'cloudy_linewl',
                                       'cloudy_linelum'])
        cloudylines = cloudylines_type(data.id, data.trial, data.time, 
                                       linelist, linewl, linelum)
        write_cluster_cloudylines(
            cloudylines,
            osp.join(outpath,
                     osp.basename(args.slug_model_name)),
            file_info['format'])
    else:
        # Deal with outputs where we have set the luminosity to zero
        nline = 1
        for i in range(len(linelum)):
            for j in range(len(linelum[i])):
                if linelum[i][j] is not None:
                    nline = linelum[i][j].shape[0]
                    break
            if nline > 1:
                break
        for i in range(len(linelum)):
            for j in range(len(linelum[i])):
                if linelum[i][j] is None:
                    linelum[i][j] = np.zeros(nline)
                    
        # Write line data
        linelum = np.array(linelum)
        cloudylines_type = namedtuple('cluster_cloudylines',
                                      ['time', 'cloudy_linelabel',
                                       'cloudy_linewl',
                                       'cloudy_linelum'])
        cloudylines = cloudylines_type(data.time, linelist,
                                       linewl, 
                                       np.transpose(linelum, (2,0,1)))
        write_integrated_cloudylines(
            cloudylines,
            osp.join(outpath,
                     osp.basename(args.slug_model_name)),
            file_info['format'])

######################################################################
# Step 12: write photometry to file
######################################################################
if compute_continuum:

    # Cluster or integrated mode
    if args.clustermode:

        # Set clusters we skipped to have photometric values of either
        # 0 (for non-magnitude systems) or +infinity (for magnitude
        # systems)
        nfilter = len(data.filter_names)
        for i in range(len(cloudyphot)):
            if cloudyphot[i] is None:
                cloudyphot[i] = np.zeros((3,nfilter))
                for j in range(nfilter):
                    if 'mag' in data.filter_units[j]:
                        cloudyphot[i][:,j] = np.inf

        # Convert to array
        cloudyphot = np.array(cloudyphot)

        # Build namedtuple
        cloudyphot_type = namedtuple('cluster_cloudyphot',
                                     ['id', 'trial', 'time', 
                                      'cloudy_filter_names', 
                                      'cloudy_filter_units',
                                      'cloudy_filter_wl_eff', 
                                      'cloudy_filter_wl',
                                      'cloudy_filter_response',
                                      'cloudy_filter_beta',
                                      'cloudy_filter_wl_c',
                                      'cloudy_phot_trans',
                                      'cloudy_phot_emit', 
                                      'cloudy_phot_trans_emit'])
        cloudyphot_data \
            = cloudyphot_type(data.id, data.trial, data.time,
                              data.filter_names, data.filter_units,
                              data.filter_wl_eff, data.filter_wl,
                              data.filter_response, data.filter_beta,
                              data.filter_wl_c, cloudyphot[:,0,:],
                              cloudyphot[:,1,:], cloudyphot[:,2,:])

        # Write
        write_cluster_cloudyphot(
            cloudyphot_data,
            osp.join(outpath,
                     osp.basename(args.slug_model_name)),
            file_info['format'])

    else:

        # Set times / trials we skipped to have photometric values of either
        # 0 (for non-magnitude systems) or +infinity (for magnitude
        # systems)
        nfilter = len(data.filter_names)
        for i in range(len(cloudyphot)):
            for j in range(len(cloudyphot[i])):
                if cloudyphot[i][j] is None:
                    cloudyphot[i][j] = np.zeros((3,nfilter))
                    for k in range(nfilter):
                        if 'mag' in data.filter_units[k]:
                            cloudyphot[i][j][:,k] = np.inf
                        
        # Convert to array
        cloudyphot = np.array(cloudyphot)

        # Build namedtuple
        cloudyphot_type = namedtuple('integrated_cloudyphot',
                                     ['time', 
                                      'cloudy_filter_names', 
                                      'cloudy_filter_units',
                                      'cloudy_filter_wl_eff', 
                                      'cloudy_filter_wl',
                                      'cloudy_filter_response',
                                      'cloudy_filter_beta',
                                      'cloudy_filter_wl_c',
                                      'cloudy_phot_trans',
                                      'cloudy_phot_emit', 
                                      'cloudy_phot_trans_emit'])
        cloudyphot_data \
            = cloudyphot_type(data.time, 
                              data.filter_names, data.filter_units,
                              data.filter_wl_eff, data.filter_wl,
                              data.filter_response, data.filter_beta,
                              data.filter_wl_c,
                              np.transpose(cloudyphot[:,:,0,:], (2,0,1)),
                              np.transpose(cloudyphot[:,:,1,:], (2,0,1)),
                              np.transpose(cloudyphot[:,:,2,:], (2,0,1)))

        # Write
        write_integrated_cloudyphot(
            cloudyphot_data,
            osp.join(outpath,
                     osp.basename(args.slug_model_name)),
            file_info['format'])


######################################################################
# Step 13: write cloudy parameters to file
######################################################################
if args.clustermode:
    # Build namedtuple
    fields = ['id', 'trial', 'time',
              'cloudy_hden', 'cloudy_r0', 'cloudy_r1',
              'cloudy_QH0', 'cloudy_covFac',
              'cloudy_U', 'cloudy_U0', 'cloudy_Omega', 'cloudy_zeta']
    values = [
        np.array([p['id'] for p in cloudy_params], dtype=int),
        np.array([p['trial'] for p in cloudy_params], dtype=int),
        np.array([p['time'] for p in cloudy_params]),
        np.array([p['hden'] for p in cloudy_params]),
        np.array([p['r0'] for p in cloudy_params]),
        np.array([p['r1'] for p in cloudy_params]),
        np.array([p['QH0'] for p in cloudy_params]),
        np.array([p['coveringfac'] for p in cloudy_params]),
        np.array([p['U'] for p in cloudy_params]),
        np.array([p['U0'] for p in cloudy_params]),
        np.array([p['Omega'] for p in cloudy_params]),
        np.array([p['zeta'] for p in cloudy_params])
        ]
    # See if any of our runs have the cloudy output parameters; this
    # needs a bit of care, because if we skipped some runs they may
    # not have these fields, but non-skipped ones may. The final
    # object we construct needs to have the fields included if any of
    # the individual runs has them
    if np.any(['r1_out' in p.keys() for p in cloudy_params]):
        fields.append('cloudy_r1_out')
        val = []
        for p in cloudy_params:
            if 'r1_out' in p.keys(): val.append(p['r1_out'])
            else: val.append(0.0)
        values.append(np.array(val))
    if np.any(['hden_out' in p.keys() for p in cloudy_params]):
        fields.append('cloudy_hden_out')
        val = []
        for p in cloudy_params:
            if 'hden_out' in p.keys(): val.append(p['hden_out'])
            else: val.append(0.0)
        values.append(np.array(val))
    if np.any(['Omega_out' in p.keys() for p in cloudy_params]):
        fields.append('cloudy_Omega_out')
        val = []
        for p in cloudy_params:
            if 'Omega_out' in p.keys(): val.append(p['Omega_out'])
            else: val.append(0.0)
        values.append(np.array(val))
    if np.any(['zeta_out' in p.keys() for p in cloudy_params]):
        fields.append('cloudy_zeta_out')
        val = []
        for p in cloudy_params:
            if 'zeta_out' in p.keys(): val.append(p['zeta_out'])
            else: val.append(0.0)
        values.append(np.array(val))
    cloudyparams_type = namedtuple('cluster_cloudyparams', fields)
    cloudyparams_data = cloudyparams_type(*values)
    
    # Write to file
    write_cluster_cloudyparams(
        cloudyparams_data, 
        osp.join(outpath,
                 osp.basename(args.slug_model_name)),
        file_info['format'])
else:
    # Build namedtuple
    fields = ['time',
              'cloudy_hden', 'cloudy_r0', 'cloudy_r1',
              'cloudy_QH0', 'cloudy_covFac',
              'cloudy_U', 'cloudy_U0', 'cloudy_Omega', 'cloudy_zeta']
    values = [
        np.array([[q['time'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['hden'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['r0'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['r1'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['QH0'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['coveringfac'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['U'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['U0'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['Omega'] for q in p]
                  for p in cloudy_params]),
        np.array([[q['zeta'] for q in p]
                  for p in cloudy_params])
        ]
    if np.any([['r1_out' in q.keys() for q in p]
               for p in cloudy_params]):
        fields.append('cloudy_r1_out')
        val1 = []
        for p in cloudy_params:
            val2 = []
            for q in p:
                if 'r1_out' in q.keys(): val2.append(q['r1_out'])
                else: val2.append(0.0)
            val1.append(val2)
        values.append(np.array(val1))
    if np.any([['hden_out' in q.keys() for q in p]
               for p in cloudy_params]):
        fields.append('cloudy_hden_out')
        val1 = []
        for p in cloudy_params:
            val2 = []
            for q in p:
                if 'hden_out' in q.keys(): val2.append(q['hden_out'])
                else: val2.append(0.0)
            val1.append(val2)
        values.append(np.array(val1))
    if np.any([['Omega_out' in q.keys() for q in p]
               for p in cloudy_params]):
        fields.append('cloudy_Omega_out')
        val1 = []
        for p in cloudy_params:
            val2 = []
            for q in p:
                if 'Omega_out' in q.keys(): val2.append(q['Omega_out'])
                else: val2.append(0.0)
            val1.append(val2)
        values.append(np.array(val1))
    if np.any([['zeta_out' in q.keys() for q in p]
               for p in cloudy_params]):
        fields.append('cloudy_zeta_out')
        val1 = []
        for p in cloudy_params:
            val2 = []
            for q in p:
                if 'zeta_out' in q.keys(): val2.append(q['zeta_out'])
                else: val2.append(0.0)
            val1.append(val2)
        values.append(np.array(val1))
    cloudyparams_type = namedtuple('cluster_cloudyparams', fields)
    cloudyparams_data = cloudyparams_type(*values)

    # Write to file
    write_integrated_cloudyparams(
        cloudyparams_data, 
        osp.join(outpath,
                 osp.basename(args.slug_model_name)),
        file_info['format'])

######################################################################
# Step 14: final cleanup
######################################################################
if not args.save:
    if args.verbose:
        print("Cleaning up temporary directory...")
    try:
        os.rmdir(osp.join(tmpdirname))
    except OSError:
        # If we fail, wait 10 seconds and try again; the file system
        # may just need time to catch up
        sleep(10)
        try:
            os.rmdir(osp.join(tmpdirname))
        except OSError:
            warnings.warn("unable to clean up temporary directory "+
                          osp.join(tmpdirname))
