"""
This is an automator script to run multiple instances of SLUG in order
to speed the computation of large numbers of trials. It requires that
the slug executable be compiled in and in the same directory as this
script. Usage is as follows:

python slug.py PARAMFILE --nproc=NPROC --batchsize=NBATCH

where PARAMFILE is the name of the parameter file for the run, 
NPROC is the number of simultaneous processes to run (defaults to the
number of cores available on the local machine), and NBATCH is the
number of trials to run per SLUG process (defaults to 20, or number of
trials / NPROC is that is smaller).
"""

import argparse
import copy
import multiprocessing
import numpy as np
import os
import os.path as osp
import subprocess
import sys
import threading
import warnings
try:
    from slugpy import *    # If slugpy is already in our path
except ImportError:
    # If import failed, try to find slugPy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugPy")


# Step 1: grab command line arguments and working directory
parser = argparse. \
         ArgumentParser(
             description="Wrapper script to run slug in parallel")
parser.add_argument('paramfile', help="name of slug parameter file")
parser.add_argument('-n', '--nproc', default=None, type=int,
                    help="number of slug processes (default: "+
                    "number of cores)")
parser.add_argument('-b', '--batchsize', default=None, type=int,
                    help="number of trials per slug process "+
                    "(default: ntrial/nproc)")
args = parser.parse_args()
cwd = osp.dirname(osp.realpath(__file__))


# Step 2: suck the slug paramter file into memory
try:
    fp = open(args.paramfile, 'r')
except IOError:
    fp = None
    if not osp.isabs(args.paramfile) and 'SLUG_DIR' in os.environ:
        fname = osp.join(os.environ['SLUG_DIR'], args.paramfile)
        try:
            fp = open(fname, 'r')
        except IOError:
            pass
if fp is None:
    raise IOError("slug error: unable to open file "+args.paramfile)
pfile = fp.read().split('\n')
fp.close()


# Step 3: parse the file to grab the data we need out of it
ntrials = 1
ntrials_line = -1
model_name = 'SLUG_DEF'
model_name_line = -1
if 'SLUG_DIR' in os.environ:
    out_dir = osp.join(os.environ['SLUG_DIR'], 'output')
else:
    out_dir = 'output'
out_dir_line = -1
output_mode = 'ascii'
verbosity = 0
sim_type = 'galaxy'
rng_offset_line = -1
for i, line in enumerate(pfile):
    linesplit = line.split()
    if len(linesplit) == 0:
        continue
    if linesplit[0].lower() == 'n_trials':
        try:
            ntrials = int(linesplit[1])
            ntrials_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'model_name':
        try:
            model_name = linesplit[1]
            model_name_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'out_dir':
        try:
            out_dir = linesplit[1]
            if not osp.isabs(out_dir) and 'SLUG_DIR' in os.environ:
                out_dir = osp.join(os.environ['SLUG_DIR'], out_dir)
            out_dir_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'output_mode':
        try:
            output_mode = linesplit[1]
            output_mode = output_mode.lower()
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'rng_offset':
        try:
            rng_offset_line = i
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'verbosity':
        try:
            verbosity = int(linesplit[1])
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)
    if linesplit[0].lower() == 'sim_type':
        try:
            sim_type = linesplit[1]
        except IndexError, ValueError:
            raise IOError("slug: error: couldn't parse the following"
                          " parameter file line:\n"+line)


# Step 3: set number of processors and batch size
if args.nproc is None:
    nproc = multiprocessing.cpu_count()
else:
    nproc = args.nproc
if args.batchsize is None:
    batchsize = int(np.ceil(float(ntrials)/nproc))
    if batchsize == 0:
        batchsize = 1
else:
    batchsize = args.batchsize


# Step 4: set up the function that will be responsible for displaying
# output; we'll call threads to run this later
def display_out(out, procnum):
    for line in iter(out.readline, ''):
        print("proc {:d}: ".format(procnum) + line.split('\n')[0])
    out.close()


# Step 5: start launching slug processes
if verbosity > 0:
    print(("Starting parallel SLUG runs with {:d} threads, " +
           "{:d} trials per thread").format(nproc, batchsize))
completed_trials = 0
proc_list = [None] * nproc
proc_ctr = [0] * nproc
io_threads = [None] * nproc
err_threads = [None] * nproc
out_names = []
trial_num = []
ON_POSIX = 'posix' in sys.builtin_module_names
try: 
    os.mkdir(osp.join(cwd, 'slug_par_tmp'))  # Temporary working directory
except OSError: pass           # Probably failed because dir exists
while completed_trials < ntrials:

    # Step 5a: see if any running processes have completed; flag those
    # that have, and increment their counters
    proc_avail_list = []
    for i, proc in enumerate(proc_list):
        if proc is None:
            # Process not started, flag as available
            proc_avail_list.append(i)
        elif proc.poll() != None:
            # Process completed, flag as available
            proc_avail_list.append(i)
            proc_ctr[i] = proc_ctr[i] + 1

    # Step 5b: if there is work remaining to be done, and we have
    # processors available, launch jobs on them
    for p in proc_avail_list:

        # Stop if all tasks have already been assigned
        if completed_trials >= ntrials:
            break

        # Generate model name and number of trials for this run
        new_model_name = model_name + "_p{:03d}_n{:03d}". \
                         format(p, proc_ctr[p])
        new_ntrials = min(batchsize, ntrials-completed_trials)

        # Modify number of trials, output name, and output directory
        # in parameter file; also add an offset to the rng
        pfile_tmp = copy.deepcopy(pfile)
        if ntrials_line != -1:
            pfile_tmp[ntrials_line] \
                = "n_trials   {:d}".format(new_ntrials)
        else:
            pfile_tmp.append("n_trials   {:d}".
                             format(new_ntrials))
        if model_name_line != -1:
            pfile_tmp[model_name_line] \
                = "model_name   "+new_model_name
        else:
            pfile_tmp.append("model_name   "+new_model_name)
        if out_dir_line != -1:
            pfile_tmp[out_dir_line] \
                = "out_dir   "+osp.join(cwd, "slug_par_tmp")
        else:
            pfile_tmp.append("out_dir   "+
                             osp.join(cwd, "slug_par_tmp"))
        if rng_offset_line != -1:
            pfile_tmp[rng_offset_line] \
                = "rng_offset   " + str(10000*p)
        else:
            pfile_tmp.append("rng_offset   " + str(10000*p))

        # Record characteristics of this run
        out_names.append(osp.join(cwd, "slug_par_tmp",
                                  new_model_name))
        trial_num.append(new_ntrials)

        # Write temporary parameter file to disk
        pfile_name = osp.join(cwd, "slug_par_tmp", 
                              "slug_par_p{:03d}".format(p))
        fp = open(pfile_name, 'w')
        for line in pfile_tmp:
            fp.write(line+'\n')
        fp.close()

        # Start new process
        cmd = osp.join(cwd, 'slug') + " " + pfile_name
        proc_list[p] \
            = subprocess.Popen(cmd, bufsize=1, shell=True,
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE,
                               close_fds=ON_POSIX)

        # Increment number of trials assigned
        completed_trials = completed_trials + new_ntrials

        # Start threads to monitor stdout and stderr from this process
        io_threads[p] \
            = threading.Thread(target=display_out,
                               args=(proc_list[p].stdout, p))
        io_threads[p].daemon = True
        io_threads[p].start()
        err_threads[p] \
            = threading.Thread(target=display_out,
                               args=(proc_list[p].stderr, p))
        err_threads[p].daemon = True
        err_threads[p].start()


# Step 6: wait for final processes to complete; we do this because we
# don't want to kill the output threads or do file cleanup before all
# processes are done
while True:
    nrunning = len(proc_list)
    for i, proc in enumerate(proc_list):
        if proc is None:
            # Process never started
            nrunning = nrunning - 1
        elif proc.poll() != None:
            # Process completed
            nrunning = nrunning - 1
    if nrunning == 0:
        break

# Step 7: consolidate output files
combined_name = osp.join(out_dir, model_name)
if verbosity > 0:
    print("Completed parallel runs, consolidating outputs to " + 
          combined_name)

# Step 7a: summary files: change the model name, output directory, and
# number of trials in the first output, and delete all the rest
fp = open(out_names[0]+'_summary.txt', 'r')
fpout = open(osp.join(out_dir, model_name+'_summary.txt'), 'w')
for line in fp:
    linesplit = line.split()
    if linesplit[0] == 'model_name':
        fpout.write("model_name           "+model_name+"\n")
    elif linesplit[0] == 'out_dir':
        fpout.write("out_dir              "+out_dir+"\n")
    elif linesplit[0] == 'n_trials':
        fpout.write("n_trials             {:d}\n".format(ntrials))
    else:
        fpout.write(line)
fp.close()
fpout.close()
for f in out_names:
    try: os.remove(f+'_summary.txt')
    except OSError:
        warnings.warn("unable to clean up temporary file "+f)

# Step 7b: integrated files: read data from all files, combine, then
# write back out
if sim_type != 'cluster':
    data = []
    for f in out_names:
        if verbosity > 1:
            print("Reading integrated data from "+f+"...")
        data.append(read_integrated(f, fmt=output_mode,
                                    nofilterdata=True))
    combined_data = combine_integrated(data)
    write_integrated(combined_data, combined_name, fmt=output_mode)

# Step 7c: cluster files; same as integrated files
data = []
for f in out_names:
    if verbosity > 1:
        print("Reading cluster data from "+f+"...")
    data.append(read_cluster(f, fmt=output_mode,
                             nofilterdata=True))
combined_data = combine_cluster(data)
write_cluster(combined_data, combined_name, fmt=output_mode)


# Step 8: clean up remaining temporary files
if verbosity > 0:
    print("Consolidation complete, cleaning up temporary files")

# Delete parameter files
for i in range(nproc):
    pfile_name = osp.join(cwd, "slug_par_tmp", 
                          "slug_par_p{:03d}".format(i))
    try:
        os.remove(pfile_name)
    except OSError:
        warnings.warn("unable to clean up temporary file "+pfile_name)

# Delete output files
if output_mode == 'ascii':
    extension = '.txt'
elif output_mode == 'binary':
    extension = '.bin'
elif output_mode == 'fits':
    extension = '.fits'
for f in out_names:
    try:
        fname = f+'_summary.txt'
        if osp.isfile(fname):
            os.remove(fname)
        fname = f+'_integrated_prop'+extension
        if osp.isfile(fname):
            os.remove(fname)
        fname = f+'_integrated_spec'+extension
        if osp.isfile(fname):
            os.remove(fname)
        fname = f+'_integrated_phot'+extension
        if osp.isfile(fname):
            os.remove(fname)
        fname = f+'_cluster_prop'+extension
        if osp.isfile(fname):
            os.remove(fname)
        fname = f+'_cluster_spec'+extension
        if osp.isfile(fname):
            os.remove(fname)
        fname = f+'_cluster_phot'+extension
        if osp.isfile(fname):
            os.remove(fname)
    except OSError:
        warnings.warn("unable to clean up temporary file "+fname)

# Remove temporary directory
try:
    os.rmdir(osp.join(cwd, "slug_par_tmp"))
except OSError:
    warnings.warn("unable to clean up temporary directory "+
                  osp.join(cwd, "slug_par_tmp"))