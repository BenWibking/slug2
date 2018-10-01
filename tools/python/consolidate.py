"""
This script takes outputs produced by an unconsolidated run of slug.py
and consolidates them into a single file, as would happen if
noconsolidate was not set.
"""

import argparse
import os
import os.path as osp
import copy
import glob
import sys
import warnings
import numpy as np
try:
    from slugpy import read_cluster, combine_cluster, \
        read_integrated, combine_integrated, write_cluster, \
        write_integrated
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import read_cluster, combine_cluster, \
            read_integrated, combine_integrated, write_cluster, \
            write_integrated
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

# Step 1: grab command line arguments
parser = argparse. \
         ArgumentParser(
             description=
             "Script to consolidate unconsolidate slug.py outputs")
parser.add_argument("name", nargs='?',
                    help="base name(s) of run to be consolidated;"
                    " default is all runs in the target directory")
parser.add_argument("-d", "--dir", default=None,
                    help="directory containing files to be "
                    "consolidated; default is the current directory")
parser.add_argument("--delete", default=False, action="store_true",
                    help="delete original files after consolidation")
parser.add_argument("-f", "--fmt", default=None,
                    help="format for consolidated output; valid "+
                    "options are ascii, bin, fits, and fits2")
parser.add_argument("-cd", "--checkdata", default=False, 
                    action="store_true",
	 	    help="check data integrity")
parser.add_argument("--halt", default=False,
                    action="store_true",
                    help="halt execution on read failure; if this"
                    "option is not set, a warning is issued instead, "
                    "and processing continues")
parser.add_argument("-v", "--verbose", action="store_true",
                    default=False, help="produce verbose output")
args = parser.parse_args()

# Step 2: construct list of runs to be consolidated
if args.dir is None:
    workdir = os.getcwd()
else:
    workdir = args.dir
extension = '_p[0-9][0-9][0-9][0-9][0-9]' \
            '_n[0-9][0-9][0-9][0-9][0-9]_summary.txt'
if args.name is not None:
    runs = args.name
else:
    runs = glob.glob(osp.join(workdir, '*'+extension))
    for i in range(len(runs)):
        runs[i] = runs[i][:-len('_p00000_n00000_summary.txt')]
    runs = list(set(runs))
if args.verbose:
    print("Runs to be consolidated:")
    for runname in runs:
        print("   "+runname)

# Step 3: loop over runs
for basename in runs:

    # Step 3a: grab list of files to consolidate
    files = glob.glob(basename+extension)
    for i in range(len(files)):
        files[i] = files[i][:-len('_summary.txt')]
    files.sort()
    if args.verbose:
        print("Files to be consolidated:")
        for f in files:
            print("   "+f)

    # Step 3b: read all cluster data
    read = []
    data = []
    for f in files:
        clfile = glob.glob(f+'_cluster_*')
        if len(clfile) == 0:
            continue
        if args.verbose:
            print("... reading cluster data for "+f)
        rinfo={}
        try:
            cldata = read_cluster(f, nofilterdata=True, read_info=rinfo) 
            if args.checkdata:
                if np.sum(cldata.target_mass == 0.0) > 0 or \
                   np.sum(cldata.time == 0.0) > 0 or \
                   np.sum(cldata.A_V == 0.0) > 0:
                    errstr = "data for "+f+" appears to be damaged"
                    if args.halt:
                        raise RuntimeError(errstr)
                    else:
                        warnings.warn(errstr+"; processing continues",
                                      RuntimeWarning)
                        continue
            data.append(cldata)
            for k in rinfo.keys():
                if k != 'format':
                    read.append(rinfo[k])
        except IOError as e:
            if args.halt:
                raise IOError(e)
            else:
                warnings.warn("IOError during processing of cluster "
                              "data for "+ f + ":" 
                              + str(e) + "; processing continues", 
                              RuntimeWarning)

    # Step 3c: combine data in memory and write
    if len(data) > 0:
        comb = combine_cluster(data)
        if args.fmt is not None:
            fmt = args.fmt
        else:
            fmt = rinfo['format']
        if args.verbose:
            print("Writing {:d} clusters for {:s}"
                  .format(comb.id.size, basename))
        write_cluster(comb, basename, fmt=fmt)

    # Step 3d: repeat 3b for integrated data
    data = []
    for f in files:
        infile = glob.glob(f+'_integrated_*')
        if len(infile) == 0:
            continue
        if args.verbose:
            print("... reading integrated data for "+f)
        rinfo={}
        try:
            indata = read_integrated(f, nofilterdata=True, read_info=rinfo)
            data.append(indata)
            for k in rinfo.keys():
                if k != 'format':
                    read.append(rinfo[k])
        except IOError as e:
            if args.halt:
                raise IOError(e)
            else:
                warnings.warn("IOError during processing of integrated "
                              "data for " + f +":" 
                              + str(e) + "; processing continues", 
                              RuntimeWarning)
                
    # Step 3e: repeat 3c for integrated data
    if len(data) > 0:
        comb = combine_integrated(data)
        if args.fmt is not None:
            fmt = args.fmt
        else:
            fmt = rinfo['format']
        if args.verbose:
            print("Writing {:d} times, {:d}} trials for {:s}"
                  .format(comb.time.size, data.target_mass.shape[-1], 
                          basename))
        write_integrated(comb, basename, fmt=fmt)

    # Step 3f: combine summary files
    ntrials = 0
    for f in files:
        try:
            fp = open(f+'_summary.txt', 'r')
            read.append(f+'_summary.txt')
            for line in fp:
                linesplit = line.split()
                if linesplit[0] == 'n_trials':
                    ntrials += int(linesplit[1])
            fp.close()
        except IOError as e:
            if args.halt:
                raise IOError(e)
            else:
                warnings.warn("IOError during processing of summary file "
                              "for " + f +":" 
                              + str(e) + "; processing continues", 
                              RuntimeWarning)
 
    print("Writing summary file for "+basename)
    fpout = open(basename+'_summary.txt', 'w')
    fp = open(files[0]+'_summary.txt', 'r')
    for line in fp:
        linesplit = line.split()
        if linesplit[0] == 'model_name':
            fpout.write("model_name           "+osp.basename(basename)+"\n")
        elif linesplit[0] == 'out_dir':
            fpout.write("out_dir              "+workdir+"\n")
        elif linesplit[0] == 'n_trials':
            fpout.write("n_trials             {:d}\n".format(ntrials))
        else:
            fpout.write(line)
    fp.close()
    fpout.close()

    # Step 3g: clean up if requested
    if args.delete:
        for r in read:
            os.remove(r)
