"""
This script downloads and untar's the entire MIST v1
database. Warning: run this where you have plenty of space, because
the MIST data are stored in plain text, which makes them BIG. This
script is to be used together with slime.py to produce the much more
svelte version of the MIST tracks that ships with slug.
"""

import subprocess
import os.path as osp

# Local destination directory
localdir = '/Users/krumholz/Data/mist'

# MIST version tag
ver = 'v1.2'

# Locations of the MIST data
urlbase = 'http://waps.cfa.harvard.edu/MIST/data/tarballs_'+ver+'/'
fnames = [ 'MIST_'+ver+'_feh_ZZZZZ_afe_p0.0_vvcrit0.4_EEPS.txz',
           'MIST_'+ver+'_feh_ZZZZZ_afe_p0.0_vvcrit0.0_EEPS.txz']
zstr = ['m4.00', 'm3.50', 'm3.00', 'm2.50', 'm2.00', 'm1.75',
        'm1.50', 'm1.25', 'm1.00', 'm0.75', 'm0.50', 'm0.25',
        'p0.00', 'p0.25', 'p0.50']

procs = []
for f in fnames:
    for z in zstr:
        fn = f.replace('ZZZZZ', z)
        cmd = "curl "+urlbase+fn+" > "+osp.join(localdir, fn)
        if not osp.isfile(osp.join(localdir, fn)):
            print("Executing "+cmd)
            procs.append(
                (subprocess.Popen(cmd, shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE),
                 osp.join(localdir, fn)))
        else:
            procs.append((None, osp.join(localdir, fn)))

while len(procs) > 0:
    pop = False
    for i, p in enumerate(procs):
        if p[0] is None:
            done = True
        elif p[0].poll() is not None:
            done = True
        else:
            done = False
        if done:
            cmd = "tar -C "+localdir+" -xzf "+p[1]
            print("Executing "+cmd)
            subprocess.Popen(cmd, shell=True)
            pop = True
            pop_idx = i
            break
    if pop:
        procs.pop(i)

