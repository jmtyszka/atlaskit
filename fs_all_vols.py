#!/usr/bin/env python
"""
Extract aseg volumes from a set of FS6 processed subjects
- Output an R-compatible CSV summary file
- Include volumes normalized to estimated total intracranial volume (eTIV)

Usage
----
label_volumes.py <FS subjects directory>
label_volumes.py -h

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2017-09-22 JMT From scratch

License
----
This file is part of atlaskit.

    atlaskit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    atlaskit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with atlaskit.  If not, see <http://www.gnu.org/licenses/>.

Copyright
----
2017 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Extract Freesurfer aseg volumes')
    parser.add_argument('-f', '--fs_dir', help="Freesurfer subjects directory")

    args = parser.parse_args()

    if args.fs_dir:
        fs_dir = args.fs_dir
    else:
        fs_dir = os.getcwd()

    # Loop over all subjects in FS directory
    for subj_dir in os.listdir(fs_dir):

        # ASEG stats filename
        aseg_fname = os.path.join(fs_dir, subj_dir, 'stats', 'aseg.stats')

        if os.path.isfile(aseg_fname):

            # Subject ID from directory name
            sid = subj_dir.replace('sub-','')

            # Load volume stats
            vol_stats = load_stats(aseg_fname)


    # Clean exit
    sys.exit(0)


def load_stats(aseg_fname):

    import numpy as np

    # Extract eTIV from comment line
    with open(aseg_fname, 'r') as fd:

        for line in fd.readlines():
            if "EstimatedTotalIntraCranialVol" in line:
                chunks = line.split(',')
                eTIV = float(chunks[3])

    print('eTIV : %0.1f ul' % eTIV)

    # Load all data from non-commented lines in aseg.stats
    stats = np.genfromtxt(aseg_fname,
                          dtype="u4,u4,u4,f8,U32,f8,f8,u4,u4,u4",
                          names=['Index','Seg_ID','N_vox','Vol_ul','Name','Imean','Istd','Imin','Imax','Irng'])

    # Parse stats, create hemisphere column, normalize to eTIV, discard intensity stats
    for label in stats:
        index, id, n, v, name, _, _, _, _, _ = label
        if 'Left' in name:
            name = name.replace('Left','')
            hemi = 'Left'

    return stats


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
