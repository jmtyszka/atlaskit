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
import numpy as np
import pandas as pd
from glob import glob


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Extract Freesurfer aseg volumes')
    parser.add_argument('-f', '--fs_dir', help="Freesurfer subjects directory")
    parser.add_argument('-o', '--out_dir', help="Output volumes directory")

    args = parser.parse_args()

    if args.fs_dir:
        fs_dir = args.fs_dir
    else:
        fs_dir = os.getcwd()

    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_dir = os.path.join(os.getcwd(), 'volumes')

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    subj_vols = []

    # Loop over all subjects in FS directory
    for subj_dir in glob(os.path.join(fs_dir, '*')):

        subj_id = os.path.basename(subj_dir)
        print('Processing subject : {}'.format(subj_id))

        lh_cort_vols = load_fs_cort_vols(subj_dir, 'lh')
        rh_cort_vols = load_fs_cort_vols(subj_dir, 'rh')
        subcort_vols = load_fs_subcort_vols(subj_dir)

        if not (lh_cort_vols.empty or rh_cort_vols.empty or subcort_vols.empty):

            # Concatenate volume dfs and transpose - regions become columns
            vols = pd.concat([lh_cort_vols, rh_cort_vols, subcort_vols], sort=False).T

            # Replace index with subject ID
            vols.index = [subj_id]

            # Append to subject volumes list
            subj_vols.append(vols)

    # Create grand volume dataframe for all subjects
    all_vols = pd.concat(subj_vols)

    # Write grand CSV to output directory
    out_fname = os.path.join(out_dir, 'fs_volumes.csv')
    print('Writing all volumes to {}'.format(out_fname))
    all_vols.to_csv(out_fname, index=False, float_format='%0.1f')

    # Clean exit
    sys.exit(0)


def load_fs_cort_vols(subj_dir, hemi='lh'):

    fname = os.path.join(subj_dir, 'stats', '{}.aparc.a2009s.stats'.format(hemi))

    if os.path.isfile(fname):

        print('Loading {}'.format(fname))

        vols = pd.read_csv(fname,
                           sep=" ",
                           header=None,
                           names=['label', 'vol_mm3'],
                           usecols=[0, 3],
                           skipinitialspace=True,
                           comment='#')

        # Prefix labels with 'Left-' or 'Right-'
        prefix = 'Left-' if hemi == 'lh' else 'Right-'
        vols['label'] = prefix + vols['label']

        # Set label as row index
        vols = vols.set_index(['label'])

    else:

        print('{} not found - skipping'.format(fname))
        vols = pd.DataFrame()

    return vols


def load_fs_subcort_vols(subj_dir):

    fname = os.path.join(subj_dir, 'stats', 'aseg.stats')

    if os.path.isfile(fname):

        print('Loading {}'.format(fname))

        vols = pd.read_csv(fname,
                           sep=" ",
                           header=None,
                           names=['vol_mm3', 'label'],
                           usecols=[3, 4],
                           skipinitialspace=True,
                           comment='#')

        # Set label as row index
        vols = vols.set_index(['label'])

        # Extract eTIV from header comments
        eTIV = np.nan
        with open(fname, 'r') as fd:
            for line in fd.readlines():
                if 'EstimatedTotalIntraCranialVol' in line:
                    chunks = line.split(',')
                    eTIV = float(chunks[3])
                    break

        # Append eTIV to dataframe
        vols.loc['eTIV'] = eTIV

    else:

        print('{} not found - skipping'.format(fname))
        vols = pd.DataFrame()

    return vols


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
