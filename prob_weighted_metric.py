#!/usr/bin/env python3
"""
Calculate probability-weighted metrics given a metric image and a probabilistic atlas

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

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
2020 California Institute of Technology.
"""

__version__ = '0.1.0'

import os
import sys
import argparse
import nibabel as nib
import numpy as np


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Probability weighted metric averages')
    parser.add_argument('-p', '--probatlas', required=True, help="4D probabilistic atlas labels")
    parser.add_argument('-m', '--metrics', required=True, nargs='+', help="Space separated list of metric images")

    # Parse command line arguments
    args = parser.parse_args()
    probatlas_fname = args.probatlas
    metric_fname_list = args.metrics

    # Load prob atlas
    try:
        probatlas_nii = nib.load(probatlas_fname)
        p = probatlas_nii.get_fdata()
    except:
        print('* Problem loading {} - exiting'.format(probatlas_fname))
        raise ()

    n_labels = p.shape[3]

    for metric_fname in metric_fname_list:

        # Force absolute path
        metric_fname = os.path.abspath(metric_fname)

        # Atlas voxel volume in mm^3 (microliters)
        atlas_vox_vol_ul = np.array(p_nii.header.get_zooms()).prod()

        # Treat probabilities as partial volumes and integrate
        if nd == 3:

            V = np.sum(p) * atlas_vox_vol_ul
            print('%0.3f' % V)

        elif nd == 4:

            nx, ny, nz, nt = p.shape

            for t in range(0, nt):
                V = np.sum(p[:, :, :, t]) * atlas_vox_vol_ul
                print('%0.3f' % V),

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
