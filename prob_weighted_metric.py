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

import sys
import argparse
import nibabel as nib
import numpy as np
import pandas as pd
from glob import glob


def main():

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Probability weighted metric averages')
    parser.add_argument('-p', '--probatlas', required=True, help="4D probabilistic atlas labels")
    parser.add_argument('-m', '--metrics', required=True, nargs='+', help="Metric image(s) (wildcards supported)")

    # Parse command line arguments
    args = parser.parse_args()
    probatlas_fname = args.probatlas

    # Compile metric filename list
    metric_fname_list = []
    for chunk in args.metrics:
        metric_fname_list += glob(chunk)

    # Load prob atlas
    try:
        print('Loading probabilistic atlas {}'.format(probatlas_fname))
        probatlas_nii = nib.load(probatlas_fname)
        p_all = probatlas_nii.get_fdata()
    except IOError as err:
        print('* Problem loading {} - exiting'.format(probatlas_fname))
        raise (err)

    # Number of atlas labels
    n_labels = p_all.shape[3]
    print('Atlas contains {} labels'.format(n_labels))

    # Atlas voxel volume in mm^3 (microliters)
    vox_vol_um = np.array(probatlas_nii.header.get_zooms()).prod()
    print('Atlas voxel volume is {:0.3f} ul'.format(vox_vol_um))

    # Init stats list
    stats = []

    # Loop over prob labels
    for lc in range(n_labels):

        print()
        print('Label {}'.format(lc))

        p = p_all[:, :, :, lc]

        # Integral of p over image space
        p_tot = np.sum(p)

        for mc, metric_fname in enumerate(metric_fname_list):

            print('  {}'.format(metric_fname))

            # Load this metric image
            try:
                metric_nii = nib.load(metric_fname)
                m = metric_nii.get_fdata()
            except IOError as err:
                print('* Problem loading {} - exiting'.format(metric_fname))
                raise(err)

            # Probability-weighted mean metric for current label
            mp = np.sum(m * p) / p_tot
            label_vol = p_tot * vox_vol_um

            stats.append([lc, label_vol, metric_fname, mp])

    # Convert stats list to dataframe
    df = pd.DataFrame(stats, columns=['Label', 'ProbLabelVolume_ul', 'MetricImageFile', 'ProbWeightedMeanMetric'])

    # Export to CSV
    csv_fname = 'prob_weighted_metrics.csv'

    print()
    print('Saving results to {}'.format(csv_fname))

    df.to_csv(csv_fname, index=False)

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
