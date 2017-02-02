#!/usr/bin/env python3
"""
Analyse mutual volume overlap of lesion labels with probabilistic atlas labels
- Outputs relative and absolute volume overlaps
- Overlap volumes relative to both lesion and atlas labels

Usage
----
lesion_atlas_analysis.py
    -l <3D lesion label image>
    -a <4D prob atlas image>
    [-lk <lesion label key>]
    [-ak <atlas label key>]

Authors
----
Mike Tyszka, Caltech Brain Imaging Center

Dates
----
2016-10-26 JMT From scratch

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
2016 California Institute of Technology.
"""

__version__ = '0.1.0'

import sys
import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def main():

    # Construct a command line argument parser
    parser = argparse.ArgumentParser(description='Atlas-based lesion volumetrics')
    parser.add_argument('-l', '--lesion', required=True, help='3D multivalued lesion label image')
    parser.add_argument('-a', '--atlas', required=True, help='4D probabilistic atlas label image')
    parser.add_argument('-lk', '--lesionkey', required=False, help='Lesion label key (ITKSNAP format)')
    parser.add_argument('-ak', '--atlaskey', required=False, help='Atlas label key (ITKSNAP format)')

    # Parse command line arguments
    args = parser.parse_args()
    lesion_fname = args.lesion
    atlas_fname = args.atlas

    if args.lesionkey:
        lesion_key = args.lesionkey
    else:
        lesion_key = []

    if args.atlaskey:
        atlas_key = args.lesionkey
    else:
        atlas_key = []

    # Load lesion label image
    try:
        lesion_nii = nib.load(lesion_fname)
        lesion = lesion_nii.get_data()

    except:
        print('* Problem loading lesion image')

    # Load probabilistic atlas
    try:
        atlas_nii = nib.load(atlas_fname)
        atlas = atlas_nii.get_data()
    except:
        print('* Problem loading atlas image')

    # Construct list of unique label values in lesion image
    # Drop first label (assumed background = 0)
    lesions = np.unique(lesion)
    lesions = lesions[1:]
    n_lesions = lesions.size

    # Number of atlas labels (4th dimension size)
    n_atlas_labels = (atlas.shape)[3]
    atlas_labels=range(0,n_atlas_labels)

    print('  Detected %d labels in lesion image' % n_lesions)
    print('  Detected %d labels in atlas image' % n_atlas_labels)

    # Init result arrays
    intersect_lesion = np.zeros([n_lesions, n_atlas_labels])
    intersect_atlas = np.zeros_like(intersect_lesion)

    for l_i, lesion_i in enumerate(lesions):

        print('Lesion label %d' % lesion_i)

        # Extract binary mask for current lesion label
        lesion_mask = lesion == lesion_i

        # Lesion mask volume
        lesion_volume = lesion_mask.sum()

        # Loop over each atlas label
        for a_i, atlas_label_i in enumerate(atlas_labels):

            print('  Atlas label %d' % atlas_label_i)

            # Extract probability field for current atlas label
            atlas_prob = atlas[:,:,:,a_i]

            # Integrated volume of prob atlas label
            atlas_volume = atlas_prob.sum()
            print('    Atlas volume : %0.1f voxels' % atlas_volume)

            # Multiply lesion label and atlas prob image (voxelwise)
            lesion_atlas = lesion_mask * atlas_prob

            # Integrate volume of lesion-masked atlas label
            intersect_volume = lesion_atlas.sum()
            print('    Lesion masked atlas volume : %0.1f voxels' % intersect_volume)

            # Lesion-masked atlas volume normalized to lesion label volume
            intersect_lesion[l_i, a_i] = intersect_volume / lesion_volume * 100.0
            print('    Intersection / lesion : %0.1f%%' % intersect_lesion[l_i, a_i])

            # Lesion-masked atlas volume normalized to atlas label volume
            intersect_atlas[l_i, a_i] = intersect_volume / atlas_volume * 100.0
            print('    Intersection / atlas : %0.1f%%' % intersect_atlas[l_i, a_i])

    # Plot overlap volume arrays
    fig = plt.figure()

    # Hardwire atlas labels and lesion names for now
    atlas_labels = np.array(['La','BL (BLD+BLI)','BM','CEN','CMN','BL (BLV)','ATA','ATA (ASTA)','AAA', 'Amy (Other)'])
    lesions = np.array(['Right Calcification', 'Left Calcification'])
    lesion_inds = range(1, n_lesions+1)

    # Grouped bar chart parameters
    space = 0.3
    width = (1 - space) / n_atlas_labels

    # Intersection / lesion volume
    ax = fig.add_subplot(111)

    for a_i, atlas_label_i in enumerate(atlas_labels):

        vals = intersect_atlas[:, a_i]
        pos = [l_i - (1 - space) / 2. + a_i * width for l_i in lesion_inds]

        ax.bar(pos, vals, width=width, label=atlas_label_i, color=cm.rainbow(float(a_i)/n_atlas_labels))

    # Set x-axis tick labels
    ax.set_xticks(lesion_inds)
    ax.set_xticklabels(lesions)

    ax.set_ylabel("Nucleus Volume Lesioned (Percent))")
    ax.set_xlabel("Lesion")

    # Add a legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper right')

    plt.show()

    # Clean exit
    sys.exit(0)


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
