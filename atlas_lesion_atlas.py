#!/usr/bin/env python3
"""
Analyse mutual volume overlap of lesion labels with probabilistic atlas labels
- Outputs relative and absolute volume overlaps
- Overlap volumes relative to both lesion and atlas labels

Usage
----
atlas_lesion_atlas.py
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


import os
import sys
import argparse
import nibabel as nib
import numpy as np

__version__ = '0.1.0'


def main():

    # Construct a command line argument parser
    parser = argparse.ArgumentParser(description='Atlas-based lesion volumetrics')
    parser.add_argument('-l', '--lesion', required=True, help='3D lesion labels')
    parser.add_argument('-a', '--atlas', required=True, help='4D bilateral probabilistic atlas labels')
    parser.add_argument('-lk', '--lesionkey', required=False, help='Lesion label key (ITKSNAP format)')
    parser.add_argument('-ak', '--atlaskey', required=False, help='Atlas label key (ITKSNAP format)')

    # Parse command line arguments
    args = parser.parse_args()
    lesion_fname = args.lesion
    atlas_fname = args.atlas

    if args.lesionkey:
        lesion_keyfname = args.lesionkey
    else:
        lesion_keyfname = []

    if args.atlaskey:
        atlas_keyfname = args.atlaskey
    else:
        atlas_keyfname = []

    # Load lesion label image
    try:
        print('  Loading lesion labels from %s' % lesion_fname)
        lesion_nii = nib.load(lesion_fname)
        lesion = lesion_nii.get_data()

    except:
        print('* Problem loading lesion image')

    # Load probabilistic atlas
    try:
        print('  Loading probabilistic atlas from %s' % atlas_fname)
        atlas_nii = nib.load(atlas_fname)
        atlas = atlas_nii.get_data()
    except:
        print('* Problem loading atlas image')

    # Check that atlas and lesion 3D dimensions match
    lesion_dims = lesion.shape
    atlas_dims = (atlas.shape)[0:3]
    if not np.array_equal(lesion_dims, atlas_dims):
        print('* Lesion and atlas image dimensions do not match - exiting')
        sys.exit(1)

    # Load ITK-SNAP label key lists
    lesion_key = load_key(lesion_keyfname)
    atlas_key = load_key(atlas_keyfname)

    # Remove first element of each key (clear label) - unused in prob atlases
    del lesion_key[0]
    del atlas_key[0]

    # Atlas voxel volume in ul
    vox_mm = np.array(atlas_nii.header.get_zooms())
    vox_ul = vox_mm.prod()

    # Split probalistic atlas into left and right hemisphere labels
    # Update atlas label key accordingly
    atlas_split, atlas_key_split = split_brain(atlas, atlas_key)

    # Init result list
    results = []

    # Iterate over lesion labels
    for l_c, l_label in enumerate(lesion_key):

        l_index, l_name = l_label[0], l_label[7]

        # Extract binary mask for current lesion label
        l_mask = lesion == l_index

        # Current lesion label volume in ul
        l_vol_ul = l_mask.sum() * vox_ul

        print('Processing %s' % l_name)

        lesion_results = []

        # Loop over each atlas label
        for a_c, a_label in enumerate(atlas_key_split):

            a_i, a_name = a_label[0], a_label[7]

            print('  Atlas label %s (%d)' % (a_name, a_i))

            # Extract probability field for current atlas label
            a_prob = atlas_split[:,:,:,a_c]

            # Integrated volume of prob atlas label
            a_vol_ul = a_prob.sum() * vox_ul
            print('    Atlas label volume : %0.1f ul' % a_vol_ul)

            # Voxel-wise multiply lesion label and atlas prob image
            intersect = l_mask * a_prob

            # Lesion-atlas intersection volume in ul
            intersect_vol_ul = intersect.sum() * vox_ul
            print('    Lesion-atlas label intersection volume : %0.1f ul' % intersect_vol_ul)

            # Intersection as a percentage of lesion volume
            l_perc = intersect_vol_ul / l_vol_ul * 100.0
            print('    Intersection as a fraction of lesion voume : %0.1f%%' % l_perc)

            # Lesion-masked atlas volume normalized to atlas label volume
            a_perc = intersect_vol_ul / a_vol_ul * 100.0
            print('    Intersection as a fraction of atlas label volume : %0.1f%%' % a_perc)

            # Populate current row of results list
            lesion_results.append([l_name, a_name, l_vol_ul, a_vol_ul, intersect_vol_ul, l_perc, a_perc])

        results.append(lesion_results)

    # Create HTML results report
    out_dir = os.path.dirname(os.path.abspath(lesion_fname))
    report_results(results, out_dir)

    # Clean exit
    sys.exit(0)


def report_results(results, out_dir):
    """
    Generate HTML plot report of absolute and relative intersection volumes

    Parameters
    ----------
    results: list of results lists
    out_dir: output directory name

    Returns
    -------

    """

    from bokeh.io import output_file, show
    from bokeh.layouts import gridplot
    from bokeh.charts import Bar, defaults
    import bokeh.palettes as bp

    # Init output HTML report page
    out_fname = os.path.join(out_dir, 'lesion_intersection_report.html')
    output_file(out_fname)

    # Init bar chart defaults
    defaults.width = 500

    # 20-element Brewer palette
    pal = bp.d3['Category20b'][20]

    # Init plot list
    plots = []

    # Loop over each lesion
    for lesion_results in results:

        # Convert results list to data table (aka dictionary)
        # Source: results.append([l_name, a_name, l_vol_ul, a_vol_ul, intersect_vol_ul, l_perc, a_perc])

        # Convert nested results list into a data dictionary for Bokeh
        l_labels = []
        a_labels = []
        i_vols = []
        l_perc = []
        a_perc = []

        for atlas_result in lesion_results:
            l_labels.append(atlas_result[0])
            a_labels.append(atlas_result[1])
            i_vols.append(atlas_result[4])
            l_perc.append(atlas_result[5])
            a_perc.append(atlas_result[6])

        lesion_name = l_labels[0]

        print(lesion_name)

        data = {
            'Atlas Label': a_labels,
            'Intersect Vol': i_vols,
            'Lesion Perc': l_perc,
            'Atlas Perc': a_perc
        }

        bar_i_vol = Bar(data,
                        values='Intersect Vol',
                        label='Atlas Label',
                        color='Atlas Label', palette=pal,
                        title="%s : Intersection Volume (ul)" % lesion_name,
                        legend=[])

        bar_l_perc = Bar(data,
                        values='Lesion Perc',
                        label='Atlas Label',
                        color='Atlas Label', palette=pal,
                        title="%s : Intersection / Lesion Volume (%%)" % lesion_name,
                        legend=[])

        bar_a_perc = Bar(data,
                        values='Atlas Perc',
                        label='Atlas Label',
                        color='Atlas Label', palette=pal,
                        title="%s : Intersection / Atlas Label Volume (%%)" % lesion_name,
                        legend=[])

        plots.append([bar_i_vol, bar_l_perc, bar_a_perc])

    show(gridplot(plots))


def load_key(key_fname):
    """
    Parse an ITK-SNAP label key file

    Parameters
    ----------
    key_fname: ITK-SNAP label key filename

    Returns
    -------
    key: List containing ITK-SNAP style label key
    """

    import pandas as pd

    # Import key as a pandas dataframe
    # Note the partially undocumented delim_whitespace flag
    key_df = pd.read_table(key_fname,
                         comment='#',
                         header=None,
                         names=['Index','R','G','B','A','Vis','Mesh','Name'],
                         delim_whitespace=True)

    # Convert dataframe to list
    key = key_df.values.tolist()

    return key


def split_brain(atlas, atlas_key):
    """
    Split bilateral prob atlas into left and right hemisphere labels
    Concatentate left and right labels into single atlas
    Update atlas key accordingly

    Parameters
    ----------
    atlas: 4D numpy float array of prob labels
    atlas_key: dictionary of label keys

    Returns
    -------

    """

    # Get atlas dimensions
    nx, ny, nz, nl = atlas.shape

    # Find sagittal midplane index
    hx = int(nx/2.0)

    # Create left and right hemisphere masks
    left_mask = np.ones([nx, ny, nz])
    left_mask[hx:,:,:] = 0
    right_mask = 1 - left_mask

    # Create copies of the atlas
    atlas_left = atlas.copy()
    atlas_right = atlas.copy()

    # Apply left and right masks to atlas copies
    for l in range(0, nl):
        atlas_left[:,:,:,l] = atlas[:,:,:,l] * left_mask
        atlas_right[:,:,:,l] = atlas[:,:,:,l] * right_mask

    # Concatenate left and right atlases
    atlas_split = np.concatenate((atlas_left, atlas_right), axis=3)

    atlas_key_left = []
    atlas_key_right = []

    for label in atlas_key:

        # Prefix left label names with "L_"
        label_left = list(label)
        label_left[7] = "L_" + label_left[7]
        atlas_key_left.append(label_left)

        # Prefix right label names with "R_"
        label_right = list(label)
        label_right[7] = "R_" + label_right[7]
        label_right[0] += len(atlas_key)
        atlas_key_right.append(label_right)

    # Concatenate left and right label lists
    atlas_key_split = atlas_key_left + atlas_key_right

    return atlas_split, atlas_key_split


# This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()
